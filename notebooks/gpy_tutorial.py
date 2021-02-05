import math
import torch
import numpy as np
import gpytorch
from matplotlib import pyplot as plt

from botorch.models.gpytorch import GPyTorchModel
from botorch.acquisition import ExpectedImprovement
from botorch.optim import optimize_acqf

from ax import ParameterType, RangeParameter, SearchSpace
from ax import SimpleExperiment
from ax.modelbridge import get_sobol
from ax.modelbridge.factory import get_botorch

# Training data is 100 points in [0,1] inclusive regularly spaced
train_x = torch.linspace(0, 1, 4)
# True function is sin(2*pi*x) with Gaussian noise
train_y = torch.sin(train_x * (2 * math.pi)) + torch.randn(train_x.size()) * math.sqrt(0.1)

# We will use the simplest form of GP model, exact inference
class ExactGPModel(gpytorch.models.ExactGP, GPyTorchModel):
    _num_outputs = 1  # to inform GPyTorchModel API
    def __init__(self, train_x, train_y, likelihood):
        super(ExactGPModel, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ConstantMean()
        self.covar_module = gpytorch.kernels.ScaleKernel(gpytorch.kernels.RBFKernel())
        self.to(train_x)  # make sure we're on the right device/dtype

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)


def _get_and_fit_simple_custom_gp(train_x, train_y, **kwargs):
    # initialize likelihood and model
    likelihood = gpytorch.likelihoods.GaussianLikelihood()
    model = ExactGPModel(train_x[0], train_y[0], likelihood)

    training_iter = 50


    # Find optimal model hyperparameters
    model.train()
    likelihood.train()

    for param_name, param in model.named_parameters():
        print(f'Parameter name: {param_name:42} value = {param.item()}')

    hypers = {
        'likelihood.noise': torch.tensor(0.01),
    }
    model.initialize(**hypers)

    # Use the adam optimizer
    # optimizer = torch.optim.Adam([
    #     {"params": model.mean_module.parameters()},
    #     {"params": model.covar_module.parameters()},
    # ], lr=0.1)  # Includes GaussianLikelihood parameters
    optimizer = torch.optim.SGD(model.parameters(),
        lr=0.1
    )
    list(model.likelihood.noise_covar.parameters())[0].requires_grad=False
    # "Loss" for GPs - the marginal log likelihood
    mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)

    for i in range(training_iter):
        # Zero gradients from previous iteration
        optimizer.zero_grad()
        # Output from model
        output = model(train_x)
        # Calc loss and backprop gradients
        loss = -mll(output, train_y)
        loss.backward()
        print('Iter %d/%d - Loss: %.3f   lengthscale: %.3f   noise: %.3f' % (
            i + 1, training_iter, loss.item(),
            model.covar_module.base_kernel.lengthscale.item(),
            model.likelihood.noise.item()
        ))
        optimizer.step()
    return model

def branin(parameterization, *args):
    x1, x2 = parameterization["x1"], parameterization["x2"]
    y = (x2 - 5.1 / (4 * np.pi ** 2) * x1 ** 2 + 5 * x1 / np.pi - 6) ** 2
    y += 10 * (1 - 1 / (8 * np.pi)) * np.cos(x1) + 10
    # let's add some synthetic observation noise
    return {"branin": (y, 0.0)}

search_space = SearchSpace(
    parameters=[
        RangeParameter(
            name="x1", parameter_type=ParameterType.FLOAT, lower=-5, upper=10
        ),
        RangeParameter(
            name="x2", parameter_type=ParameterType.FLOAT, lower=0, upper=15
        ),
    ]
)

exp = SimpleExperiment(
    name="test_branin",
    search_space=search_space,
    evaluation_function=branin,
    objective_name="branin",
    minimize=True,
)

sobol = get_sobol(exp.search_space)
exp.new_batch_trial(generator_run=sobol.gen(5))

for i in range(5):
    print(f"Running optimization batch {i+1}/5...")
    model = get_botorch(
        experiment=exp,
        data=exp.eval(),
        search_space=exp.search_space,
        model_constructor=_get_and_fit_simple_custom_gp,
    )
    batch = exp.new_trial(generator_run=model.gen(1))

# Get into evaluation (predictive posterior) mode
model.eval()
likelihood.eval()

best_value = train_y.max()
EI = ExpectedImprovement(model=model, best_f=best_value)
new_point_analytic, _ = optimize_acqf(
    acq_function=EI,
    bounds=torch.tensor([[0.0] * 1, [1.0] * 1]),
    q=1,
    num_restarts=20,
    raw_samples=100,
    options={},
)
print('new point', new_point_analytic)
# Test points are regularly spaced along [0,1]
# Make predictions by feeding model through likelihood
with torch.no_grad(), gpytorch.settings.fast_pred_var():
    test_x = torch.linspace(0, 1, 51)
    observed_pred = likelihood(model(test_x))

with torch.no_grad():
    # Initialize plot
    f, ax = plt.subplots(1, 1, figsize=(4, 3))

    # Get upper and lower confidence bounds
    lower, upper = observed_pred.confidence_region()
    # Plot training data as black stars
    ax.plot(train_x.numpy(), train_y.numpy(), 'k*')
    # Plot predictive means as blue line
    ax.plot(test_x.numpy(), observed_pred.mean.numpy(), 'b')
    # Shade between the lower and upper confidence bounds
    ax.fill_between(test_x.numpy(), lower.numpy(), upper.numpy(), alpha=0.5)
    ax.set_ylim([-3, 3])
    ax.legend(['Observed Data', 'Mean', 'Confidence'])
    plt.show()

f_preds = model(test_x)
y_preds = likelihood(model(test_x))

f_mean = f_preds.mean
f_var = f_preds.variance
f_covar = f_preds.covariance_matrix
f_samples = f_preds.sample(sample_shape=torch.Size(1000,))