import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pymc3 as pm
import theano.tensor as tt

from pymc3 import Model, Normal, Slice, sample, traceplot
from pymc3.distributions import Interpolated
from scipy import stats
from theano import as_op
import corner
from fastkde import fastKDE
import pylab as PP
plt.style.use("seaborn-darkgrid")
print(f"Running on PyMC3 v{pm.__version__}")

# Initialize random number generator
np.random.seed(93457)


#Generate two random variables dataset (representing 100000 pairs of datapoints)
N = int(2e5)
var1 = 50*np.random.normal(size=N) + 0.1
var2 = 0.01*np.random.normal(size=N) - 300

#Do the self-consistent density estimate
myPDF,axes = fastKDE.pdf(var1,var2)

#Extract the axes from the axis list
v1,v2 = axes

#Plot contours of the PDF should be a set of concentric ellipsoids centered on
#(0.1, -300) Comparitively, the y axis range should be tiny and the x axis range
#should be large
PP.contour(v1,v2,myPDF)
PP.show()

def rbf_kernel(x1, x2, variance = 1):
    return np.exp(-1 * ((x1-x2) ** 2) / (2*variance))

def gram_matrix(xs):
    return [[rbf_kernel(x1,x2) for x2 in xs] for x1 in xs]

xs = np.arange(-1, 1, 0.01)
mean = [0 for x in xs]
gram = gram_matrix(xs)

plt_vals = []
for i in range(0, 50):
    ys = np.random.multivariate_normal(mean, gram)
    plt_vals.extend([xs, ys, "k"])
plt.plot(*plt_vals)
plt.show()
# # A one dimensional column vector of inputs.
# X = np.linspace(0, 1, 10)[:,None]

# # set the seed
# np.random.seed(1)

# n = 100  # The number of data points
# X = np.linspace(0, 10, n)[:, None]  # The inputs to the GP, they must be arranged as a column vector

# # Define the true covariance function and its parameters
# ℓ_true = 1.0
# η_true = 3.0
# cov_func = η_true ** 2 * pm.gp.cov.Matern52(1, ℓ_true)

# # A mean function that is zero everywhere
# mean_func = pm.gp.mean.Zero()

# # The latent function values are one sample from a multivariate normal
# # Note that we have to call `eval()` because PyMC3 built on top of Theano
# f_true = np.random.multivariate_normal(
#     mean_func(X).eval(), cov_func(X).eval() + 1e-8 * np.eye(n), 1
# ).flatten()

# # The observed data is the latent function plus a small amount of IID Gaussian noise
# # The standard deviation of the noise is `sigma`
# σ_true = 2.0
# y = f_true + σ_true * np.random.randn(n)

# ## Plot the data and the unobserved latent function
# fig = plt.figure(figsize=(12, 5))
# ax = fig.gca()
# ax.plot(X, f_true, "dodgerblue", lw=3, label="True f")
# ax.plot(X, y, "ok", ms=3, alpha=0.5, label="Data")
# ax.set_xlabel("X")
# ax.set_ylabel("The true f(x)")
# plt.legend()
# plt.show()

# with pm.Model() as model:
#     ℓ = pm.Gamma("ℓ", alpha=2, beta=1)
#     η = pm.HalfCauchy("η", beta=5)

#     cov = η ** 2 * pm.gp.cov.Matern52(1, ℓ)
#     gp = pm.gp.Marginal(cov_func=cov)

#     σ = pm.HalfCauchy("σ", beta=5)
#     y_ = gp.marginal_likelihood("y", X=X, y=y, noise=σ)

#     mp = pm.find_MAP()

# # new values from x=0 to x=20
# X_new = np.linspace(0, 20, 600)[:, None]

# # add the GP conditional to the model, given the new X values
# with model:
#     f_pred = gp.conditional("f_pred", X_new)

# # To use the MAP values, you can just replace the trace with a length-1 list with `mp`
# with model:
#     pred_samples = pm.sample_posterior_predictive([mp], vars=[f_pred], samples=2000)

#     # plot the results
# fig = plt.figure(figsize=(12, 5))
# ax = fig.gca()

# # plot the samples from the gp posterior with samples and shading
# from pymc3.gp.util import plot_gp_dist

# plot_gp_dist(ax, pred_samples["f_pred"], X_new)

# # plot the data and the true latent function
# plt.plot(X, f_true, "dodgerblue", lw=3, label="True f")
# plt.plot(X, y, "ok", ms=3, alpha=0.5, label="Observed data")

# # axis labels and title
# plt.xlabel("X")
# plt.ylim([-13, 13])
# plt.title("Posterior distribution over $f(x)$ at the observed values")
# plt.legend()

# with model:
#     values = gp.prior('values', X=X_new)

# with pm.Model() as model:
#     # Model definition
#     pass
# with pm.Model() as model:
#     mu = pm.Normal("mu", mu=0, sd=1)
#     obs = pm.Normal("obs", mu=mu, sd=1, observed=np.random.randn(100))



# True parameter values
alpha_true = 5
beta0_true = 7
beta1_true = 13

# Size of dataset
size = 100

# Predictor variable
X1 = np.random.randn(size)
X2 = np.random.randn(size) * 0.2

# Simulate outcome variable
Y = alpha_true + beta0_true * X1 + beta1_true * X2 + np.random.randn(size)

basic_model = Model()

with basic_model:
    # Priors for unknown model parameters
    alpha = Normal("alpha", mu=0, sd=1)
    beta0 = Normal("beta0", mu=12, sd=1)
    beta1 = Normal("beta1", mu=18, sd=1)

    # Expected value of outcome
    mu = alpha + beta0 * X1 + beta1 * X2

    # Likelihood (sampling distribution) of observations
    Y_obs = Normal("Y_obs", mu=mu, sd=1, observed=Y)

    # draw 1000 posterior samples
    trace = sample(1000)

traceplot(trace)

def from_posterior(param, samples):
    smin, smax = np.min(samples), np.max(samples)
    width = smax - smin
    x = np.linspace(smin, smax, 100)
    y = stats.gaussian_kde(samples)(x)

    # what was never sampled should have a small probability but not 0,
    # so we'll extend the domain and use linear approximation of density on it
    x = np.concatenate([[x[0] - 3 * width], x, [x[-1] + 3 * width]])
    y = np.concatenate([[0], y, [0]])
    return Interpolated(param, x, y)

traces = [trace]

for _ in range(2):

    # generate more data
    X1 = np.random.randn(size)
    X2 = np.random.randn(size) * 0.2
    Y = alpha_true + beta0_true * X1 + beta1_true * X2 + np.random.randn(size)

    model = Model()
    with model:
        # Priors are posteriors from previous iteration
        alpha, beta0 = from_posterior("alpha", trace["alpha"])
        alpha = from_posterior("alpha", trace["alpha"])
        beta0 = from_posterior("beta0", trace["beta0"])
        beta1 = from_posterior("beta1", trace["beta1"])

        # Expected value of outcome
        mu = alpha + beta0 * X1 + beta1 * X2

        # Likelihood (sampling distribution) of observations
        Y_obs = Normal("Y_obs", mu=mu, sd=1, observed=Y)

        # draw 10000 posterior samples
        trace = sample(1000)
        traces.append(trace)

print("Posterior distributions after " + str(len(traces)) + " iterations.")
cmap = mpl.cm.autumn
for param in ["alpha", "beta0", "beta1"]:
    plt.figure(figsize=(8, 2))
    for update_i, trace in enumerate(traces):
        samples = trace[param]
        smin, smax = np.min(samples), np.max(samples)
        x = np.linspace(smin, smax, 100)
        y = stats.gaussian_kde(samples)(x)
        y = stats.gaussian_kde([trace['beta0'],trace['beta1']])(x)
        plt.plot(x, y, color=cmap(1 - update_i / len(traces)))
    plt.axvline({"alpha": alpha_true, "beta0": beta0_true, "beta1": beta1_true}[param], c="k")
    plt.ylabel("Frequency")
    plt.title(param)

plt.tight_layout()
plt.show()
parameters = ["alpha", "beta0", "beta1"]
samples2 = np.zeros((len(samples),3))
for update_i, parameter in enumerate(["alpha", "beta0", "beta1"]):
    for trace in [traces[-1]]:
        samples2[:,update_i] = trace[parameters[update_i]]
fig =  corner.corner(samples2,  bins=40, hist_kwargs={'density':True, 'lw':3}, plot_datapoints=False, fill_contours=False,  show_titles=True, \
                        color='#348ABD',  truth_color='k', use_math_test=True,\
                         levels=[0.9], title_kwargs={"fontsize": 12})
plt.show()
print('s')