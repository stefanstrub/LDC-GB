#%%
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rcParams
import scipy
from scipy import misc
import numpy as np
import xarray as xr
from astropy import units as u
import pandas as pd
import time
from copy import deepcopy
import multiprocessing as mp
from tqdm import tqdm

import torch
import gpytorch
from sklearn.metrics import mean_squared_error
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel, Matern, RationalQuadratic, ExpSineSquared, RBF
from sklearn.kernel_approximation import Nystroem
from sklearn import linear_model, pipeline

from pure_sklearn.map import convert_estimator

from botorch.models.gpytorch import GPyTorchModel
from botorch.acquisition.monte_carlo import qExpectedImprovement, qUpperConfidenceBound
from botorch.optim import optimize_acqf

import optuna

import ldc.io.hdf5 as hdfio
from ldc.lisa.noise import get_noise_model
from ldc.lisa.orbits import Orbits
from ldc.lisa.projection import ProjectedStrain
from ldc.common.series import TimeSeries, FrequencySeries, window
import ldc.waveform.fastGB as fastGB
from ldc.common.tools import compute_tdi_snr
from ldc.waveform.waveform import HpHc

from LISAhdf5 import LISAhdf5, ParsUnits
import tdi

import tensorflow.compat.v2 as tf
import tensorflow_probability as tfp
import pymc3 as pm

import corner

# customized settings
plot_parameter = {  # 'backend': 'ps',
    "font.family": "serif",
    "font.serif": "times",
    "font.size": 16,
    "axes.labelsize": "medium",
    "axes.titlesize": "medium",
    "legend.fontsize": "medium",
    "xtick.labelsize": "medium",
    "ytick.labelsize": "medium",
    "grid.color": "k",
    "grid.linestyle": ":",
    "grid.linewidth": 0.5,
    "savefig.dpi": 150,
}

# tell matplotlib about your param_plots
rcParams.update(plot_parameter)
# set nice figure sizes
fig_width_pt = 464.0  # Get this from LaTeX using \showthe\columnwidth
golden_mean = (np.sqrt(5.0) - 1.0) / 2.0  # Aesthetic ratio
ratio = golden_mean
inches_per_pt = 1.0 / 72.27  # Convert pt to inches
fig_width = fig_width_pt * inches_per_pt  # width in inches
fig_height = fig_width * ratio  # height in inches
fig_size = [fig_width, fig_height]
fig_size_squared = [fig_width, fig_width]
rcParams.update({"figure.figsize": fig_size})

# use a GPU if available
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
dtype = torch.float


def Window(tm, offs=1000.0):
    xl = offs
    ind_r = np.argwhere(tm[-1] - tm <= offs)[0][0]
    xr = tm[ind_r]
    kap = 0.005
    winl = 0.5 * (1.0 + np.tanh(kap * (tm - xl)))
    winr = 0.5 * (1.0 - np.tanh(kap * (tm - xr)))
    return winl * winr

np.random.seed(40)

parameters = [
    "Amplitude",
    "EclipticLatitude",
    "EclipticLongitude",
    "Frequency",
    "FrequencyDerivative",
    "Inclination",
    "InitialPhase",
    "Polarization",
]

DATAPATH = "/home/stefan/LDC/Radler/data"
sangria_fn = DATAPATH + "/dgb-tdi.h5"
sangria_fn = DATAPATH + "/LDC1-5_SOBBH_FD_v2.hdf5"
# sangria_fn = DATAPATH + "/LDC1-3_VGB_v2_FD_noiseless.hdf5"
FD5 = LISAhdf5(sangria_fn)
Nsrc = FD5.getSourcesNum()
GWs = FD5.getSourcesName()
print("Found %d GW sources: " % Nsrc, GWs)
### TOD make sure GalBin is there
if GWs[0] != "SOBBH":
    raise NotImplementedError
p = FD5.getSourceParameters(GWs[0])
td = FD5.getPreProcessTDI()
del_t = float(p.get("Cadence"))
Tobs = float(p.get("ObservationDuration"))

dt = del_t

# Build timeseries and frequencyseries object for X,Y,Z
tdi_ts = xr.Dataset(dict([(k, TimeSeries(td[:, n], dt=dt)) for k, n in [["X", 1], ["Y", 2], ["Z", 3]]]))
# tdi_ts = xr.Dataset(dict([(k,TimeSeries(tdi_ts[k][:,1], dt=dt)) for k in ["X", "Y", "Z"]]))
tdi_fs = xr.Dataset(dict([(k, tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))

noise_model = "MRDv1"
Nmodel = get_noise_model(noise_model, np.logspace(-5, -1, 100))
Npsd = Nmodel.psd()

def loglikelihood(pGBs):
    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs, oversample=4, simulator="synthlisa")
    index_low = np.searchsorted(Xs.f, dataX.f[0])
    Xs = Xs[index_low : index_low + len(dataX)]
    Ys = Ys[index_low : index_low + len(dataY)]
    Zs = Zs[index_low : index_low + len(dataZ)]
    diff = np.abs(dataX - Xs.values) ** 2 + np.abs(dataY - Ys.values) ** 2 + np.abs(dataZ - Zs.values) ** 2
    # p1 = -float(np.sum(diff / Sn)*Xs.attrs['df'])/2.0
    p1 = float(np.sum(diff / Sn) * Xs.df) / 2.0
    # p1 = np.exp(p1)
    return -p1


def loglikelihoodratio(pGBs):
    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs, oversample=4, simulator="synthlisa")
    index_low = np.searchsorted(Xs.f, dataX.f[0])
    Xs = Xs[index_low : index_low + len(dataX)]
    Ys = Ys[index_low : index_low + len(dataY)]
    Zs = Zs[index_low : index_low + len(dataZ)]
    Xs.values = Xs - dataX.values
    Ys.values = Ys - dataY.values
    Zs.values = Zs - dataZ.values

    source = dict({"X": Xs, "Y": Ys, "Z": Zs})
    data = dict({"X": dataX, "Y": dataY, "Z": dataZ})
    # SNRd = compute_tdi_snr(source, Nmodel, data=data)["tot2"]
    SNRs = compute_tdi_snr(source, Nmodel)["tot2"]
    # p1 = SNRd - 0.5*SNRs
    # diff = np.abs(dataX - Xs.values)**2 + np.abs(dataY - Ys.values)**2 + np.abs(dataZ - Zs.values)**2
    # # p1 = -float(np.sum(diff / Sn)*Xs.attrs['df'])/2.0
    # p1 = float(np.sum(diff / (Sn+noise))*Xs.df)/2.0
    # p1 = np.exp(p1)
    return SNRs

def sampler(
    number_of_samples,
    parameters,
    pGB,
    boundaries,
    p1,
    uniform=False,
    MCMC=False,
    only=False,
    onlyparameter="Frequency",
    twoD=False,
    secondparameter="Amplitude",
    calculate_loglikelihood=True,
    gaussian=False,
    pGB_centered=True,
    std=0.2
):
    samples = {}
    pGBs01 = {}
    for parameter in parameters:
        samples[parameter] = []
        if parameter in ["EclipticLatitude"]:
            samples[parameter].append((np.sin(pGB[parameter]) - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0]))
        elif parameter in ["Inclination"]:
            samples[parameter].append((np.cos(pGB[parameter]) - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0]))
        elif parameter in ['Amplitude',"FrequencyDerivative"]:
            samples[parameter].append((np.log10(pGB[parameter]) - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0]))
        else:
            samples[parameter].append((pGB[parameter] - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0]))
        pGBs01[parameter] = samples[parameter][0]
    samples["Likelihood"] = []
    psample = loglikelihood(pGB)
    samples["Likelihood"].append(psample)
    j = 0
    number_of_sampels_sqrt = np.sqrt(number_of_samples)
    for i in range(1, number_of_samples):
        if only:
            parameter = onlyparameter
            if uniform:
                pGBs01[parameter] = i / number_of_samples
            else:
                pGBs01[parameter] = np.random.rand()
            if parameter != "InitialPhase":
                pGBs01["InitialPhase"] = samples["InitialPhase"][0]
            if parameter != "Polarization":
                pGBs01["Polarization"] = samples["Polarization"][0]
        elif twoD:
            parameter = onlyparameter
            parameter2 = secondparameter
            if uniform:
                if i % number_of_sampels_sqrt == 0:
                    j += 1
                pGBs01[parameter] = ((i - 1) % number_of_sampels_sqrt + 1) / (number_of_sampels_sqrt + 2)
                pGBs01[parameter2] = (j + 1) / (number_of_sampels_sqrt + 2)
            else:
                pGBs01[parameter] = np.random.rand()
                pGBs01[parameter2] = np.random.rand()
            # if parameter != 'InitialPhase':
            #     pGBs01['InitialPhase'] = samples['InitialPhase'][0]
        else:
            for parameter in parameters:
                pGBs01[parameter] = np.random.rand()
                if gaussian:
                    center = 0.5
                    if pGB_centered:
                        center = samples[parameter][0]
                    pGBs01[parameter] = scipy.stats.truncnorm.rvs((0-center)/std,(1-center)/std,loc=center,scale=std,size=1)[0]
                # if parameter in ['Amplitude']:#,'FrequencyDerivative','Amplitude','EclipticLongitude']:
                #     pGBs01[parameter] = i/number_of_samples
                # elif parameter in ['FrequencyDerivative']:
                #     pass
                # if parameter in ['Inclination']:
                #     pGBs01[parameter] = np.random.rand()
                # elif parameter in ['Inclination']:
                #     pGBs01[parameter] = np.random.rand()
                # else:
                #     pGBs01[parameter] = np.random.rand()
        for parameter in parameters:
            # if parameter in ['FrequencyDerivative']:
            #     pass
            if parameter in ["EclipticLatitude"]:
                pGBs[parameter] = np.arcsin((pGBs01[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
            elif parameter in ["Inclination"]:
                pGBs[parameter] = np.arccos((pGBs01[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
            elif parameter in ['Amplitude',"FrequencyDerivative"]:
                pGBs[parameter] = 10**((pGBs01[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
            else:
                pGBs[parameter] = (pGBs01[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]
        # pGBs01array = np.zeros(len(parametersfd))
        # i = 0
        # for name in parametersfd:
        #     pGBs01array[i] = pGBs01[name]
        #     i +=1
        # pGBscaled = scaletooriginal(pGBs01array ,boundaries)
        # print(loglikelihood(pGBs), loglikelihood(pGBscaled))
        if calculate_loglikelihood:
            psample = loglikelihood(pGBs)
            samples["Likelihood"].append(psample)
        for parameter in parameters:
            samples[parameter].append(pGBs01[parameter])
    for parameter in parameters:
        samples[parameter] = np.asarray(samples[parameter])
    samples["Likelihood"] = np.asarray(samples["Likelihood"])
    return samples
#%%
class ExactGPModel(gpytorch.models.ExactGP, GPyTorchModel):
    _num_outputs = 1  # to inform GPyTorchModel API

    def __init__(self, train_x, train_y, likelihood, kernel):
        super(ExactGPModel, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ZeroMean()
        kernelA = gpytorch.kernels.RBFKernel(
            active_dims=(0),
            lengthscale_constraint=gpytorch.constraints.GreaterThan(0.2),
        )
        kernelLat = gpytorch.kernels.PeriodicKernel(
            active_dims=torch.tensor([1]),
            lengthscale_constraint=gpytorch.constraints.GreaterThan(0.8),
        )
        kernelLong = gpytorch.kernels.PeriodicKernel(
            active_dims=torch.tensor([2]),
            lengthscale_constraint=gpytorch.constraints.GreaterThan(0.8),
        )
        kernelF = gpytorch.kernels.RBFKernel(
            active_dims=(3),
            lengthscale_constraint=gpytorch.constraints.Interval(0.06, 0.1),
        )
        # kernelFD = gpytorch.kernels.RBFKernel(active_dims=(4))
        kernelI = gpytorch.kernels.PolynomialKernel(power=2, active_dims=(4))
        # kernelIP = gpytorch.kernels.PeriodicKernel(active_dims=torch.tensor([6]))
        kernelIP = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.CosineKernel(
                active_dims=torch.tensor([5]),
                period_length_constraint=gpytorch.constraints.Interval(0.499, 0.501),
            ),
            outputscale_constraint=gpytorch.constraints.Interval(0.9, 1.1),
        )
        # kernelP = gpytorch.kernels.PeriodicKernel(active_dims=torch.tensor([7]))
        # kernelP = gpytorch.kernels.CosineKernel(active_dims=torch.tensor([6]))*gpytorch.kernels.CosineKernel(active_dims=torch.tensor([6]))
        kernelP = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.CosineKernel(
                active_dims=torch.tensor([6]),
                period_length_constraint=gpytorch.constraints.Interval(0.249, 0.251),
            ),
            outputscale_constraint=gpytorch.constraints.Interval(0.9, 1.1),
        )
        # kernel = kernelA + kernelLat + kernelLong + kernelF + kernelFD + kernelI + kernelIP + kernelP
        kernel = kernelA * kernelLat * kernelLong * kernelF * kernelI * kernelIP * kernelP
        # self.covar_module = gpytorch.kernels.ScaleKernel(gpytorch.kernels.RBFKernel(ard_num_dims=8))
        self.covar_module = gpytorch.kernels.RBFKernel(ard_num_dims=8)
        self.to(train_x)  # make sure we're on the right device/dtype

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)

def objective(trial,changeableparameters, maxpGB2):
    for parameter in changeableparameters:
        parametervalue = trial.suggest_uniform(parameter, boundaries[parameter][0], boundaries[parameter][1])
        maxpGB2[parameter] = parametervalue
        if parameter in ["EclipticLatitude"]:
            maxpGB2[parameter] = np.arcsin(parametervalue)
        elif parameter in ["Inclination"]:
            maxpGB2[parameter] = np.arccos(parametervalue)
        elif parameter in ['Amplitude',"FrequencyDerivative"]:
            maxpGB2[parameter] = 10**(parametervalue)
    p = loglikelihood(maxpGB2)
    return p


def model(optimize_parameters):
    i = 0
    for parameter in changeableparameters:
        parametervalue = optimize_parameters[i]
        maxpGB2[parameter] = (parametervalue * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]
        if parameter in ["EclipticLatitude"]:
            maxpGB2[parameter] = troch.arcsin(parametervalue)
        elif parameter in ["Inclination"]:
            maxpGB2[parameter] = torch.arccos(parametervalue)
    p = loglikelihood(maxpGB2)
    i += 1
    return p


def scaletooriginal(previous_max, boundaries):
    maxpGB = deepcopy(pGBs)
    for parameter in parametersfd:
        if parameter in ["EclipticLatitude"]:
            maxpGB[parameter] = np.arcsin((previous_max[parametersfd.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
        elif parameter in ["Inclination"]:
            maxpGB[parameter] = np.arccos((previous_max[parametersfd.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
        elif parameter in ['Amplitude',"FrequencyDerivative"]:
            maxpGB[parameter] = 10**((previous_max[parametersfd.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
        else:
            maxpGB[parameter] = (previous_max[parametersfd.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]
    return maxpGB


def plotplanes(parameterstocheck, parameter2, plot_x, plot_y):
    fig, ax = plt.subplots(2, 4, figsize=(15, 15))
    plt.suptitle("loglikelihood true")
    i = 0
    for parameter in parameterstocheck:
        if parameter != parameter2:
            j = 0
            if i > 3:
                j = 1
            with torch.no_grad():
                ax[j, i % 4].axvline(x=previous_max[parametersfd.index(parameter)], color="k")
                ax[j, i % 4].axhline(y=previous_max[parametersfd.index(parameter2)], color="k")
                im = ax[j, i % 4].scatter(
                    plot_x[parameter + parameter2].numpy()[:, parametersfd.index(parameter)],
                    plot_x[parameter + parameter2].numpy()[:, parametersfd.index(parameter2)],
                    c=plot_y[parameter + parameter2].numpy()[:],
                )
                ax[j, i % 4].set_xlabel(parameter)
                ax[j, i % 4].set_ylabel(parameter2)
                fig.colorbar(im, ax=ax[j, i % 4])
        i += 1
    plt.show()

def traingpmodelsk(train_x, train_y, kernel):
    # train_x = train_x.numpy()
    # train_y = train_y.numpy()
    gpr = GaussianProcessRegressor(kernel=kernel,
            random_state=0).fit(train_x, train_y)
    gpr.score(train_x, train_y)
    return gpr

def traingpmodel(train_x, train_y, kernel, sigma, nu):
    bo_iterations = 1
    for bo_iter in range(bo_iterations):
        # initialize likelihood and model
        likelihood = gpytorch.likelihoods.GaussianLikelihood()
        model = ExactGPModel(train_x, train_y, likelihood, kernel)

        training_iter = 15

        hypers = {
            "likelihood.noise": torch.tensor(0.0001),
        }
        model.initialize(**hypers)
        # Find optimal model hyperparameters
        model.train()
        likelihood.train()

        # Use the adam optimizer
        optimizer = torch.optim.Adam(
            [
                # {"params": model.mean_module.parameters()},
                {"params": model.covar_module.parameters()},
            ],
            lr=0.1,
        )  # Includes GaussianLikelihood parameters
        # optimizer = torch.optim.Adam(model.parameters(), lr=0.1)

        # "Loss" for GPs - the marginal log likelihood
        mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)

        for i in range(training_iter):
            torch.cuda.empty_cache()
            # Zero gradients from previous iteration
            optimizer.zero_grad()
            # Output from model
            output = model(train_x)
            # Calc loss and backprop gradients
            loss = -mll(output, train_y)
            loss.backward()
            print(
                "Iter %d/%d - Loss: %.3f  Len: %.3f Len: %.3f Len: %.3f Len: %.3f Len: %.3f Len: %.3f Len: %.3f Len: %.3f"
                % (
                    i + 1,
                    training_iter,
                    loss.item(),
                    model.covar_module.lengthscale[0][0],
                    model.covar_module.lengthscale[0][1],
                    model.covar_module.lengthscale[0][2],
                    model.covar_module.lengthscale[0][3],
                    model.covar_module.lengthscale[0][4],
                    model.covar_module.lengthscale[0][5],
                    model.covar_module.lengthscale[0][6],
                    model.covar_module.lengthscale[0][7],
                )
            )
            optimizer.step()

        best_value = train_y.max()
        print("best value", bo_iter, (best_value * sigma) + nu, train_x[train_y.argmax()])
        best_value = best_value
        if bo_iter < bo_iterations - 1:
            qEI = qExpectedImprovement(model=model, best_f=best_value)
            qUCB = qUpperConfidenceBound(model=model, beta=0.001)
            UCB = UpperConfidenceBound(model=model, beta=0.001, maximize=True)
            candidates = 1
            new_point_analytic, _ = optimize_acqf(
                acq_function=UCB,
                bounds=torch.tensor([[0.0] * 8, [1.0] * 8]),
                q=candidates,
                num_restarts=1,
                raw_samples=100,
                options={},
            )
            train_x = torch.cat((train_x, new_point_analytic), 0)
            new_point_analytic = new_point_analytic.cpu().numpy()
            for candi_n in range(candidates):
                para_n = 0
                for parameter in parameters:
                    if parameter in ["EclipticLatitude"]:
                        pGBs[parameter] = np.arcsin((new_point_analytic[candi_n][para_n] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
                    elif parameter in ["Inclination"]:
                        pGBs[parameter] = np.arccos((new_point_analytic[candi_n][para_n] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
                    else:
                        pGBs[parameter] = (new_point_analytic[candi_n][para_n] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]
                    para_n += 1
                loglike = loglikelihood(pGBs)
                print("new point", i, new_point_analytic[candi_n], loglike)

                loglike = (loglike - nu) / sigma
                loglike = torch.tensor([loglike]).float()
                train_y = torch.cat((train_y, loglike), 0)
    return model, likelihood


def CoordinateMC(n):
    maxpGB = deepcopy(pGBs)
    parameters_recorded1 = {}
    no_improvement_counter = 0
    np.random.seed(n)
    for parameter in parametersfd:
        parametervalue = np.random.uniform(boundaries[parameter][0], boundaries[parameter][1])
        maxpGB[parameter] = parametervalue
        if parameter in ["EclipticLatitude"]:
            maxpGB[parameter] = np.arcsin(parametervalue)
        elif parameter in ["Inclination"]:
            maxpGB[parameter] = np.arccos(parametervalue)
        elif parameter in ['Amplitude',"FrequencyDerivative"]:
            maxpGB[parameter] = 10**(parametervalue)
        parameters_recorded1[parameter] = []
        parameters_recorded1[parameter].append(maxpGB[parameter])
    parameters_recorded1["Loglikelihood"] = []
    if n == 0:
        maxpGB = deepcopy(pGBs)
    try:
        best_value
    except:
        best_value = p1
    previous_best = loglikelihood(maxpGB)
    maxpGB2 = deepcopy(maxpGB)
    for i in range(100):
        parameter1 = parametersfd[i % 7]
        parameter2 = parametersfd[np.random.randint(0, 6)]
        parameter3 = parametersfd[np.random.randint(0, 6)]
        # parameter2 = 'InitialPhase'
        # parameter1 = 'Inclination'
        while parameter2 == parameter1:
            parameter2 = parametersfd[np.random.randint(0, 6)]
        while parameter3 == parameter1 or parameter3 == parameter2:
            parameter3 = parametersfd[np.random.randint(0, 6)]
        # if parameter1 == 'Frequency':
        #     parameter2 = 'Polarization'
        parametersreduced = [parameter1]
        changeableparameters = [parameter1, parameter2]#, parameter3]
        params = np.zeros(len(changeableparameters))

        optuna.logging.set_verbosity(optuna.logging.WARNING)
        start = time.time()
        study = optuna.create_study(sampler=optuna.samplers.RandomSampler(), direction="maximize")
        study.optimize(lambda trial: objective(trial, changeableparameters, maxpGB2), n_trials=50)
        # print("optuna time", time.time() - start)
        if study.best_value > previous_best:
            no_improvement_counter = 0
            previous_best = study.best_value
            for parameter in changeableparameters:
                maxpGB[parameter] = study.best_params[parameter]
                if parameter in ["EclipticLatitude"]:
                    maxpGB[parameter] = np.arcsin(study.best_params[parameter])
                elif parameter in ["Inclination"]:
                    maxpGB[parameter] = np.arccos(study.best_params[parameter])
                elif parameter in ['Amplitude',"FrequencyDerivative"]:
                    maxpGB[parameter] = 10**(study.best_params[parameter])
        else:
            no_improvement_counter += 1
        if no_improvement_counter > 16:
            past_mean = 0
            sum_count = 0
            for l in range(n):
                s = 0
                try:
                    past_mean += parameters_recorded[l]["Loglikelihood"][i]
                    sum_count += 1
                except:
                    pass
            try:
                past_mean = past_mean / sum_count
                if previous_best > past_mean:
                    no_improvement_counter = 0
                else:
                    # print("no improvement")
                    break
            except:
                # no_improvement_counter = 0
                break

        if i in [20,30,40,50,60,70]:
            past_mean = 0
            sum_count = 0
            for l in range(n):
                s = 0
                try:
                    past_mean += parameters_recorded[l]["Loglikelihood"][i]
                    sum_count += 1
                except:
                    pass
            try:
                past_mean = past_mean / sum_count
                if previous_best > past_mean:
                    pass
                else:
                    break
            except:
                pass

        start = time.time()
        # print(n, i, previous_best, loglikelihood(maxpGB), maxpGB)
        parameters_recorded1["Loglikelihood"].append(loglikelihood(maxpGB))
        for parameter in parametersfd:
            parameters_recorded1[parameter].append(maxpGB[parameter])

        maxpGB2 = deepcopy(maxpGB)
        x = 1
    if previous_best > best_value:
        best_run = n
        best_value = previous_best
        best_params = deepcopy(maxpGB)
    parameters_recorded[n] = parameters_recorded1
    pbar.update(1)
    return parameters_recorded1

def Reduce_boundaries(maxpGB, boundaries, ratio=0.1):
    boundaries_reduced = deepcopy(boundaries)
    for parameter in parameters:
        length = boundaries[parameter][1] - boundaries[parameter][0]
        if parameter == "EclipticLatitude":
            boundaries_reduced[parameter] = [
                np.sin(maxpGB[parameter]) - length * ratio / 2,
                np.sin(maxpGB[parameter]) + length * ratio / 2,
            ]
        # elif parameter == "Inclination":
        #     boundaries_reduced[parameter] = [
        #         np.cos(maxpGB[parameter]) - length * ratio / 2*10,
        #         np.cos(maxpGB[parameter]) + length * ratio / 2*10,
        #     ]
        elif parameter == "Frequency":
            boundaries_reduced[parameter] = [
                maxpGB[parameter] - length * ratio / 16,
                maxpGB[parameter] + length * ratio / 16,
            ]
        elif parameter == "FrequencyDerivative":
            boundaries_reduced[parameter] = [boundaries[parameter][0],-16.5]
            boundaries_reduced[parameter] = [
                np.log10(maxpGB[parameter]) - length * ratio / 2*8,
                np.log10(maxpGB[parameter]) + length * ratio / 2*8,
            ]
        # elif parameter == "Amplitude":
        #     boundaries_reduced[parameter] = [
        #         np.log10(maxpGB[parameter]) - length * ratio / 2*20,
        #         np.log10(maxpGB[parameter]) + length * ratio / 2*20,
        #     ]
        elif parameter in ["InitialPhase",'Polarization']:
            boundaries_reduced[parameter] = [
                maxpGB[parameter] - length * ratio / 8,
                maxpGB[parameter] + length * ratio / 8,
            ]
        elif parameter == "EclipticLongitude":
            boundaries_reduced[parameter] = [
                maxpGB[parameter] - length * ratio / 6,
                maxpGB[parameter] + length * ratio / 6,
            ]
        if boundaries_reduced[parameter][0] < boundaries[parameter][0]:
            boundaries_reduced[parameter][0] = boundaries[parameter][0]
        if boundaries_reduced[parameter][1] > boundaries[parameter][1]:
            boundaries_reduced[parameter][1] = boundaries[parameter][1]
    return boundaries_reduced

def function(pGBs01, boundaries_reduced):
    i = 0
    for parameter in parameters:
        if parameter in ["EclipticLatitude"]:
            pGBs[parameter] = np.arcsin((pGBs01[i] * (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])) + boundaries_reduced[parameter][0])
        elif parameter in ["Inclination"]:
            pGBs[parameter] = np.arccos((pGBs01[i] * (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])) + boundaries_reduced[parameter][0])
        elif parameter in ['Amplitude',"FrequencyDerivative"]:
            pGBs[parameter] = 10**((pGBs01[i] * (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])) + boundaries_reduced[parameter][0])
        else:
            pGBs[parameter] = (pGBs01[i] * (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])) + boundaries_reduced[parameter][0]
        i += 1
    p = -loglikelihood(pGBs)
    return p

class Distribution(object):
    """
    draws samples from a one dimensional probability distribution,
    by means of inversion of a discrete inverstion of a cumulative density function

    the pdf can be sorted first to prevent numerical error in the cumulative sum
    this is set as default; for big density functions with high contrast,
    it is absolutely necessary, and for small density functions,
    the overhead is minimal

    a call to this distibution object returns indices into density array
    """
    def __init__(self, pdf, sort = True, interpolation = True, transform = lambda x: x):
        self.shape          = pdf.shape
        self.pdf            = pdf.ravel()
        self.sort           = sort
        self.interpolation  = interpolation
        self.transform      = transform

        #a pdf can not be negative
        assert(np.all(pdf>=0))

        #sort the pdf by magnitude
        if self.sort:
            self.sortindex = np.argsort(self.pdf, axis=None)
            self.pdf = self.pdf[self.sortindex]
        #construct the cumulative distribution function
        self.cdf = np.cumsum(self.pdf)
    @property
    def ndim(self):
        return len(self.shape)
    @property
    def sum(self):
        """cached sum of all pdf values; the pdf need not sum to one, and is imlpicitly normalized"""
        return self.cdf[-1]
    def __call__(self, N):
        """draw """
        #pick numbers which are uniformly random over the cumulative distribution function
        choice = np.random.uniform(high = self.sum, size = N)
        #find the indices corresponding to this point on the CDF
        index = np.searchsorted(self.cdf, choice)
        pdfs = self.pdf[index]
        #if necessary, map the indices back to their original ordering
        if self.sort:
            index = self.sortindex[index]
        #map back to multi-dimensional indexing
        index = np.unravel_index(index, self.shape)
        index = np.vstack(index)
        #is this a discrete or piecewise continuous distribution?
        if self.interpolation:
            index = index + np.random.uniform(size=index.shape)
        return self.transform(index), pdfs

# Neil2005
# pGB['Amplitude'] = 1.4e-22
# pGB['EclipticLatitude'] = 0.399
# pGB['EclipticLongitude'] = 5.71
# pGB['Frequency'] = 1.0020802e-3
# pGB['FrequencyDerivative'] = 1e-19
# pGB['Inclination'] = 0.96
# pGB['InitialPhase'] = 1.0
# pGB['Polarization'] = 1.3

pGBadded = {}
pGBadded['Amplitude'] = 4*10**-23
pGBadded['EclipticLatitude'] = -0.2
pGBadded['EclipticLongitude'] = 1.4
pGBadded['Frequency'] = 0.00459687
pGBadded['FrequencyDerivative'] = 3*10**-17
pGBadded['Inclination'] = 0.5
pGBadded['InitialPhase'] = 3
pGBadded['Polarization'] = 2
pGBadded2 = {}
pGBadded2['Amplitude'] = 1*10**-22
pGBadded2['EclipticLatitude'] = 0.4
pGBadded2['EclipticLongitude'] = 0.4
pGBadded2['Frequency'] = 0.00459687
pGBadded2['FrequencyDerivative'] = 1*10**-17
pGBadded2['Inclination'] = 0.9
pGBadded2['InitialPhase'] = 1
pGBadded2['Polarization'] = 3


GB = fastGB.FastGB(delta_t=dt, T=Tobs)  # in seconds
Xs_added, Ys_added, Zs_added = GB.get_fd_tdixyz(template=pGBadded, oversample=4, simulator="synthlisa")
source_added = dict({"X": Xs_added, "Y": Ys_added, "Z": Zs_added})
index_low = np.searchsorted(tdi_fs["X"].f, Xs_added.f[0])
index_high = index_low+len(Xs_added)
# tdi_fs['X'] = tdi_fs['X'] #+ Xs_added
for k in ["X", "Y", "Z"]:
    tdi_fs[k].data[index_low:index_high] = tdi_fs[k].data[index_low:index_high] + source_added[k].data

Xs_added, Ys_added, Zs_added = GB.get_fd_tdixyz(template=pGBadded2, oversample=4, simulator="synthlisa")
source_added = dict({"X": Xs_added, "Y": Ys_added, "Z": Zs_added})
index_low = np.searchsorted(tdi_fs["X"].f, Xs_added.f[0])
index_high = index_low+len(Xs_added)
# tdi_fs['X'] = tdi_fs['X'] #+ Xs_added
for k in ["X", "Y", "Z"]:
    tdi_fs[k].data[index_low:index_high] = tdi_fs[k].data[index_low:index_high] + source_added[k].data

tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft(dt=dt)) for k, n in [["X", 1], ["Y", 2], ["Z", 3]]]))

pGB = {}
ind = 0

first_start = time.time()
# for ind in range(1,len(p.get('Frequency'))):
for ind in [8]: #[3,8,9]
    print('index',ind)
    for parameter in parameters:
        pGB[parameter] = p.get(parameter)[ind]
    pGB = pGBadded
    print('pGB', pGB)

    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB, oversample=4, simulator="synthlisa")
    fmin, fmax = float(Xs.f[0]), float(Xs.f[-1] + Xs.attrs["df"])
    source = dict({"X": Xs, "Y": Ys, "Z": Zs})
    f_noise = np.logspace(-5, -1, 100)
    Nmodel = get_noise_model(noise_model, f_noise)
    freq = np.array(source["X"].sel(f=slice(fmin, fmax)).f)
    Sn = Nmodel.psd(freq=freq, option="X")

    parametersfd = [
        "Amplitude",
        "EclipticLatitude",
        "EclipticLongitude",
        "Frequency",
        "FrequencyDerivative",
        "Inclination",
        "InitialPhase",
        "Polarization",
    ]
    boundaries = {
        "Amplitude": [-23.0,-21.4],
        "EclipticLatitude": [-1.0, 1.0],
        "EclipticLongitude": [-np.pi, np.pi],
        # "Frequency": [pGB["Frequency"] * 0.99995, pGB["Frequency"] * 1.00015],
        "Frequency": [pGB["Frequency"] - 3*10e-8, pGB["Frequency"] + 3*10e-8],
        "FrequencyDerivative": [-20.0,-15.0],
        "Inclination": [-1.0, 1.0],
        "InitialPhase": [0.0, 1.0 * np.pi],
        "Polarization": [0.0, 1.0 * np.pi],
    }

    previous_max = np.random.rand(8)
    previous_max[0] = np.random.rand(1)*0.1 +0.9
    previous_max[3] = np.random.rand(1)*0.1 +0.5
    i = 0
    pGBs = deepcopy(pGB)
    for parameter in parameters:
        if parameter in ["FrequencyDerivative"]:
            i -= 1
        elif parameter in ["EclipticLatitude"]:
            pGBs[parameter] = np.arcsin((previous_max[i] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
        elif parameter in ["Inclination"]:
            pGBs[parameter] = np.arccos((previous_max[i] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
        elif parameter in ['Amplitude',"FrequencyDerivative"]:
            pGBs[parameter] = 10**((previous_max[i] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
        else:
            pGBs[parameter] = (previous_max[i] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]
        i += 1


    number_of_samples = 8 * 10 ** 1
    cutoff_ratio = 1000
    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs, oversample=4, simulator="synthlisa")
    psd_signal = np.abs(Xs.values) ** 2 + np.abs(Ys.values) ** 2 + np.abs(Zs.values) ** 2
    highSNR = psd_signal > np.max(psd_signal) / cutoff_ratio
    lowerindex = np.where(highSNR)[0][0] - 30
    higherindex = np.where(highSNR)[0][-1] + 30
    dataX = tdi_fs["X"].isel(f=slice(Xs.kmin, Xs.kmin + len(Xs)))[lowerindex:higherindex]
    dataY = tdi_fs["Y"].isel(f=slice(Ys.kmin, Ys.kmin + len(Ys)))[lowerindex:higherindex]
    dataZ = tdi_fs["Z"].isel(f=slice(Zs.kmin, Zs.kmin + len(Zs)))[lowerindex:higherindex]
    Xs, Ys, Zs = (
        Xs[lowerindex:higherindex],
        Ys[lowerindex:higherindex],
        Zs[lowerindex:higherindex],
    )
    spd_data = np.abs(dataX) ** 2 + np.abs(dataY) ** 2 + np.abs(dataZ) ** 2
    noise = (np.mean(spd_data[:2]) + np.mean(spd_data[-2:])).values / 2
    noise = 0  # (np.mean(spd_data).values)/2
    fmin, fmax = float(Xs.f[0]), float(Xs.f[-1] + Xs.attrs["df"])
    freq = np.array(Xs.sel(f=slice(fmin, fmax)).f)
    Nmodel = get_noise_model(noise_model, freq)
    Sn = Nmodel.psd(freq=freq, option="X")
    diff = np.abs(dataX - Xs.values) ** 2 + np.abs(dataY - Ys.values) ** 2 + np.abs(dataZ - Zs.values) ** 2
    p1 = float(np.sum(diff / (Sn + noise)) * Xs.df) / 2.0
    p1 = -p1

    # index_low = np.searchsorted(Xs_added.f, dataX.f[0])
    # Xs_added = Xs_added[index_low : index_low + len(dataX)]
    # Ys_added = Ys_added[index_low : index_low + len(dataY)]
    # Zs_added = Zs_added[index_low : index_low + len(dataZ)]
    # plt.figure(figsize=fig_size)
    # ax1 = plt.subplot(111)
    # # plt.plot(dataX_training.f*1000,dataX_training.values, label='data')
    # ax1.plot(dataX.f * 1000, dataX.values.real, label="data", marker="o", zorder=5)
    # ax1.plot(Xs.f * 1000, Xs.values.real, label="VGB", marker=".", zorder=5)
    # Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs, oversample=4, simulator="synthlisa")
    # index_low = np.searchsorted(Xs.f, dataX.f[0])
    # Xs = Xs[index_low : index_low + len(dataX)]
    # Ys = Ys[index_low : index_low + len(dataY)]
    # Zs = Zs[index_low : index_low + len(dataZ)]
    # ax1.plot(Xs.f * 1000, Xs.values.real, label="start", marker=".", zorder=5)
    # # ax1.plot(Xs_added2.f * 1000, Xs_added2.values.real, label="VGB2", marker=".", zorder=5)
    # ax1.axvline(boundaries['Frequency'][0]* 1000, color= 'red')
    # ax1.axvline(boundaries['Frequency'][1]* 1000, color= 'red')
    # # ax1.plot(Xs.f * 1000, dataX.values.real - Xs.values.real, label="residual", alpha=0.8, color="red", marker=".")
    # plt.legend()
    # plt.show()

    # f, psdX =  scipy.signal.welch(tdi_ts["X"], fs=1.0/dt, window='hanning', nperseg=len(tdi_ts["X"])/10)
    # peaks, properties = scipy.signal.find_peaks(np.log(np.sqrt(psdX)/np.mean(np.sqrt(psdX))),width=2,  prominence=(0.93,None))
    # peaks = peaks[np.logical_and(f[peaks]>10**-4, f[peaks]<2*10**-2)]
    # plt.figure(figsize=fig_size)
    # ax1 = plt.subplot(111)
    # ax1.loglog(f*1000, np.sqrt(psdX), label="TDI X")
    # ax1.loglog(f[peaks]*1000, np.sqrt(psdX[peaks]),'.', label="peaks")
    # # ax1.plot(Xs_added2.f * 1000, Xs_added2.values.real, label="VGB2", marker=".", zorder=5)
    # # ax1.axvline(boundaries['Frequency'][0]* 1000, color= 'red')
    # # ax1.axvline(boundaries['Frequency'][1]* 1000, color= 'red')
    # ax1.set_xlim(1,10)
    # plt.legend()
    # plt.show()

    print("p start", p1, "p true", loglikelihood(pGB))

    boundaries_reduced = deepcopy(boundaries)
    best_params = deepcopy(pGBs)

    parameters_recorded = [None] * 64
    pbar = tqdm(total=len(parameters_recorded))
    pool = mp.Pool(mp.cpu_count())
    start = time.time()
    parameters_recorded = pool.map(CoordinateMC, [n for n in range(len(parameters_recorded))])
    pool.close()
    pool.join()
    pbar.close()
    print('parallel time', time.time()-start)

    best_value = parameters_recorded[0]['Loglikelihood'][-1]
    best_run = 0
    loglikelihoodofruns = np.zeros(len(parameters_recorded))
    for i in range(len(parameters_recorded)):
        loglikelihoodofruns[i] = parameters_recorded[i]['Loglikelihood'][-1]
    best_value = np.max(loglikelihoodofruns)
    best_run = np.argmax(loglikelihoodofruns)
    good_runs = loglikelihoodofruns > best_value*1.8
    indices = (-loglikelihoodofruns).argsort()[:len(good_runs)]
    pGBmodes = []
    for i in range(len(good_runs)):
        if good_runs[i]:
            pGBmodes.append({})
    indices = (-loglikelihoodofruns).argsort()[:len(pGBmodes)]
    pGBmodes = []
    for i in indices:
        pGBmodes.append({})
        for parameter in parameters + ["Loglikelihood"]:
            pGBmodes[-1][parameter] = parameters_recorded[i][parameter][-1]


    #vgb0
    # maxpGB = {'Amplitude': 1.080513721009205e-22, 'EclipticLatitude': 0.35146042068944394, 'EclipticLongitude': -2.7339382470202227, 'Frequency': 0.0013596195825758522, 'FrequencyDerivative': 1.149294847024252e-17, 'Inclination': 0.5085588022277326, 'InitialPhase': 2.08481566510721, 'Polarization': 3.1035650828185584}
    # maxpGB = {'Amplitude': 1.1387744606299033e-22, 'EclipticLatitude': 0.34065690302493784, 'EclipticLongitude': -2.7525936350006894, 'Frequency': 0.0013596198275430923, 'FrequencyDerivative': 7.043152913703897e-18, 'Inclination': 0.6115962268957064, 'InitialPhase': 2.987037385921129, 'Polarization': 0.3988720064203186}
    # vgb1
    # maxpGB = {'Amplitude': 1.8878124482548156e-22, 'EclipticLatitude': -0.4945665539262134, 'EclipticLongitude': -2.556425644869919, 'Frequency': 0.0012531299256783774, 'FrequencyDerivative': 2.612138275081537e-18, 'Inclination': 0.8378469295715183, 'InitialPhase': 1.05572105044955, 'Polarization': 1.4407699699324252}
    #vgb2
    # maxpGB = {'Amplitude': 1.8197822740805281e-22, 'EclipticLatitude': -0.0227171575038093, 'EclipticLongitude': -2.1852674475987226, 'Frequency': 0.0018132382498430247, 'FrequencyDerivative': 5.618399870615021e-17, 'Inclination': 0.555087061273414, 'InitialPhase': 0.40617952161336923, 'Polarization': 2.389277698273662}
    #vgb3
    # maxpGB = {'Amplitude': 5.972895098376508e-23, 'EclipticLatitude': 0.4910386553462666, 'EclipticLongitude': 2.2789424229984716, 'Frequency': 0.001666671504560602, 'FrequencyDerivative': 1.3550038826189221e-19, 'Inclination': 1.1243835377244706, 'InitialPhase': 0.656426573704772, 'Polarization': 1.599445906356936}
    #vgb4
    # maxpGB ={'Amplitude': 1.2896328615704332e-22, 'EclipticLatitude': 0.6470805386714888, 'EclipticLongitude': 2.969991327893853, 'Frequency': 0.001944139730784408, 'FrequencyDerivative': 2.0004671277000396e-17, 'Inclination': 0.5250921144847437, 'InitialPhase': 2.778230376173407, 'Polarization': 1.9774563736644553}
    #vgb5 second is shifted
    # maxpGB = {'Amplitude': 4.759583244583299e-23, 'EclipticLatitude': -0.30228233054788095, 'EclipticLongitude': 0.4217098156672088, 'Frequency': 0.003220610165188476, 'FrequencyDerivative': 1.8762496834485188e-17, 'Inclination': 1.5910831797571503, 'InitialPhase': 3.0022326099685652, 'Polarization': 2.533989150735752}
    # maxpGB = {'Amplitude': 5.2338558008604886e-23, 'EclipticLatitude': -0.28836200951264557, 'EclipticLongitude': 0.4410248829693293, 'Frequency': 0.003220609930596796, 'FrequencyDerivative': 7.662392857337062e-17, 'Inclination': 1.599482971586445, 'InitialPhase': 0.09406951705562117, 'Polarization': 0.8310838280696345}
    #vgb6
    # maxpGB = {'Amplitude': 3.30762576172529e-23, 'EclipticLatitude': 0.7996139227673765, 'EclipticLongitude': -1.1289788491784982, 'Frequency': 0.003512502674691513, 'FrequencyDerivative': 7.545928043483127e-19, 'Inclination': 1.5883745476051856, 'InitialPhase': 0.3089308933545025, 'Polarization': 1.3201741777568354}
    #vgb7
    # maxpGB = {'Amplitude': 6.712522478939933e-23, 'EclipticLatitude': 0.43272984046449575, 'EclipticLongitude': 2.3119366800112173, 'Frequency': 0.0016834931242863315, 'FrequencyDerivative': 1.607079720227523e-16, 'Inclination': 1.6380039800624906, 'InitialPhase': 2.7981190242628093, 'Polarization': 3.1035650828185584}
    #vgb8
    # maxpGB = {'Amplitude': 5.916085522914226e-23, 'EclipticLatitude': -0.08794877047292846, 'EclipticLongitude': 2.102955853197464, 'Frequency': 0.0062202800609208134, 'FrequencyDerivative': 7.497316758880453e-16, 'Inclination': 0.5298895406862852, 'InitialPhase': 0.6952592815551397, 'Polarization': 1.664743447592667}
    # maxpGB = {'Amplitude': 5.958188508736942e-23, 'EclipticLatitude': -0.20876046473326043, 'EclipticLongitude': 2.100206957505942, 'Frequency': 0.006220292548425058, 'FrequencyDerivative': 2.6882778683832e-16, 'Inclination': 0.7227920918062414, 'InitialPhase': 3.0927903442057665, 'Polarization': 2.6779311950651894}
    # maxpGB = {'Amplitude': 8.126656596001908e-23, 'EclipticLatitude': -0.04484340674690194, 'EclipticLongitude': 2.1005306165381326, 'Frequency': 0.00622027634608126, 'FrequencyDerivative': 8.500086694814612e-16, 'Inclination': 0.9365397412117447, 'InitialPhase': 1.6731189287589059, 'Polarization': 2.2592867175778575}
    #vgb9
    # maxpGB = {'Amplitude': 1.6569883812641871e-22, 'EclipticLatitude': 0.19959021519888565, 'EclipticLongitude': 1.7866188257135631, 'Frequency': 0.0026130062714148075, 'FrequencyDerivative': 1.1463323442476781e-16, 'Inclination': 1.536127072699318, 'InitialPhase': 2.652971140157701, 'Polarization': 0.6738994920582831}
    # pGBmodes = [maxpGB]
    if len(pGBmodes) > 10:
        pGBmodes = pGBmodes[:10]
    for i in range(len(pGBmodes)):
        maxpGB = {}
        for parameter in parameters:
            maxpGB[parameter] = pGBmodes[i][parameter]
        # print(maxpGB)
        boundaries_reduced = Reduce_boundaries(maxpGB, boundaries,ratio=0.1)
        boundaries_reduced = deepcopy(boundaries)
        x = []
        pGBs01 = {}
        for parameter in parameters:
            if parameter in ["EclipticLatitude"]:
                pGBs01[parameter] = (np.sin(maxpGB[parameter]) - boundaries_reduced[parameter][0]) / (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])
            elif parameter in ["Inclination"]:
                pGBs01[parameter] = (np.cos(maxpGB[parameter]) - boundaries_reduced[parameter][0]) / (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])
            elif parameter in ['Amplitude',"FrequencyDerivative"]:
                pGBs01[parameter] = (np.log10(maxpGB[parameter]) - boundaries_reduced[parameter][0]) / (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])
            else:
                pGBs01[parameter] = (maxpGB[parameter] - boundaries_reduced[parameter][0]) / (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])
        for parameter in parameters:
            x.append(pGBs01[parameter])
        # print(loglikelihood(maxpGB))
        res = scipy.optimize.minimize(function, x, args=boundaries_reduced, method='SLSQP', tol=1e-6, bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)))
        maxpGB = scaletooriginal(res.x,boundaries_reduced)
        # print('optimized loglikelihood', loglikelihood(maxpGB))
        for n in range(2):
            boundaries_reduced = Reduce_boundaries(maxpGB, boundaries,ratio=0.1)
            x = []
            pGBs01 = {}
            for parameter in parameters:
                if parameter in ["EclipticLatitude"]:
                    pGBs01[parameter] = (np.sin(maxpGB[parameter]) - boundaries_reduced[parameter][0]) / (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])
                elif parameter in ["Inclination"]:
                    pGBs01[parameter] = (np.cos(maxpGB[parameter]) - boundaries_reduced[parameter][0]) / (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])
                elif parameter in ['Amplitude',"FrequencyDerivative"]:
                    pGBs01[parameter] = (np.log10(maxpGB[parameter]) - boundaries_reduced[parameter][0]) / (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])
                else:
                    pGBs01[parameter] = (maxpGB[parameter] - boundaries_reduced[parameter][0]) / (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])
            for parameter in parameters:
                x.append(pGBs01[parameter])
            # print(loglikelihood(maxpGB))
            # print(pGBs01)
            res = scipy.optimize.minimize(function, x, args=boundaries_reduced, method='SLSQP', tol=1e-6, bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)))
            maxpGB = scaletooriginal(res.x,boundaries_reduced)
            best_value = loglikelihood(maxpGB)
        if i == 0:
            current_best_value = best_value
            current_maxpGB = maxpGB
        try:
            if current_best_value < best_value:
                current_best_value = best_value
                current_maxpGB = maxpGB
        except:
            pass
        print('optimized loglikelihood', loglikelihood(maxpGB),maxpGB)
    
    maxpGB = current_maxpGB
    best_value = loglikelihood(maxpGB)
    boundaries_reduced = Reduce_boundaries(maxpGB, boundaries,ratio=0.1)
    for n in range(2):
        resolution = 2000
        pGB_centered = True
        std=0.2
        if n == 1:
            resolution = 1000
            pGB_centered = True
            std=0.2
        resolution_reduced = int(20 ** 2)
        resolution_reduced = resolution
        start = time.time()
        parameter = "Frequency"
        train_samples = sampler(resolution_reduced,parameters,maxpGB, boundaries_reduced, p1, uniform=False, twoD=False, gaussian=True, pGB_centered=pGB_centered, std=std)
        print('sample time of', resolution, 'samples ',time.time() - start)
        train_x = np.zeros((resolution_reduced, len(parameters)))
        i = 0
        for parameter in parametersfd:
            if parameter == 'FrequencyDerivative':
                train_x[:, i] = train_samples[parameter]
            else:
                train_x[:, i] = train_samples[parameter]
            i += 1
        train_y = train_samples["Likelihood"]
        normalizer = np.mean(train_y - train_y.max())
        train_y2 = (train_y - train_y.max())/normalizer
        train_y = train_samples["Likelihood"]
        nu = np.mean(train_y)
        sigma = np.std(train_y)
        train_y = (train_y - nu) / sigma
        cutoff = 2.5
        if n == 1:
            cutoff = 1.15
        for parameter in parameters:
            ranges = train_x[:,parameters.index(parameter)][train_samples["Likelihood"]>best_value*cutoff]
            lowerbound = ranges.min()*(boundaries_reduced[parameter][1]-boundaries_reduced[parameter][0])+boundaries_reduced[parameter][0]
            upperbound = ranges.max()*(boundaries_reduced[parameter][1]-boundaries_reduced[parameter][0])+boundaries_reduced[parameter][0]
            boundaries_reduced[parameter] = [lowerbound, upperbound]
        print(boundaries_reduced)
        print(len(ranges))

    split_fd = -17
    if boundaries_reduced['FrequencyDerivative'][1] < split_fd:
        boundaries_reduced['FrequencyDerivative'][0] = -19
        boundaries_reduced['FrequencyDerivative'][1] = -16
    if boundaries_reduced['FrequencyDerivative'][0] < split_fd and boundaries_reduced['FrequencyDerivative'][1] > split_fd:
        boundaries_reduced1 = deepcopy(boundaries_reduced)
        boundaries_reduced1['FrequencyDerivative'] = [boundaries_reduced['FrequencyDerivative'][0], split_fd]
        boundaries_reduced2 = deepcopy(boundaries_reduced)
        boundaries_reduced2['FrequencyDerivative'] = [split_fd, boundaries_reduced['FrequencyDerivative'][1]]
    else:
        halfway = boundaries_reduced['FrequencyDerivative'][0]+(boundaries_reduced['FrequencyDerivative'][1]-boundaries_reduced['FrequencyDerivative'][0])/2
        boundaries_reduced1 = deepcopy(boundaries_reduced)
        boundaries_reduced1['FrequencyDerivative'] = [boundaries_reduced['FrequencyDerivative'][0], halfway]
        boundaries_reduced2 = deepcopy(boundaries_reduced)
        boundaries_reduced2['FrequencyDerivative'] = [halfway, boundaries_reduced['FrequencyDerivative'][1]]
    print(boundaries_reduced)
    resolution = 2000
    resolution_reduced = int(20 ** 2)
    resolution_reduced = resolution
    start = time.time()
    parameter = "Frequency"
    train_samples = sampler(resolution_reduced,parameters,maxpGB, boundaries_reduced1, p1, uniform=False, twoD=False)
    print('sample time of', resolution, 'samples ',time.time() - start)
    train_x = np.zeros((resolution_reduced, len(parameters)))
    i = 0
    boundary_ratio = (boundaries_reduced1['FrequencyDerivative'][1]-boundaries_reduced1['FrequencyDerivative'][0])/(boundaries_reduced2['FrequencyDerivative'][1]-boundaries_reduced1['FrequencyDerivative'][0])
    for parameter in parametersfd:
        if parameter == 'FrequencyDerivative':
            train_x[:, i] = train_samples[parameter]*boundary_ratio
        else:
            train_x[:, i] = train_samples[parameter]
        i += 1
    train_y = train_samples["Likelihood"]
    resolution = 500
    test_samples = sampler(resolution, parameters, maxpGB, boundaries_reduced1, p1, uniform=False, twoD=False, calculate_loglikelihood=True)
    test_x = np.zeros((resolution, len(parameters)))
    i = 0
    for parameter in parametersfd:
        if parameter == 'FrequencyDerivative':
            test_x[:, i] = test_samples[parameter]*boundary_ratio
        else:
            test_x[:, i] = test_samples[parameter]
        i += 1
    test_y = test_samples["Likelihood"]

    nu = np.mean(train_y)
    sigma = np.std(train_y)
    train_y = (train_y - nu) / sigma
    kernel = RBF(length_scale=[0.5,2,1,1,1,1,1,1],length_scale_bounds=[(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,10)])
    start = time.time()
    gpr = GaussianProcessRegressor(kernel=kernel, random_state=0).fit(train_x, train_y)
    # gpr = traingpmodelsk(train_x, train_y, kernel)
    print('train',time.time() - start)
    start = time.time()
    observed_pred_sk = gpr.predict(test_x)
    print('eval',time.time() - start)
    observed_pred_sk_scaled = observed_pred_sk*sigma +nu
    print("RMSE ",np.sqrt(mean_squared_error(test_y,observed_pred_sk_scaled)))

    # cov_func = pm.gp.cov.ExpQuad(input_dim = 8, ls=[0.5,2,1,1,1,1,1,1])
    # with pm.Model() as lgcp_model:
    #     rho = pm.Uniform("rho", lower=0.1, upper=2, shape=8)
    #     cov_func = pm.gp.cov.ExpQuad(input_dim = 8, ls=rho)
        
    # with lgcp_model:
    #     gp = pm.gp.Marginal(cov_func=cov_func)

    #     # Since the normal noise model and the GP are conjugates, we use `Marginal` with the `.marginal_likelihood` method
    #     y_ = gp.marginal_likelihood("y", X=train_x, y=train_y, noise=0)

    #     # this line calls an optimizer to find the MAP
    #     mp = pm.sample(include_transformed=True)
    #     gp = pm.gp.Marginal(cov_func=cov_func)
    # sorted([name + ":" + str(mp[name]) for name in mp.keys() if not name.endswith("_")])

    resolution = 2000
    maxpGBmod = deepcopy(maxpGB)
    maxpGBmod['FrequencyDerivative'] = 10**(np.random.rand()*(boundaries_reduced2['FrequencyDerivative'][1]-boundaries_reduced2['FrequencyDerivative'][0])+boundaries_reduced2['FrequencyDerivative'][0])
    start = time.time()
    parameter = "Frequency"
    train_samples2 = sampler(resolution,parameters,maxpGBmod, boundaries_reduced2, p1, uniform=False, twoD=False)
    print('sample time of', resolution, 'samples ',time.time() - start)
    train_x = np.zeros((resolution, len(parameters)))
    i = 0
    for parameter in parametersfd:
        if parameter == 'FrequencyDerivative':
            train_x[:, i] = train_samples2[parameter]*(1-boundary_ratio)+boundary_ratio
        else:
            train_x[:, i] = train_samples2[parameter]
        i += 1
    train_y = train_samples2["Likelihood"]
    resolution = 500
    test_samples2 = sampler(resolution, parameters, maxpGBmod, boundaries_reduced2, p1, uniform=False, twoD=False, calculate_loglikelihood=True,)
    test_x = np.zeros((resolution, len(parameters)))
    i = 0
    for parameter in parametersfd:
        if parameter == 'FrequencyDerivative':
            test_x[:, i] = test_samples2[parameter]*(1-boundary_ratio)+boundary_ratio
        else:
            test_x[:, i] = test_samples2[parameter]
        i += 1
    test_y = test_samples2["Likelihood"]

    nu2 = np.mean(train_y)
    sigma2 = np.std(train_y)
    train_y = (train_y - nu2) / sigma2
    kernel = RBF(length_scale=[0.5,2,1,1,1,1,1,1],length_scale_bounds=[(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,20),(0.1,20)])
    start = time.time()
    gpr2 = traingpmodelsk(train_x, train_y, kernel)
    print('train',time.time() - start)
    start = time.time()
    observed_pred_sk = gpr2.predict(test_x)
    print('eval',time.time() - start)
    observed_pred_sk_scaled2 = observed_pred_sk*sigma2 +nu2
    print("RMSE ",np.sqrt(mean_squared_error(test_y,observed_pred_sk_scaled2)))


    resolution = 1*10 ** 6
    start = time.time()
    test_x_m = np.random.uniform(size=(resolution,len(parameters)))
    # test_x_m = np.random.normal(loc= 0.5,scale=0.2,size=(resolution,len(parameters)))
    # std_scale = 0.4
    # test_x_m = scipy.stats.truncnorm.rvs((0-0.5)/std_scale,(1-0.5)/std_scale,loc=0.5,scale=std_scale,size=(resolution,len(parameters)))
    test_x_m = test_x_m[test_x_m[:,4].argsort()]
    def next(arr, target):
        start = 0
        end = len(arr) - 1
        ans = -1
        while (start <= end):
            mid = (start + end) // 2
            # Move to right side if target is
            # greater.
            if (arr[mid] <= target):
                start = mid + 1
            # Move left side.
            else:
                ans = mid
                end = mid - 1
        return ans
    numberinlowerbatch = next(test_x_m[:,4],boundary_ratio)
    numberinupperbatch = resolution-numberinlowerbatch
    test_x_m1 = test_x_m[:numberinlowerbatch]
    test_x_m2 = test_x_m[numberinlowerbatch:]
    print('sample time', time.time()-start)

    partial_length = 1*10**3
    start = time.time()
    observed_pred_mean = np.zeros(resolution)
    observed_pred_sk = np.zeros(resolution)
    def Evaluate(i):
        prediction = gpr.predict(test_x_m1[(i)*partial_length:(i+1)*partial_length])
        return prediction
    def Evaluate2(i):
        prediction = gpr2.predict(test_x_m2[(i)*partial_length:(i+1)*partial_length])
        return prediction
    pool = mp.Pool(mp.cpu_count()-1)
    observed_pred_sk = pool.map(Evaluate, [n for n in range(int(numberinlowerbatch/partial_length))])
    pool.close()
    pool.join()
    print('eval time', time.time()-start)
    try:
        observed_pred_sk = np.append(observed_pred_sk, gpr.predict(test_x_m1[int(numberinlowerbatch/partial_length)*partial_length:]))
    except:
        pass
    observed_pred_sk = np.asarray(observed_pred_sk)
    observed_pred_sk = observed_pred_sk.reshape(numberinlowerbatch)
    observed_pred_mean[:numberinlowerbatch] = observed_pred_sk[:numberinlowerbatch]*sigma +nu
    start = time.time()
    pool = mp.Pool(mp.cpu_count()-1)
    observed_pred_sk2 = pool.map(Evaluate2, [n for n in range(int(numberinupperbatch/partial_length))])
    pool.close()
    pool.join()
    print('eval time', time.time()-start)
    try:
        observed_pred_sk2 = np.append(observed_pred_sk2, gpr2.predict(test_x_m2[int(numberinupperbatch/partial_length)*partial_length:]))
    except:
        pass
    observed_pred_sk2 = np.asarray(observed_pred_sk2)
    observed_pred_sk2 = observed_pred_sk2.reshape((numberinupperbatch))
    observed_pred_mean[-numberinupperbatch:] = observed_pred_sk2*sigma2 +nu2

    flatsamples = np.zeros(resolution)
    flatsamplesparameters = np.zeros((resolution,len(parameters)+1))
    i = 0
    flatsamples[:] = observed_pred_mean
    flatsamplesparameters[:,1:] = test_x_m
    flatsamplesparameters[:,0] = observed_pred_mean

    maxindx = np.unravel_index(flatsamplesparameters[:,0].argmax(), flatsamplesparameters[:,0].shape)
    max_parameters = flatsamplesparameters[maxindx[0],1:]
    max_loglike = flatsamplesparameters[:,0].max()
    maxpGBpredicted = scaletooriginal(max_parameters, boundaries_reduced)
    if loglikelihood(maxpGBpredicted) > loglikelihood(maxpGB):
        maxpGB = maxpGBpredicted
    best_value = loglikelihood(maxpGB)
    print("pred", max_loglike, "true", loglikelihood(scaletooriginal(max_parameters, boundaries_reduced)), "max", loglikelihood(maxpGB), maxpGB)

    np.random.shuffle(flatsamplesparameters)
    start = time.time()
    normalizer = sum(np.exp(flatsamplesparameters[:,0]-best_value))
    flatsamples_normalized = np.exp(flatsamplesparameters[:,0]-best_value)/normalizer
    mcmc_samples = []
    mcmc_samples.append(flatsamplesparameters[0,1:])
    previous_p = flatsamples_normalized[0]
    if previous_p == 0:
        previous_p == 1e-300
    rejection_count = 0
    current_parameters = flatsamplesparameters[0,1:]
    # probability = scipy.stats.multivariate_normal.pdf(flatsamplesparameters[:,1:], mean=[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5], cov=[std_scale,std_scale,std_scale,std_scale,std_scale,std_scale,std_scale,std_scale])
    # previous_probability = probability[0]
    accepted = 0
    for i in range(len(flatsamples_normalized)-1):
        # if (flatsamples_normalized[i+1] / previous_p * previous_probability/probability[i+1]) > np.random.uniform():
        #     previous_p = flatsamples_normalized[i+1]
        #     previous_probability = probability[i+1]
        #     current_parameters = flatsamplesparameters[i+1,1:]
        #     mcmc_samples.append(current_parameters)
        # else:
        #     mcmc_samples.append(current_parameters)
        if (flatsamples_normalized[i+1] / previous_p) > np.random.uniform():
            previous_p = flatsamples_normalized[i+1]
            # previous_probability = probability[i+1]
            current_parameters = flatsamplesparameters[i+1,1:]
            rejection_count = 0
            accepted += 1
        else:
            rejection_count += 1
            if rejection_count == 10:
                mcmc_samples.append(current_parameters)
                rejection_count = 0
    mcmc_samples = np.asarray(mcmc_samples)
    print('time MHMC', time.time()-start)
    print('acceptance rate %',accepted/resolution*100)

    # plt.figure()
    # plt.scatter(mcmc_samples[:10**4,3],mcmc_samples[:10**4,4])
    for round in range(2):
        resolution = 2*10**6
        numPoints = 33
        if round == 1:
            resolution = 1*10**5
            numPoints = 2**6+1
        test_x_m = np.zeros((resolution,len(parameters)))
        ax = np.linspace(-0.15,1.15,numPoints)
        from fastkde import fastKDE
        mypdf,axes = fastKDE.pdf(mcmc_samples[:,1],mcmc_samples[:,2], axes=[ax,ax])
        dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
        data, pdfs = dist(resolution)
        test_x_m[:,1] = data[1]
        test_x_m[:,2] = data[0]
        probability = pdfs
        mypdf,axes = fastKDE.pdf(mcmc_samples[:,0],mcmc_samples[:,5], axes=[ax,ax])
        dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
        data, pdfs = dist(resolution)
        test_x_m[:,0] = data[1]
        test_x_m[:,5] = data[0]
        # plt.figure()
        # plt.scatter(np.log10(test_x_m[:10000,0]),np.arccos(test_x_m[:10000,5]* (boundaries_reduced['Inclination'][1] - boundaries_reduced['Inclination'][0]) + boundaries_reduced['Inclination'][0]), c=pdfs[:10000])
        probability *= pdfs
        mypdf,axes = fastKDE.pdf(mcmc_samples[:,3],mcmc_samples[:,4], axes=[ax,ax])
        dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
        data, pdfs = dist(resolution)
        test_x_m[:,3] = data[1]
        test_x_m[:,4] = data[0]
        probability *= pdfs
        # plt.figure()
        # plt.scatter(test_x_m[:,3],test_x_m[:,4], c=pdfs)
        # plt.figure()
        # plt.scatter(test_x_m[:1000,3],test_x_m[:1000,4], c=probability[:1000])
        # plt.show()
        # mypdf,axes = fastKDE.pdf(mcmc_samples[:,6],mcmc_samples[:,7],numPoints= 33, axisExpansionFactor=0)
        # dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
        # data, pdfs = dist(resolution)
        test_x_m[:,6] = np.random.uniform(size=resolution)
        test_x_m[:,7] = np.random.uniform(size=resolution)
        # probability *= pdfs
        for n in range(6):
            index = np.where(test_x_m[:,n] > 0)
            test_x_m = test_x_m[index]
            probability = probability[index]
            index = np.where(test_x_m[:,n] < 1)
            test_x_m = test_x_m[index]
            probability = probability[index]
        
        # fig =  corner.corner(test_x_m,  bins=40, hist_kwargs={'density':True, 'lw':3}, plot_datapoints=False, fill_contours=False,  show_titles=True, \
        #                 color='#348ABD',  truth_color='k', use_math_test=True, \
        #                  levels=[0.9], title_kwargs={"fontsize": 12})
        # plt.show()
        if round == 0:
            start = time.time()
            
            # test_x_m = np.random.normal(loc= 0.5,scale=0.2,size=(resolution,len(parameters)))
            # std_scale = 0.4
            # test_x_m = scipy.stats.truncnorm.rvs((0-0.5)/std_scale,(1-0.5)/std_scale,loc=0.5,scale=std_scale,size=(resolution,len(parameters)))
            new_index = test_x_m[:,4].argsort()
            test_x_m = test_x_m[new_index]
            probability = probability[new_index]
            numberinlowerbatch = next(test_x_m[:,4],boundary_ratio)
            numberinupperbatch = len(probability)-numberinlowerbatch
            test_x_m1 = test_x_m[:numberinlowerbatch]
            test_x_m2 = test_x_m[numberinlowerbatch:]
            print('sample time', time.time()-start)
            partial_length = 1*10**3
            start = time.time()
            observed_pred_mean = np.zeros(len(probability))
            observed_pred_sk = np.zeros(len(probability))
            def Evaluate(i):
                prediction = gpr.predict(test_x_m1[(i)*partial_length:(i+1)*partial_length])
                return prediction
            def Evaluate2(i):
                prediction = gpr2.predict(test_x_m2[(i)*partial_length:(i+1)*partial_length])
                return prediction
            pool = mp.Pool(mp.cpu_count()-1)
            observed_pred_sk = pool.map(Evaluate, [n for n in range(int(numberinlowerbatch/partial_length))])
            pool.close()
            pool.join()
            print('eval time', time.time()-start)
            try:
                observed_pred_sk = np.append(observed_pred_sk, gpr.predict(test_x_m1[int(numberinlowerbatch/partial_length)*partial_length:]))
            except:
                pass
            observed_pred_sk = np.asarray(observed_pred_sk)
            observed_pred_sk = observed_pred_sk.reshape(numberinlowerbatch)
            observed_pred_mean[:numberinlowerbatch] = observed_pred_sk[:numberinlowerbatch]*sigma +nu
            start = time.time()
            pool = mp.Pool(mp.cpu_count()-1)
            observed_pred_sk2 = pool.map(Evaluate2, [n for n in range(int(numberinupperbatch/partial_length))])
            pool.close()
            pool.join()
            print('eval time', time.time()-start)
            try:
                observed_pred_sk2 = np.append(observed_pred_sk2, gpr2.predict(test_x_m2[int(numberinupperbatch/partial_length)*partial_length:]))
            except:
                pass
            observed_pred_sk2 = np.asarray(observed_pred_sk2)
            observed_pred_sk2 = observed_pred_sk2.reshape((numberinupperbatch))
            observed_pred_mean[-numberinupperbatch:] = observed_pred_sk2*sigma2 +nu2

            flatsamples = np.zeros(len(probability))
            flatsamplesparameters = np.zeros((len(probability),len(parameters)+1))
            flatsamples[:] = observed_pred_mean
            flatsamplesparameters[:,1:] = test_x_m
            flatsamplesparameters[:,0] = observed_pred_mean

            maxindx = np.unravel_index(flatsamplesparameters[:,0].argmax(), flatsamplesparameters[:,0].shape)
            max_parameters = flatsamplesparameters[maxindx[0],1:]
            max_loglike = flatsamplesparameters[:,0].max()
            maxpGBpredicted = scaletooriginal(max_parameters, boundaries_reduced)
            if loglikelihood(maxpGBpredicted) > loglikelihood(maxpGB):
                maxpGB = maxpGBpredicted
            best_value = loglikelihood(maxpGB)
            print("pred", max_loglike, "true", loglikelihood(scaletooriginal(max_parameters, boundaries_reduced)), "max", loglikelihood(maxpGB), maxpGB)

            indexes = np.arange(len(probability))
            np.random.shuffle(indexes)
            flatsamplesparameters = flatsamplesparameters[indexes]
            probability = probability[indexes]
            start = time.time()
            normalizer = sum(np.exp(flatsamplesparameters[:,0]-best_value))
            flatsamples_normalized = np.exp(flatsamplesparameters[:,0]-best_value)/normalizer
            mcmc_samples = []
            mcmc_samples.append(flatsamplesparameters[0,1:])
            previous_p = flatsamples_normalized[0]
            if previous_p == 0:
                previous_p == 1e-300
            rejection_count = 0
            current_parameters = flatsamplesparameters[0,1:]
            accepted = 0
            previous_probability = probability[0]
            for i in range(len(flatsamples_normalized)-1):
                if (flatsamples_normalized[i+1] / previous_p) * (previous_probability/probability[i+1]) > np.random.uniform():
                    previous_p = flatsamples_normalized[i+1]
                    previous_probability = probability[i+1]
                    current_parameters = flatsamplesparameters[i+1,1:]
                    mcmc_samples.append(current_parameters)
                    accepted += 1
                else:
                    mcmc_samples.append(current_parameters)
            mcmc_samples = np.asarray(mcmc_samples)
            print('time MHMC', time.time()-start)
            print('acceptance rate %',accepted/len(probability)*100)
            start = time.time()
            mcmc_samples_rescaled = np.zeros(np.shape(mcmc_samples))
            i = 0
            for parameter in parametersfd:
                if parameter in ["EclipticLatitude"]:
                    mcmc_samples_rescaled[:,i] = np.arcsin((mcmc_samples[:,parametersfd.index(parameter)] * (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])) + boundaries_reduced[parameter][0])
                elif parameter in ["Inclination"]:
                    mcmc_samples_rescaled[:,i] = np.arccos((mcmc_samples[:,parametersfd.index(parameter)] * (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])) + boundaries_reduced[parameter][0])
                elif parameter in ['Amplitude',"FrequencyDerivative"]:
                    mcmc_samples_rescaled[:,i] = 10**((mcmc_samples[:,parametersfd.index(parameter)] * (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])) + boundaries_reduced[parameter][0])
                else:
                    mcmc_samples_rescaled[:,i] = (mcmc_samples[:,parametersfd.index(parameter)] * (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])) + boundaries_reduced[parameter][0]
                i += 1
            print('time rescale', time.time()-start)
            start = time.time()
            df = pd.DataFrame(data=mcmc_samples_rescaled, columns=parametersfd)
            df.to_csv('/home/stefan/Repositories/ldc1_evaluation_data/submission/Stefan_full/GW'+str(int(np.round(maxpGB['Frequency']*10**8)))+'.csv',index=False)
            print('saving time', time.time()-start)

    # mcmc_samples_reshaped = mcmc_samples.reshape((len(parameters),len(mcmc_samples[:,0])))
    # kernel = scipy.stats.gaussian_kde(mcmc_samples_reshaped)
    # new_samples = kernel.resample(100000)
    # new_samples = np.transpose(new_samples)
    # filter_arr = []
    # for element in new_samples:
    #     accepted = True
    #     for value in element:
    #         if value < 0 or value > 1 :
    #             accepted = False
    #             pass
    #     filter_arr.append(accepted)
    # new_samples = new_samples[filter_arr]
    # mcmc_samples = new_samples
    mcmc_samples = test_x_m
    start = time.time()
    mcmc_samples_rescaled = np.zeros(np.shape(mcmc_samples))
    i = 0
    for parameter in parametersfd:
        if parameter in ["EclipticLatitude"]:
            mcmc_samples_rescaled[:,i] = np.arcsin((mcmc_samples[:,parametersfd.index(parameter)] * (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])) + boundaries_reduced[parameter][0])
        elif parameter in ["Inclination"]:
            mcmc_samples_rescaled[:,i] = np.arccos((mcmc_samples[:,parametersfd.index(parameter)] * (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])) + boundaries_reduced[parameter][0])
        elif parameter in ['Amplitude',"FrequencyDerivative"]:
            mcmc_samples_rescaled[:,i] = 10**((mcmc_samples[:,parametersfd.index(parameter)] * (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])) + boundaries_reduced[parameter][0])
        else:
            mcmc_samples_rescaled[:,i] = (mcmc_samples[:,parametersfd.index(parameter)] * (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])) + boundaries_reduced[parameter][0]
        i += 1
    print('time rescale', time.time()-start)
    start = time.time()
    df = pd.DataFrame(data=mcmc_samples_rescaled, columns=parametersfd)
    df.to_csv('/home/stefan/Repositories/ldc1_evaluation_data/submission/Stefan/GW'+str(int(np.round(maxpGB['Frequency']*10**8)))+'.csv',index=False)
    print('saving time', time.time()-start)

print('full time', time.time()-first_start)
length = len(probability)
datS = np.zeros(np.shape(mcmc_samples))
datS[:,0] = mcmc_samples_rescaled[:,2]
datS[:,1] = mcmc_samples_rescaled[:,1]
datS[:,2] = mcmc_samples_rescaled[:,3]
datS[:,3] = np.log10(mcmc_samples_rescaled[:,4])
datS[:,4] = mcmc_samples_rescaled[:,5]
datS[:,5] = np.log10(mcmc_samples_rescaled[:,0])
datS[:,6] = mcmc_samples_rescaled[:,6]
datS[:,7] = mcmc_samples_rescaled[:,7]
rng = []
lbls = [r'$\lambda$', r'$\beta$', 'f(Hz)', r'$\dot{f}$', r'$\iota$', r'$A$', r'$\phi$', r'$\Phi$']
for i in range(len(lbls)):
    oner = ( datS[:,i].min(), datS[:, i].max())
    rng.append(oner)
tr_s = np.zeros(len(parameters))
maxvalues = np.zeros(len(parameters))
i = 0
for parameter in ['EclipticLongitude','EclipticLatitude','Frequency','FrequencyDerivative','Inclination','Amplitude','InitialPhase','Polarization']:
    if parameter in ['Amplitude','FrequencyDerivative']:
        tr_s[i] = np.log10(pGB[parameter])
        maxvalues[i] = np.log10(maxpGB[parameter])
    else:
        tr_s[i] = pGB[parameter]
        maxvalues[i] = maxpGB[parameter]
    i += 1
fig =  corner.corner(datS,  bins=40, hist_kwargs={'density':True, 'lw':3}, plot_datapoints=False, fill_contours=False,  show_titles=True, \
                        color='#348ABD', truths= tr_s, truth_color='k', use_math_test=True, labels= lbls,\
                         levels=[0.9], title_kwargs={"fontsize": 12})
# Extract the axes
ndim = len(parameters)
axes = np.array(fig.axes).reshape((ndim, ndim))
# Loop over the diagonal
for i in range(ndim):
    ax = axes[i, i]
    ax.axvline(maxvalues[i], color="g")
# Loop over the histograms
for yi in range(ndim):
    for xi in range(yi):
        ax = axes[yi, xi]
        ax.axvline(maxvalues[xi], color="g")
        ax.axhline(maxvalues[yi], color="g")
        ax.plot(maxvalues[xi], maxvalues[yi], "sg")
# corner.ove(fig, maxvalues[None], marker="s", color="C1")
plt.show()

parameters_recorded[best_run]["Loglikelihood"].append(loglikelihood(maxpGB))
for parameter in parametersfd:
    parameters_recorded[best_run][parameter].append(maxpGB[parameter])
fig, ax = plt.subplots(2, 5, figsize=np.asarray(fig_size) * 2)
# plt.suptitle("Intermediate Parameters")
for n in range(len(parameters_recorded)):
    i = 0
    for parameter in parametersfd + ["Loglikelihood"]:
        j = 0
        if i > 3:
            j = 1
        if parameter in ["Frequency"]:
            ax[j, i % 5].plot(
                np.asarray(parameters_recorded[n][parameter]) * 10 ** 3,
                np.arange(0, len(parameters_recorded[n][parameter])),
            )
        else:
            ax[j, i % 5].plot(
                parameters_recorded[n][parameter],
                np.arange(0, len(parameters_recorded[n][parameter])),
            )
        i += 1
i = 0
for parameter in parametersfd + ["Log-likelihood"]:
    j = 0
    if i > 3:
        j = 1
    if parameter in ["Log-likelihood"]:
        ax[j, i % 5].axvline(x=loglikelihood(pGB), color="k", label="True")
    elif parameter in ["Frequency"]:
        ax[j, i % 5].axvline(x=pGB[parameter] * 10 ** 3, color="k", label="True")
    else:
        ax[j, i % 5].axvline(x=pGB[parameter], color="k", label="True")
    if parameter in ["Amplitude"]:
        ax[j, i % 5].set_xlabel(parameter, ha="right", va="top")
    elif parameter in ["Frequency"]:
        ax[j, i % 5].xaxis.set_major_locator(plt.MaxNLocator(3))
        ax[j, i % 5].set_xlabel(parameter + " (mHz)", ha="right", va="top")
    else:
        ax[j, i % 5].xaxis.set_major_locator(plt.MaxNLocator(3))
        ax[j, i % 5].set_xlabel(parameter)
    # if parameter in ['Log-likelihood']:
    #     ax[j,i%5].set_xlabel('$$')
    ax[j, i % 1].set_ylabel("Step")
    ax[j, i % 5].set_ylim(0, 150)
    if parameter in ["Log-likelihood"]:
        pass
    elif parameter == "EclipticLatitude":
        ax[j, i % 5].set_xlim(np.arcsin(boundaries[parameter][0]), np.arcsin(boundaries[parameter][1]))
    elif parameter == "Inclination":
        ax[j, i % 5].set_xlim(np.arccos(boundaries[parameter][1]), np.arccos(boundaries[parameter][0]))
    elif parameter in ["Frequency"]:
        ax[j, i % 5].set_xlim(boundaries[parameter][0] * 10 ** 3, boundaries[parameter][1] * 10 ** 3)
    else:
        ax[j, i % 5].set_xlim(boundaries[parameter][0], boundaries[parameter][1])
    i += 1
ax[0, 0].legend()
# plt.tight_layout()
plt.show()

fig, ax = plt.subplots(2, 5, figsize=np.asarray(fig_size) * 2)
# plt.suptitle("Intermediate Parameters")
for n in range(len(pGBmodes)):
    i = 0
    for parameter in parametersfd + ["Loglikelihood"]:
        j = 0
        if i > 3:
            j = 1
        if parameter in ["Frequency"]:
            ax[j, i % 5].plot(
                np.asarray(pGBmodes[n][parameter]) * 10 ** 3,
                0.5,'.'
            )
        else:
            ax[j, i % 5].plot(
                pGBmodes[n][parameter],
                0.5,'.'
            )
        i += 1
i = 0
for parameter in parametersfd + ["Log-likelihood"]:
    j = 0
    if i > 3:
        j = 1
    if parameter in ["Log-likelihood"]:
        ax[j, i % 5].axvline(x=loglikelihood(pGB), color="k", label="True")
    elif parameter in ["Frequency"]:
        ax[j, i % 5].axvline(x=pGB[parameter] * 10 ** 3, color="k", label="True")
    else:
        ax[j, i % 5].axvline(x=pGB[parameter], color="k", label="True")
    if parameter in ["Amplitude"]:
        ax[j, i % 5].set_xlabel(parameter, ha="right", va="top")
    elif parameter in ["Frequency"]:
        ax[j, i % 5].xaxis.set_major_locator(plt.MaxNLocator(3))
        ax[j, i % 5].set_xlabel(parameter + " (mHz)", ha="right", va="top")
    else:
        ax[j, i % 5].xaxis.set_major_locator(plt.MaxNLocator(3))
        ax[j, i % 5].set_xlabel(parameter)
    # if parameter in ['Log-likelihood']:
    #     ax[j,i%5].set_xlabel('$$')
    ax[j, i % 1].set_ylabel("Step")
    ax[j, i % 5].set_ylim(0, 1)
    if parameter in ["Log-likelihood"]:
        pass
    elif parameter == "EclipticLatitude":
        ax[j, i % 5].set_xlim(np.arcsin(boundaries[parameter][0]), np.arcsin(boundaries[parameter][1]))
    elif parameter == "Inclination":
        ax[j, i % 5].set_xlim(np.arccos(boundaries[parameter][1]), np.arccos(boundaries[parameter][0]))
    elif parameter in ["Frequency"]:
        ax[j, i % 5].set_xlim(boundaries[parameter][0] * 10 ** 3, boundaries[parameter][1] * 10 ** 3)
    else:
        ax[j, i % 5].set_xlim(boundaries[parameter][0], boundaries[parameter][1])
    i += 1
ax[0, 0].legend()
# plt.tight_layout()
plt.show()

Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB, oversample=4, simulator="synthlisa")
index_low = np.searchsorted(Xs.f, dataX.f[0])
Xs = Xs[index_low : index_low + len(dataX)]
Ys = Ys[index_low : index_low + len(dataY)]
Zs = Zs[index_low : index_low + len(dataZ)]

plt.figure(figsize=np.asarray(fig_size) * 2)
ax1 = plt.subplot(231)
# plt.plot(dataX_training.f*1000,dataX_training.values, label='data')
ax1.plot(dataX.f * 1000, dataX.values.real, label="Data")
ax1.plot(Xs.f * 1000, Xs.values.real, label="GB")
ax2 = plt.subplot(232)
# plt.plot(dataY_training.f*1000,dataY_training.values, label='data')
ax2.plot(dataY.f * 1000, dataY.values.real, label="Data")
ax2.plot(Ys.f * 1000, Ys.values.real, label="True Response")
ax3 = plt.subplot(233)
# plt.plot(dataZ_training.f*1000,dataZ_training.values, label='data')
ax3.plot(dataZ.f * 1000, dataZ.values.real, label="Data")
ax3.plot(Zs.f * 1000, Zs.values.real, label="True Response")
ax4 = plt.subplot(234)
# plt.plot(dataX_training.f*1000,dataX_training.values.imag, label='data')
ax4.plot(dataX.f * 1000, dataX.values.imag, label="Data")
ax4.plot(Xs.f * 1000, Xs.values.imag, label="True Response")
ax5 = plt.subplot(235)
# plt.plot(dataY_training.f*1000,dataY_training.values.imag, label='data')
ax5.plot(dataY.f * 1000, dataY.values.imag, label="Data")
ax5.plot(Ys.f * 1000, Ys.values.imag, label="True Response")
ax6 = plt.subplot(236)
# plt.plot(dataZ_training.f*1000,dataZ_training.values.imag, label='data')
ax6.plot(dataZ.f * 1000, dataZ.values.imag, label="Data")
ax6.plot(Zs.f * 1000, Zs.values.imag, label="True Response")
Xs, Ys, Zs = GB.get_fd_tdixyz(template=maxpGB, oversample=4, simulator="synthlisa")
index_low = np.searchsorted(Xs.f, dataX.f[0])
Xs = Xs[index_low : index_low + len(dataX)]
Ys = Ys[index_low : index_low + len(dataY)]
Zs = Zs[index_low : index_low + len(dataZ)]
ax1.plot(Xs.f * 1000, Xs, label="Found GB")
ax1.xaxis.set_major_locator(plt.MaxNLocator(3))
ax1.set_xlabel("f (mHz)")
ax1.set_ylabel("X-TDI real (1/Hz)")
ax1.legend(loc="upper left", numpoints=0.5)
ax2.plot(Ys.f * 1000, Ys, label="Found Response")
ax2.xaxis.set_major_locator(plt.MaxNLocator(3))
ax2.set_xlabel("f (mHz)")
ax2.set_ylabel("Y-TDI real (1/Hz)")
# ax2.legend()
ax3.plot(Zs.f * 1000, Zs, label="Found Response")
ax3.xaxis.set_major_locator(plt.MaxNLocator(3))
ax3.set_xlabel("f (mHz)")
ax3.set_ylabel("Z-TDI real (1/Hz)")
# ax3.legend()
ax4.plot(Xs.f * 1000, Xs.imag, label="Found Response")
ax4.xaxis.set_major_locator(plt.MaxNLocator(3))
ax4.set_xlabel("f (mHz)")
ax4.set_ylabel("X-TDI imag (1/Hz)")
# ax4.legend()
ax5.plot(Ys.f * 1000, Ys.imag, label="Found Response")
ax5.xaxis.set_major_locator(plt.MaxNLocator(3))
ax5.set_xlabel("f (mHz)")
ax5.set_ylabel("Y-TDI imag (1/Hz)")
# ax5.legend()
ax6.plot(Zs.f * 1000, Zs.imag, label="Found Response")
ax6.xaxis.set_major_locator(plt.MaxNLocator(3))
ax6.set_xlabel("f (mHz)")
ax6.set_ylabel("Z-TDI imag (1/Hz)")
# ax6.legend()
plt.tight_layout()
plt.show()


observed_pred = {}

number_of_test_samples = 100
test_x = {}
test_y = {}
for parameter in parametersfd:
    test_samples = sampler(
        number_of_test_samples,
        parameters,
        maxpGB,
        boundaries_reduced,
        p1,
        uniform=True,
        only=True,
        onlyparameter=parameter,
    )
    test_x[parameter] = np.zeros((number_of_test_samples, len(parameters)))
    i = 0
    for name in parametersfd:
        test_x[parameter][:, i] = test_samples[name]
        i += 1
    test_y[parameter] = test_samples["Likelihood"]
    test_x[parameter] = torch.from_numpy(test_x[parameter]).float()
    test_x[parameter] = test_x[parameter]
    test_y[parameter] = torch.from_numpy(test_y[parameter]).float()
number_of_test_samples = 500
test_x2 = {}
test_y2 = {}

for parameter in parametersfd:
    test_samples = sampler(
        number_of_test_samples,
        parameters,
        maxpGB,
        boundaries,
        p1,
        uniform=True,
        only=True,
        onlyparameter=parameter,
    )
    test_x2[parameter] = np.zeros((number_of_test_samples, len(parameters)))
    i = 0
    for name in parametersfd:
        test_x2[parameter][:, i] = test_samples[name]
        i += 1
    test_y2[parameter] = test_samples["Likelihood"]
    test_x2[parameter] = torch.from_numpy(test_x2[parameter]).float()
    test_x2[parameter] = test_x2[parameter]
    test_y2[parameter] = torch.from_numpy(test_y2[parameter]).float()
number_of_test_samples2d = 50 ** 2

parameter = "EclipticLongitude"
parameters_reduced = [
    "EclipticLatitude",
    "Frequency",
    "Inclination",
    "Amplitude",
    'InitialPhase',
    'Polarization'
]
observed_pred_mean = {}
pred_sigma = {}
for n in range(len(parameters_reduced)):
    for parameter2 in parameters_reduced:
        print(parameter+parameter2)
        if parameter != parameter2:
            test_samples = sampler(
                number_of_test_samples2d,
                parameters,
                maxpGB,
                boundaries_reduced,
                p1,
                uniform=True,
                twoD=True,
                onlyparameter=parameter,
                secondparameter=parameter2,
                calculate_loglikelihood=False,
            )
            test_x[parameter + parameter2] = np.zeros((number_of_test_samples2d, len(parameters)))
            i = 0
            for name in parametersfd:
                test_x[parameter + parameter2][:, i] = test_samples[name]
                i += 1
            test_x[parameter + parameter2] = torch.from_numpy(test_x[parameter + parameter2]).float()

            observed_pred_sk, pred_sigma[parameter + parameter2] = gpr.predict(test_x[parameter + parameter2],return_std=True)
            observed_pred_mean[parameter + parameter2] = observed_pred_sk*sigma +nu
            pred_sigma[parameter + parameter2] = pred_sigma[parameter + parameter2]*sigma +nu
            # with torch.no_grad(), gpytorch.settings.fast_pred_var():
            #     observed_pred[parameter + parameter2] = likelihood(model(test_x[parameter + parameter2]))
    parameter = parameters_reduced[0]
    parameters_reduced.pop(0)


for parameter in parametersfd:
    # Make predictions by feeding model through likelihood
    observed_pred_sk, pred_sigma[parameter] = gpr.predict(test_x[parameter],return_std=True)
    observed_pred_mean[parameter] = observed_pred_sk*sigma +nu
    pred_sigma[parameter] = pred_sigma[parameter]*sigma +nu
    # with torch.no_grad(), gpytorch.settings.fast_pred_var():
    #     observed_pred[parameter] = likelihood(model(test_x[parameter]))
    #     observed_pred_print = (observed_pred[parameter].mean.numpy() * sigma) + nu

    print(
        "sqrt(MSE) ",
        parameter,
        np.sqrt(mean_squared_error(test_y[parameter].numpy(), observed_pred_mean[parameter])),
    )
    # if parameter != parameter2:
    #     with torch.no_grad(), gpytorch.settings.fast_pred_var():
    #         observed_pred[parameter + parameter2] = likelihood(model(test_x[parameter + parameter2]))

# parameter = 'random'
# with torch.no_grad(), gpytorch.settings.fast_pred_var():
#     observed_pred[parameter] = likelihood(model(test_x[parameter]))
# observed_pred_print = (observed_pred[parameter]*sigma)+nu
# print('sqrt(MSE) ','random', np.sqrt(mean_squared_error(test_y[parameter].numpy(),observed_pred_print.mean.cpu().numpy())))


for parameter in parameters:
    if parameter in ["EclipticLatitude"]:
        samples[parameter] = np.arcsin((samples[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
        test_samples[parameter] = np.arcsin((test_samples[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
    elif parameter in ["Inclination"]:
        samples[parameter] = np.arccos((samples[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
        test_samples[parameter] = np.arccos((test_samples[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
    else:
        samples[parameter] = (samples[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]
        test_samples[parameter] = (test_samples[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]


train_x = train_x.cpu()
train_y = train_y.cpu()
# train_y = (train_y*sigma)+nu

pGB01 = {}
for parameter in parameters:
    if parameter in ["EclipticLatitude"]:
        pGB01[parameter] = (np.sin(pGB[parameter]) - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0])
    elif parameter in ["Inclination"]:
        pGB01[parameter] = (np.cos(pGB[parameter]) - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0])
    else:
        pGB01[parameter] = (pGB[parameter] - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0])

fig, ax = plt.subplots(2, 4, figsize=np.asarray(fig_size) * 2.5)
# plt.suptitle("loglikelihood")
i = 0
mean = {}
for parameter in parametersfd:
    j = 0
    if i > 3:
        j = 1
    test_x_rescaled = []
    for k in range(len(test_x[parameter])):
        test_x_rescaled.append(scaletooriginal(test_x[parameter][k].numpy(), boundaries_reduced)[parameter])
    test_x2_rescaled = []
    for k in range(len(test_x2[parameter])):
        test_x2_rescaled.append(scaletooriginal(test_x2[parameter][k].numpy(), boundaries)[parameter])
    # Plot training data as black stars
    if parameter == "Frequency":
        ax[j, i % 4].axvline(x=pGB[parameter] * 10 ** 3, color="k", label="GB")
        ax[j, i % 4].plot(
            np.asarray(test_x_rescaled[1:]) * 10 ** 3,
            np.exp(test_y[parameter].numpy()[1:]-best_value)/normalizer,
            "g",linewidth=3, alpha=0.5,
            label="True",
        )
        ax[j, i % 4].plot(
            np.asarray(test_x_rescaled[1:]) * 10 ** 3,
            np.exp(observed_pred_mean[parameter][1:]-best_value)/normalizer,
            "b",
            label="Mean",
        )
        # ax[j, i % 4].fill_between(
        #     np.asarray(test_x_rescaled[1:]) * 10 ** 3,
        #     np.exp(observed_pred_mean[parameter][1:]+pred_sigma[parameter][1:]-best_value)/normalizer,
        #     np.exp(observed_pred_mean[parameter][1:]-pred_sigma[parameter][1:]-best_value)/normalizer,
        #     alpha=0.5,
        #     label="Confidence",
        # )
    else:
        ax[j, i % 4].axvline(x=pGB[parameter], color="k", label="GB")
        ax[j, i % 4].plot(test_x_rescaled[1:], np.exp(test_y[parameter].numpy()[1:]-best_value)/normalizer, "g",linewidth=5, alpha=0.5, label="True")
        ax[j, i % 4].plot(test_x_rescaled[1:], np.exp(observed_pred_mean[parameter][1:]-best_value)/normalizer, "b", label="Mean")
        # ax[j, i % 4].fill_between(
        #     test_x_rescaled[1:],
        #     np.exp(observed_pred_mean[parameter][1:]+pred_sigma[parameter][1:]-best_value)/normalizer,
        #     np.exp(observed_pred_mean[parameter][1:]-pred_sigma[parameter][1:]-best_value)/normalizer,
        #     alpha=0.5,
        #     label="Confidence",
        # )
    ax[0, 0].legend()
    if parameter in ["Amplitude"]:
        ax[j, i % 4].set_xlabel(parameter, ha="right", va="top")
    elif parameter in ["Frequency"]:
        ax[j, i % 4].set_xlabel(parameter + " (mHz)", ha="right", va="top")
    else:
        ax[j, i % 4].set_xlabel(parameter)
    ax[j, i % 4].set_ylabel("Likelihood")
    if parameter == "EclipticLatitude":
        ax[j, i % 4].set_xlim(np.arcsin(boundaries_reduced[parameter][0]), np.arcsin(boundaries_reduced[parameter][1]))
    elif parameter == "Inclination":
        ax[j, i % 4].set_xlim(np.arccos(boundaries_reduced[parameter][1]), np.arccos(boundaries_reduced[parameter][0]))
    elif parameter in ["Frequency"]:
        ax[j, i % 4].set_xlim(boundaries_reduced[parameter][0] * 10 ** 3, boundaries_reduced[parameter][1] * 10 ** 3)
    else:
        ax[j, i % 4].set_xlim(boundaries_reduced[parameter][0], boundaries_reduced[parameter][1])
    i += 1
plt.tight_layout()
plt.show()

fig, ax = plt.subplots(2, 4, figsize=np.asarray(fig_size) * 2.5)
# plt.suptitle("loglikelihood")
i = 0
for parameter in parametersfd:
    j = 0
    if i > 3:
        j = 1
    with torch.no_grad():
        # Get upper and lower confidence bounds
        test_x_rescaled = []
        for k in range(len(test_x[parameter])):
            test_x_rescaled.append(scaletooriginal(test_x[parameter][k].numpy(), boundaries_reduced)[parameter])
        test_x2_rescaled = []
        for k in range(len(test_x2[parameter])):
            test_x2_rescaled.append(scaletooriginal(test_x2[parameter][k].numpy(), boundaries)[parameter])
        # Plot training data as black stars
        if parameter == "Frequency":
            ax[j, i % 4].axvline(x=pGB[parameter] * 10 ** 3, color="k", label="GB")
            ax[j, i % 4].plot(
                np.asarray(test_x2_rescaled[1:]) * 10 ** 3,
                test_y2[parameter].numpy()[1:],
                "g",
                label="True",
            )
            ax[j, i % 4].plot(
                np.asarray(test_x_rescaled[1:]) * 10 ** 3,
                observed_pred_mean[parameter][1:],
                "b",
                label="Mean",
            )
            # ax[j, i % 4].fill_between(
            #     np.asarray(test_x_rescaled[1:]) * 10 ** 3,
            #     observed_pred_mean[parameter][1:]+pred_sigma[parameter][1:],
            #     observed_pred_mean[parameter][1:]-pred_sigma[parameter][1:],
            #     alpha=0.5,
            #     label="Confidence",
            # )
        else:
            ax[j, i % 4].axvline(x=pGB[parameter], color="k", label="GB")
            ax[j, i % 4].plot(test_x2_rescaled[1:], test_y2[parameter].numpy()[1:], "g", label="True")
            ax[j, i % 4].plot(test_x_rescaled[1:], observed_pred_mean[parameter][1:], "b", label="Mean")
            # ax[j, i % 4].fill_between(
            #     test_x_rescaled[1:],
            #     observed_pred_mean[parameter][1:]+pred_sigma[parameter][1:],
            #     observed_pred_mean[parameter][1:]-pred_sigma[parameter][1:],
            #     alpha=0.5,
            #     label="Confidence",
            # )
        ax[0, 0].legend()
    if parameter in ["Amplitude"]:
        ax[j, i % 4].set_xlabel(parameter, ha="right", va="top")
    elif parameter in ["Frequency"]:
        ax[j, i % 4].set_xlabel(parameter + " (mHz)", ha="right", va="top")
    else:
        ax[j, i % 4].set_xlabel(parameter)
    ax[j, i % 4].set_ylabel("Log-likelihood")
    # if parameter == "EclipticLatitude":
    #     ax[j, i % 4].set_xlim(np.arcsin(boundaries_reduced[parameter][0]), np.arcsin(boundaries_reduced[parameter][1]))
    # elif parameter == "Inclination":
    #     ax[j, i % 4].set_xlim(np.arccos(boundaries_reduced[parameter][1]), np.arccos(boundaries_reduced[parameter][0]))
    # elif parameter in ["Frequency"]:
    #     ax[j, i % 4].set_xlim(boundaries_reduced[parameter][0] * 10 ** 3, boundaries_reduced[parameter][1] * 10 ** 3)
    # else:
    #     ax[j, i % 4].set_xlim(boundaries_reduced[parameter][0], boundaries_reduced[parameter][1])
    i += 1
plt.tight_layout()
plt.show()
# fig, ax = plt.subplots(2, 4,figsize=(15,15))
# plt.suptitle("loglikelihood true")
# i = 0
# for parameter in parametersfd:
#     if parameter != parameter2:
#         j = 0
#         if i > 3:
#             j = 1
#         with torch.no_grad():
#             ax[j,i%4].axvline(x=pGB01[parameter], color='k')
#             ax[j,i%4].axhline(y=pGB01[parameter2], color='k')
#             im = ax[j,i%4].scatter(test_x[parameter+parameter2].numpy()[:,i],test_x[parameter+parameter2].numpy()[:,parametersfd.index(parameter2)],c=test_y[parameter+parameter2].numpy()[:])
#             ax[j,i%4].set_xlabel(parameter)
#             ax[j,i%4].set_ylabel(parameter2)
#             fig.colorbar(im, ax=ax[j,i%4])
#     else:
#         ax[j,i%4].plot(test_x[parameter].numpy()[1:,i], test_y[parameter].numpy()[1:], 'g.')
#         ax[j,i%4].set_xlabel(parameter)
#         ax[j,i%4].set_ylabel('loglikelihood')
#         ax[j,i%4].legend(['True'])
#     i += 1

fig, ax = plt.subplots(2, 4, figsize=(15, 15))
plt.suptitle("loglikelihood predicted mean")
i = 0
for parameter in parametersfd:
    if parameter != parameter2:
        j = 0
        if i > 3:
            j = 1
        with torch.no_grad():
            # Get upper and lower confidence bounds
            test_x_rescaled1 = []
            test_x_rescaled2 = []
            for k in range(len(test_x[parameter + parameter2])):
                test_x_rescaled1.append(scaletooriginal(test_x[parameter + parameter2][k].numpy(), boundaries_reduced)[parameter])
                test_x_rescaled2.append(scaletooriginal(test_x[parameter + parameter2][k].numpy(), boundaries_reduced)[parameter2])
            ax[j, i % 4].axvline(x=pGB[parameter], color="k")
            ax[j, i % 4].axhline(y=pGB[parameter2], color="k")
            im = ax[j, i % 4].scatter(test_x_rescaled1, test_x_rescaled2, c=observed_pred_mean[parameter + parameter2][:])
            ax[j, i % 4].set_xlabel(parameter)
            ax[j, i % 4].set_ylabel(parameter2)
            fig.colorbar(im, ax=ax[j, i % 4])

            if parameter == "EclipticLatitude":
                ax[j, i % 4].set_xlim(
                    np.arcsin(boundaries_reduced[parameter][0]),
                    np.arcsin(boundaries_reduced[parameter][1]),
                )
            elif parameter == "Inclination":
                ax[j, i % 4].set_xlim(
                    np.arccos(boundaries_reduced[parameter][1]),
                    np.arccos(boundaries_reduced[parameter][0]),
                )
            elif parameter in ["Frequency"]:
                ax[j, i % 4].set_xlim(
                    boundaries_reduced[parameter][0] * 10 ** 3,
                    boundaries_reduced[parameter][1] * 10 ** 3,
                )
            else:
                ax[j, i % 4].set_xlim(boundaries_reduced[parameter][0], boundaries_reduced[parameter][1])
            if parameter2 == "EclipticLatitude":
                ax[j, i % 4].set_ylim(
                    np.arcsin(boundaries_reduced[parameter2][0]),
                    np.arcsin(boundaries_reduced[parameter2][1]),
                )
            elif parameter2 == "Inclination":
                ax[j, i % 4].set_ylim(
                    np.arccos(boundaries_reduced[parameter2][1]),
                    np.arccos(boundaries_reduced[parameter2][0]),
                )
            elif parameter2 in ["Frequency"]:
                ax[j, i % 4].set_ylim(
                    boundaries_reduced[parameter2][0] * 10 ** 3,
                    boundaries_reduced[parameter2][1] * 10 ** 3,
                )
            else:
                ax[j, i % 4].set_ylim(boundaries_reduced[parameter2][0], boundaries_reduced[parameter2][1])
    else:
        with torch.no_grad():
            lower, upper = observed_pred[parameter].confidence_region()
            lower = lower.cpu()
            upper = upper.cpu()
            lower = (lower * sigma) + nu
            upper = (upper * sigma) + nu
            ax[j, i % 4].plot(test_x[parameter].numpy()[1:, i], test_y[parameter].numpy()[1:], "g.")
            ax[j, i % 4].plot(test_x[parameter].numpy()[1:, i], observed_pred_mean[parameter].numpy()[1:], "b.")
            ax[j, i % 4].fill_between(
                test_x[parameter].numpy()[1:, i],
                lower.numpy()[1:],
                upper.numpy()[1:],
                alpha=0.5,
            )
            ax[j, i % 4].set_xlabel(parameter)
            ax[j, i % 4].set_ylabel("loglikelihood")
            ax[j, i % 4].legend(["True", "Mean", "Confidence"])
    i += 1
plt.show()


fig, ax = plt.subplots(7, 7, figsize=np.asarray(fig_size)*3)
ax[0,1].axis('off')
ax[0,2].axis('off')
ax[0,3].axis('off')
ax[0,4].axis('off')
ax[1,2].axis('off')
ax[1,3].axis('off')
ax[1,4].axis('off')
ax[2,3].axis('off')
ax[2,4].axis('off')
ax[3,4].axis('off')
i = 0
parameter = "EclipticLatitude"
parameters_reduced = [
    "EclipticLongitude",
    "EclipticLatitude",
    "Frequency",
    "Inclination",
    "Amplitude",
    'InitialPhase',
    'Polarization'
]
for parameter in parameters_reduced:
    test_x_rescaled = []
    for k in range(len(test_x[parameter])):
        test_x_rescaled.append(scaletooriginal(test_x[parameter][k].numpy(), boundaries_reduced)[parameter])
    test_x2_rescaled = []
    for k in range(len(test_x2[parameter])):
        test_x2_rescaled.append(scaletooriginal(test_x2[parameter][k].numpy(), boundaries)[parameter])
    # Plot training data as black stars
    if parameter == "Frequency":
        ax[i,i].axvline(x=pGB[parameter] * 10 ** 3, color="k", label="GB")
        ax[i,i].plot(
            np.asarray(test_x_rescaled[1:]) * 10 ** 3,
            np.exp(test_y[parameter].numpy()[1:]-best_value)/normalizer,
            "g",
            label="Target",
        )
        ax[i,i].plot(
            np.asarray(test_x_rescaled[1:]) * 10 ** 3,
            np.exp(observed_pred_mean[parameter][1:]-best_value)/normalizer,
            "b",
            label="Mean",
        )
        # ax[i,i].fill_between(
        #     np.asarray(test_x_rescaled[1:]) * 10 ** 3,
        #     np.exp(lower.numpy()[1:]-best_value)/normalizer,
        #     np.exp(upper.numpy()[1:]-best_value)/normalizer,
        #     alpha=0.5,
        #     label="Confidence",
        # )
    else:
        ax[i,i].axvline(x=pGB[parameter], color="r", label="GB")
        ax[i,i].plot(test_x_rescaled[1:], np.exp(test_y[parameter].numpy()[1:]-best_value)/normalizer, "g", label="Target")
        ax[i,i].plot(test_x_rescaled[1:], np.exp(observed_pred_mean[parameter][1:]-best_value)/normalizer, "b", label="Mean")
        # ax[i,i].fill_between(
        #     test_x_rescaled[1:],
        #     np.exp(lower.numpy()[1:]-best_value)/normalizer,
        #     np.exp(upper.numpy()[1:]-best_value)/normalizer,
        #     alpha=0.5,
        #     label="Confidence",
        # )
    ax[4, 4].legend()
    ax[4,4].set_xlabel("Amplitude", ha="right", va="top")
    ax[0,0].set_ylabel('Likelihood')
    if i != 4:
        ax[i,i].xaxis.set_ticklabels([])
    if i != 0:
        ax[i,i].yaxis.set_ticklabels([])
    # ax[i,i].set_ylabel("Likelihood")
    if parameter == "EclipticLatitude":
        ax[i,i].set_xlim(np.arcsin(boundaries_reduced[parameter][0]), np.arcsin(boundaries_reduced[parameter][1]))
    elif parameter == "Inclination":
        ax[i,i].set_xlim(np.arccos(boundaries_reduced[parameter][1]), np.arccos(boundaries_reduced[parameter][0]))
    elif parameter in ["Frequency"]:
        ax[i,i].set_xlim(boundaries_reduced[parameter][0] * 10 ** 3, boundaries_reduced[parameter][1] * 10 ** 3)
    else:
        ax[i,i].set_xlim(boundaries_reduced[parameter][0], boundaries_reduced[parameter][1])
    i += 1
parameters_reduced = [
    "EclipticLatitude",
    "Frequency",
    "Inclination",
    "Amplitude",
    'InitialPhase',
    'Polarization'
]
parameter = "EclipticLongitude"
for n in range(len(parameters_reduced)):
    print(n)
    print(parameter2)
    i = n+1
    for parameter2 in parameters_reduced:
        test_x_rescaled1 = []
        test_x_rescaled2 = []
        for k in range(len(test_x[parameter + parameter2])):
            test_x_rescaled1.append(scaletooriginal(test_x[parameter + parameter2][k].numpy(), boundaries_reduced)[parameter])
            test_x_rescaled2.append(scaletooriginal(test_x[parameter + parameter2][k].numpy(), boundaries_reduced)[parameter2])
        ax[i,n].axvline(x=pGB[parameter], color="r")
        ax[i,n].axhline(y=pGB[parameter2], color="r")
        ax[i,n].scatter(test_x_rescaled1, test_x_rescaled2, c=np.exp(observed_pred_mean[parameter + parameter2][:]-best_value)/normalizer)
        if parameter == 'Frequency':
            ax[i,n].scatter(np.asarray(test_x_rescaled1)*1e3, test_x_rescaled2, c=np.exp(observed_pred_mean[parameter + parameter2][:]-best_value)/normalizer)
            ax[i,n].axvline(x=pGB[parameter]*1e3, color="r")
        if parameter2 == 'Frequency':
            ax[i,n].scatter(test_x_rescaled1, np.asarray(test_x_rescaled2)*1e3, c=np.exp(observed_pred_mean[parameter + parameter2][:]-best_value)/normalizer)
        ax[i,n].axhline(y=pGB[parameter2]*1e3, color="r")
        if i == 6:
            ax[i,n].set_xlabel(parameter)
            if parameter in ["Frequency"]:
                ax[i,n].set_xlabel(parameter + " (mHz)", ha="right", va="top" )
                ax[i,n].xaxis.set_major_locator(plt.MaxNLocator(2))
        else:
            ax[i,n].xaxis.set_ticklabels([])
        if n == 0:
            ax[i,n].set_ylabel(parameter2)
            if parameter in ["Frequency"]:
                ax[i,n].set_ylabel(parameter + " (mHz)")
        else:
            ax[i,n].yaxis.set_ticklabels([])
        # fig.colorbar(im, ax=ax[i,n])

        if parameter == "EclipticLatitude":
            ax[i,n].set_xlim(
                np.arcsin(boundaries_reduced[parameter][0]),
                np.arcsin(boundaries_reduced[parameter][1]),
            )
        elif parameter == "Inclination":
            ax[i,n].set_xlim(
                np.arccos(boundaries_reduced[parameter][1]),
                np.arccos(boundaries_reduced[parameter][0]),
            )
        elif parameter in ["Frequency"]:
            ax[i,n].set_xlim(
                boundaries_reduced[parameter][0] * 10 ** 3,
                boundaries_reduced[parameter][1] * 10 ** 3,
            )
        else:
            ax[i,n].set_xlim(boundaries_reduced[parameter][0], boundaries_reduced[parameter][1])
        if parameter2 == "EclipticLatitude":
            ax[i,n].set_ylim(
                np.arcsin(boundaries_reduced[parameter2][0]),
                np.arcsin(boundaries_reduced[parameter2][1]),
            )
        elif parameter2 == "Inclination":
            ax[i,n].set_ylim(
                np.arccos(boundaries_reduced[parameter2][1]),
                np.arccos(boundaries_reduced[parameter2][0]),
            )
        elif parameter2 in ["Frequency"]:
            ax[i,n].set_ylim(
                boundaries_reduced[parameter2][0] * 10 ** 3,
                boundaries_reduced[parameter2][1] * 10 ** 3,
            )
        else:
            ax[i,n].set_ylim(boundaries_reduced[parameter2][0], boundaries_reduced[parameter2][1])
        i += 1
    parameter = parameters_reduced[0]
    parameters_reduced.pop(0)
plt.show()

prediction = observed_pred["Frequency"].mean.cpu().numpy()
n_bin = 50
i = 0
fig, axes = plt.subplots(2, 4, figsize=(15, 15))
plt.suptitle("sampled posterior")
for parameter in parameters:
    j = 0
    if i > 3:
        j = 1
    axes[j, i % 4].axvline(x=pGB[parameter], color="r")
    axes[j, i % 4].plot(samples[parameter], samples["Likelihood"], ".", label="train")
    axes[j, i % 4].plot(test_samples[parameter], test_samples["Likelihood"], ".", label="true")
    axes[j, i % 4].errorbar(
        test_samples[parameter][1:],
        prediction[1:],
        prediction[1:] - lower.numpy()[1:],
        fmt=".",
        label="prediction",
    )
    # plt.fill_between(test_samples[parameter][1:],prediction[1:]-p_std[1:],prediction[1:]+p_std[1:], color='g', alpha= 0.4)
    if i == 0:
        axes[j, i % 4].set_ylabel("log-Likelihood")
    # plt.subplot(3,number_of_parameters,i+number_of_parameters)
    # plt.axvline(x=pGB[parameter], color='r')
    # n, bins, patches = plt.hist(samples[parameter], n_bin, density=True, facecolor='k', alpha=0.5)

    # plt.subplot(3,number_of_parameters,i+2*number_of_parameters)
    # if parameter == 'Frequency':
    #     plt.axvline(x=pGB[parameter]*1000, color='r')
    #     plt.plot(samples[parameter]*1000, range(number_of_samples), 'k')
    # else:
    #     plt.axvline(x=pGB[parameter], color='r')
    #     plt.plot(samples[parameter], range(number_of_samples), 'k')
    axes[j, i % 4].set_xlabel(parameter)
    i += 1
plt.legend()

name = "Frequency"
plt.figure(figsize=(10, 8))
plt.suptitle("sampled posterior")
plt.subplot(1, 1, 1)
plt.axvline(x=pGB[name], color="r")
plt.plot(samples[name], samples["Likelihood"], ".", label="train", zorder=1)
plt.plot(test_samples[name], test_samples["Likelihood"], ".", label="true")
plt.plot(test_samples[name][1:], prediction[1:], label="prediction")
plt.fill_between(test_samples[name][1:], lower.numpy()[1:], upper.numpy()[1:], color="g", alpha=0.4)
plt.ylabel("log-Likelihood")
plt.xlabel(name)
plt.legend()
# plt.figure()
# plt.scatter(samples['Amplitude'],samples['Frequency'],c=samples['Likelihood'])
# plt.xlim(max(samples['Amplitude']),max(samples['Amplitude']))
# plt.ylim(max(samples['Frequency']),max(samples['Frequency']))
# plt.colorbar()
# plt.figure()
# # plt.title('sqrt(MSE)',np.round(np.sqrt(mean_squared_error(test_y,prediction)),2))
# plt.scatter(samples['EclipticLatitude'],samples['EclipticLongitude'],c=samples['Likelihood'])
# plt.xlim(max(samples['EclipticLatitude']),max(samples['EclipticLatitude']))
# plt.ylim(max(samples['EclipticLongitude']),max(samples['EclipticLongitude']))
# plt.colorbar()
plt.show()


# %%
