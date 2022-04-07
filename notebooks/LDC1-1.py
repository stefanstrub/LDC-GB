#%%
from re import escape
from gpytorch.likelihoods import likelihood
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rcParams
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
import scipy
from scipy import misc
from scipy.stats import multivariate_normal
from scipy.optimize import rosen, differential_evolution
import numpy as np
from six import b
import xarray as xr
from astropy import units as u
import pandas as pd
import time
from copy import deepcopy
import multiprocessing as mp
from functools import partial
import itertools

import Cosmology

import torch
import gpytorch
from sklearn.metrics import mean_squared_error
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel, Matern, RationalQuadratic, ExpSineSquared, RBF
from sklearn.kernel_approximation import Nystroem
from sklearn import linear_model, pipeline
from sklearn import preprocessing

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

from fastkde import fastKDE
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
fig_width_pt = 1.5*464.0  # Get this from LaTeX using \showthe\columnwidth
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


def semi_fast_tdi(config, pMBHBi, t_min, t_max, dt):
    hphc = HpHc.type("MBHB-%d"%s_index, "MBHB", "IMRPhenomD")
    hphc.set_param(pMBHBi)
    orbits = Orbits.type(config)
    P = ProjectedStrain(orbits)    
    yArm = P.arm_response(t_min, t_max, dt, [hphc], tt_order=1)
    X = P.compute_tdi_x(np.arange(t_min, t_max, dt))
    return TimeSeries(X, dt=dt)

def Window(tm, offs=1000.0):
    xl = offs
    ind_r = np.argwhere(tm[-1] - tm <= offs)[0][0]
    xr = tm[ind_r]
    kap = 0.005
    winl = 0.5 * (1.0 + np.tanh(kap * (tm - xl)))
    winr = 0.5 * (1.0 - np.tanh(kap * (tm - xr)))
    return winl * winr

def loglikelihoodflat(pGBs):
    for i in range(len(pGBs)):
        Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs[i], oversample=4, simulator="synthlisa")
        index_low = np.searchsorted(Xs.f, dataX.f[0])
        if i == 0:
            Xs_total = Xs[index_low : index_low + len(dataX)]
            Ys_total = Ys[index_low : index_low + len(dataY)]
            Zs_total = Zs[index_low : index_low + len(dataZ)]
        else:
            Xs_total += Xs[index_low : index_low + len(dataX)]
            Ys_total += Ys[index_low : index_low + len(dataY)]
            Zs_total += Zs[index_low : index_low + len(dataZ)]
        
    Af = (Zs_total - Xs_total)/np.sqrt(2.0)
    Ef = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)
    SNR2 = np.sum( np.real(DAf * np.conjugate(Af.data) + DEf * np.conjugate(Ef.data))/SA )
    hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /SA)
    dd = np.sum((np.absolute(DAf.data)**2 + np.absolute(DEf.data)**2) /SA)
    plotIt = False
    if plotIt:
        fig, ax = plt.subplots(nrows=2, sharex=True) 
        ax[0].plot(Af.f, np.abs(DAf))
        ax[0].plot(Af.f, np.abs(Af.data))
        
        ax[1].plot(Af.f, np.abs(DEf))
        ax[1].plot(Af.f, np.abs(Ef.data))
        plt.show()
        
    p2 = np.sum((np.absolute(DAf - Af.data)**2 + np.absolute(DEf - Ef.data)**2) /SA) * Xs.df *2
    diff = np.abs(DAf - Af.data) ** 2 + np.abs(DEf - Ef.data) ** 2
    p1 = float(np.sum(diff / SA) * Xs.df) / 2.0
    loglik = 4.0*Xs.df*( SNR2 - 0.5 * hh - 0.5 * dd)
    print(p2, loglik)
    return (loglik) 

def loglikelihood10(pGBs):
    for i in range(len(pGBs)):
        Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs[i], oversample=4, simulator="synthlisa")
        index_low = np.searchsorted(Xs.f, dataX.f[0])
        if i == 0:
            Xs_total = Xs[index_low : index_low + len(dataX)]
            Ys_total = Ys[index_low : index_low + len(dataY)]
            Zs_total = Zs[index_low : index_low + len(dataZ)]
        else:
            Xs_total += Xs[index_low : index_low + len(dataX)]
            Ys_total += Ys[index_low : index_low + len(dataY)]
            Zs_total += Zs[index_low : index_low + len(dataZ)]

    diff = np.abs(dataX - Xs_total.values*100) ** 2 + np.abs(dataY - Ys_total.values*100) ** 2 + np.abs(dataZ - Zs_total.values*100) ** 2
    # p1 = -float(np.sum(diff / Sn)*Xs.attrs['df'])/2.0
    p1 = float(np.sum(diff / Sn) * Xs_total.df) / 2.0
    # p1 = np.exp(p1)
    return -p1/10000



def sampler(number_of_samples, parameters, pGB, boundaries, uniform=False, MCMC=False, only=False, onlyparameter="Frequency", twoD=False, secondparameter="Amplitude",
    calculate_loglikelihood=True, gaussian=False, pGB_centered=True, std=0.2):
    samples = {}
    pGBs01 = {}
    pGBs = {}
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
    psample = loglikelihood([pGB])
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
            psample = loglikelihood([pGBs])
            samples["Likelihood"].append(psample)
        for parameter in parameters:
            samples[parameter].append(pGBs01[parameter])
    for parameter in parameters:
        samples[parameter] = np.asarray(samples[parameter])
    samples["Likelihood"] = np.asarray(samples["Likelihood"])
    return samples
#%%
def objective(trial,changeableparameters, maxpGB2, signal):
    for parameter in changeableparameters:
        parametervalue = trial.suggest_uniform(parameter, boundaries[parameter][0], boundaries[parameter][1])
        maxpGB2[signal][parameter] = parametervalue
        if parameter in ["EclipticLatitude"]:
            maxpGB2[signal][parameter] = np.arcsin(parametervalue)
        elif parameter in ["Inclination"]:
            maxpGB2[signal][parameter] = np.arccos(parametervalue)
        elif parameter in ['Amplitude',"FrequencyDerivative"]:
            maxpGB2[signal][parameter] = 10**(parametervalue)
        # elif parameter in ['Amplitude']:
        #     maxpGB2[signal][parameter] = 10**(parametervalue)/100
    p = loglikelihood(maxpGB2)
    return p


def changeparameterstoinput(pMBHB01, boundaries):
    pMBHBnew = deepcopy(pMBHBs)
    for parameter in parametersboundary:
        if parameter in parameters:
            pMBHBnew[parameter] = (pMBHB01[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
        elif parameter in ['ChirpMass']:
            ChirpMass = (pMBHB01[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
        elif parameter in ['MassRatio']:
            MassRatio = (pMBHB01[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
        elif parameter in ['sinEclipticLatitude']:
            pMBHBnew['EclipticLatitude'] = np.arcsin((pMBHB01[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
        else:
            print('parameter missed',parameter)
    pMBHBnew['Mass1'], pMBHBnew['Mass2'] = funcm1m2ofMchirpq(ChirpMass,MassRatio)
    # pMBHBnew['Distance'] = 10**pMBHB01['Distance']
    # v = hubbleconstant*(10**3*pMBHBnew['Distance'])
    # pMBHBnew['Redshift'] = np.sqrt((1+v/speedoflight)/(1-v/speedoflight)) - 1
    return pMBHBnew   

def scaletooriginal(previous_max, boundaries):
    maxpGB = {}
    for parameter in parameters:
        if parameter in ["EclipticLatitude"]:
            maxpGB[parameter] = np.arcsin((previous_max[parameters.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
        elif parameter in ["Inclination"]:
            maxpGB[parameter] = np.arccos((previous_max[parameters.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
        elif parameter in ['Amplitude',"FrequencyDerivative"]:
            maxpGB[parameter] = 10**((previous_max[parameters.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
        else:
            maxpGB[parameter] = (previous_max[parameters.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]
    return maxpGB

def scaletooriginalparameter(previous_max, boundaries):
    maxpGB = {}
    for parameter in parameters:
        if parameter in ["EclipticLatitude"]:
            maxpGB[parameter] = np.arcsin((previous_max[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
        elif parameter in ["Inclination"]:
            maxpGB[parameter] = np.arccos((previous_max[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
        elif parameter in ['Amplitude',"FrequencyDerivative"]:
            maxpGB[parameter] = 10**((previous_max[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
        else:
            maxpGB[parameter] = (previous_max[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]
    return maxpGB

def scaleto01(previous_max, boundaries):
    maxpGB = {}
    for parameter in parameters:
        if parameter in ["EclipticLatitude"]:
            maxpGB[parameter] = (np.sin(previous_max[parameter]) - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0])
        elif parameter in ["Inclination"]:
            maxpGB[parameter] = (np.cos(previous_max[parameter]) - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0])
        elif parameter in ['Amplitude',"FrequencyDerivative"]:
            maxpGB[parameter] = (np.log10(previous_max[parameter]) - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0])
        else:
            maxpGB[parameter] = (previous_max[parameter] - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0])
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
    # gpr.score(train_x, train_y)
    return gpr


def CoordinateMC(n, pGBs, boundaries, parameters_recorded, loglikelihood):

    maxpGB = []
    parameters_recorded1 = []
    no_improvement_counter = 0
    for i in range(number_of_signals):
        maxpGB.append(deepcopy(pGBs))
        parameters_recorded1.append({})
        for parameter in parameters:
            parametervalue = np.random.uniform(boundaries[parameter][0], boundaries[parameter][1])
            maxpGB[i][parameter] = parametervalue
            if parameter in ["EclipticLatitude"]:
                maxpGB[i][parameter] = np.arcsin(parametervalue)
            elif parameter in ["Inclination"]:
                maxpGB[i][parameter] = np.arccos(parametervalue)
            elif parameter in ['Amplitude',"FrequencyDerivative"]:
                maxpGB[i][parameter] = 10**(parametervalue)
            parameters_recorded1[i][parameter] = []
            parameters_recorded1[i][parameter].append(maxpGB[i][parameter])
        parameters_recorded1[i]["Loglikelihood"] = []
    # if n == 0:
    #     maxpGB = deepcopy(pGBs)
    try:
        best_value
    except:
        best_value = loglikelihood(maxpGB)
    previous_best = deepcopy(best_value)
    maxpGB2 = deepcopy(maxpGB)
    number_per_row = 5
    n_trials = number_per_row**2
    n_trials = 50
    for i in range(100):
        # if i > 48:
        #     n_trials = 50
        parameter1 = parameters[i % 8]
        parameter2 = parameters[np.random.randint(0, 7)]
        parameter3 = parameters[np.random.randint(0, 7)]
        # parameter2 = 'InitialPhase'
        # parameter1 = 'Inclination'
        while parameter2 == parameter1:
            parameter2 = parameters[np.random.randint(0, 7)]
        while parameter3 == parameter1 or parameter3 == parameter2:
            parameter3 = parameters[np.random.randint(0, 7)]
        # if parameter1 == 'Frequency':
        #     parameter2 = 'Polarization'
        changeableparameters = [parameter1, parameter2, parameter3]
        signal = i % number_of_signals
        for j in range(n_trials):
            k = 0
            for parameter in changeableparameters:
                parametervalue = np.random.uniform(boundaries[parameter][0], boundaries[parameter][1],1)[0]
                # if k == 0:
                #     parametervalue = np.random.uniform(j%number_per_row * 1/number_per_row, (j%number_per_row+1)*1/number_per_row)
                # else:
                #     parametervalue = np.random.uniform(int(j/number_per_row) * 1/number_per_row, (int(j/number_per_row)+1)*1/number_per_row)
                # parametervalue = parametervalue * (boundaries[parameter][1] - boundaries[parameter][0]) + boundaries[parameter][0]
                maxpGB2[signal][parameter] = parametervalue
                if parameter in ["EclipticLatitude"]:
                    maxpGB2[signal][parameter] = np.arcsin(parametervalue)
                elif parameter in ["Inclination"]:
                    maxpGB2[signal][parameter] = np.arccos(parametervalue)
                elif parameter in ['Amplitude',"FrequencyDerivative"]:
                    maxpGB2[signal][parameter] = 10**(parametervalue)
                k += 1
            suggestion = loglikelihood(maxpGB2)
            if suggestion > previous_best:
                no_improvement_counter = 0
                previous_best = suggestion
                for parameter in changeableparameters:
                    maxpGB[signal][parameter] = maxpGB2[signal][parameter]
        if no_improvement_counter != 0:
            no_improvement_counter += 1
        if no_improvement_counter > 16:
            past_mean = 0
            sum_count = 0
            for l in range(n):
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

        if i in [30,40,50,60,70]:
            past_mean = 0
            sum_count = 0
            for l in range(n):
                try:
                    past_mean += parameters_recorded[l][0]["Loglikelihood"][i]
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

        # start = time.time()
        # print(n, i, previous_best, loglikelihood(maxpGB), maxpGB)
        parameters_recorded1[0]["Loglikelihood"].append(loglikelihood(maxpGB))
        for i in range(number_of_signals):
            for parameter in parameters:
                parameters_recorded1[i][parameter].append(maxpGB[i][parameter])

        maxpGB2 = deepcopy(maxpGB)
    if previous_best > best_value:
        best_value = previous_best
    parameters_recorded[n] = parameters_recorded1
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
                np.log10(maxpGB[parameter]) - length * ratio / 2*2,
                np.log10(maxpGB[parameter]) + length * ratio / 2*2,
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

class Search():
    def __init__(self,tdi_fs,Tobs) -> None:
        self.tdi_fs = tdi_fs
        self.GB = fastGB.FastGB(delta_t=dt, T=Tobs)
        self.reduced_frequency_boundaries = None

        amplitude = np.sqrt(np.max(psd))
        frequencyrange =  [lower_frequency-padding, upper_frequency+padding]
        indexes = np.argsort(p.get('Frequency'))
        index_low = np.searchsorted(p.get('Frequency')[indexes], frequencyrange[0])
        index_high = np.searchsorted(p.get('Frequency')[indexes], frequencyrange[1])
        strongest_source_in_window = np.argmax(p.get('Amplitude')[indexes][index_low:index_high])
        self.pGB = {}
        for parameter in parameters:
            self.pGB[parameter] = p.get(parameter)[indexes][index_low:index_high][strongest_source_in_window]

        # pGB = deepcopy(pGBadded)
        # self.pGB = {'Amplitude': 3.676495e-22, 'EclipticLatitude': 0.018181, 'EclipticLongitude': 1.268061, 'Frequency': 0.01158392, 'FrequencyDerivative': 8.009579e-15, 'Inclination': 0.686485, 'InitialPhase': 4.201455, 'Polarization': 2.288223}
        

        self.boundaries = {
            "Amplitude": [np.log10(amplitude)-3,np.log10(amplitude)-1],
            # "Amplitude": [-23.5,-21],
            # "Amplitude": [np.log10(self.pGB['Amplitude'])-2,np.log10(self.pGB['Amplitude'])+1],
            "EclipticLatitude": [-1.0, 1.0],
            "EclipticLongitude": [-np.pi, np.pi],
            # "Frequency": [self.pGB["Frequency"] * 0.99995, self.pGB["Frequency"] * 1.00015],
            # "Frequency": [self.pGB["Frequency"] - 3e-7, self.pGB["Frequency"] + 3e-7],
            "Frequency": frequencyrange,
            "FrequencyDerivative": [-20.0,-13.0],
            # "FrequencyDerivative": [np.log10(5e-6*self.pGB['Frequency']**(13/3)),np.log10(8e-8*self.pGB['Frequency']**(11/3))],
            "Inclination": [-1.0, 1.0],
            "InitialPhase": [0.0, 2.0 * np.pi],
            "Polarization": [0.0, 1.0 * np.pi],
        }
        if self.boundaries['FrequencyDerivative'][0] > self.boundaries['FrequencyDerivative'][1]:
            c = self.boundaries['FrequencyDerivative'][0]
            self.boundaries['FrequencyDerivative'][0] = self.boundaries['FrequencyDerivative'][1]
            self.boundaries['FrequencyDerivative'][1] = c

        print('amplitude boundaries',amplitude, 10**(self.boundaries['Amplitude'][0]), 10**(self.boundaries['Amplitude'][1]))
        previous_max = np.random.rand(8)
        # previous_max[0] = np.random.rand(1)*0.1 +0.5
        # previous_max[3] = np.random.rand(1)*0.1 +0.5
        i = 0
        self.pGBs = deepcopy(self.pGB)
        for parameter in parameters:
            # if parameter in ["FrequencyDerivative"]:
            #     i -= 1
            if parameter in ["EclipticLatitude"]:
                self.pGBs[parameter] = np.arcsin((previous_max[i] * (self.boundaries[parameter][1] - self.boundaries[parameter][0])) + self.boundaries[parameter][0])
            elif parameter in ["Inclination"]:
                self.pGBs[parameter] = np.arccos((previous_max[i] * (self.boundaries[parameter][1] - self.boundaries[parameter][0])) + self.boundaries[parameter][0])
            elif parameter in ['Amplitude',"FrequencyDerivative"]:
                self.pGBs[parameter] = 10**((previous_max[i] * (self.boundaries[parameter][1] - self.boundaries[parameter][0])) + self.boundaries[parameter][0])
            else:
                self.pGBs[parameter] = (previous_max[i] * (self.boundaries[parameter][1] - self.boundaries[parameter][0])) + self.boundaries[parameter][0]
            i += 1
        # self.pGBs = {'Amplitude': 4.0900673126042746e-22, 'EclipticLatitude': 0.8718477251317046, 'EclipticLongitude': 0.48599945403230693, 'Frequency': 0.003995220986111426, 'FrequencyDerivative': 1.0985841703423861e-16, 'Inclination': 1.0262955111380103, 'InitialPhase': 5.453865686076588, 'Polarization': 1.089057196561609}

        # cutoff_ratio = 1000
        # Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=self.pGBs, oversample=4, simulator="synthlisa")
        # psd_signal = np.abs(Xs.values) ** 2 + np.abs(Ys.values) ** 2 + np.abs(Zs.values) ** 2
        # highSNR = psd_signal > np.max(psd_signal) / cutoff_ratio
        # lowerindex = np.where(highSNR)[0][0] - 30
        # higherindex = np.where(highSNR)[0][-1] + 30
        # self.dataX = tdi_fs["X"].isel(f=slice(Xs.kmin, Xs.kmin + len(Xs)))[lowerindex:higherindex]
        # self.dataY = tdi_fs["Y"].isel(f=slice(Ys.kmin, Ys.kmin + len(Ys)))[lowerindex:higherindex]
        # self.dataZ = tdi_fs["Z"].isel(f=slice(Zs.kmin, Zs.kmin + len(Zs)))[lowerindex:higherindex]
        indexes = np.logical_and(tdi_fs['X'].f > frequencyrange[0]-padding, tdi_fs['X'].f < frequencyrange[1]+padding) 
        self.dataX = tdi_fs["X"][indexes]
        self.dataY = tdi_fs["Y"][indexes]
        self.dataZ = tdi_fs["Z"][indexes]

        self.DAf = (self.dataZ - self.dataX)/np.sqrt(2.0)
        self.DEf = (self.dataZ - 2.0*self.dataY + self.dataX)/np.sqrt(6.0)
        # self.DAf = (2/3*self.dataX - self.dataY - self.dataZ)/3.0
        # self.DEf = (self.dataZ - self.dataY)/np.sqrt(3.0)

        # Xs, Ys, Zs = (
        #     Xs[lowerindex:higherindex],
        #     Ys[lowerindex:higherindex],
        #     Zs[lowerindex:higherindex],
        # )
        spd_data = np.abs(self.dataX) ** 2 + np.abs(self.dataY) ** 2 + np.abs(self.dataZ) ** 2
        noise = (np.mean(spd_data[:2]) + np.mean(spd_data[-2:])).values / 2
        noise = 0  # (np.mean(spd_data).values)/2
        fmin, fmax = float(self.dataX.f[0]), float(self.dataX.f[-1] + self.dataX.attrs["df"])
        freq = np.array(self.dataX.sel(f=slice(fmin, fmax)).f)
        Nmodel = get_noise_model(noise_model, freq)
        self.Sn = Nmodel.psd(freq=freq, option="X")
        self.SA = Nmodel.psd(freq=freq, option="A")
        self.SE = Nmodel.psd(freq=freq, option="E")
        # diff = np.abs(self.dataX - Xs.values) ** 2 + np.abs(self.dataY - Ys.values) ** 2 + np.abs(self.dataZ - Zs.values) ** 2
        # p1 = float(np.sum(diff / (self.Sn + noise)) * Xs.df) / 2.0
        # p1 = -p1
        # diff = np.abs(self.dataX) ** 2 + np.abs(self.dataY) ** 2 + np.abs(self.dataZ) ** 2
        # null_pGBs = deepcopy(self.pGBs)
        # null_pGBs['Amplitude'] = 4*10**-25
        print('pGB', self.pGB, signal_peak, self.loglikelihood([self.pGB]))

    def f_statistic(self, lower_frequency, upper_frequency, N_frequency, N_sky):
        # We construct a global proposal density using the single
        # source F statistic to compute the individual likelihoods
        # ln pðdjλ i Þ maximized over the extrinsic parameters A, φ 0 , ι,
        # ψ.
        F_stat = []
        frequency = []
        eclipticlatitude = []
        eclipticlongitude = []
        pGBf = {}
        for parameter in parameters:
            pGBf[parameter] = 0
        pGBf['Amplitude'] = 1e-24
        pGBf['FrequencyDerivative'] = 0
        frequency_boundaries = [lower_frequency,upper_frequency]
        for n in range(N_sky):
            eclipticlatitude.append(self.boundaries['EclipticLatitude'][0]+(self.boundaries['EclipticLatitude'][1]-self.boundaries['EclipticLatitude'][0])*n/N_sky)
            eclipticlongitude.append(self.boundaries['EclipticLongitude'][0]+(self.boundaries['EclipticLongitude'][1]-self.boundaries['EclipticLongitude'][0])*n/N_sky)
        for k in range(N_frequency):
            F_stat.append([])
            frequency.append(frequency_boundaries[0] + (frequency_boundaries[1]-frequency_boundaries[0])*k/(N_frequency-1))
            for l in range(N_sky):
                F_stat[-1].append([])
                for m in range(N_sky):
                    F_stat[-1][-1].append(self.F(frequency[-1],eclipticlatitude[l],eclipticlongitude[m],pGBf))
        F_stat = np.asarray(F_stat)
        return F_stat, frequency, eclipticlatitude, eclipticlongitude

    def F(self, f0, theta, phi, pGBs):
        g = []
        pGBs['Frequency'] = f0
        pGBs['EclipticLongitude'] = theta
        pGBs['EclipticLatitude'] = phi
        pGBs['InitialPhase'] = 0
        pGBs['Inclination'] = np.pi/2
        pGBs['Polarization'] = 0
        g.append(self.GB.get_fd_tdixyz(template=pGBs, oversample=4, simulator="synthlisa"))
        pGBs['InitialPhase'] = np.pi
        pGBs['Inclination'] = np.pi/2
        pGBs['Polarization'] = np.pi/4
        g.append(self.GB.get_fd_tdixyz(template=pGBs, oversample=4, simulator="synthlisa"))
        pGBs['InitialPhase'] = 3*np.pi/2
        pGBs['Inclination'] = np.pi/2
        pGBs['Polarization'] = 0
        g.append(self.GB.get_fd_tdixyz(template=pGBs, oversample=4, simulator="synthlisa"))
        pGBs['InitialPhase'] = np.pi/2
        pGBs['Inclination'] = np.pi/2
        pGBs['Polarization'] = np.pi/4
        g.append(self.GB.get_fd_tdixyz(template=pGBs, oversample=4, simulator="synthlisa"))
        g2 = []
        for i in range(4):
            g2.append([])
            for j in range(3):
                g2[i].append(xr.align(self.dataX, g[i][j], join='left',fill_value=0)[1])
        g = g2
        data = [self.dataX,self.dataY,self.dataZ]
        f = 0
        for i in range(4):
            for j in range(4):
                if i != j:
                    f += 1/2* self.scalarproduct(g[i],g[j])**(-1)*self.scalarproduct(data,g[j])*self.scalarproduct(data,g[j])
        return f

    def scalarproduct(self, a, b):
        diff = np.real(a[0].values * np.conjugate(b[0].values)) ** 2 + np.real(a[1].values * np.conjugate(b[1].values)) ** 2 + np.real(a[2].values * np.conjugate(b[2].values)) ** 2
        res = 4*float(np.sum(diff / self.Sn) * self.dataX.df)
        return res


    def plot(self, maxpGBs=None):
        plt.figure(figsize=fig_size)
        ax1 = plt.subplot(111)
        # plt.plot(dataX_training.f*1000,dataX_training.values, label='data')
        # ax1.plot(self.dataX.f * 1000, self.dataX.values.real, label="data", marker="o", zorder=5)
        Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=self.pGB, oversample=4, simulator="synthlisa")
        index_low = np.searchsorted(Xs.f, self.dataX.f[0])
        Xs = Xs[index_low : index_low + len(self.dataX)]
        Ys = Ys[index_low : index_low + len(self.dataY)]
        Zs = Zs[index_low : index_low + len(self.dataZ)]
        # Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=self.pGB, oversample=8, simulator="synthlisa")
        # index_low = np.searchsorted(Xs.f, self.dataX.f[0])
        # Xs = Xs[index_low : index_low + len(self.dataX)]
        # Ys = Ys[index_low : index_low + len(self.dataY)]
        # Zs = Zs[index_low : index_low + len(self.dataZ)]
        # ax1.plot(Xs.f * 1000, Xs.values.real, label="VGB2", marker=".", zorder=5)
                    
        Af = (Zs - Xs)/np.sqrt(2.0)

        ax1.plot(Af.f* 1000, np.abs(self.DAf), label='Data')
        ax1.plot(Af.f* 1000, np.abs(Af.data), marker='.', label='Injection')

        # ax1.plot(Xs_added2.f * 1000, Xs_added2.values.real, label="VGB2", marker=".", zorder=5)
        ax1.axvline(self.boundaries['Frequency'][0]* 1000, color= 'red', label='Boundaries')
        ax1.axvline(self.boundaries['Frequency'][1]* 1000, color= 'red')
        if self.reduced_frequency_boundaries != None:
            ax1.axvline(self.reduced_frequency_boundaries[0]* 1000, color= 'green', label='Reduced Boundaries')
            ax1.axvline(self.reduced_frequency_boundaries[1]* 1000, color= 'green')

        # ax1.plot(Xs.f * 1000, dataX.values.real - Xs.values.real, label="residual", alpha=0.8, color="red", marker=".")
        plt.xlabel('f [mHz]')
        plt.ylabel('|A|')
        plt.legend()
        plt.show()
        # print("p true", self.loglikelihood([pGB]), "null hypothesis", self.loglikelihood([null_pGBs]))


    def SNR(self, pGBs):
        for i in range(len(pGBs)):
            Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=pGBs[i], oversample=4, simulator="synthlisa")
            index_low = np.searchsorted(Xs.f, self.dataX.f[0])
            if i == 0:
                Xs_total = Xs[index_low : index_low + len(self.dataX)]
                Ys_total = Ys[index_low : index_low + len(self.dataY)]
                Zs_total = Zs[index_low : index_low + len(self.dataZ)]
            else:
                Xs_total += Xs[index_low : index_low + len(self.dataX)]
                Ys_total += Ys[index_low : index_low + len(self.dataY)]
                Zs_total += Zs[index_low : index_low + len(self.dataZ)]
            if len(Xs_total) < len(self.dataX):
                a,Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)
                a,Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)
                a,Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)

        # Af = (Zs_total - Xs_total)/np.sqrt(2.0)
        # Ef = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)
        # diff = np.abs(self.DAf - Af.values) ** 2 + np.abs(self.DEf - Ef.values) ** 2
        # diff = np.abs(Xs_total.values) ** 2 + np.abs( Ys_total.values) ** 2 + np.abs(Zs_total.values) ** 2
        # p1 = float(np.sum(diff / self.Sn) * Xs_total.df)
        source = dict({"X": Xs, "Y": Ys, "Z": Zs})
        SNR = compute_tdi_snr(source=source, noise=Nmodel)
        # p1 = np.exp(p1)
        sum_SNR = SNR['X2'] + SNR['Y2'] + SNR['Z2']
        return sum_SNR#/10000

    def loglikelihoodXYZ(self, pGBs):
        for i in range(len(pGBs)):
            Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=pGBs[i], oversample=4, simulator="synthlisa")
            index_low = np.searchsorted(Xs.f, self.dataX.f[0])
            if i == 0:
                Xs_total = Xs[index_low : index_low + len(self.dataX)]
                Ys_total = Ys[index_low : index_low + len(self.dataY)]
                Zs_total = Zs[index_low : index_low + len(self.dataZ)]
            else:
                Xs_total += Xs[index_low : index_low + len(self.dataX)]
                Ys_total += Ys[index_low : index_low + len(self.dataY)]
                Zs_total += Zs[index_low : index_low + len(self.dataZ)]
            if len(Xs_total) < len(self.dataX):
                a,Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)
                a,Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)
                a,Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)

        # Af = (Zs_total - Xs_total)/np.sqrt(2.0)
        # Ef = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)
        # diff = np.abs(self.DAf - Af.values) ** 2 + np.abs(self.DEf - Ef.values) ** 2
        diff = np.abs(self.dataX - Xs_total.values) ** 2 + np.abs(self.dataY - Ys_total.values) ** 2 + np.abs(self.dataZ - Zs_total.values) ** 2
        # p1 = -float(np.sum(diff / Sn)*Xs.attrs['df'])/2.0
        p1 = float(np.sum(diff / self.Sn) * Xs_total.df) / 2.0
        # p1 = np.exp(p1)
        return -p1#/10000

    def loglikelihood(self,pGBs, plotIt = False):
        for i in range(len(pGBs)):
            Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=pGBs[i], oversample=4, simulator="synthlisa")
            index_low = np.searchsorted(Xs.f, self.dataX.f[0])
            if i == 0:
                Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
                Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
                Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            else:
                Xs_total += xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
                Ys_total += xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
                Zs_total += xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            
        Af = (Zs_total - Xs_total)/np.sqrt(2.0)
        Ef = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)
        # Af = (2/3*Xs_total - Ys_total - Zs_total)/3.0
        # Ef = (Zs_total - Ys_total)/np.sqrt(3.0)
        if plotIt:
            fig, ax = plt.subplots(nrows=2, sharex=True) 
            ax[0].plot(Af.f, np.abs(self.DAf))
            ax[0].plot(Af.f, np.abs(Af.data))
            
            ax[1].plot(Af.f, np.abs(self.DEf))
            ax[1].plot(Af.f, np.abs(Ef.data))
            plt.show()

        diff = np.abs(self.DAf - Af.values) ** 2 + np.abs(self.DEf - Ef.values) ** 2
        loglik = -float(np.sum(diff / self.SA) * self.dataX.df) /2

        # scalarproduct_signal = 4*np.real(np.sum((Af*np.conjugate(Af) / self.SA).values) * self.dataX.df)
        # scalarproduct_signal += 4*np.real(np.sum((Ef*np.conjugate(Ef) / self.SE).values) * self.dataX.df)
        # scalarproduct_data_signal = 4*np.real(np.sum((self.DAf*np.conjugate(Af) / self.SA).values) * self.dataX.df)
        # scalarproduct_data_signal += 4*np.real(np.sum((self.DEf*np.conjugate(Ef) / self.SE).values) * self.dataX.df) 
        # loglik = scalarproduct_data_signal - scalarproduct_signal/2   

        return loglik

    def loglikelihood3(self, pGBs):
        for i in range(len(pGBs)):
            Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=pGBs[i], oversample=4, simulator="synthlisa")

            a,Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)
            a,Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)
            a,Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)


        diff = np.abs(self.dataX - Xs_total) ** 2 + np.abs(self.dataY - Ys_total) ** 2 + np.abs(self.dataZ - Zs_total) ** 2
        # p1 = -float(np.sum(diff / Sn)*Xs.attrs['df'])/2.0
        p1 = float(np.sum(diff / self.Sn) * Xs_total.df) / 2.0
        # p1 = np.exp(p1)
        return -p1#/10000

    def loglikelihood2(self, pGBs):
        for i in range(len(pGBs)):
            Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=pGBs[i], oversample=4, simulator="synthlisa")

            _ , Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0, copy=False)
            _ , Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0, copy=False)
            _ , Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0, copy=False)


        diff = np.abs(self.dataX - Xs_total.values) ** 2 + np.abs(self.dataY - Ys_total.values) ** 2 + np.abs(self.dataZ - Zs_total.values) ** 2
        # p1 = -float(np.sum(diff / Sn)*Xs.attrs['df'])/2.0
        p1 = float(np.sum(diff / self.Sn) * Xs_total.df) / 2.0
        # p1 = np.exp(p1)
        return -p1#/10000

    def differential_evolution_search(self, frequency_boundaries):
        bounds = []
        for signal in range(number_of_signals):
            for i in range(8):
                bounds.append((0,1))

        maxpGB = []
        self.boundaries_reduced = deepcopy(self.boundaries)
        self.boundaries_reduced['Frequency'] = frequency_boundaries

        start = time.time()
        res = differential_evolution(self.function_evolution, bounds=bounds, disp=True, strategy='best1exp', popsize=15,tol= 1e-4, recombination=0.75, mutation=(0.5,1))
        print('time',time.time()-start)
        for signal in range(number_of_signals):
            maxpGB.append(scaletooriginal(res.x[signal*8:signal*8+8],self.boundaries_reduced))
        print(res)
        print(maxpGB)
        print(self.loglikelihood(maxpGB))
        # print(pGB)
        return [maxpGB]

    def search(self):
        # np.random.seed(42)

        parameters_recorded = [None] * 10
        for n in range(len(parameters_recorded)):
            start = time.time()
            parameters_recorded[n] = CoordinateMC(n, self.pGBs, self.boundaries, parameters_recorded, self.loglikelihood)

            print('n',n+1,'time', int(time.time()-start), np.round(parameters_recorded[n][0]['Loglikelihood'][-1],2),len(parameters_recorded[n][0]['Loglikelihood']),np.round(self.loglikelihood([self.pGB]),3),parameters_recorded[n][0]['Amplitude'][-1])
        # pbar = tqdm(total=len(parameters_recorded))
        # pool = mp.Pool(mp.cpu_count())
        # parameters_recorded = pool.map(CoordinateMC, [n for n in range(len(parameters_recorded))])
        # pool.close()
        # pool.join()
        # pbar.close()

        best_run = 0
        loglikelihoodofruns = np.zeros(len(parameters_recorded))
        for i in range(len(parameters_recorded)):
            loglikelihoodofruns[i] = parameters_recorded[i][0]['Loglikelihood'][-1]
        best_value = np.max(loglikelihoodofruns)
        best_run = np.argmax(loglikelihoodofruns)
        good_runs = loglikelihoodofruns > best_value*3
        indices = (-loglikelihoodofruns).argsort()[:len(good_runs)]
        pGBmodes = []
        for i in range(len(good_runs)):
            if good_runs[i]:
                pGBmodes.append({})
        indices = (-loglikelihoodofruns).argsort()[:len(pGBmodes)]
        pGBmodes = []
        for n in indices:
            pGBmodes.append([])
            for i in range(number_of_signals):
                pGBmodes[-1].append({})
                for parameter in parameters + ['Loglikelihood']:
                    if parameter == 'Loglikelihood':
                        pGBmodes[-1][i][parameter] = parameters_recorded[n][0][parameter][-1]
                    else:
                        pGBmodes[-1][i][parameter] = parameters_recorded[n][i][parameter][-1]
        if len(pGBmodes) > 5:
            for signal in range(number_of_signals):
                pGBmodes = pGBmodes[:5]
        return pGBmodes

    def optimize(self, pGBmodes):
        bounds = ()
        for signal in range(number_of_signals):
            bounds += ((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1))
        for i in range(len(pGBmodes)):
            maxpGB = []
            boundaries_reduced = []
            pGBs01 = []

            for j in range(5):
                x = []
                for signal in range(number_of_signals):
                    if j == 0:
                        maxpGB.append({})
                        boundaries_reduced.append({})
                        for parameter in parameters:
                            maxpGB[signal][parameter] = pGBmodes[i][signal][parameter]
                    # print(maxpGB)
                    boundaries_reduced[signal] = Reduce_boundaries(maxpGB[signal], self.boundaries,ratio=0.1)
                    # if j == 0:
                    #     boundaries_reduced[signal] = deepcopy(self.boundaries)
                    pGBs01.append({})
                    for parameter in parameters:
                        if parameter in ["EclipticLatitude"]:
                            pGBs01[signal][parameter] = (np.sin(maxpGB[signal][parameter]) - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
                        elif parameter in ["Inclination"]:
                            pGBs01[signal][parameter] = (np.cos(maxpGB[signal][parameter]) - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
                        elif parameter in ['Amplitude',"FrequencyDerivative"]:
                            pGBs01[signal][parameter] = (np.log10(maxpGB[signal][parameter]) - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
                        else:
                            pGBs01[signal][parameter] = (maxpGB[signal][parameter] - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
                    for parameter in parameters:
                        x.append(pGBs01[signal][parameter])
                # print(loglikelihood(maxpGB))
                res = scipy.optimize.minimize(self.function, x, args=boundaries_reduced, method='SLSQP', bounds=bounds, tol=1e-10)
                for signal in range(number_of_signals):
                    maxpGB[signal] = scaletooriginal(res.x[signal*8:signal*8+8],boundaries_reduced[signal])
                # print('optimized loglikelihood', loglikelihood(maxpGB),maxpGB)
                # print('boundaries reduced', boundaries_reduced)

            best_value = self.loglikelihood(maxpGB)
            if i == 0:
                current_best_value = best_value
                current_maxpGB = maxpGB
            try:
                if current_best_value < best_value:
                    current_best_value = best_value
                    current_maxpGB = maxpGB
            except:
                pass
            print('optimized loglikelihood', self.loglikelihood(maxpGB),self.loglikelihood([self.pGB]))
        maxpGB = current_maxpGB
        print('final optimized loglikelihood', self.loglikelihood(maxpGB),self.loglikelihood([self.pGB]),maxpGB[0]['Frequency'])
        return maxpGB, self.pGB

    def function(self, pGBs01, boundaries_reduced):
        pGBs = []
        for signal in range(number_of_signals):
            pGBs.append({})
            i = 0
            for parameter in parameters:
                if parameter in ["EclipticLatitude"]:
                    pGBs[signal][parameter] = np.arcsin((pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
                elif parameter in ["Inclination"]:
                    pGBs[signal][parameter] = np.arccos((pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
                elif parameter in ['Amplitude',"FrequencyDerivative"]:
                    pGBs[signal][parameter] = 10**((pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
                else:
                    pGBs[signal][parameter] = (pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0]
                i += 1
        p = -self.loglikelihood(pGBs)
        return p#/10**4

    def function_evolution(self, pGBs01):
        pGBs = []
        for signal in range(number_of_signals):
            pGBs.append({})
            i = 0
            for parameter in parameters:
                if parameter in ["EclipticLatitude"]:
                    pGBs[signal][parameter] = np.arcsin((pGBs01[signal*8:signal*8+8][i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
                elif parameter in ["Inclination"]:
                    pGBs[signal][parameter] = np.arccos((pGBs01[signal*8:signal*8+8][i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
                elif parameter in ['Amplitude',"FrequencyDerivative"]:
                    pGBs[signal][parameter] = 10**((pGBs01[signal*8:signal*8+8][i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
                else:
                    pGBs[signal][parameter] = (pGBs01[signal*8:signal*8+8][i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0]
                i += 1
        p = -self.loglikelihood(pGBs)
        return p#/10**4

    def fisher_information(self, maxpGB):
        maxpGB_changed = deepcopy(maxpGB)
        maxpGB01 = scaleto01(maxpGB, self.boundaries)
        maxpGB01_changed = deepcopy(maxpGB01)
        step_size = {}
        pGB_low = {}
        pGB_high = {}
        derivativeAf = {}
        derivativeEf = {}
        inner_product = {}
        for i in range(5):
            for parameter in parameters:
                if i == 0:
                    step_size[parameter] = 1e-2
                    # if parameter == 'Frequency':
                    #     step_size[parameter] = 0.00001
                else:
                    step_size[parameter] = 0.001/np.sqrt(inner_product[parameter][parameter])
                
                pGB_low = maxpGB01[parameter] - step_size[parameter]/2
                pGB_high = maxpGB01[parameter] + step_size[parameter]/2
                if pGB_low < 0:
                    pGB_low = 0
                if pGB_high > 1:
                    pGB_high = 1
                maxpGB01_changed[parameter] = pGB_low
                maxpGB_changed = scaletooriginalparameter(maxpGB01_changed,self.boundaries)
                # print(maxpGB_changed)
                Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=maxpGB_changed, oversample=4, simulator="synthlisa")
                index_low = np.searchsorted(Xs.f, self.dataX.f[0])
                Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
                Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
                Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
                Af_low = (Zs_total - Xs_total)/np.sqrt(2.0)
                Ef_low = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)

                maxpGB01_changed[parameter] = pGB_high
                maxpGB_changed = scaletooriginalparameter(maxpGB01_changed,self.boundaries)
                Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=maxpGB_changed, oversample=4, simulator="synthlisa")
                index_low = np.searchsorted(Xs.f, self.dataX.f[0])
                Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
                Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
                Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
                Af_high = (Zs_total - Xs_total)/np.sqrt(2.0)
                Ef_high = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)

                derivativeAf[parameter] = (Af_high - Af_low)/step_size[parameter]
                derivativeEf[parameter] = (Ef_high - Ef_low)/step_size[parameter]

                maxpGB01_changed[parameter] = maxpGB01[parameter]

            for parameter1 in parameters:
                inner_product[parameter1] = {}
                for parameter2 in parameters:
                    AE = derivativeAf[parameter1]*np.conjugate(derivativeAf[parameter2]) + derivativeAf[parameter1]*np.conjugate(derivativeAf[parameter2])
                    inner_product[parameter1][parameter2] = 4*float(np.real(np.sum(AE / self.SA) * self.dataX.df))
            print(step_size['Amplitude'],inner_product['Amplitude']['Amplitude'],step_size['Frequency'],inner_product['Frequency']['Frequency'])
        return inner_product

def objective(n,tdi_fs,Tobs):
    print('index',n)
    search = Search(n,tdi_fs,Tobs)
    # search.plot(search.pGBs)
    pGBmodes =  search.search()
    maxpGB, pGB =  search.optimize(pGBmodes)
    return maxpGB, pGB

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
pGBadded['Amplitude'] = 4*10**-21
pGBadded['EclipticLatitude'] = -0.2
pGBadded['EclipticLongitude'] = 1.4
pGBadded['Frequency'] = 0.00459687
pGBadded['FrequencyDerivative'] = 3*10**-15
pGBadded['Inclination'] = 0.5
pGBadded['InitialPhase'] = 3
pGBadded['Polarization'] = 2
pGBadded2 = {}
pGBadded2['Amplitude'] = 4*10**-24
pGBadded2['EclipticLatitude'] = 0.4
pGBadded2['EclipticLongitude'] = 0.4
pGBadded2['Frequency'] = 0.00459687
pGBadded2['FrequencyDerivative'] = 1*10**-17
pGBadded2['Inclination'] = 0.5
pGBadded2['InitialPhase'] = 1
pGBadded2['Polarization'] = 3
pGBadded3 = {}
pGBadded3['Amplitude'] = 4*10**-24
pGBadded3['EclipticLatitude'] = 0.6
pGBadded3['EclipticLongitude'] = 3
pGBadded3['Frequency'] = 0.00459687
pGBadded3['FrequencyDerivative'] = 1*10**-17
pGBadded3['Inclination'] = 0.5
pGBadded3['InitialPhase'] = 1
pGBadded3['Polarization'] = 3
pGBadded4 = {}
pGBadded4['Amplitude'] = 4*10**-24
pGBadded4['EclipticLatitude'] = 0.1
pGBadded4['EclipticLongitude'] = 2
pGBadded4['Frequency'] = 0.00459687
pGBadded4['FrequencyDerivative'] = 1*10**-17
pGBadded4['Inclination'] = 0.7
pGBadded4['InitialPhase'] = 1
pGBadded4['Polarization'] = 3
pGBadded5 = {}
pGBadded5['Amplitude'] = 6.37823e-23
pGBadded5['EclipticLatitude'] = -0.2
pGBadded5['EclipticLongitude'] = 1.4
pGBadded5['Frequency'] = 0.00622028
pGBadded5['FrequencyDerivative'] = 3*10**-15
pGBadded5['Inclination'] = 0.5
pGBadded5['InitialPhase'] = 3
pGBadded5['Polarization'] = 2
# pGBadded = {'Amplitude': 6.486747e-24, 'EclipticLatitude': -0.318064, 'EclipticLongitude': 4.369758, 'Frequency': 0.01620542, 'FrequencyDerivative': 8.508619e-14, 'Inclination': 0.866814, 'InitialPhase': 6.28136, 'Polarization': 5.978979}


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
# sangria_fn = DATAPATH + "/dgb-tdi.h5"
sangria_fn = DATAPATH + "/LDC1-1_MBHB_v2_TD_noiseless.hdf5"
# sangria_fn = DATAPATH + "/LDC1-4_GB_v2.hdf5"
# sangria_fn = DATAPATH + "/LDC1-3_VGB_v2_FD_noiseless.hdf5"
FD5 = LISAhdf5(sangria_fn)
Nsrc = FD5.getSourcesNum()
GWs = FD5.getSourcesName()
print("Found %d GW sources: " % Nsrc, GWs)
### TOD make sure MBHB is there
if GWs[0] != "MBHB-0":
    raise NotImplementedError
p = FD5.getSourceParameters(GWs[0])
td = FD5.getPreProcessTDI()
del_t = float(p.get("Cadence"))
reduction = 1
Tobs = float(int(p.get("ObservationDuration")/reduction))

dt = del_t

# Build timeseries and frequencyseries object for X,Y,Z
tdi_ts = xr.Dataset(dict([(k, TimeSeries(td[:int(len(td[:,1])/reduction), n], dt=dt)) for k, n in [["X", 1], ["Y", 2], ["Z", 3]]]))
# tdi_ts = xr.Dataset(dict([(k,TimeSeries(tdi_ts[k][:,1], dt=dt)) for k in ["X", "Y", "Z"]]))
tdi_fs = xr.Dataset(dict([(k, tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
GB = fastGB.FastGB(delta_t=dt, T=Tobs)  # in seconds

noise_model = "MRDv1"
Nmodel = get_noise_model(noise_model, np.logspace(-5, -1, 100))
Npsd = Nmodel.psd()

default_units = {'EclipticLatitude':'rad','EclipticLongitude':'rad',
         'PolarAngleOfSpin1':'rad','PolarAngleOfSpin2':'rad',
         'Spin1': '1','Spin2':'1',
         'Mass1':'Msun','Mass2':'Msun',
         'CoalescenceTime': 's','PhaseAtCoalescence':'rad',
         'InitialPolarAngleL':'rad','InitialAzimuthalAngleL':'rad',
         'Cadence': 's','Redshift': '1','Distance': 'Gpc',
         'ObservationDuration':'s'}


parameters = ['EclipticLatitude','EclipticLongitude',
         'PolarAngleOfSpin1','PolarAngleOfSpin2',
         'Spin1','Spin2',
         'Mass1','Mass2',
         'CoalescenceTime','PhaseAtCoalescence',
         'InitialPolarAngleL','InitialAzimuthalAngleL',
         'Distance']

pMBHB = {}
for parameter in parameters+['Redshift']:
    pMBHB[parameter] = p.get(parameter)
secondsperyear = 60*60*24*365.25
dt = 5 # waveform sampling
def funcMchirpofm1m2(m1, m2):
    return pow(m1*m2, 3./5) / pow(m1+m2, 1./5)
def funcm1m2ofMchirpq(Mc, q):
    m2 = (1+q)**(1/5)*Mc/q**(3/5)
    return m2*q, m2
Mc = funcMchirpofm1m2(pMBHB['Mass1'], pMBHB['Mass2'])
q = pMBHB['Mass1']/pMBHB['Mass2']
m1, m2 = funcm1m2ofMchirpq(Mc,q)

pGB = {}
ind = 0
found_sources = []
target_sources = []
first_start = time.time()
np.random.seed(42) #40
# for ind in range(1,len(p.get('Frequency'))):
number_of_signals = 1
signals_per_subtraction = 1



pMBHBs = deepcopy(pMBHB)
# pMBHBs['PolarAngleOfSpin1'] = 0
# pMBHBs = scaletooriginal(psample,boundaries)

t_max = pMBHB["CoalescenceTime"]+3000#*1.0001#60*60*24*365 # time of observation = 1yr
t_min = pMBHB["CoalescenceTime"]-3000#*0.9998
shift = t_min#int(pMBHB["CoalescenceTime"]*0.9998)
t_max -= shift
t_min -= shift
s_index = 0
# t_min = 0
coalescencetime = pMBHB['CoalescenceTime']
pMBHB['CoalescenceTime'] -= shift
pMBHBs['CoalescenceTime'] -= shift
initial_position = 2*np.pi*(((shift)/secondsperyear)%1)
config = {"initial_position": initial_position, "initial_rotation": 0, 
          "nominal_arm_length": 2500000000, "orbit_type": 'analytic'}
lisa_orbits = Orbits.type(config)
pMBHB["ObservationDuration"] = t_max
pMBHB["Cadence"] = dt
pMBHBs["ObservationDuration"] = t_max
pMBHBs["Cadence"] = dt


boundaries = {'ChirpMass': [Mc*(1.0 - 0.01), Mc*(1.0 + 0.01)], ### broadening just lead to longer run time
'MassRatio': [1.0, 3.0],
'CoalescenceTime': [pMBHB['CoalescenceTime'] - 50., pMBHB['CoalescenceTime'] + 50.], ### could be reduced to this range by computing the slide
'Spin1': [-0.99, 0.99],
'Spin2': [-0.99, 0.99],
'PolarAngleOfSpin1': [0, np.pi],
'PolarAngleOfSpin2': [0, np.pi],
# 'Distance': [np.log10(pMBHB['Distance'])*(0.9), np.log10(pMBHB['Distance'])*(1.1)],
'sinEclipticLatitude': [-1.0, 1.0],
'EclipticLongitude': [0.0, 2.0*np.pi],
'InitialAzimuthalAngleL': [0.0, 2.0*np.pi],
'InitialPolarAngleL': [0.0, np.pi],
'PhaseAtCoalescence': [0.0, 1.0*np.pi]}
parametersboundary = ['ChirpMass','MassRatio','CoalescenceTime','Spin1','Spin2','PolarAngleOfSpin1','PolarAngleOfSpin2','InitialPolarAngleL',
'sinEclipticLatitude','EclipticLongitude','InitialAzimuthalAngleL','PhaseAtCoalescence']
psample = np.random.rand(13)
samples1 = {}
for parameter in parametersboundary:
    samples1[parameter] = psample[parametersboundary.index(parameter)]
pMBHBs = changeparameterstoinput(samples1, boundaries)

# redshift = Cosmology.zofDl(pMBHB['Distance'],w=0, tolerance=0.00001)
Distance = Cosmology.DL(pMBHB['Redshift'],w=0)[0]
pMBHB['Distance'] *= 1000
pMBHBs['Distance'] = pMBHB['Distance']
trangeshift = np.arange(t_min, t_max, dt)
start = time.time()
tdi_X = semi_fast_tdi(config, pMBHB, t_min, t_max, dt)
print('pMBHB time',time.time()- start)
start = time.time()
Xs = semi_fast_tdi(config, pMBHBs, t_min, t_max, dt)
print('pMBHBs time',time.time()- start)

index_low = np.searchsorted(tdi_ts.t,  t_min+shift)

secondsperday = 24*3600
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
plt.figure(figsize=fig_size)
plt.plot(tdi_ts.t[index_low:index_low+len(tdi_X)]/secondsperday, tdi_ts['X'][index_low:index_low+len(tdi_X)], label="Data")
plt.plot((trangeshift+shift)/secondsperday, tdi_X, label="MBHB", alpha=1, color= colors[2])
plt.plot((trangeshift+shift)/secondsperday, Xs, label="sample", alpha=0.5)
# plt.plot(trangeshift-t_min, Xs, label="strain to TDI", alpha=0.5)
plt.xlabel("Time (days)")
# plt.axis([coalescencetime-1000, coalescencetime+600, None, None])
plt.locator_params(nbins=4)
plt.ylabel("TDI X strain")
plt.tight_layout()
plt.legend()
tdi_Xfd = tdi_X.ts.fft(win=window)
Xsfd = Xs.ts.fft(win=window)
tdi_fdX = tdi_ts['X'][index_low:index_low+len(tdi_X)].ts.fft(win=window)
plt.figure()
plt.semilogx(tdi_fdX.f, tdi_fdX, label="Data")
plt.semilogx(tdi_Xfd.f, tdi_Xfd, label="MBHB", alpha=1, color=colors[2])
plt.semilogx(Xsfd.f, Xsfd, label="sample", alpha=0.5)
plt.xlabel("f (Hz)")
plt.ylabel("TDI X real (1/Hz)")
plt.legend()
plt.tight_layout()
plt.show()


# maxpGB = [[{'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}]]
# # maxpGB = [[{'Amplitude': 4.083357969533303e-22, 'EclipticLatitude': 0.8719966900463507, 'EclipticLongitude': 0.4861128797587986, 'Frequency': 0.00399522108238035, 'FrequencyDerivative': 1.0720262111754569e-16, 'Inclination': 1.0241926728648307, 'InitialPhase': 2.319788782634353, 'Polarization': 2.6588421907673028}]]
# for k in range(1):
#     if k == 1:
#         maxpGB = [[{'Amplitude': 1.1706831455114382e-22, 'EclipticLatitude': -1.182657374135248, 'EclipticLongitude': -2.671010079711571, 'Frequency': 0.0039946199549690566, 'FrequencyDerivative': 9.547621993738103e-17, 'Inclination': 1.9399086433607453, 'InitialPhase': 5.612220707908651, 'Polarization': 0.9418521680342067}]]
#     for j in range(signals_per_subtraction):
#         for i in range(1):
#             found_sources.append(maxpGB[j][i])
#             Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=maxpGB[j][i], oversample=4, simulator="synthlisa")
#             source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
#             index_low = np.searchsorted(tdi_fs["X"].f, Xs_subtracted.f[0])
#             index_high = index_low+len(Xs_subtracted)
#             for k in ["X", "Y", "Z"]:
#                 tdi_fs[k].data[index_low:index_high] = tdi_fs[k].data[index_low:index_high] - source_subtracted[k].data
#             tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft()) for k in ["X", "Y", "Z"]]))
#             Xs_subtracted, Ys_subtracted, Zs_subtracted = GB_long.get_fd_tdixyz(template=maxpGB[j][i], oversample=4, simulator="synthlisa")
#             source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
#             index_low = np.searchsorted(tdi_fs["X"].f, Xs_subtracted.f[0])
#             index_high = index_low+len(Xs_subtracted)
#             for k in ["X", "Y", "Z"]:
#                 tdi_fs_long[k].data[index_low:index_high] = tdi_fs[k].data[index_low:index_high] - source_subtracted[k].data
#             tdi_ts_long = xr.Dataset(dict([(k, tdi_fs_long[k].ts.ifft()) for k in ["X", "Y", "Z"]]))

for ind in range(1): #[3,8,9]
    signal_peak = -1
    f, psdX =  scipy.signal.welch(tdi_ts["X"], fs=1.0/dt, window='hanning', nperseg=len(tdi_ts["X"])/1)
    f, psdY =  scipy.signal.welch(tdi_ts["Y"], fs=1.0/dt, window='hanning', nperseg=len(tdi_ts["X"])/1)
    f, psdZ =  scipy.signal.welch(tdi_ts["Z"], fs=1.0/dt, window='hanning', nperseg=len(tdi_ts["X"])/1)
    psd = psdX + psdY + psdZ
    # indexes = np.logical_and(f[peaks]>10**-4, f[peaks]<2*10**-2)
    indexes = np.logical_and(f>lower_frequency-padding, f<upper_frequency+padding)
    psd = psd[indexes]
    f = f[indexes]
    peaks, properties = scipy.signal.find_peaks(np.log(np.sqrt(psd)/np.mean(np.sqrt(psd))),width=1,  prominence=(0.22,None))
    peak_index = 2
    # plt.figure()
    # plt.semilogy(f,psd)
    # # plt.scatter(f[peaks],psd[peaks])
    # # plt.scatter(f[peaks[peak_index]],psd[peaks[peak_index]])
    # # plt.scatter(f[properties['left_bases'][peak_index]],psd[properties['left_bases'][peak_index]])
    # # plt.scatter(f[properties['right_bases'][peak_index]],psd[properties['right_bases'][peak_index]])
    # for k in range(len(pGB_injected)):
    #     plt.axvline(pGB_injected[k]['Frequency'])
    # plt.show()

    indexes = np.logical_and(f[peaks]>lower_frequency-padding, f[peaks]<upper_frequency+padding)
    peaks = peaks[indexes]
    # if len(peaks) < signals_per_subtraction:
    #     break
    prominences = properties['prominences'][indexes] 
    indexes_peaks = np.argsort(prominences)


    range_index = np.logical_and(tdi_fs.f > lower_frequency-padding, tdi_fs.f < upper_frequency+padding)


    # pGB = {'Amplitude': 2.360780938411229e-22, 'EclipticLatitude': 0.8726817299401121, 'EclipticLongitude': 0.4857930058990844, 'Frequency': 0.0039952210295307105, 'FrequencyDerivative': 1.0762922030519564e-16, 'Inclination': 0.008495799066858739, 'InitialPhase': 0.013026742125237894, 'Polarization': 1.507449983207914}
    search1 = Search(tdi_fs,Tobs)
    search1.plot()
    # F_stat = search1.F(pGB['Frequency'],pGB['EclipticLatitude'],pGB['EclipticLongitude'],pGB)
    # maxpGBsearch =  search1.differential_evolution_search()


    start = time.time()
    # F_stat = []
    # frequency = []
    # eclipticlatitude = []
    # eclipticlongitude = []
    # pGB['Amplitude'] = 1e-22
    # pGB['FrequencyDeviation'] = 0
    # N = 10
    # Nsky = 10
    # frequency_boundaries = [lower_frequency,upper_frequency]
    # for n in range(Nlat):
    #     eclipticlatitude.append(search1.boundaries['EclipticLatitude'][0]+(search1.boundaries['EclipticLatitude'][1]-search1.boundaries['EclipticLatitude'][0])*n/Nlat)
    #     eclipticlongitude.append(search1.boundaries['EclipticLongitude'][0]+(search1.boundaries['EclipticLongitude'][1]-search1.boundaries['EclipticLongitude'][0])*n/Nlat)
    # for k in range(N):
    #     F_stat.append([])
    #     frequency.append(frequency_boundaries[0] + (frequency_boundaries[1]-frequency_boundaries[0])*k/N)
    #     for l in range(Nlat):
    #         F_stat[-1].append([])
    #         for m in range(Nlat):
    #             F_stat[-1][-1].append(search1.F(frequency[-1],eclipticlatitude[l],eclipticlongitude[m],pGB))
    # F_stat = np.asarray(F_stat)

    # F-statistic
    N_frequency = 10
    N_sky = 10
    F_stat, frequency, eclipticlatitude, eclipticlongitude = search1.f_statistic(lower_frequency-padding, upper_frequency+padding, N_frequency, N_sky)
    ind = np.unravel_index(np.argmax(F_stat, axis=None), F_stat.shape)
    if ind[0]>0:
        lower_index = ind[0]-1
    else:
        lower_index = ind[0]
    if ind[0]<N_frequency-1:
        upper_index = ind[0]+1
    else:
        upper_index = ind[0]
    search1.reduced_frequency_boundaries = [frequency[lower_index],frequency[upper_index]]
    max = np.max(F_stat, axis=(1,2))
    print('time',time.time()-start)
    Eclipticlatitude,Eclipticlongitude = np.meshgrid(np.arccos(eclipticlatitude),eclipticlongitude)

    # plt.figure()
    # plt.plot(np.asarray(frequency),max)
    # for k in range(len(pGB_injected)):
    #     plt.axvline(pGB_injected[k]['Frequency'])
    # plt.xlabel('f [Hz]')
    # plt.ylabel('F statistic')
    # plt.show()
    # plt.figure()
    # plt.scatter(Eclipticlongitude,Eclipticlatitude,c=F_stat[0])
    # plt.colorbar()
    # plt.show()

    # start = time.time()
    # N_frequency = 10
    # F_stat, frequency, eclipticlatitude, eclipticlongitude = search1.f_statistic(search1.reduced_frequency_boundaries[0], search1.reduced_frequency_boundaries[1], N_frequency, N_sky)
    # max = np.max(F_stat, axis=(1,2))
    # print('time',time.time()-start)
    # Eclipticlatitude,Eclipticlongitude = np.meshgrid(np.arccos(eclipticlatitude),eclipticlongitude)
    # plt.figure()
    # plt.plot(np.asarray(frequency)*1e3,max)
    # plt.xlabel('f [mHz]')
    # plt.ylabel('F statistic')
    # plt.show()
    # ind = np.unravel_index(np.argmax(F_stat, axis=None), F_stat.shape)
    # if ind[0]>1:
    #     lower_index = ind[0]-2
    # else:
    #     lower_index = ind[0]
    # if ind[0]<N_frequency-2:
    #     upper_index = ind[0]+2
    # else:
    #     upper_index = ind[0]
    # search1.reduced_frequency_boundaries = [frequency[lower_index],frequency[upper_index]]
    # plt.figure()
    # plt.scatter(Eclipticlongitude,Eclipticlatitude,c=F_stat[ind[0],:,:])
    # plt.xlabel('Eclipticlongitude')
    # plt.ylabel('Eclipticlatitude')
    # plt.colorbar()
    # plt.show()
    # search1.reduced_frequency_boundaries = [lower_frequency,upper_frequency]

    search1.plot()
    print(search1.reduced_frequency_boundaries)
    maxpGBsearch =  search1.differential_evolution_search(search1.reduced_frequency_boundaries)
    # pGBmodes =  search1.search()
    # maxpGBsearch, pGB =  search1.optimize(pGBmodes)

    # plt.figure(figsize=fig_size)
    # ax1 = plt.subplot(111)
    # ax1.semilogy(tdi_fs.f[range_index],np.abs(tdi_fs['X'][range_index])**2,'k',zorder= 2)
    # for i in range(len(p.get('Amplitude')[indexes][index_low:index_high])):
    #     pGB = {}
    #     for parameter in parameters:
    #         pGB[parameter] = p.get(parameter)[indexes][index_low:index_high][i]
    #     Xs, Ys, Zs = GB.get_fd_tdixyz(template= pGB, oversample=4, simulator="synthlisa")
    #     ax1.semilogy(Xs.f,np.abs(Xs)**2)
    # plt.show()

    # plt.figure(figsize=fig_size)
    # ax1 = plt.subplot(111)
    # ax1.loglog(f*1000, np.sqrt(psd), label="TDI X")
    # ax1.loglog(f[peaks]*1000, np.sqrt(psd[peaks]),'.', label="peaks")
    # ax1.loglog(f[peaks[indexes_peaks[-1]]]*1000, np.sqrt(psd[peaks[indexes_peaks[-1]]]),'.', label="most prominent peak")
    # # ax1.plot(Xs_added2.f * 1000, Xs_added2.values.real, label="VGB2", marker=".", zorder=5)
    # # ax1.axvline(boundaries['Frequency'][0]* 1000, color= 'red')
    # # ax1.axvline(boundaries['Frequency'][1]* 1000, color= 'red')
    # ax1.set_xlim(0.01,100)
    # plt.legend()
    # plt.show()

    # search1 = Search(0,tdi_fs,Tobs)
    # search1.plot(search1.pGBs)
    # pGBmodes =  search1.search()
    # maxpGB, pGB =  search1.optimize(pGBmodes)
    # maxpGB, pGB = objective(0,tdi_fs,Tobs)

    # reduction = 1
    # Tobs = float(int(p.get("ObservationDuration")/reduction))
    # # Build timeseries and frequencyseries object for X,Y,Z
    # # tdi_ts = xr.Dataset(dict([(k, TimeSeries(td[:int(len(td[:,1])/reduction), n], dt=dt)) for k, n in [["X", 1], ["Y", 2], ["Z", 3]]]))
    # # tdi_fs = xr.Dataset(dict([(k, tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
    # tdi_ts = xr.Dataset(dict([["X", TimeSeries(td[:int(len(td[:,1])/reduction), 1], dt=dt)],["Y", TimeSeries(td[:int(len(td[:,1])/reduction), 2], dt=dt)],
    # ["Z", TimeSeries(td[:int(len(td[:,1])/reduction), 3], dt=dt)]]))
    # tdi_fs = xr.Dataset(dict([["X", tdi_ts['X'].ts.fft(win=window)],["Y", tdi_ts['Y'].ts.fft(win=window)],["Z", tdi_ts['Z'].ts.fft(win=window)]]))
    # print(len(tdi_fs['X']),len(tdi_ts['X']),len(td[:int(len(td[:,1])/reduction),1]),'timeseries',TimeSeries(td[:int(len(td[:,1])/reduction), 1], dt=dt),'xarray',xr.Dataset(dict([(k, TimeSeries(td[:int(len(td[:,1])/reduction), n], dt=dt)) for k, n in [["X", 1]]])),
    # 'next',dict([["X", TimeSeries(td[:int(len(td[:,1])/reduction), 1], dt=dt)]]),TimeSeries(td[:int(len(td[:,1])/reduction), 1], dt=dt))
    # GB = fastGB.FastGB(delta_t=dt, T=Tobs)  # in seconds
    # search0 = Search(0,tdi_fs,GB)
    # search0.plot(maxpGB)
    # maxpGB2, pGB = search0.optimize([[maxpGB2]])
    
    # pool = mp.Pool(mp.cpu_count())
    # start_parallel = time.time()
    # maxpGBsearch = pool.map(partial(objective, tdi_fs=tdi_fs, Tobs=Tobs), [n for n in range(signals_per_subtraction)] )
    # pool.close()
    # pool.join()
    # print('time parallel', time.time()-start_parallel)

    # maxpGB = []
    # for j in range(signals_per_subtraction):
    #     maxpGB.append(maxpGBsearch[j])
    #     print(maxpGB[-1])
    # for j in range(signals_per_subtraction):
    #     for i in range(number_of_signals):
    #         found_sources.append(maxpGB[j][i])
    #         # print(found_sources[-1])
    #         # target_sources.append(pGB[j])
    #         Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=maxpGB[j][i], oversample=4, simulator="synthlisa")
    #         source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
    #         index_low = np.searchsorted(tdi_fs["X"].f, Xs_subtracted.f[0])
    #         index_high = index_low+len(Xs_subtracted)
    #         for k in ["X", "Y", "Z"]:
    #             tdi_fs[k].data[index_low:index_high] = tdi_fs[k].data[index_low:index_high] - source_subtracted[k].data
    #         tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft()) for k in ["X", "Y", "Z"]]))

    #         Xs_subtracted, Ys_subtracted, Zs_subtracted = GB_long.get_fd_tdixyz(template=maxpGB[j][i], oversample=4, simulator="synthlisa")
    #         source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
    #         index_low = np.searchsorted(tdi_fs["X"].f, Xs_subtracted.f[0])
    #         index_high = index_low+len(Xs_subtracted)
    #         for k in ["X", "Y", "Z"]:
    #             tdi_fs_long[k].data[index_low:index_high] = tdi_fs[k].data[index_low:index_high] - source_subtracted[k].data
    #         tdi_ts_long = xr.Dataset(dict([(k, tdi_fs_long[k].ts.ifft()) for k in ["X", "Y", "Z"]]))

# # found_sourcesrun = [{'Amplitude': 1.5364465535602838e-21, 'EclipticLatitude': 0.22842923568790388, 'EclipticLongitude': 3.9876628102916634, 'Frequency': 0.0034068399355733194, 'FrequencyDerivative': 1.7265154770131203e-16, 'Inclination': 1.140895523332838, 'InitialPhase': 1.3231144225092324, 'Polarization': 2.7082821014391674}, {'Amplitude': 1.8299785587003754e-21, 'EclipticLatitude': -0.5717437471770699, 'EclipticLongitude': 5.013423186770661, 'Frequency': 0.0021468057546732903, 'FrequencyDerivative': 2.5679108781904502e-17, 'Inclination': 2.381977466585251, 'InitialPhase': 1.7795097897215049, 'Polarization': 1.361815899661127}, {'Amplitude': 1.3740545095358953e-21, 'EclipticLatitude': 0.4586355029484734, 'EclipticLongitude': 2.326341135168954, 'Frequency': 0.002171903889305508, 'FrequencyDerivative': 2.656711011190014e-17, 'Inclination': 2.6231811067672464, 'InitialPhase': 2.246734702442237, 'Polarization': 0.4505374916175284}, {'Amplitude': 4.85984439130295e-22, 'EclipticLatitude': -0.6103079146464152, 'EclipticLongitude': 3.8165981027013838, 'Frequency': 0.007219358007407393, 'FrequencyDerivative': 3.215874944181044e-15, 'Inclination': 0.3611841295412159, 'InitialPhase': 2.8258044166152523, 'Polarization': 2.426471942484511}, {'Amplitude': 6.653710374755356e-22, 'EclipticLatitude': -0.359628257716563, 'EclipticLongitude': 4.9343500673177365, 'Frequency': 0.009211757029070448, 'FrequencyDerivative': 4.897974394062554e-15, 'Inclination': 0.8418300854577668, 'InitialPhase': 2.120331638811925, 'Polarization': 1.2267534109224667}, {'Amplitude': 3.312821037152804e-22, 'EclipticLatitude': 0.7326377959505177, 'EclipticLongitude': 6.056532678360872, 'Frequency': 0.004022512317404639, 'FrequencyDerivative': 7.928660484261939e-17, 'Inclination': 2.687294151927051, 'InitialPhase': 1.6080815997044122, 'Polarization': 1.9588214370682089}, {'Amplitude': 4.686184942845765e-22, 'EclipticLatitude': 0.20352748572849222, 'EclipticLongitude': 5.007749923410212, 'Frequency': 0.009211757978548653, 'FrequencyDerivative': 4.895550870611016e-15, 'Inclination': 1.1415963061394274, 'InitialPhase': 3.14159265358979, 'Polarization': 2.2644816658566564}, {'Amplitude': 5.18155313331795e-22, 'EclipticLatitude': 0.19478251263514385, 'EclipticLongitude': 5.219095549872743, 'Frequency': 0.009445485983463868, 'FrequencyDerivative': 7.023800557316083e-15, 'Inclination': 0.5777619297133901, 'InitialPhase': 3.141592653589793, 'Polarization': 2.0928411308295223}, {'Amplitude': 5.876055535475725e-22, 'EclipticLatitude': -0.46618198434204294, 'EclipticLongitude': 4.289702767312378, 'Frequency': 0.003077973937546665, 'FrequencyDerivative': 1.025303603104508e-16, 'Inclination': 0.7494287470138635, 'InitialPhase': 0.056843419329868985, 'Polarization': 2.8300541786136124}, {'Amplitude': 4.0839934840674097e-22, 'EclipticLatitude': 0.8720661480722399, 'EclipticLongitude': 0.4861465045465657, 'Frequency': 0.003995221091227311, 'FrequencyDerivative': 1.0688851160394643e-16, 'Inclination': 1.0242692798927258, 'InitialPhase': 2.319609919910819, 'Polarization': 2.6582674105848603}]
# # target_sourcesrun = [{'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}, {'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}, {'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}, {'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}, {'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}, {'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}, {'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}, {'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}, {'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}, {'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}]
# # found_sources19 = [{'Amplitude': 7.158007512978195e-23, 'EclipticLatitude': -0.17879074789424193, 'EclipticLongitude': 4.603147784376865, 'Frequency': 0.00399097103940054, 'FrequencyDerivative': 2.29199470942113e-16, 'Inclination': 2.3395967762010357, 'InitialPhase': 1.6414925267236022, 'Polarization': 1.409625714095948}, {'Amplitude': 5.417316688438823e-23, 'EclipticLatitude': 0.09299982963031442, 'EclipticLongitude': 4.709278461472791, 'Frequency': 0.003991793230604934, 'FrequencyDerivative': 1.85302117460904e-16, 'Inclination': 1.8243956188626644, 'InitialPhase': 0.3119979003693603, 'Polarization': 1.9954273407108052}, {'Amplitude': 2.6586721928556783e-23, 'EclipticLatitude': -0.23347086835824327, 'EclipticLongitude': 4.667459964579811, 'Frequency': 0.003991718511582032, 'FrequencyDerivative': 1e-20, 'Inclination': 2.052394297268333, 'InitialPhase': 2.713778167855026, 'Polarization': 2.6371523097400362}, {'Amplitude': 1.4603651612868757e-23, 'EclipticLatitude': 0.5088538078457533, 'EclipticLongitude': 5.061705367435939, 'Frequency': 0.003991299832643548, 'FrequencyDerivative': 1.422149445394665e-16, 'Inclination': 3.1293675414452404, 'InitialPhase': 0.9692647620129494, 'Polarization': 1.6681641616303722}, {'Amplitude': 3.059691174442375e-23, 'EclipticLatitude': 0.4317147381501876, 'EclipticLongitude': 4.997451644492614, 'Frequency': 0.003989147910089541, 'FrequencyDerivative': 1.8726122581411511e-16, 'Inclination': 1.1753859885192686, 'InitialPhase': 1.6649392950095585e-10, 'Polarization': 1.7744204264185908}, {'Amplitude': 1.0744473573747681e-23, 'EclipticLatitude': -0.04537250792371965, 'EclipticLongitude': 4.75590442505031, 'Frequency': 0.003991905382027591, 'FrequencyDerivative': 1.1687840960558461e-16, 'Inclination': 3.1085189766581087, 'InitialPhase': 2.4798219832899684, 'Polarization': 0.23152794448077565}, {'Amplitude': 6.3758855191988285e-24, 'EclipticLatitude': -0.5711395779756441, 'EclipticLongitude': 4.719785340796603, 'Frequency': 0.003991184090209936, 'FrequencyDerivative': 5.197211730593386e-16, 'Inclination': 0.03112334053551802, 'InitialPhase': 0.1840427721876706, 'Polarization': 2.665465621636476}, {'Amplitude': 2.695963984061525e-23, 'EclipticLatitude': 0.698020789183839, 'EclipticLongitude': 4.544599641191273, 'Frequency': 0.00398921600077084, 'FrequencyDerivative': 2.748523299944693e-17, 'Inclination': 1.566152596953124, 'InitialPhase': 1.6976257328944115, 'Polarization': 0.17969203457850272}, {'Amplitude': 1.7140045975600926e-23, 'EclipticLatitude': 0.5353309026749867, 'EclipticLongitude': 4.661338785721523, 'Frequency': 0.003990068561124402, 'FrequencyDerivative': 1.448478884234823e-20, 'Inclination': 1.0325571150373603, 'InitialPhase': 2.8574124537257295, 'Polarization': 1.2248048144455297}, {'Amplitude': 8.600732555086791e-24, 'EclipticLatitude': -0.464957284959853, 'EclipticLongitude': 4.537557676317077, 'Frequency': 0.0039912439454365765, 'FrequencyDerivative': 1.811792639171847e-16, 'Inclination': 1.009201467440333, 'InitialPhase': 1.7030607900434433e-19, 'Polarization': 3.141592653589793}]
# maxpGBsearch = [[{'Amplitude': 5.857579179174051e-23, 'EclipticLatitude': -0.08784800634433576, 'EclipticLongitude': 2.102933791009969, 'Frequency': 0.006220280040884535, 'FrequencyDerivative': 7.503323987637861e-16, 'Inclination': 0.5104961907537644, 'InitialPhase': 3.8035280736734642, 'Polarization': 0.07816858489674128}]]
# maxpGBsearch = [[{'Amplitude': 2.1808846772946646e-22, 'EclipticLatitude': -0.5283847702958109, 'EclipticLongitude': -2.5389708751966173, 'Frequency': 0.0012531301631903925, 'FrequencyDerivative': 3.6807471100533655e-19, 'Inclination': 0.9737933771974954, 'InitialPhase': 3.95558129985279, 'Polarization': 2.8786856073448313}]]
# maxpGBsearch = [[{'Amplitude': 1.5454899498120303e-22, 'EclipticLatitude': 0.3239828534237335, 'EclipticLongitude': -2.7552505031608576, 'Frequency': 0.0013596198681856582, 'FrequencyDerivative': 6.7913371080273325e-18, 'Inclination': 0.9607561157105087, 'InitialPhase': 0.8835396431378469, 'Polarization': 2.498067711496641}]]
# maxpGBsearch = [[{'Amplitude': 5.79066183343146e-23, 'EclipticLatitude': -0.0878417628451765, 'EclipticLongitude': 2.102946967840941, 'Frequency': 0.006220280042028942, 'FrequencyDerivative': 7.503066211114356e-16, 'Inclination': 0.4881731730267816, 'InitialPhase': 0.7708964347916855, 'Polarization': 1.703361386089743}]]
tdi_fs_long_subtracted = deepcopy(tdi_fs_long)
tdi_fs_subtracted = deepcopy(tdi_fs)
reduction = 1
Tobs = float(int(p.get("ObservationDuration")/reduction))
# Build timeseries and frequencyseries object for X,Y,Z
tdi_ts = {}
tdi_fs = {}
n = 0
for k in ["X", "Y", "Z"]:
    n += 1
    tdi_ts[k] = TimeSeries(td[:int(len(td[:,1])/reduction),n], dt=dt)
    tdi_fs[k] = tdi_ts[k].ts.fft(win=window)
tdi_ts = xr.Dataset(tdi_ts)
tdi_fs = xr.Dataset(tdi_fs)
# tdi_ts = xr.Dataset(dict([(k,TimeSeries(tdi_ts[k][:,1], dt=dt)) for k in ["X", "Y", "Z"]]))
GB = fastGB.FastGB(delta_t=dt, T=Tobs)  # in seconds
# maxpGB = [[{'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}]]
# maxpGB = [[{'Amplitude': 4.083357969533303e-22, 'EclipticLatitude': 0.8719966900463507, 'EclipticLongitude': 0.4861128797587986, 'Frequency': 0.00399522108238035, 'FrequencyDerivative': 1.0720262111754569e-16, 'Inclination': 1.0241926728648307, 'InitialPhase': 2.319788782634353, 'Polarization': 2.6588421907673028}]]
# # # maxpGB = [[found_sources[0]]]
# for k in range(1):
#     # if k == 1:
#     #     maxpGB = [[{'Amplitude': 1.1706831455114382e-22, 'EclipticLatitude': -1.182657374135248, 'EclipticLongitude': -2.671010079711571, 'Frequency': 0.0039946199549690566, 'FrequencyDerivative': 9.547621993738103e-17, 'Inclination': 1.9399086433607453, 'InitialPhase': 5.612220707908651, 'Polarization': 0.9418521680342067}]]
#     for j in range(signals_per_subtraction):
#         for i in range(number_of_signals):
#             Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=maxpGB[j][i], oversample=4, simulator="synthlisa")
#             source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
#             index_low = np.searchsorted(tdi_fs["X"].f, Xs_subtracted.f[0])
#             index_high = index_low+len(Xs_subtracted)
#             for k in ["X", "Y", "Z"]:
#                 tdi_fs[k].data[index_low:index_high] = tdi_fs[k].data[index_low:index_high] - source_subtracted[k].data
#             tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft()) for k in ["X", "Y", "Z"]]))
            # Xs_subtracted, Ys_subtracted, Zs_subtracted = GB_long.get_fd_tdixyz(template=maxpGB[j][i], oversample=4, simulator="synthlisa")
            # source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
            # index_low = np.searchsorted(tdi_fs["X"].f, Xs_subtracted.f[0])
            # index_high = index_low+len(Xs_subtracted)
            # for k in ["X", "Y", "Z"]:
            #     tdi_fs_long[k].data[index_low:index_high] = tdi_fs[k].data[index_low:index_high] - source_subtracted[k].data
            # tdi_ts_long = xr.Dataset(dict([(k, tdi_fs_long[k].ts.ifft()) for k in ["X", "Y", "Z"]]))

indexes = np.argsort(p.get('Frequency'))
index_low = np.searchsorted(p.get('Frequency')[indexes], lower_frequency)
index_high = np.searchsorted(p.get('Frequency')[indexes], upper_frequency)
range_index = np.logical_and(tdi_fs.f > lower_frequency-padding, tdi_fs.f < upper_frequency+padding)


# plt.figure(figsize=fig_size)
# ax1 = plt.subplot(111)
# for i in range(len(found_sources)):
#     Xs, Ys, Zs = GB.get_fd_tdixyz(template= found_sources[i], oversample=4, simulator="synthlisa")
#     ax1.plot(found_sources[i]['Frequency'],np.log10(found_sources[i]['Amplitude']), marker = 'o', markersize = 3)
# for i in range(len(p.get('Amplitude')[indexes][index_low:index_high])):
#     pGB = {}
#     for parameter in parameters:
#         pGB[parameter] = p.get(parameter)[indexes][index_low:index_high][i]
#     ax1.plot(pGB['Frequency'],np.log10(pGB['Amplitude']),color='grey', marker = 'o', markersize = 5)
# plt.show()

# maxpGB = {'Amplitude': 2.3606573655101142e-22, 'EclipticLatitude': 0.8726571396465939, 'EclipticLongitude': 0.4857771559989381, 'Frequency': 0.003995221016631018, 'FrequencyDerivative': 1.0802976360457847e-16, 'Inclination': 0.0006744803859335543, 'InitialPhase': 4.122296471788966, 'Polarization': 0.42105304476782063}
# maxpGB = {'Amplitude': 2.360780938411229e-22, 'EclipticLatitude': 0.8726817299401121, 'EclipticLongitude': 0.4857930058990844, 'Frequency': 0.0039952210295307105, 'FrequencyDerivative': 1.0762922030519564e-16, 'Inclination': 0.008495799066858739, 'InitialPhase': 0.013026742125237894, 'Polarization': 1.507449983207914}
# maxpGB = found_sources[1]
# maxpGB = {'Amplitude': 1.1704311160697462e-22, 'EclipticLatitude': -1.1826967434100408, 'EclipticLongitude': -2.6708875378846444, 'Frequency': 0.003994619931753887, 'FrequencyDerivative': 9.619861969590325e-17, 'Inclination': 1.9399691543655102, 'InitialPhase': 2.4682379038762075, 'Polarization': 2.512857354623118}
#second highest found
# maxpGB = {'Amplitude': 1.1684041146866079e-22, 'EclipticLatitude': -1.181910906694611, 'EclipticLongitude': -2.6694499096862354, 'Frequency': 0.003994618496148221, 'FrequencyDerivative': 1.3469993511889025e-16, 'Inclination': 1.9408782301800371, 'InitialPhase': 2.295524550856767, 'Polarization': 2.5328368226079614}
# maxpGB = {'Amplitude': 1.2117026142890994e-22, 'EclipticLatitude': -1.2206202633308574, 'EclipticLongitude': -2.6796337238866297, 'Frequency': 0.003994617791176147, 'FrequencyDerivative': 1.2782188418443533e-15, 'Inclination': 1.895745380989841, 'InitialPhase': 2.8823071421535413, 'Polarization': 2.4804104605148782}
# maxpGB = {'Amplitude': 1.165244087502911e-22, 'EclipticLatitude': -1.1808445115412942, 'EclipticLongitude': -2.6631741296680835, 'Frequency': 0.003994616096440512, 'FrequencyDerivative': 1.9820859862264132e-16, 'Inclination': 1.940250722025412, 'InitialPhase': 2.099175010008629, 'Polarization': 2.5093628935404544}
# maxpGB = {'Amplitude': 1.248193e-22, 'EclipticLatitude': -1.185356, 'EclipticLongitude': 3.593803, 'Frequency': 0.003994621, 'FrequencyDerivative': 6.709408e-17, 'Inclination': 1.906596, 'InitialPhase': 5.663538, 'Polarization': 0.96414}
# maxpGB = {'Amplitude': 2.5050879071884518e-23, 'EclipticLatitude': 0.26226502465165796, 'EclipticLongitude': 1.55571265787166, 'Frequency': 0.003995000207460324, 'FrequencyDerivative': 8.501729703518769e-17, 'Inclination': 1.960387154895422, 'InitialPhase': 5.722461960271562, 'Polarization': 2.015593509835858}
# maxpGB = target_sources[0]
# if maxpGB['EclipticLongitude'] > np.pi:
#     maxpGB['EclipticLongitude'] -= 2*np.pi
maxpGB = maxpGBsearch[0][0]
for i in range(len(p.get('Amplitude')[indexes][index_low:index_high])):
    pGB_small = {}
    for parameter in parameters:
        pGB_small[parameter] = p.get(parameter)[indexes][index_low:index_high][i]
# pGB = {'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}
# pGB = {'Amplitude': 1.248193e-22, 'EclipticLatitude': -1.185356, 'EclipticLongitude': 3.593803, 'Frequency': 0.003994621, 'FrequencyDerivative': 6.709408e-17, 'Inclination': 1.906596, 'InitialPhase': 5.663538, 'Polarization': 0.96414}
# maxpGB = deepcopy(pGB)
search1 = Search(tdi_fs,Tobs)
pGB = search1.pGB
# maxpGB, pGB =  search1.optimize([[maxpGB]])
# maxpGB = maxpGB[0]
# search1.plot(maxpGB)
loglikelihood = search1.loglikelihood
SNR = search1.SNR
boundaries = search1.boundaries
best_value = loglikelihood([maxpGB])
boundaries_reduced = Reduce_boundaries(maxpGB, boundaries,ratio=0.1)

colors = cm.jet
# plt.figure(figsize=fig_size)
# ax1 = plt.subplot(111)
# ax1.semilogy(tdi_fs.f[range_index],np.abs(tdi_fs['X'][range_index])**2,'k',zorder= 5)
# # ax1.semilogy(tdi_fs_long_subtracted.f[range_index],np.abs(tdi_fs_long_subtracted['X'][range_index])**2,'b',zorder= 5)
# for i in range(len(found_sources)):
#     Xs, Ys, Zs = GB.get_fd_tdixyz(template= found_sources[i], oversample=4, simulator="synthlisa")
#     ax1.semilogy(Xs.f,np.abs(Xs)**2,'--')
# for i in range(len(pGB_injected)):
#     Xs, Ys, Zs = GB.get_fd_tdixyz(template= pGB_injected[i], oversample=4, simulator="synthlisa")
#     a,Xs = xr.align(search1.dataX, Xs, join='left',fill_value=0)
#     ax1.semilogy(Xs.f,np.abs(Xs)**2,label= str(np.round(pGB_injected[i]['Frequency'],0)))
# # Xs, Ys, Zs = GB.get_fd_tdixyz(template= pGB, oversample=4, simulator="synthlisa")
# # a,Xs = xr.align(search1.dataX, Xs, join='left',fill_value=0)
# # ax1.semilogy(Xs.f,np.abs(tdi_fs['X'][range_index]-Xs)**2,label= 'residual')
# ax1.axhline(y=0)
# plt.legend()
# plt.show()

# plt.figure(figsize=fig_size)
# ax1 = plt.subplot(121)
# for i in range(len(pGB_injected)):
#     ax1.scatter(pGB_injected[i]['Frequency'],pGB_injected[i]['Amplitude'], color='gray', marker= 'o', s=100)
# for i in range(len(found_sources)):
#     ax1.scatter(found_sources[i]['Frequency'],found_sources[i]['Amplitude'], marker='o', alpha=0.8, zorder = 4)
# ax1.set_yscale('log')
# ax1.set_ylim(1e-24,1e-21)
# ax1.set_xlim(lower_frequency-padding,upper_frequency+padding)
# ax2 = plt.subplot(122)
# for i in range(len(pGB_injected)):
#     if pGB_injected[i]['EclipticLongitude'] > np.pi:
#         pGB_injected[i]['EclipticLongitude'] -= np.pi*2
#     ax2.scatter(pGB_injected[i]['EclipticLongitude'],pGB_injected[i]['EclipticLatitude'], color='gray', marker= 'o', s=100)
# for i in range(len(found_sources)):
#     ax2.scatter(found_sources[i]['EclipticLongitude'],found_sources[i]['EclipticLatitude'], marker='o', alpha=0.8, zorder = 4)
# plt.show()


# maxpGB = pGB

# Xs_td, Ys_td, Zs_td = GB.get_td_tdixyz(template=pGB, oversample=4, simulator="synthlisa")
# from ldc.lisa.orbits import Orbits
# config = {"initial_position": 0, "initial_rotation": 0, 
#           "nominal_arm_length": 2500000000, "orbit_type": 'analytic'}
# lisa_orbits = Orbits.type(config)

# t_max = 60*60*24*365*2
# t_min = 0
# dt = 15
# from ldc.lisa.projection import ProjectedStrain
# from ldc.waveform.waveform import HpHc

# GB_ldc = HpHc.type('GW1111', "GB", "TD_fdot")
# GB_ldc.set_param(pGB)
# projector = ProjectedStrain(lisa_orbits)
# yArm = projector.arm_response(t_min, t_max, dt, [GB_ldc])
# projector.links
# tdi_X = projector.compute_tdi_x(np.arange(t_min, t_max, dt))
# tdi_X = projector.compute_tdi_x(Xs_td.t.values)
# tdi_X_time = TimeSeries(tdi_X[:int(len(tdi_X))], dt=dt)
# tdi_X_fd = tdi_X_time.ts.fft(win=window)
# index_low = np.searchsorted(tdi_X_fd.f, search1.dataX.f[0])
# tdi_X_fd = tdi_X_fd[index_low : index_low + len(search1.dataX)]
# if len(tdi_X_fd) < len(search1.dataX):
#     a,tdi_X_fd = xr.align(search1.dataX, tdi_X_fd, join='left',fill_value=0)

# plt.figure()
# plt.plot(tdi_X_time.t,tdi_X_time)
# plt.show()

# Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB, oversample=4, simulator="synthlisa")
# index_low = np.searchsorted(Xs.f, search1.dataX.f[0])
# Xs = Xs[index_low : index_low + len(search1.dataX)]
# Ys = Ys[index_low : index_low + len(search1.dataY)]
# Zs = Zs[index_low : index_low + len(search1.dataZ)]
# if len(Xs) < len(search1.dataX):
#     a,Xs = xr.align(search1.dataX, Xs, join='left',fill_value=0)
#     a,Ys = xr.align(search1.dataY, Ys, join='left',fill_value=0)
#     a,Zs = xr.align(search1.dataZ, Zs, join='left',fill_value=0)
# sum = np.abs(Xs.values) ** 2 + np.abs( Ys.values) ** 2 + np.abs(Zs.values) ** 2
# SNR_pGB = float(np.sum(sum / search1.Sn) * Xs.df)
# diffX = tdi_fs['X'][range_index]-Xs
# diffY = np.abs(tdi_fs['Y'][range_index])**2-np.abs(Ys)**2
# diffZ = np.abs(tdi_fs['Z'][range_index])**2-np.abs(Zs)**2
# diff_sum = (np.abs(diffX) + np.abs(diffY) + np.abs(diffZ))/3
# SNR_diff = float(np.sum(diff_sum / search1.Sn) * Xs.df)
# SNR_noise = float(np.sum(search1.Sn / search1.Sn) * Xs.df)
# ratio = SNR_diff/SNR_noise
# # pGB_small = deepcopy(pGB)
# # pGB_small['Amplitude'] = pGB['Amplitude']/100
# Xs_small, Ys_small, Zs_small = GB.get_fd_tdixyz(template=pGB_small, oversample=4, simulator="synthlisa")
# index_low = np.searchsorted(Xs_small.f, search1.dataX.f[0])
# Xs_small = Xs_small[index_low : index_low + len(search1.dataX)]
# Ys_small = Ys_small[index_low : index_low + len(search1.dataY)]
# Zs_small = Zs_small[index_low : index_low + len(search1.dataZ)]
# if len(Xs_small) < len(search1.dataX):
#     a,Xs_small = xr.align(search1.dataX, Xs_small, join='left',fill_value=0)
#     a,Ys_small = xr.align(search1.dataY, Ys_small, join='left',fill_value=0)
#     a,Zs_small = xr.align(search1.dataZ, Zs_small, join='left',fill_value=0)
# plt.figure(figsize=fig_size)
# ax1 = plt.subplot(111)
# ax1.semilogy(Xs.f*1000,np.abs(tdi_fs['X'][range_index])**2,label='Data', color='black')
# ax1.semilogy(Xs.f*1000,np.abs(Xs)**2,label='Signal')
# # ax1.semilogy(Xs.f*1000,tdi_X_fd)**2,label='tdi_X')
# ax1.semilogy(Xs.f*1000,np.abs(diffX)**2,label='Residual')
# ax1.semilogy(Xs.f*1000,np.abs(Xs_small)**2, label= 'Weaker Signal')
# ax1.semilogy(Xs.f*1000,search1.Sn/Xs.df, label='Noise')
# plt.xlabel('$f$ mHz')
# plt.ylabel('Strain $|X|^2$')
# ax1.xaxis.set_major_locator(plt.MaxNLocator(5))
# plt.legend()
# plt.tight_layout()
# plt.show()


import matplotlib.transforms as transforms

def confidence_ellipse(cov, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of `x` and `y`

    Parameters
    ----------
    x, y : array_like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    Returns
    -------
    matplotlib.patches.Ellipse

    Other parameters
    ----------------
    kwargs : `~matplotlib.patches.Patch` properties
    """

    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0),
        width=ell_radius_x * 2,
        height=ell_radius_y * 2,
        facecolor=facecolor,
        **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)

fisher_information_matrix = search1.fisher_information(maxpGB)
FIM = np.zeros((len(parameters),len(parameters)))
for i,parameter1 in enumerate(parameters):
    for j,parameter2 in enumerate(parameters):
        FIM[i,j] = fisher_information_matrix[parameter1][parameter2]
covariance_matrix = scipy.linalg.inv(FIM)
maxpGB01 = scaleto01(maxpGB, search1.boundaries)
cov2d = np.zeros((2,2))
index_parameter1 = parameters.index('Frequency')
index_parameter2 = parameters.index('FrequencyDerivative')
index_parameter1 = parameters.index('EclipticLatitude')
index_parameter2 = parameters.index('EclipticLongitude')
cov2d[0,0] = covariance_matrix[index_parameter1,index_parameter1]
cov2d[0,1] = covariance_matrix[index_parameter1,index_parameter2]
cov2d[1,0] = covariance_matrix[index_parameter2,index_parameter1]
cov2d[1,1] = covariance_matrix[index_parameter2,index_parameter2]
# cov2d = scipy.linalg.inv([[2,1],[1,1]])
x = maxpGB01[parameters[index_parameter1]]
y = maxpGB01[parameters[index_parameter2]]
samples = np.random.multivariate_normal([x,y],cov2d,1000)
samples_x = samples[:,0]
samples_y = samples[:,1]
def eigsorted(cov):
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:,order]
vals, vecs = eigsorted(cov2d)
theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
w, h = 2 * np.sqrt(vals)
lambda_, v = np.linalg.eig(cov2d)
lambda_ = np.sqrt(lambda_)
plt.figure()
ax = plt.subplot(111, )
for j in range(1, 7):
    ell = Ellipse(xy=(x, y),
                width=w*j, height=h*j,
                angle=theta, facecolor='none', edgecolor='k')
    ax.add_artist(ell)
# confidence_ellipse(cov2d, ax, edgecolor='red')
# ax.scatter(samples_x,samples_y)
# ax.quiver(x,y,v[0,0], v[0,1],color='r',scale = 1/lambda_[0], scale_units='xy')
# ax.quiver(x,y,v[1,0], v[1,1],scale = 1/lambda_[1], scale_units='xy')
# ax.scatter(x+v[0,0]*lambda_[0]*5, y+v[1,0]*lambda_[0]*5)
# ax.scatter(x+v[0,1]*lambda_[1]*5, y+v[1,1]*lambda_[1]*5)
ax.scatter(x, y, color='g',label='Found parameters')
plt.xlabel(parameters[index_parameter1])
plt.ylabel(parameters[index_parameter2])
xlim = np.max(v[:,0]*lambda_[0])
ylim = np.max(v[:,1]*lambda_[1])
ax.vlines([x-xlim*5,x+xlim*5],y-ylim*5,y+ylim*5,'r',label='5$\sigma$ Boundary')
ax.hlines([y-ylim*5,y+ylim*5],x-xlim*5,x+xlim*5,'r')
plt.xlim(x-xlim*8,x+xlim*8)
plt.ylim(y-ylim*8,y+ylim*8)
# plt.axis('scaled')
plt.legend()
plt.show()

lambda_, v = scipy.linalg.eig(covariance_matrix)
transf = np.zeros((len(lambda_),len(lambda_)))
for i in range(len(lambda_)):
    transf[i] = v[i] * np.sqrt(lambda_[i])
covariance2 = v @ np.diag(lambda_) @ scipy.linalg.inv(v)
transf = v @ np.diag(np.sqrt(lambda_))
scalematrix = np.max(np.abs(transf), axis=1)

# for parameter in parameters:
#     scalematrix[parameters.index(parameter)] = np.sqrt(covariance_matrix)[parameters.index(parameter)][parameters.index(parameter)]

maxpGB01_low = deepcopy(maxpGB01)
maxpGB01_high = deepcopy(maxpGB01)
boundaries_reduced_fisher = {}
sigma_multiplyer = 5
for parameter in parameters:
    maxpGB01_low[parameter] = maxpGB01[parameter] - scalematrix[parameters.index(parameter)] * sigma_multiplyer 
    maxpGB01_high[parameter] = maxpGB01[parameter] + scalematrix[parameters.index(parameter)] * sigma_multiplyer 
    if parameter in [ 'InitialPhase', 'Polarization']:
        maxpGB01_low[parameter] = maxpGB01[parameter] - scalematrix[parameters.index(parameter)] * 0.01
        maxpGB01_high[parameter] = maxpGB01[parameter] + scalematrix[parameters.index(parameter)] * 0.01
    if parameter == 'FrequencyDerivative':
        if scalematrix[parameters.index(parameter)] > 1:
            maxpGB01_low[parameter] = maxpGB01[parameter] - 0.1
            maxpGB01_high[parameter] = maxpGB01[parameter] + 0.1
    if maxpGB01_low[parameter] < 0:
        maxpGB01_low[parameter] = 0
    if maxpGB01_high[parameter] > 1:
        maxpGB01_high[parameter] = 1
    maxpGB_fisher_low = (maxpGB01_low[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]
    maxpGB_fisher_high = (maxpGB01_high[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]
    boundaries_reduced_fisher[parameter] = [maxpGB_fisher_low, maxpGB_fisher_high]
print(boundaries_reduced_fisher)

print("max SNR", SNR([maxpGB]))
# for n in range(3):
#     resolution = 2000
#     pGB_centered = True
#     std=0.2
#     if n == 1:
#         resolution = 2000
#         pGB_centered = True
#         std=0.3
#     if n == 2:
#         resolution = 2000
#         pGB_centered = True
#         std=0.2
#     resolution_reduced = int(20 ** 2)
#     resolution_reduced = resolution
#     start = time.time()
#     parameter = "Frequency"
#     train_samples = sampler(resolution_reduced,parameters,maxpGB, boundaries_reduced, uniform=False, twoD=False, gaussian=True, pGB_centered=pGB_centered, std=std)
#     print('sample time of', resolution, 'samples ',time.time() - start)
#     train_x = np.zeros((resolution_reduced, len(parameters)))
#     i = 0
#     for parameter in parameters:
#         if parameter == 'FrequencyDerivative':
#             train_x[:, i] = train_samples[parameter]
#         else:
#             train_x[:, i] = train_samples[parameter]
#         i += 1
#     # train_y = train_samples["Likelihood"]
#     # normalizer = np.mean(train_y - train_y.max())
#     # train_y2 = (train_y - train_y.max())/normalizer
#     # train_y = train_samples["Likelihood"]
#     # nu = np.mean(train_y)
#     # sigma = np.std(train_y)
#     # train_y = (train_y - nu) / sigma
#     best_value2 = np.max(train_samples["Likelihood"])
#     cutoff = 2.5
#     if n == 1:
#         cutoff = 1.2
#     if n == 2:
#         cutoff = 1+np.abs(50/best_value)
#     for parameter in parameters:
#         ranges = train_x[:,parameters.index(parameter)][train_samples["Likelihood"]>best_value*cutoff]
#         length_boundary_reduced = boundaries_reduced[parameter][1]-boundaries_reduced[parameter][0]
#         lowerbound = ranges.min()*length_boundary_reduced+boundaries_reduced[parameter][0]
#         upperbound = ranges.max()*length_boundary_reduced+boundaries_reduced[parameter][0]
#         length_bound = upperbound - lowerbound
#         length_bound = boundaries_reduced[parameter][1]- boundaries_reduced[parameter][0]
#         boundaries_reduced[parameter] = [lowerbound-0.1*length_bound, upperbound+0.1*length_bound]
#         if boundaries_reduced[parameter][0] < boundaries[parameter][0]:
#             boundaries_reduced[parameter][0] = boundaries[parameter][0]
#         if boundaries_reduced[parameter][1] > boundaries[parameter][1]:
#             boundaries_reduced[parameter][1] = boundaries[parameter][1]
#     print(len(ranges))
#     print(boundaries_reduced)
# boundaries_reduced = {'Amplitude': [-22.09199095650623, -21.754896777575016], 'EclipticLatitude': [0.23342965233168822, 0.41564173360417933], 'EclipticLongitude': [-2.8350052298807924, -2.6815358667660765], 'Frequency': [0.0013596183498943753, 0.0013596205738444077], 'FrequencyDerivative': [-18.327331554013387, -16.17261318323499], 'Inclination': [0.5693017818133534, 0.9999985157689428], 'InitialPhase': [3.1233384853789494, 3.141436732941647], 'Polarization': [0.4574549865923378, 0.5258289654969943]}
# boundaries_reduced = {'Amplitude': [-22.041512485441654, -21.66085838829452], 'EclipticLatitude': [0.25060684012980305, 0.4090230660047114], 'EclipticLongitude': [-2.8162326961746995, -2.6979109639619887], 'Frequency': [0.0013596183309468168, 0.0013596206907834918], 'FrequencyDerivative': [-18.904545279263328, -16.086135962520995], 'Inclination': [0.3241884292841432, 0.9651799058749854], 'InitialPhase': [0.8480897009220866, 0.9087738353980884], 'Polarization': [2.468374989050734, 2.5358247512505083]}
boundaries_reduced = deepcopy(boundaries_reduced_fisher)
split_fd = -17
if boundaries_reduced['FrequencyDerivative'][1] < split_fd:
    boundaries_reduced['FrequencyDerivative'][0] = -18
    boundaries_reduced['FrequencyDerivative'][1] = -16
# boundaries_reduced['FrequencyDerivative'][0] = -18
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
train_size = 1000
resolution = train_size + 500
start = time.time()
parameter = "Frequency"
samples = sampler(resolution,parameters,maxpGB, boundaries_reduced1, uniform=False, twoD=False)
print('sample time of', resolution, 'samples ',time.time() - start)
train_x = np.zeros((resolution, len(parameters)))
i = 0
boundary_ratio = (boundaries_reduced1['FrequencyDerivative'][1]-boundaries_reduced1['FrequencyDerivative'][0])/(boundaries_reduced2['FrequencyDerivative'][1]-boundaries_reduced1['FrequencyDerivative'][0])
for parameter in parameters:
    if parameter == 'FrequencyDerivative':
        train_x[:, i] = samples[parameter]*boundary_ratio
    else:
        train_x[:, i] = samples[parameter]
    i += 1
train_y = samples["Likelihood"][:train_size]
test_y = samples["Likelihood"][-500:]
test_x = train_x[-500:]
train_x = train_x[:train_size]

nu = np.mean(train_y)
sigma = np.std(train_y)
train_y = (train_y - nu) / sigma
kernel = RBF(length_scale=[1,2,5,1,1,1,1,1],length_scale_bounds=[(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,30),(0.1,30)])
start = time.time()
gpr = GaussianProcessRegressor(kernel=kernel, random_state=0).fit(train_x, train_y)
# gpr = traingpmodelsk(train_x, train_y, kernel)
print('train',time.time() - start)
start = time.time()
observed_pred_sk = gpr.predict(test_x)
print('eval',time.time() - start)
observed_pred_sk_scaled = observed_pred_sk*sigma +nu
print("RMSE ",np.sqrt(mean_squared_error(test_y,observed_pred_sk_scaled)))

maxpGBmod = deepcopy(maxpGB)
maxpGBmod['FrequencyDerivative'] = 10**(np.random.rand()*(boundaries_reduced2['FrequencyDerivative'][1]-boundaries_reduced2['FrequencyDerivative'][0])+boundaries_reduced2['FrequencyDerivative'][0])
start = time.time()
parameter = "Frequency"
samples2 = sampler(resolution,parameters,maxpGBmod, boundaries_reduced2, uniform=False, twoD=False)
print('sample time of', resolution, 'samples ',time.time() - start)
train_x = np.zeros((resolution, len(parameters)))
i = 0
for parameter in parameters:
    if parameter == 'FrequencyDerivative':
        train_x[:, i] = samples2[parameter]*(1-boundary_ratio)+boundary_ratio
    else:
        train_x[:, i] = samples2[parameter]
    i += 1
train_y = samples2["Likelihood"][:train_size]
test_y = samples2["Likelihood"][-500:]
test_x = train_x[-500:]
train_x = train_x[:train_size]

nu2 = np.mean(train_y)
sigma2 = np.std(train_y)
train_y_scaled = (train_y - nu2) / sigma2
test_y_scaled = (test_y - nu2) / sigma2
kernel = RBF(length_scale=[1,2,5,1,1,1,1,1],length_scale_bounds=[(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,30),(0.1,30)])
start = time.time()
# gpr2 = traingpmodelsk(train_x, train_y_scaled, kernel)
gpr2 = GaussianProcessRegressor(kernel=kernel, random_state=0).fit(train_x, train_y_scaled)
print('train',time.time() - start)
start = time.time()
observed_pred_sk = gpr2.predict(test_x)
print('eval',time.time() - start)
observed_pred_sk_scaled2 = observed_pred_sk*sigma2 +nu2
print("RMSE ",np.sqrt(mean_squared_error(test_y,observed_pred_sk_scaled2)))

# first sampling for calculating prior
resolution = 2*10 ** 6
start = time.time()
test_x_m = np.random.uniform(size=(resolution,len(parameters)))
# test_x_m = np.random.normal(loc= 0.5,scale=0.2,size=(resolution,len(parameters)))
# std_scale = 0.4
# test_x_m = scipy.stats.truncnorm.rvs((0-0.5)/std_scale,(1-0.5)/std_scale,loc=0.5,scale=std_scale,size=(resolution,len(parameters)))
test_x_m = test_x_m[test_x_m[:,4].argsort()]

start = time.time()
for i in range(10**4):
    prediction = gpr.predict([test_x_m[0]])
print('time', time.time()-start)
start = time.time()
for i in range(1000):
    prediction = gpr.predict(test_x_m[:10**1])
print('time', time.time()-start)


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

# def proposal(x,x_prime):
#     p = scipy.stats.truncnorm.pdf(x[0],(0-x_prime[0])/std,(1-x_prime[0])/std,loc=x_prime[0],scale=std)
#     q = 1
#     for i in range(len(p)):
#         q *= p[i]
#     return q
# start = time.time()
# x = np.random.uniform(size= 8)
# current_parameters = x
# if x[4] < boundary_ratio:
#     previous_prediction = gpr.predict([x])*sigma +nu
# else:
#     previous_prediction = gpr2.predict([x])*sigma2 +nu2
# mcmc_samples = []
# accepted = 0
# center = 0.5
# std = 0.01
# for i in range(int(resolution/100)):
#     x = np.random.random(size=8)
#     x = scipy.stats.truncnorm.rvs((0-current_parameters)/std,(1-current_parameters)/std,loc=current_parameters,scale=std,size=8)
#     if x[4] < boundary_ratio:
#         prediction = gpr.predict([x])*sigma +nu
#     else:
#         prediction = gpr2.predict([x])*sigma2 +nu2
#     if (np.exp((prediction - previous_prediction)))**(1/10) > np.random.uniform():
#         previous_prediction = prediction
#         # previous_probability = probability[i+1]
#         current_parameters = x
#         rejection_count = 0
#         accepted += 1
#     mcmc_samples.append(current_parameters)    
# mcmc_samples = np.asarray(mcmc_samples)
# print('time MHMC', time.time()-start)
# print('acceptance rate %',accepted/resolution/100*100)

# center = 0.1
# x = scipy.stats.truncnorm.rvs((0-center)/std,(1-center)/std,loc=center,scale=std,size=800)
# plt.figure()
# plt.hist(x, density=True, histtype='stepfilled', alpha=0.2)
# plt.legend(loc='best', frameon=False)
# plt.show()

def plot_corner(mcmc_samples, pGB):
    start = time.time()
    mcmc_samples_rescaled = np.zeros(np.shape(mcmc_samples))
    i = 0
    for parameter in parameters:
        if parameter in ["EclipticLatitude"]:
            mcmc_samples_rescaled[:,i] = np.arcsin((mcmc_samples[:,parameters.index(parameter)] * (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])) + boundaries_reduced[parameter][0])
        elif parameter in ["Inclination"]:
            mcmc_samples_rescaled[:,i] = np.arccos((mcmc_samples[:,parameters.index(parameter)] * (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])) + boundaries_reduced[parameter][0])
        elif parameter in ['Amplitude',"FrequencyDerivative"]:
            mcmc_samples_rescaled[:,i] = 10**((mcmc_samples[:,parameters.index(parameter)] * (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])) + boundaries_reduced[parameter][0])
        else:
            mcmc_samples_rescaled[:,i] = (mcmc_samples[:,parameters.index(parameter)] * (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])) + boundaries_reduced[parameter][0]
        i += 1
    print('time rescale', time.time()-start)
    # start = time.time()
    # df = pd.DataFrame(data=mcmc_samples_rescaled, columns=parameters)
    # df.to_csv('/home/stefan/Repositories/ldc1_evaluation_data/submission/Stefan_LDC14/GW'+str(int(np.round(maxpGB['Frequency']*10**8)))+'.csv',index=False)
    # print('saving time', time.time()-start)

    print('full time', time.time()-first_start)
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
    # pGB = search1.pGB
    for parameter in ['EclipticLongitude','EclipticLatitude','Frequency','FrequencyDerivative','Inclination','Amplitude','InitialPhase','Polarization']:
        if parameter in ['Amplitude','FrequencyDerivative']:
            tr_s[i] = np.log10(pGB[parameter])
            maxvalues[i] = np.log10(maxpGB[parameter])
        else:
            tr_s[i] = pGB[parameter]
            maxvalues[i] = maxpGB[parameter]
        i += 1
    if tr_s[0] > np.pi:
        tr_s[0] -= 2*np.pi
    rng = [0.999, 0.999, 0.999, 0.999, (0, np.pi), (-22,-21.4), 0.999, 0.999]
    fig =  corner.corner(datS[:,:6],  bins=40, hist_kwargs={'density':True, 'lw':3}, plot_datapoints=False, fill_contours=False,  show_titles=True, \
                        color='#348ABD', truths= tr_s[:6], truth_color='k', use_math_test=True, labels= lbls[:6],\
                        levels=[0.9,0.63], title_kwargs={"fontsize": 12})
    # Extract the axes
    ndim = len(parameters)-2
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
    legend_elements = [Line2D([0], [0], color='k', ls='-',lw=2, label='Injected Parameters'),
                     Line2D([0], [0], color='g', ls='-',lw=2, label='Found Parameters')]
    axes[0,3].legend(handles=legend_elements, loc='upper left')
    plt.show()

# plot_corner(mcmc_samples, pGB)

# center = 0.2
# r = scipy.stats.truncnorm.rvs((0-center)/std,(1-center)/std,loc=[0.1,0.2],scale=std,size=2)
# plt.hist(r, density=True, histtype='stepfilled', alpha=0.2)
# plt.legend(loc='best', frameon=False)
# plt.show()

partial_length = 1*10**3
start = time.time()
observed_pred_mean = np.zeros(resolution)
observed_pred_sk = np.zeros(numberinlowerbatch)
observed_pred_sk2 = np.zeros(numberinupperbatch)
def Evaluate(i):
    prediction = gpr.predict(test_x_m1[(i)*partial_length:(i+1)*partial_length])
    return prediction
def Evaluate2(i):
    prediction = gpr2.predict(test_x_m2[(i)*partial_length:(i+1)*partial_length])
    return prediction
for n in range(int(numberinlowerbatch/partial_length)):
    observed_pred_sk[n*partial_length:(n+1)*partial_length] = gpr.predict(test_x_m1[(n)*partial_length:(n+1)*partial_length])
# pool = mp.Pool(mp.cpu_count())
# observed_pred_sk = pool.map(Evaluate, [n for n in range(int(numberinlowerbatch/partial_length))])
# pool.close()
# pool.join()
# print('eval time', time.time()-start)
try:
    observed_pred_sk[int(numberinlowerbatch/partial_length)*partial_length:] = gpr.predict(test_x_m1[int(numberinlowerbatch/partial_length)*partial_length:])
except:
    pass
observed_pred_sk = np.asarray(observed_pred_sk)
observed_pred_sk = observed_pred_sk.reshape(numberinlowerbatch)
observed_pred_mean[:numberinlowerbatch] = observed_pred_sk[:numberinlowerbatch]*sigma +nu
print('eval time', time.time()-start)

start = time.time()
for n in range(int(numberinupperbatch/partial_length)):
    observed_pred_sk2[n*partial_length:(n+1)*partial_length] = gpr2.predict(test_x_m2[(n)*partial_length:(n+1)*partial_length])
try:
    observed_pred_sk2[int(numberinupperbatch/partial_length)*partial_length:] = gpr2.predict(test_x_m2[int(numberinupperbatch/partial_length)*partial_length:])
except:
    pass
observed_pred_sk2 = np.asarray(observed_pred_sk2)
observed_pred_sk2 = observed_pred_sk2.reshape(numberinupperbatch)
observed_pred_mean[-numberinupperbatch:] = observed_pred_sk2*sigma2 +nu2
print('eval time', time.time()-start)

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
if loglikelihood([maxpGBpredicted]) > loglikelihood([maxpGB]):
    maxpGB = maxpGBpredicted
best_value = loglikelihood([maxpGB])
print("pred", max_loglike, "true", loglikelihood([scaletooriginal(max_parameters, boundaries_reduced)]), "max", loglikelihood([maxpGB]), maxpGB)
print("true", SNR([scaletooriginal(max_parameters, boundaries_reduced)]), "max", SNR([maxpGB]))

np.random.shuffle(flatsamplesparameters)
start = time.time()
normalizer = np.sum(np.exp(flatsamplesparameters[:,0]-best_value))
flatsamples_normalized = np.exp(flatsamplesparameters[:,0]-best_value)/normalizer
flatsamples_normalized = flatsamplesparameters[:,0]
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
    if np.exp((flatsamples_normalized[i+1] - previous_p))**(1/10) > np.random.uniform():
        previous_p = flatsamples_normalized[i+1]
        # previous_probability = probability[i+1]
        current_parameters = flatsamplesparameters[i+1,1:]
        rejection_count = 0
        accepted += 1
    mcmc_samples.append(current_parameters)
    # else:
    #     rejection_count += 1
    #     if rejection_count == 10:
    #         mcmc_samples.append(current_parameters)
    #         rejection_count = 0
mcmc_samples = np.asarray(mcmc_samples)
print('time MHMC', time.time()-start)
print('acceptance rate %',accepted/resolution*100)

plot_corner(mcmc_samples, pGB)

# plt.figure()
# plt.scatter(mcmc_samples[:10**4,3],mcmc_samples[:10**4,4])
for round in range(2):

    # covariance_matrix = np.cov(mcmc_samples.T)
    # mean = np.mean(mcmc_samples, axis=0)
    resolution = 2*10**6
    numPoints = 2**6+1
    if round == 1:
        resolution = 1*10**6
        numPoints = 2**6+1
    test_x_m = np.zeros((resolution,len(parameters)))
    ax = np.linspace(-0.15,1.15,numPoints)
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

    # rv = multivariate_normal(mean, covariance_matrix)
    # test_x_m = rv.rvs(resolution)
    # probability = rv.pdf(test_x_m)

    # probability *= pdfs
    for n in range(8):
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
    observed_pred_sk = np.zeros(numberinlowerbatch)
    observed_pred_sk2 = np.zeros(numberinupperbatch)
    for n in range(int(numberinlowerbatch/partial_length)):
        observed_pred_sk[n*partial_length:(n+1)*partial_length] = gpr.predict(test_x_m1[(n)*partial_length:(n+1)*partial_length])
    try:
        observed_pred_sk[int(numberinlowerbatch/partial_length)*partial_length:] = gpr.predict(test_x_m1[int(numberinlowerbatch/partial_length)*partial_length:])
    except:
        pass
    observed_pred_sk = np.asarray(observed_pred_sk)
    observed_pred_sk = observed_pred_sk.reshape(numberinlowerbatch)
    observed_pred_mean[:numberinlowerbatch] = observed_pred_sk[:numberinlowerbatch]*sigma +nu
    print('eval time', time.time()-start)

    start = time.time()
    for n in range(int(numberinupperbatch/partial_length)):
        observed_pred_sk2[n*partial_length:(n+1)*partial_length] = gpr2.predict(test_x_m2[(n)*partial_length:(n+1)*partial_length])
    try:
        observed_pred_sk2[int(numberinupperbatch/partial_length)*partial_length:] = gpr2.predict(test_x_m2[int(numberinupperbatch/partial_length)*partial_length:])
    except:
        pass
    observed_pred_sk2 = np.asarray(observed_pred_sk2)
    observed_pred_sk2 = observed_pred_sk2.reshape(numberinupperbatch)
    observed_pred_mean[-numberinupperbatch:] = observed_pred_sk2*sigma2 +nu2
    print('eval time', time.time()-start)


    flatsamples = np.zeros(len(probability))
    flatsamplesparameters = np.zeros((len(probability),len(parameters)+1))
    flatsamples[:] = observed_pred_mean
    flatsamplesparameters[:,1:] = test_x_m
    flatsamplesparameters[:,0] = observed_pred_mean

    maxindx = np.unravel_index(flatsamplesparameters[:,0].argmax(), flatsamplesparameters[:,0].shape)
    max_parameters = flatsamplesparameters[maxindx[0],1:]
    max_loglike = flatsamplesparameters[:,0].max()
    maxpGBpredicted = scaletooriginal(max_parameters, boundaries_reduced)
    if loglikelihood([maxpGBpredicted]) > loglikelihood([maxpGB]):
        maxpGB = maxpGBpredicted
    best_value = loglikelihood([maxpGB])
    print("pred", max_loglike, "true", loglikelihood([scaletooriginal(max_parameters, boundaries_reduced)]), "max", loglikelihood([maxpGB]), maxpGB)
    print("true", SNR([scaletooriginal(max_parameters, boundaries_reduced)]), "max", SNR([maxpGB]))

    indexes = np.arange(len(probability))
    np.random.shuffle(indexes)
    flatsamplesparameters = flatsamplesparameters[indexes]
    probability = probability[indexes]
    start = time.time()
    normalizer = np.sum(np.exp(flatsamplesparameters[:,0]-best_value))
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
    if round == 0:
        temperature = 2
    else:
        temperature = 1
    for i in range(len(flatsamples_normalized)-1):
        if ((flatsamples_normalized[i+1] / previous_p) * (previous_probability/probability[i+1]))**(1/temperature) > np.random.uniform():
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

    plot_corner(mcmc_samples, pGB)

start = time.time()
mcmc_samples_rescaled = np.zeros(np.shape(mcmc_samples))
i = 0
for parameter in parameters:
    if parameter in ["EclipticLatitude"]:
        mcmc_samples_rescaled[:,i] = np.arcsin((mcmc_samples[:,parameters.index(parameter)] * (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])) + boundaries_reduced[parameter][0])
    elif parameter in ["Inclination"]:
        mcmc_samples_rescaled[:,i] = np.arccos((mcmc_samples[:,parameters.index(parameter)] * (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])) + boundaries_reduced[parameter][0])
    elif parameter in ['Amplitude',"FrequencyDerivative"]:
        mcmc_samples_rescaled[:,i] = 10**((mcmc_samples[:,parameters.index(parameter)] * (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])) + boundaries_reduced[parameter][0])
    else:
        mcmc_samples_rescaled[:,i] = (mcmc_samples[:,parameters.index(parameter)] * (boundaries_reduced[parameter][1] - boundaries_reduced[parameter][0])) + boundaries_reduced[parameter][0]
    i += 1
print('time rescale', time.time()-start)
# start = time.time()
# df = pd.DataFrame(data=mcmc_samples_rescaled, columns=parameters)
# df.to_csv('/home/stefan/Repositories/ldc1_evaluation_data/submission/Stefan_LDC14/GW'+str(int(np.round(maxpGB['Frequency']*10**8)))+'.csv',index=False)
# print('saving time', time.time()-start)

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
# pGB = search1.pGB
for parameter in ['EclipticLongitude','EclipticLatitude','Frequency','FrequencyDerivative','Inclination','Amplitude','InitialPhase','Polarization']:
    if parameter in ['Amplitude','FrequencyDerivative']:
        tr_s[i] = np.log10(pGB[parameter])
        maxvalues[i] = np.log10(maxpGB[parameter])
    else:
        tr_s[i] = pGB[parameter]
        maxvalues[i] = maxpGB[parameter]
    i += 1
if tr_s[0] > np.pi:
    tr_s[0] -= 2*np.pi
rng = [0.999, 0.999, 0.999, 0.999, (0, np.pi), (-22,-21.4), 0.999, 0.999]
fig =  corner.corner(datS[:,:6],  bins=40, hist_kwargs={'density':True, 'lw':3}, plot_datapoints=False, fill_contours=False,  show_titles=True, \
                        color='#348ABD', truths= tr_s[:6], truth_color='k', use_math_test=True, labels= lbls[:6],\
                        levels=[0.9,0.63], title_kwargs={"fontsize": 12})
# Extract the axes
ndim = len(parameters)-2
axes = np.array(fig.axes).reshape((ndim, ndim))
# Loop over the diagonal
for i in range(ndim):
    ax = axes[i, i]
    ax.axvline(maxvalues[i], color="g", label= 'Found parameters')
# Loop over the histograms
for yi in range(ndim):
    for xi in range(yi):
        ax = axes[yi, xi]
        ax.axvline(maxvalues[xi], color="g")
        ax.axhline(maxvalues[yi], color="g")
        ax.scatter(maxvalues[xi], maxvalues[yi],s=130, color="g")
        # ax.xaxis.set_tick_params(labelsize=12)
        # ax.yaxis.set_tick_params(labelsize=12)
# corner.ove(fig, maxvalues[None], marker="s", color="C1")

maxvalues_previous = np.zeros(len(parameters))
maxpGBsearch_previous = {'Amplitude': 5.857579179174051e-23, 'EclipticLatitude': -0.08784800634433576, 'EclipticLongitude': 2.102933791009969, 'Frequency': 0.006220280040884535, 'FrequencyDerivative': 7.503323987637861e-16, 'Inclination': 0.5104961907537644, 'InitialPhase': 3.8035280736734642, 'Polarization': 0.07816858489674128}
i = 0
for parameter in ['EclipticLongitude','EclipticLatitude','Frequency','FrequencyDerivative','Inclination','Amplitude','InitialPhase','Polarization']:
    if parameter in ['Amplitude','FrequencyDerivative']:
        maxvalues_previous[i] = np.log10(maxpGBsearch_previous[parameter])
    else:
        maxvalues_previous[i] = maxpGBsearch_previous[parameter]
    i += 1
for i in range(ndim):
    ax = axes[i, i]
    ax.axvline(maxvalues_previous[i], color="r", label= 'Found parameters')
# Loop over the histograms
for yi in range(ndim):
    for xi in range(yi):
        ax = axes[yi, xi]
        ax.axvline(maxvalues_previous[xi], color="r")
        ax.axhline(maxvalues_previous[yi], color="r")
        ax.plot(maxvalues_previous[xi], maxvalues_previous[yi],'sr')

legend_elements = [Line2D([0], [0], color='k', ls='-',lw=2, label='Injection'),
                   Line2D([0], [0], color='g', ls='-',lw=2, label='Search1'),
                   Line2D([0], [0], color='r', ls='-',lw=2, label='Search2')]

axes[0,3].legend(handles=legend_elements, loc='upper left')
# plt.legend(loc='upper right')
plt.show()

def gauss2d(mu, sigma, to_plot=False):
    w, h = 100, 100
    print(sigma)
    std = [np.sqrt(sigma[0, 0]), np.sqrt(sigma[1, 1])]
    x = np.linspace(mu[0] - 3 * std[0], mu[0] + 3 * std[0], w)
    y = np.linspace(mu[1] - 3 * std[1], mu[1] + 3 * std[1], h)

    x, y = np.meshgrid(x, y)

    x_ = x.flatten()
    y_ = y.flatten()
    xy = np.vstack((x_, y_)).T

    print('s')
    z = multivariate_normal.pdf(xy,mu, sigma)
    print('s')
    # z = normal_rv.pdf(xy)
    z = z.reshape(w, h, order='F')

    if to_plot:
        plt.contourf(x, y, z.T)
        plt.show()

    return z


parameters_recorded[best_run]["Loglikelihood"].append(loglikelihood([maxpGB]))
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
                np.asarray(parameters_recorded[n][0][parameter]) * 10 ** 3,
                np.arange(0, len(parameters_recorded[n][0][parameter])),
            )
        else:
            ax[j, i % 5].plot(
                parameters_recorded[n][0][parameter],
                np.arange(0, len(parameters_recorded[n][0][parameter])),
            )
        i += 1
i = 0
for parameter in parametersfd + ["Log-likelihood"]:
    j = 0
    if i > 3:
        j = 1
    if parameter in ["Log-likelihood"]:
        ax[j, i % 5].axvline(x=loglikelihood([pGB]), color="k", label="True")
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
for n in range(len(pGBmodes[0])):
    i = 0
    for parameter in parametersfd + ["Loglikelihood"]:
        j = 0
        if i > 3:
            j = 1
        if parameter in ["Frequency"]:
            ax[j, i % 5].plot(
                np.asarray(pGBmodes[0][n][parameter]) * 10 ** 3,
                0.5,'.'
            )
        else:
            ax[j, i % 5].plot(
                pGBmodes[0][n][parameter],
                0.5,'.'
            )
        i += 1
i = 0
for parameter in parametersfd + ["Log-likelihood"]:
    j = 0
    if i > 3:
        j = 1
    if parameter in ["Log-likelihood"]:
        ax[j, i % 5].axvline(x=loglikelihood([pGB]), color="k", label="True")
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
