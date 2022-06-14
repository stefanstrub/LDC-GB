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
from functools import partial, total_ordering
import itertools

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
    SNR2 = np.sum( np.real(self.DAf * np.conjugate(Af.data) + self.DEf * np.conjugate(Ef.data))/SA )
    hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /SA)
    dd = np.sum((np.absolute(self.DAf.data)**2 + np.absolute(self.DEf.data)**2) /SA)
    plotIt = False
    if plotIt:
        fig, ax = plt.subplots(nrows=2, sharex=True) 
        ax[0].plot(Af.f, np.abs(self.DAf))
        ax[0].plot(Af.f, np.abs(Af.data))
        
        ax[1].plot(Af.f, np.abs(self.DEf))
        ax[1].plot(Af.f, np.abs(Ef.data))
        plt.show()
        
    p2 = np.sum((np.absolute(self.DAf - Af.data)**2 + np.absolute(self.DEf - Ef.data)**2) /SA) * Xs.df *2
    diff = np.abs(self.DAf - Af.data) ** 2 + np.abs(self.DEf - Ef.data) ** 2
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

        tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft()) for k in ["X", "Y", "Z"]]))
        f, psdX =  scipy.signal.welch(tdi_ts["X"], fs=1.0/dt, window='hanning', nperseg=len(tdi_ts["X"])/1)
        f, psdY =  scipy.signal.welch(tdi_ts["Y"], fs=1.0/dt, window='hanning', nperseg=len(tdi_ts["X"])/1)
        f, psdZ =  scipy.signal.welch(tdi_ts["Z"], fs=1.0/dt, window='hanning', nperseg=len(tdi_ts["X"])/1)
        psd = psdX + psdY + psdZ
        indexes = np.logical_and(f>lower_frequency-padding, f<upper_frequency+padding)
        psd = psd[indexes]

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
            "FrequencyDerivative": [-19.0,-14.5],
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
        indexes = np.logical_and(tdi_fs['X'].f > frequencyrange[0], tdi_fs['X'].f < frequencyrange[1]) 
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
        print('pGB', self.pGB, self.loglikelihood([self.pGB]))

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


    def plot(self, maxpGBs=None, pGBadded=None, added_label='Injection2'):
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

        ax1.semilogy(Af.f* 1000, np.abs(self.DAf), label='Data')
        ax1.semilogy(Af.f* 1000, np.abs(Af.data), marker='.', label='Injection')

        if pGBadded != None:
            Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=pGBadded, oversample=4, simulator="synthlisa")
            index_low = np.searchsorted(Xs.f, self.dataX.f[0])
            Xs = Xs[index_low : index_low + len(self.dataX)]
            Ys = Ys[index_low : index_low + len(self.dataY)]
            Zs = Zs[index_low : index_low + len(self.dataZ)]
            Af = (Zs - Xs)/np.sqrt(2.0)
            ax1.semilogy(Af.f* 1000, np.abs(Af.data), marker='.', label=added_label)

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
        data = dict({"X": self.dataX, "Y": self.dataY, "Z": self.dataZ})
        SNR = compute_tdi_snr(source=source, noise=Nmodel, data=data)
        # p1 = np.exp(p1)
        sum_SNR = SNR['X2'] + SNR['Y2'] + SNR['Z2']
        return sum_SNR#/10000

    def SNR2(self, pGBs):
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
        Af = (Zs_total - Xs_total)/np.sqrt(2.0)
        Ef = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)

        diff = np.abs(Af.values) ** 2 + np.abs(Ef.values) ** 2
        SNR = -float(np.sum(diff / self.SA) * self.dataX.df) /2

        return SNR#/10000

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

    def loglikelihoodsdf(self,pGBs, plotIt = False):
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
        loglik2 = -float(np.sum(diff / self.SA) * self.dataX.df) /2

        scalarproduct_signal_subtracted_data = 4*np.real(np.sum(((self.DAf-Af)*np.conjugate((self.DAf-Af)) / self.SA).values) * self.dataX.df)
        scalarproduct_signal_subtracted_data += 4*np.real(np.sum(((self.DEf-Ef)*np.conjugate(self.DEf-Ef) / self.SE).values) * self.dataX.df)

        scalarproduct_signal = 4*np.real(np.sum((Af*np.conjugate(Af) / self.SA).values) * self.dataX.df)
        scalarproduct_signal += 4*np.real(np.sum((Ef*np.conjugate(Ef) / self.SE).values) * self.dataX.df)
        scalarproduct_data_signal = 4*np.real(np.sum((self.DAf*np.conjugate(Af) / self.SA).values) * self.dataX.df)
        scalarproduct_data_signal += 4*np.real(np.sum((self.DEf*np.conjugate(Ef) / self.SE).values) * self.dataX.df) 
        scalarproduct_data = 4*np.real(np.sum((self.DAf*np.conjugate(self.DAf) / self.SA).values) * self.dataX.df)
        scalarproduct_data += 4*np.real(np.sum((self.DEf*np.conjugate(self.DEf) / self.SE).values) * self.dataX.df) 
        loglik = scalarproduct_data_signal - scalarproduct_signal/2   - scalarproduct_data/2

        return loglik2*4, loglik, -scalarproduct_signal_subtracted_data/2

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

    def loglikelihood(self, pGBs):
        for i in range(len(pGBs)):
            Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs[i], oversample=4, simulator="synthlisa")
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
        SNR2 = np.sum( np.real(self.DAf * np.conjugate(Af.data) + self.DEf * np.conjugate(Ef.data))/self.SA )
        hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /self.SA)
        # dd = np.sum((np.absolute(self.DAf.data)**2 + np.absolute(self.DEf.data)**2) /self.SA)
        plotIt = False
        if plotIt:
            fig, ax = plt.subplots(nrows=2, sharex=True) 
            ax[0].plot(Af.f, np.abs(self.DAf))
            ax[0].plot(Af.f, np.abs(Af.data))
            
            ax[1].plot(Af.f, np.abs(self.DEf))
            ax[1].plot(Af.f, np.abs(Ef.data))
            plt.show()
            
        # p2 = np.sum((np.absolute(self.DAf - Af.data)**2 + np.absolute(self.DEf - Ef.data)**2) /self.SA) * Xs.df *2
        # diff = np.abs(self.DAf - Af.data) ** 2 + np.abs(self.DEf - Ef.data) ** 2
        # p1 = float(np.sum(diff / self.SA) * Xs.df) / 2.0
        # loglik = 4.0*Xs.df*( SNR2 - 0.5 * hh - 0.5 * dd)
        # print(p2, loglik)
        logliks = 4.0*Xs.df*( SNR2 - 0.5 * hh )
        return logliks.values

    def SNRm(self, pGBs):
        for i in range(len(pGBs)):
            Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs[i], oversample=4, simulator="synthlisa")
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
        SNR2 = np.sum( np.real(self.DAf * np.conjugate(Af.data) + self.DEf * np.conjugate(Ef.data))/self.SA )
        hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /self.SA)
        # dd = np.sum((np.absolute(self.DAf.data)**2 + np.absolute(self.DEf.data)**2) /self.SA)
        plotIt = False
        if plotIt:
            fig, ax = plt.subplots(nrows=2, sharex=True) 
            ax[0].plot(Af.f, np.abs(self.DAf))
            ax[0].plot(Af.f, np.abs(Af.data))
            
            ax[1].plot(Af.f, np.abs(self.DEf))
            ax[1].plot(Af.f, np.abs(Ef.data))
            plt.show()
            
        # p2 = np.sum((np.absolute(self.DAf - Af.data)**2 + np.absolute(self.DEf - Ef.data)**2) /self.SA) * Xs.df *2
        # diff = np.abs(self.DAf - Af.data) ** 2 + np.abs(self.DEf - Ef.data) ** 2
        # p1 = float(np.sum(diff / self.SA) * Xs.df) / 2.0
        # loglik = 4.0*Xs.df*( SNR2 - 0.5 * hh - 0.5 * dd)
        # print(p2, loglik)
        SNR = 4.0*Xs.df* hh
        SNR2 = 4.0*Xs.df* SNR2
        return np.sqrt(SNR), np.sqrt(SNR2)

    def differential_evolution_search(self, frequency_boundaries, initial_guess = None):
        bounds = []
        for signal in range(number_of_signals):
            for i in range(8):
                bounds.append((0,1))

        maxpGB = []
        self.boundaries_reduced = deepcopy(self.boundaries)
        self.boundaries_reduced['Frequency'] = frequency_boundaries
        if initial_guess != None:
            initial_guess01 = np.zeros(len(parameters)*number_of_signals)
            for signal in range(number_of_signals):
                pGBstart01 = scaleto01(initial_guess[signal], search1.boundaries_reduced)
                for count, parameter in enumerate(parameters):
                    initial_guess01[count+len(parameters)*signal] = pGBstart01[parameter]
            start = time.time()
            res, energies = differential_evolution(self.function_evolution, bounds=bounds, disp=True, strategy='best1exp', popsize=8,tol= 1e-5, maxiter=1000, recombination=0.75, mutation=(0.5,1), x0=initial_guess01)
            print('time',time.time()-start)
        else:
            start = time.time()
            res, energies = differential_evolution(self.function_evolution, bounds=bounds, disp=True, strategy='best1exp', popsize=8,tol= 1e-5, maxiter=1000, recombination=0.75, mutation=(0.5,1))
            print('time',time.time()-start)
        for signal in range(number_of_signals):
            maxpGB.append(scaletooriginal(res.x[signal*8:signal*8+8],self.boundaries_reduced))
        print(res)
        print(maxpGB)
        print(self.loglikelihood(maxpGB))
        # print(pGB)
        return [maxpGB], energies

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

    def optimize(self, pGBmodes, boundaries = None):
        if boundaries == None:
            boundaries = self.boundaries
        bounds = ()
        for signal in range(number_of_signals_optimize):
            bounds += ((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1))
        for i in range(len(pGBmodes)):
            maxpGB = []
            boundaries_reduced = []
            pGBs01 = []

            for j in range(2):
                x = []
                for signal in range(number_of_signals_optimize):
                    if j == 0:
                        maxpGB.append({})
                        boundaries_reduced.append({})
                        for parameter in parameters:
                            maxpGB[signal][parameter] = pGBmodes[i][signal][parameter]
                    # print(maxpGB)
                    boundaries_reduced[signal] = Reduce_boundaries(maxpGB[signal], boundaries,ratio=0.1)
                    if j == 2:
                        boundaries_reduced[signal] = Reduce_boundaries(maxpGB[signal], boundaries,ratio=0.1)
                    if j in [0]:
                        boundaries_reduced[signal] = deepcopy(boundaries)
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
                # res = scipy.optimize.minimize(self.function, x, args=boundaries_reduced, method='Nelder-Mead', tol=1e-10)
                # res = scipy.optimize.least_squares(self.function, x, args=boundaries_reduced, bounds=bounds)
                for signal in range(number_of_signals_optimize):
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
        for signal in range(number_of_signals_optimize):
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
pGBadded5['Amplitude'] = 6.37823e-23*0.8
pGBadded5['EclipticLatitude'] = -0.2
pGBadded5['EclipticLongitude'] = 1.4
pGBadded5['Frequency'] = 0.0062204
pGBadded5['FrequencyDerivative'] = 3*10**-15
pGBadded5['Inclination'] = 0.5
pGBadded5['InitialPhase'] = 3
pGBadded5['Polarization'] = 2
pGBadded6 = {}
pGBadded6['Amplitude'] = 1.36368e-22
pGBadded6['EclipticLatitude'] = -0.2
pGBadded6['EclipticLongitude'] = 1.4
pGBadded6['Frequency'] = 0.00125313
pGBadded6['FrequencyDerivative'] = 9e-19
pGBadded6['Inclination'] = 0.5
pGBadded6['InitialPhase'] = 3
pGBadded6['Polarization'] = 2
pGBadded7 = {}
pGBadded7['Amplitude'] = 1.36368e-22*0.25
pGBadded7['EclipticLatitude'] = -0.2
pGBadded7['EclipticLongitude'] = 1.4
pGBadded7['Frequency'] = 0.00110457
pGBadded7['FrequencyDerivative'] = 1e-19
pGBadded7['Inclination'] = 0.5
pGBadded7['InitialPhase'] = 3
pGBadded7['Polarization'] = 2
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
sangria_fn = DATAPATH + "/LDC1-3_VGB_v2.hdf5"
# sangria_fn = DATAPATH + "/LDC1-4_GB_v2.hdf5"
# sangria_fn = DATAPATH + "/LDC1-3_VGB_v2_FD_noiseless.hdf5"
FD5 = LISAhdf5(sangria_fn)
Nsrc = FD5.getSourcesNum()
GWs = FD5.getSourcesName()
print("Found %d GW sources: " % Nsrc, GWs)
### TOD make sure GalBin is there
if GWs[0] != "GalBinaries":
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

# reduction = 1
# Tobs_long = float(int(p.get("ObservationDuration")/reduction))
# tdi_ts_long = xr.Dataset(dict([["X", TimeSeries(td[:int(len(td[:,1])/reduction), 1], dt=dt)],["Y", TimeSeries(td[:int(len(td[:,1])/reduction), 2], dt=dt)],
# ["Z", TimeSeries(td[:int(len(td[:,1])/reduction), 3], dt=dt)]]))
# tdi_fs_long = xr.Dataset(dict([["X", tdi_ts_long['X'].ts.fft(win=window)],["Y", tdi_ts_long['Y'].ts.fft(win=window)],["Z", tdi_ts_long['Z'].ts.fft(win=window)]]))
# GB_long = fastGB.FastGB(delta_t=dt, T=Tobs_long)  # in seconds
# for pGBadding in [pGBadded7]:#, pGBadded2, pGBadded3, pGBadded4]:
#     Xs_added, Ys_added, Zs_added = GB.get_fd_tdixyz(template=pGBadding, oversample=4, simulator="synthlisa")
#     source_added = dict({"X": Xs_added, "Y": Ys_added, "Z": Zs_added})
#     index_low = np.searchsorted(tdi_fs["X"].f, Xs_added.f[0])
#     index_high = index_low+len(Xs_added)
#     # tdi_fs['X'] = tdi_fs['X'] #+ Xs_added
#     for k in ["X", "Y", "Z"]:
#         tdi_fs[k].data[index_low:index_high] = tdi_fs[k].data[index_low:index_high] + source_added[k].data
# tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft(dt=dt)) for k, n in [["X", 1], ["Y", 2], ["Z", 3]]]))

for parameter in parameters:
    values = p.get(parameter)
    values = np.append(values, pGBadded7[parameter])
    unit = p.units[parameter]
    p.addPar(parameter,values,unit)


pGB = {}
ind = 0
found_sources = []
target_sources = []
first_start = time.time()
np.random.seed(40) #40
# for ind in range(1,len(p.get('Frequency'))):
number_of_signals = 1
signals_per_subtraction = 1

# pool = mp.Pool(mp.cpu_count())
# observed_pred_sk = pool.map(Evaluate, [n for n in range(int(numberinlowerbatch/partial_length))])
# pool.close()
# pool.join()

lower_frequency = 1.359*10**-3
upper_frequency = 1.360*10**-3
lower_frequency = 1.253*10**-3
upper_frequency = 1.254*10**-3
# lower_frequency = 6.220*10**-3
# upper_frequency = 6.221*10**-3
# lower_frequency = 1.104*10**-3
# upper_frequency = 1.105*10**-3
# lower_frequency = 0.0039945
# upper_frequency = 0.0039955
# lower_frequency = 0.0039955
# upper_frequency = 0.0039965
# lower_frequency = 0.0039965
# upper_frequency = 0.0039975
# lower_frequency = 0.003993
# upper_frequency = 0.003997
# lower_frequency = 0.018311
# upper_frequency = 0.018316
padding = 0.5e-6


# print(5*10**(-6)*pGB['Frequency']**(13/3),8*10**(-8)*pGB['Frequency']**(11/3))

indexes = np.argsort(p.get('Frequency'))
index_low = np.searchsorted(p.get('Frequency')[indexes], lower_frequency)
index_high = np.searchsorted(p.get('Frequency')[indexes], upper_frequency)
pGB_injected = []
for i in range(len(p.get('Amplitude')[indexes][index_low:index_high])):
    pGBs = {}
    for parameter in parameters:
        pGBs[parameter] = p.get(parameter)[indexes][index_low:index_high][i]
    pGB_injected.append(pGBs)
for i in range(len(pGB_injected)):
    print(pGB_injected[i]['Frequency'],pGB_injected[i]['Amplitude'])

maxpGB = [[{'Amplitude': 4.08091270139556e-22, 'EclipticLatitude': 0.8720294908527731, 'EclipticLongitude': 0.48611245457890506, 'Frequency': 0.003995221037609946, 'FrequencyDerivative': 1.0853165500503126e-16, 'Inclination': 1.0235918972353715, 'InitialPhase': 5.45676679259367, 'Polarization': 1.0874299088412087}]]
# maxpGB = [[{'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}]]
# # # maxpGB = [[{'Amplitude': 4.083357969533303e-22, 'EclipticLatitude': 0.8719966900463507, 'EclipticLongitude': 0.4861128797587986, 'Frequency': 0.00399522108238035, 'FrequencyDerivative': 1.0720262111754569e-16, 'Inclination': 1.0241926728648307, 'InitialPhase': 2.319788782634353, 'Polarization': 2.6588421907673028}]]
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
            # Xs_subtracted, Ys_subtracted, Zs_subtracted = GB_long.get_fd_tdixyz(template=maxpGB[j][i], oversample=4, simulator="synthlisa")
            # source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
            # index_low = np.searchsorted(tdi_fs["X"].f, Xs_subtracted.f[0])
            # index_high = index_low+len(Xs_subtracted)
            # for k in ["X", "Y", "Z"]:
            #     tdi_fs_long[k].data[index_low:index_high] = tdi_fs[k].data[index_low:index_high] - source_subtracted[k].data
            # tdi_ts_long = xr.Dataset(dict([(k, tdi_fs_long[k].ts.ifft()) for k in ["X", "Y", "Z"]]))

# found_sources = [{'Amplitude': 4.079729023951356e-22, 'EclipticLatitude': 0.8720028645326899, 'EclipticLongitude': 0.4861184116681905, 'Frequency': 0.0039952210800729485, 'FrequencyDerivative': 1.0729362427985777e-16, 'Inclination': 1.0234054245950561, 'InitialPhase': 5.460969734629043, 'Polarization': 1.0878440488246315}, {'Amplitude': 1.1705162052088404e-22, 'EclipticLatitude': -1.1826812683131485, 'EclipticLongitude': -2.6709159052066718, 'Frequency': 0.003994619913706621, 'FrequencyDerivative': 9.679235392481017e-17, 'Inclination': 1.9399892694512724, 'InitialPhase': 5.609087244892229, 'Polarization': 0.9419296285733922}, {'Amplitude': 1.2885439983397928e-23, 'EclipticLatitude': -0.2160888045559383, 'EclipticLongitude': -1.4636553007260287, 'Frequency': 0.0039953822920595515, 'FrequencyDerivative': 1.344649486483458e-19, 'Inclination': 0.2752472958012907, 'InitialPhase': 1.6595119937657414, 'Polarization': 1.762114328420676}, {'Amplitude': 1.8521127184966895e-23, 'EclipticLatitude': 0.3812186492244738, 'EclipticLongitude': 1.256068131079525, 'Frequency': 0.003994855309629712, 'FrequencyDerivative': 1.353386568754485e-15, 'Inclination': 1.8797565282129438, 'InitialPhase': 0.0, 'Polarization': 1.9493864228996722}, {'Amplitude': 3.0843284877290966e-22, 'EclipticLatitude': -0.2826645193353009, 'EclipticLongitude': 1.4808405521022285, 'Frequency': 0.003995983891137299, 'FrequencyDerivative': 3.998367893417067e-14, 'Inclination': 1.2773900117985248, 'InitialPhase': 4.427814115747999, 'Polarization': 0.9895902856309813}, {'Amplitude': 1.3941165323321472e-23, 'EclipticLatitude': 1.4578858842164413, 'EclipticLongitude': -1.4667050425253096, 'Frequency': 0.003994019592492943, 'FrequencyDerivative': 5.98283827336233e-17, 'Inclination': 2.895997278971493, 'InitialPhase': 4.199368021290813, 'Polarization': 2.407379947374118}]
# found_sources = [{'Amplitude': 4.079729023951356e-22, 'EclipticLatitude': 0.8720028645326899, 'EclipticLongitude': 0.4861184116681905, 'Frequency': 0.0039952210800729485, 'FrequencyDerivative': 1.0729362427985777e-16, 'Inclination': 1.0234054245950561, 'InitialPhase': 5.460969734629043, 'Polarization': 1.0878440488246315}, {'Amplitude': 1.1705162052088404e-22, 'EclipticLatitude': -1.1826812683131485, 'EclipticLongitude': -2.6709159052066718, 'Frequency': 0.003994619913706621, 'FrequencyDerivative': 9.679235392481017e-17, 'Inclination': 1.9399892694512724, 'InitialPhase': 5.609087244892229, 'Polarization': 0.9419296285733922}]
# global fit
# found_sources = [{'Amplitude': 3.976741322779087e-22, 'EclipticLatitude': 0.8707577355672277, 'EclipticLongitude': 0.4867723122067888, 'Frequency': 0.003995220942486532, 'FrequencyDerivative': 1.115280997337434e-16, 'Inclination': 1.0049763281778163, 'InitialPhase': 5.439165687088879, 'Polarization': 1.0869367502767089}, {'Amplitude': 1.2564114323138714e-22, 'EclipticLatitude': -1.1863122399080437, 'EclipticLongitude': -2.684386862483276, 'Frequency': 0.003994620572904299, 'FrequencyDerivative': 7.962182694912732e-17, 'Inclination': 1.9019548827453407, 'InitialPhase': 5.660451596533728, 'Polarization': 0.9494066987794656}, {'Amplitude': 1.987456195914956e-23, 'EclipticLatitude': -0.03622133684600345, 'EclipticLongitude': -1.466599785992938, 'Frequency': 0.003995381912559783, 'FrequencyDerivative': 1.6936719757479413e-17, 'Inclination': 0.8516800024186747, 'InitialPhase': 1.8735152965543351, 'Polarization': 1.8799757192130693}, {'Amplitude': 1.032573964002457e-23, 'EclipticLatitude': 0.11181636107993316, 'EclipticLongitude': 1.1724821666513028, 'Frequency': 0.003994877857204425, 'FrequencyDerivative': 6.111423477668618e-16, 'Inclination': 1.6342104935165682, 'InitialPhase': 0.23562078237472003, 'Polarization': 1.9727482842838349}, {'Amplitude': 2.4324775718877206e-23, 'EclipticLatitude': 0.09888448338517121, 'EclipticLongitude': 1.5787286386765342, 'Frequency': 0.003994999337259406, 'FrequencyDerivative': 1.2069151405362311e-16, 'Inclination': 2.0180223041727, 'InitialPhase': 2.5902331012494084, 'Polarization': 0.4430554181801801}, {'Amplitude': 6.306435210748783e-24, 'EclipticLatitude': 0.5179668223938828, 'EclipticLongitude': -0.9044559666452656, 'Frequency': 0.003995387204398293, 'FrequencyDerivative': 8.7617904575353635e-16, 'Inclination': 2.4397936833696585, 'InitialPhase': 3.955635469291709, 'Polarization': 1.6099375716201092}]
# found_sources = found_sources + [{'Amplitude': 3.0843284877290966e-22, 'EclipticLatitude': -0.2826645193353009, 'EclipticLongitude': 1.4808405521022285, 'Frequency': 0.003995983891137299, 'FrequencyDerivative': 3.998367893417067e-14, 'Inclination': 1.2773900117985248, 'InitialPhase': 4.427814115747999, 'Polarization': 0.9895902856309813}, {'Amplitude': 1.3941165323321472e-23, 'EclipticLatitude': 1.4578858842164413, 'EclipticLongitude': -1.4667050425253096, 'Frequency': 0.003994019592492943, 'FrequencyDerivative': 5.98283827336233e-17, 'Inclination': 2.895997278971493, 'InitialPhase': 4.199368021290813, 'Polarization': 2.407379947374118}]
# for i in range(len(found_sources)):
#     # print(found_sources[-1])
#     # target_sources.append(pGB[j])
#     Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources[i], oversample=4, simulator="synthlisa")
#     source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
#     index_low = np.searchsorted(tdi_fs["X"].f, Xs_subtracted.f[0])
#     index_high = index_low+len(Xs_subtracted)
#     for k in ["X", "Y", "Z"]:
#         tdi_fs[k].data[index_low:index_high] = tdi_fs[k].data[index_low:index_high] - source_subtracted[k].data
#     tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft()) for k in ["X", "Y", "Z"]]))


# initial_guess = [{'Amplitude': 3.4340841609241628e-22, 'EclipticLatitude': -0.5642312761041255, 'EclipticLongitude': -2.540164501994462, 'Frequency': 0.0012531300751455023, 'FrequencyDerivative': 5.762233382334956e-20, 'Inclination': 1.299562538559146, 'InitialPhase': 1.69817397754988, 'Polarization': 1.781921395023023}, {'Amplitude': 1.1851659270590206e-22, 'EclipticLatitude': -0.25537481681007557, 'EclipticLongitude': 1.401133480829003, 'Frequency': 0.0012531299798796004, 'FrequencyDerivative': 3.3962492707384217e-20, 'Inclination': 0.6215859432165541, 'InitialPhase': 2.584930634802693, 'Polarization': 1.808144181137842}]
# initial_guess = [{'Amplitude': 3.4340841609241628e-22, 'EclipticLatitude': -0.5642312761041255, 'EclipticLongitude': -2.540164501994462, 'Frequency': 0.0012531300751455023, 'FrequencyDerivative': 5.762233382334956e-20, 'Inclination': 1.299562538559146, 'InitialPhase': 1.69817397754988, 'Polarization': 1.781921395023023}]

# f, psdX =  scipy.signal.welch(tdi_ts["X"], fs=1.0/dt, window='hanning', nperseg=len(tdi_ts["X"])/1)
# f, psdY =  scipy.signal.welch(tdi_ts["Y"], fs=1.0/dt, window='hanning', nperseg=len(tdi_ts["X"])/1)
# f, psdZ =  scipy.signal.welch(tdi_ts["Z"], fs=1.0/dt, window='hanning', nperseg=len(tdi_ts["X"])/1)
# psd = psdX + psdY + psdZ
# # indexes = np.logical_and(f[peaks]>10**-4, f[peaks]<2*10**-2)
# indexes = np.logical_and(f>lower_frequency-padding, f<upper_frequency+padding)
# psd = psd[indexes]

tdi_fs_input = deepcopy(tdi_fs)
search_input = Search(tdi_fs,Tobs)
# previous_found_sources = [{'Amplitude': 4.079729023951356e-22, 'EclipticLatitude': 0.8720028645326899, 'EclipticLongitude': 0.4861184116681905, 'Frequency': 0.0039952210800729485, 'FrequencyDerivative': 1.0729362427985777e-16, 'Inclination': 1.0234054245950561, 'InitialPhase': 5.460969734629043, 'Polarization': 1.0878440488246315}, {'Amplitude': 1.1705162052088404e-22, 'EclipticLatitude': -1.1826812683131485, 'EclipticLongitude': -2.6709159052066718, 'Frequency': 0.003994619913706621, 'FrequencyDerivative': 9.679235392481017e-17, 'Inclination': 1.9399892694512724, 'InitialPhase': 5.609087244892229, 'Polarization': 0.9419296285733922}, {'Amplitude': 2.3648307475712127e-23, 'EclipticLatitude': 0.24973495547255922, 'EclipticLongitude': 1.5520627725385836, 'Frequency': 0.003995000325578178, 'FrequencyDerivative': 8.210112456087282e-17, 'Inclination': 1.9857061043535384, 'InitialPhase': 2.5061784908179114, 'Polarization': 0.4668885557463254}, {'Amplitude': 2.078513556387228e-23, 'EclipticLatitude': 0.009302285249180517, 'EclipticLongitude': -1.4678510005512306, 'Frequency': 0.00399537934744601, 'FrequencyDerivative': 8.922723172067281e-17, 'Inclination': 0.9997108283183515, 'InitialPhase': 2.20684418221725, 'Polarization': 2.1510593971971494}, {'Amplitude': 2.3009137999076302e-23, 'EclipticLatitude': -1.4620742621532563, 'EclipticLongitude': -2.0177735980924876, 'Frequency': 0.003994020727203907, 'FrequencyDerivative': 2.9326534966910615e-17, 'Inclination': 1.01452098829612, 'InitialPhase': 5.842707916403813, 'Polarization': 1.8952814865950574}, {'Amplitude': 8.809448653203343e-23, 'EclipticLatitude': -0.2091459628466528, 'EclipticLongitude': 1.482052276495815, 'Frequency': 0.003995809866982243, 'FrequencyDerivative': 9.268455145460366e-14, 'Inclination': 0.5396658347048873, 'InitialPhase': 3.0449279758986716, 'Polarization': 0.38958482085232776}, {'Amplitude': 1.3846071189503718e-23, 'EclipticLatitude': -0.015062600723075503, 'EclipticLongitude': 1.9977522235127187, 'Frequency': 0.0039950661341102865, 'FrequencyDerivative': 3.201401195478979e-14, 'Inclination': 1.3755184410838133, 'InitialPhase': 2.297042501340329, 'Polarization': 0.9686693671594404}]
previous_found_sources = [{'Amplitude': 4.084935966774485e-22, 'EclipticLatitude': 0.8719934546490874, 'EclipticLongitude': 0.48611009683797857, 'Frequency': 0.003995221087430858, 'FrequencyDerivative': 1.0704703957490903e-16, 'Inclination': 1.0245091695238984, 'InitialPhase': 2.320136113624083, 'Polarization': 2.65883774239409}, {'Amplitude': 1.170377953453263e-22, 'EclipticLatitude': -1.1827019140449202, 'EclipticLongitude': -2.6708716710257203, 'Frequency': 0.003994619937260686, 'FrequencyDerivative': 9.604827167870394e-17, 'Inclination': 1.9399867466326164, 'InitialPhase': 2.468693959968005, 'Polarization': 2.5128702009090644}]

# for ind in range(1): #[3,8,9]
current_SNR = 100
ind = 0

while current_SNR > 10 and ind < 1:
    ind += 1

    # indexes = np.logical_and(f[peaks]>10**-4, f[peaks]<2*10**-2)

    # f = f[indexes]
    # peaks, properties = scipy.signal.find_peaks(np.log(np.sqrt(psd)/np.mean(np.sqrt(psd))),width=1,  prominence=(0.22,None))
    # peak_index = 2

    # plt.figure()
    # plt.semilogy(f,psd)
    # # plt.scatter(f[peaks],psd[peaks])
    # # plt.scatter(f[peaks[peak_index]],psd[peaks[peak_index]])
    # # plt.scatter(f[properties['left_bases'][peak_index]],psd[properties['left_bases'][peak_index]])
    # # plt.scatter(f[properties['right_bases'][peak_index]],psd[properties['right_bases'][peak_index]])
    # for k in range(len(pGB_injected)):
    #     plt.axvline(pGB_injected[k]['Frequency'])
    # plt.show()

    # indexes = np.logical_and(f[peaks]>lower_frequency-padding, f[peaks]<upper_frequency+padding)
    # peaks = peaks[indexes]
    # if len(peaks) < signals_per_subtraction:
    #     break
    # prominences = properties['prominences'][indexes] 
    # indexes_peaks = np.argsort(prominences)


    # range_index = np.logical_and(tdi_fs.f > lower_frequency-padding, tdi_fs.f < upper_frequency+padding)


    # pGB = {'Amplitude': 2.360780938411229e-22, 'EclipticLatitude': 0.8726817299401121, 'EclipticLongitude': 0.4857930058990844, 'Frequency': 0.0039952210295307105, 'FrequencyDerivative': 1.0762922030519564e-16, 'Inclination': 0.008495799066858739, 'InitialPhase': 0.013026742125237894, 'Polarization': 1.507449983207914}
    search1 = Search(tdi_fs,Tobs)
    # search1.plot()
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
    # N_frequency = 10
    # N_sky = 15
    # F_stat, frequency, eclipticlatitude, eclipticlongitude = search1.f_statistic(lower_frequency-padding, upper_frequency+padding, N_frequency, N_sky)
    # ind = np.unravel_index(np.argmax(F_stat, axis=None), F_stat.shape)
    # if ind[0]>0:
    #     lower_index = ind[0]-1
    # else:
    #     lower_index = ind[0]
    # if ind[0]<N_frequency-1:
    #     upper_index = ind[0]+1
    # else:
    #     upper_index = ind[0]
    # search1.reduced_frequency_boundaries = [frequency[lower_index],frequency[upper_index]]
    # max = np.max(F_stat, axis=(1,2))
    # print('time F-statistic',time.time()-start)
    # Eclipticlatitude,Eclipticlongitude = np.meshgrid(np.arccos(eclipticlatitude),eclipticlongitude)

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

    print('SNR ',np.round(search1.SNR([search1.pGB])))
    print('SNR2', np.round(search1.SNR2([search1.pGB])))
    print('SNR2', np.round(search1.loglikelihood([search1.pGB])))
    print('SNR2', np.round(search1.loglikelihoodsdf([search1.pGB])))
    print('SNRm', np.round(search1.SNRm([search1.pGB]),3))
    # print('SNRflat', np.round(search1.loglikelihoodflat([search1.pGB])))
    search1.plot()#pGBadded=pGBadded5)
    print(pGBadded7["FrequencyDerivative"] * Tobs)
    print('smear f', 300*pGBadded7["Frequency"] * 10**3 / 10**9)
    print(search1.reduced_frequency_boundaries)
    if ind in [1,2]:
        maxpGBsearch = [[previous_found_sources[ind-1]]]
    else:
        maxpGBsearch, energies =  search1.differential_evolution_search(search1.boundaries['Frequency'])#search1.reduced_frequency_boundaries)
    # maxpGBsearch = [[previous_found_sources[ind-1]]]
    print(search1.loglikelihood(maxpGBsearch[0]), ', ', ind, '. Source')
    print('SNRm of found signal', np.round(search1.SNRm(maxpGBsearch[0]),3))
    # if maxpGBsearch[0][0]['Frequency'] > lower_frequency and maxpGBsearch[0][0]['Frequency'] < upper_frequency:
    current_SNR = search1.SNRm(maxpGBsearch[0])[0]
    print('current SNR', current_SNR)
    if current_SNR < 10:
        break

    maxpGB = []
    for j in range(signals_per_subtraction):
        maxpGB.append(maxpGBsearch[j])
        print(maxpGB[-1])
    for j in range(signals_per_subtraction):
        for i in range(number_of_signals):
            found_sources.append(maxpGB[j][i])
    # create two sets of found sources. found_sources_in with signals inside the boundary and founce_sources_out with outside sources
    found_sources_in = []
    found_sources_out = []
    for i in range(len(found_sources)):
        if found_sources[i]['Frequency'] > lower_frequency and found_sources[i]['Frequency'] < upper_frequency:
            found_sources_in.append(found_sources[i])
        else:
            found_sources_out.append(found_sources[i])

    #global optimization
    if len(found_sources_in) > 1:
        if maxpGBsearch[0][0]['Frequency'] > lower_frequency and maxpGBsearch[0][0]['Frequency'] < upper_frequency:
            tdi_fs_subtracted = deepcopy(tdi_fs_input)
            for i in range(len(found_sources_out)):
                Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_out[i], oversample=4, simulator="synthlisa")
                source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
                index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
                index_high = index_low+len(Xs_subtracted)
                for k in ["X", "Y", "Z"]:
                    tdi_fs_subtracted[k].data[index_low:index_high] = tdi_fs_subtracted[k].data[index_low:index_high] - source_subtracted[k].data
                tdi_ts_subtracted = xr.Dataset(dict([(k, tdi_fs_subtracted[k].ts.ifft()) for k in ["X", "Y", "Z"]]))

            search_out_subtracted = Search(tdi_fs_subtracted,Tobs)

            total_boundaries = deepcopy(search1.boundaries)
            amplitudes = []
            for i in range(len(found_sources_in)):
                amplitudes.append(found_sources_in[i]['Amplitude'])
            total_boundaries['Amplitude'] = [np.min(amplitudes),np.max(amplitudes)]
            amplitudes_length = np.log10(total_boundaries['Amplitude'][1]) - np.log10(total_boundaries['Amplitude'][0])
            total_boundaries['Amplitude'] = [np.log10(total_boundaries['Amplitude'][0]) - amplitudes_length/5, np.log10(total_boundaries['Amplitude'][1]) + amplitudes_length/5]

            number_of_signals_optimize = len(found_sources_in)
            start = time.time()
            found_sources_in, pGB = search_out_subtracted.optimize([found_sources_in], boundaries= total_boundaries)
            print(time.time()-start)

    found_sources = found_sources_in + found_sources_out

    #subtract the found sources from original
    tdi_fs = deepcopy(tdi_fs_input)
    for i in range(len(found_sources)):
        Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources[i], oversample=4, simulator="synthlisa")
        source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
        index_low = np.searchsorted(tdi_fs["X"].f, Xs_subtracted.f[0])
        index_high = index_low+len(Xs_subtracted)
        for k in ["X", "Y", "Z"]:
            tdi_fs[k].data[index_low:index_high] = tdi_fs[k].data[index_low:index_high] - source_subtracted[k].data
        tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft()) for k in ["X", "Y", "Z"]]))

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

    # for j in range(signals_per_subtraction):
    #     for i in range(number_of_signals):
    #         Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources[-1], oversample=4, simulator="synthlisa")
    #         source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
    #         index_low = np.searchsorted(tdi_fs["X"].f, Xs_subtracted.f[0])
    #         index_high = index_low+len(Xs_subtracted)
    #         for k in ["X", "Y", "Z"]:
    #             tdi_fs[k].data[index_low:index_high] = tdi_fs[k].data[index_low:index_high] - source_subtracted[k].data
    #         tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft()) for k in ["X", "Y", "Z"]]))

# found_sources = [{'Amplitude': 4.084518685817812e-22, 'EclipticLatitude': 0.8720544682694538, 'EclipticLongitude': 0.48616143577241866, 'Frequency': 0.003995221114824773, 'FrequencyDerivative': 1.0627279371996937e-16, 'Inclination': 1.0243912595492148, 'InitialPhase': 2.3275458253363115, 'Polarization': 2.661177574668723}, {'Amplitude': 1.2784116164619416e-22, 'EclipticLatitude': -1.0178871180619764, 'EclipticLongitude': -2.35290831749591, 'Frequency': 0.003996675772373996, 'FrequencyDerivative': 1.384710135509938e-16, 'Inclination': 0.5820004755676056, 'InitialPhase': 2.1816103995659355, 'Polarization': 2.61397391849868}, {'Amplitude': 1.1694525363544978e-22, 'EclipticLatitude': -1.1828904555328474, 'EclipticLongitude': -2.670951629160976, 'Frequency': 0.003994620638904597, 'FrequencyDerivative': 7.500058756238762e-17, 'Inclination': 1.9398413554462486, 'InitialPhase': 2.5324878804646227, 'Polarization': 2.5074909279157476}, {'Amplitude': 3.9688477398130994e-23, 'EclipticLatitude': -0.506364743710094, 'EclipticLongitude': -1.8889417663188044, 'Frequency': 0.003993697453743798, 'FrequencyDerivative': 2.8619361551095533e-16, 'Inclination': 0.5374525365060974, 'InitialPhase': 4.406419981959724, 'Polarization': 0.30918063627529413}]
# 13 searches failed padding
# found_sources = [{'Amplitude': 4.084518685817812e-22, 'EclipticLatitude': 0.8720544682694538, 'EclipticLongitude': 0.48616143577241866, 'Frequency': 0.003995221114824773, 'FrequencyDerivative': 1.0627279371996937e-16, 'Inclination': 1.0243912595492148, 'InitialPhase': 2.3275458253363115, 'Polarization': 2.661177574668723}, {'Amplitude': 1.2784116164619416e-22, 'EclipticLatitude': -1.0178871180619764, 'EclipticLongitude': -2.35290831749591, 'Frequency': 0.003996675772373996, 'FrequencyDerivative': 1.384710135509938e-16, 'Inclination': 0.5820004755676056, 'InitialPhase': 2.1816103995659355, 'Polarization': 2.61397391849868}, {'Amplitude': 1.1694525363544978e-22, 'EclipticLatitude': -1.1828904555328474, 'EclipticLongitude': -2.670951629160976, 'Frequency': 0.003994620638904597, 'FrequencyDerivative': 7.500058756238762e-17, 'Inclination': 1.9398413554462486, 'InitialPhase': 2.5324878804646227, 'Polarization': 2.5074909279157476}, {'Amplitude': 3.9688477398130994e-23, 'EclipticLatitude': -0.506364743710094, 'EclipticLongitude': -1.8889417663188044, 'Frequency': 0.003993697453743798, 'FrequencyDerivative': 2.8619361551095533e-16, 'Inclination': 0.5374525365060974, 'InitialPhase': 4.406419981959724, 'Polarization': 0.30918063627529413}, {'Amplitude': 2.680742041326846e-23, 'EclipticLatitude': -0.3237623760038203, 'EclipticLongitude': -0.9242656952841735, 'Frequency': 0.003997448902577133, 'FrequencyDerivative': 3.190849044907435e-19, 'Inclination': 1.141783511650644, 'InitialPhase': 5.321118979762403, 'Polarization': 2.9166143914795515}, {'Amplitude': 5.024350074118341e-23, 'EclipticLatitude': 1.1591094125821528, 'EclipticLongitude': -1.293432879111187, 'Frequency': 0.003997300847264504, 'FrequencyDerivative': 3.938743173172861e-14, 'Inclination': 1.2958483135510792, 'InitialPhase': 2.4427399939840218, 'Polarization': 0.00502302016006889}, {'Amplitude': 4.5447987860620805e-23, 'EclipticLatitude': -1.1799816348621956, 'EclipticLongitude': 2.708128011777524, 'Frequency': 0.0039968265895491085, 'FrequencyDerivative': 7.09149889495916e-17, 'Inclination': 1.9538247167592169, 'InitialPhase': 0.13771667304306426, 'Polarization': 3.008866341884971}, {'Amplitude': 1.2826409765758305e-23, 'EclipticLatitude': -0.2120385113643434, 'EclipticLongitude': -1.4676393536271393, 'Frequency': 0.003995379182969041, 'FrequencyDerivative': 9.422367915180947e-17, 'Inclination': 0.28645971756182975, 'InitialPhase': 0.14464225259938196, 'Polarization': 1.116661586733977}, {'Amplitude': 2.1660360999680058e-23, 'EclipticLatitude': 0.10911109863454742, 'EclipticLongitude': -1.0253347180980428, 'Frequency': 0.0039972571752918545, 'FrequencyDerivative': 1.3089094594434518e-14, 'Inclination': 1.1337808217353589, 'InitialPhase': 3.5879602155909414, 'Polarization': 1.2840276118726086}, {'Amplitude': 2.344855148047592e-23, 'EclipticLatitude': 0.1652297587040812, 'EclipticLongitude': 1.5617830271050357, 'Frequency': 0.003995001778895358, 'FrequencyDerivative': 2.001375061373712e-17, 'Inclination': 1.9730358903153253, 'InitialPhase': 5.668087863382861, 'Polarization': 2.040728348975599}, {'Amplitude': 8.977838448662348e-24, 'EclipticLatitude': -0.39871167926116313, 'EclipticLongitude': -1.4118616259196186, 'Frequency': 0.003997244327672964, 'FrequencyDerivative': 1.9947190067115023e-16, 'Inclination': 0.5399769842808929, 'InitialPhase': 4.802472107215573, 'Polarization': 1.8989330074655548}, {'Amplitude': 9.75371944244002e-24, 'EclipticLatitude': -0.27021772245323755, 'EclipticLongitude': -1.7075938963529784, 'Frequency': 0.003997465195487432, 'FrequencyDerivative': 4.987425475178701e-15, 'Inclination': 0.6419861070003106, 'InitialPhase': 1.2364059061328099, 'Polarization': 1.5258660903628385}, {'Amplitude': 1.5833159376414385e-23, 'EclipticLatitude': -0.11751471499147761, 'EclipticLongitude': -1.224993701314089, 'Frequency': 0.003997196499030572, 'FrequencyDerivative': 1.587086987274e-14, 'Inclination': 1.3264421354252152, 'InitialPhase': 1.7434858842768923, 'Polarization': 1.1993078386595686}]
# 13 searches first fail
# found_sources = [{'Amplitude': 2.3966413804344143e-22, 'EclipticLatitude': 0.8729393809823007, 'EclipticLongitude': 0.4859553526167626, 'Frequency': 0.003995220884562515, 'FrequencyDerivative': 1.1256773818293738e-16, 'Inclination': 0.1726427997042454, 'InitialPhase': 5.466428229870653, 'Polarization': 1.095927370758357}, {'Amplitude': 1.2684472310601388e-22, 'EclipticLatitude': -1.0180726547586445, 'EclipticLongitude': -2.3527958049295856, 'Frequency': 0.003996676082229203, 'FrequencyDerivative': 1.2915668254064134e-16, 'Inclination': 0.5693023857628196, 'InitialPhase': 2.4462267089385397, 'Polarization': 2.734663814931631}, {'Amplitude': 1.1643828696816095e-22, 'EclipticLatitude': -1.1834868987372307, 'EclipticLongitude': -2.672627672971815, 'Frequency': 0.003994620592609486, 'FrequencyDerivative': 7.495940317857837e-17, 'Inclination': 1.9388456868032802, 'InitialPhase': 2.4999357390897456, 'Polarization': 2.517229388177081}, {'Amplitude': 6.196784186818641e-23, 'EclipticLatitude': -0.5066510473279098, 'EclipticLongitude': -1.8896629885573255, 'Frequency': 0.003993698031198208, 'FrequencyDerivative': 2.698175984709252e-16, 'Inclination': 1.061123898612959, 'InitialPhase': 0.8067057525376324, 'Polarization': 1.6264986340046887}, {'Amplitude': 5.17466379658133e-23, 'EclipticLatitude': -0.08159725366121216, 'EclipticLongitude': -1.493963102080858, 'Frequency': 0.0039976376833268015, 'FrequencyDerivative': 1.4050300804248634e-20, 'Inclination': 1.2185242286849642, 'InitialPhase': 1.0448542259997833, 'Polarization': 1.3853435322206158}, {'Amplitude': 2.3871503980249664e-23, 'EclipticLatitude': 0.87897091188117, 'EclipticLongitude': 0.49008758543728925, 'Frequency': 0.0039952239063624434, 'FrequencyDerivative': 3.902732737780443e-17, 'Inclination': 3.0572883569389324, 'InitialPhase': 3.8381680825501743, 'Polarization': 2.0557792399866597}, {'Amplitude': 4.356678331738106e-23, 'EclipticLatitude': -1.1788683324441642, 'EclipticLongitude': 2.7094625700576884, 'Frequency': 0.003996826043195161, 'FrequencyDerivative': 8.519795232191783e-17, 'Inclination': 1.9937421739636965, 'InitialPhase': 3.243316941135514, 'Polarization': 1.4374245613436925}, {'Amplitude': 5.334910228399745e-23, 'EclipticLatitude': 0.13364280146766827, 'EclipticLongitude': -1.5729097593216856, 'Frequency': 0.003991794584298278, 'FrequencyDerivative': 1.5425196330033024e-16, 'Inclination': 1.8312717703435197, 'InitialPhase': 0.480897089050449, 'Polarization': 1.9815289385679795}, {'Amplitude': 1.0339250800168157e-23, 'EclipticLatitude': 0.6573957182363996, 'EclipticLongitude': -1.8142755983075673, 'Frequency': 0.00399718802223753, 'FrequencyDerivative': 9.016794956029511e-17, 'Inclination': 1.1934718351017866, 'InitialPhase': 6.283185307179586, 'Polarization': 0.22457475028681004}, {'Amplitude': 3.4499298095663157e-23, 'EclipticLatitude': -0.18624777343886156, 'EclipticLongitude': 2.2815987308068477, 'Frequency': 0.003991605646890133, 'FrequencyDerivative': 3.0815241596532174e-16, 'Inclination': 1.434846080889702, 'InitialPhase': 5.788768807570675, 'Polarization': 0.4839711910015663}, {'Amplitude': 1.8487366046597503e-23, 'EclipticLatitude': -0.30044714452278054, 'EclipticLongitude': -1.7773492502693378, 'Frequency': 0.003998405207937929, 'FrequencyDerivative': 2.241503862830622e-16, 'Inclination': 2.333031391303488, 'InitialPhase': 1.3531966149976113, 'Polarization': 0.0740933768668217}, {'Amplitude': 1.2222497200509406e-23, 'EclipticLatitude': -0.2176251856768252, 'EclipticLongitude': -1.4642754597906145, 'Frequency': 0.0039953819000715, 'FrequencyDerivative': 9.35576606973484e-18, 'Inclination': 0.1559845814326035, 'InitialPhase': 4.32060023347354, 'Polarization': 3.11365966662263}, {'Amplitude': 3.0946975414417e-23, 'EclipticLatitude': -0.5104119773380339, 'EclipticLongitude': 2.7831527046620295, 'Frequency': 0.003991516069922154, 'FrequencyDerivative': 4.483485668265565e-18, 'Inclination': 1.8239644645754165, 'InitialPhase': 2.196777936823289, 'Polarization': 1.7656508073929902}]
# 13 searches etol 1e-4
# found_sources = [{'Amplitude': 4.075196489550978e-22, 'EclipticLatitude': 0.8720314614120313, 'EclipticLongitude': 0.48614308713923204, 'Frequency': 0.003995221171442674, 'FrequencyDerivative': 1.0471090797337425e-16, 'Inclination': 1.0220560390957336, 'InitialPhase': 5.46733994120351, 'Polarization': 1.0868224083759281}, {'Amplitude': 1.345683385380931e-22, 'EclipticLatitude': -1.0179673314774975, 'EclipticLongitude': -2.3527780653351247, 'Frequency': 0.003996675838952437, 'FrequencyDerivative': 1.3663573198093563e-16, 'Inclination': 0.6628712430743936, 'InitialPhase': 2.2861973671726283, 'Polarization': 2.6633897920414764}, {'Amplitude': 1.1665007374305327e-22, 'EclipticLatitude': -1.182687405357272, 'EclipticLongitude': -2.670730631099633, 'Frequency': 0.003994620326281378, 'FrequencyDerivative': 8.356951542815989e-17, 'Inclination': 1.9397835783608024, 'InitialPhase': 2.496234407368622, 'Polarization': 2.5120989956796698}, {'Amplitude': 3.443999955287389e-23, 'EclipticLatitude': -0.5060413496094915, 'EclipticLongitude': -1.8889043773213683, 'Frequency': 0.003993697233294337, 'FrequencyDerivative': 2.9196876150506567e-16, 'Inclination': 0.08478935520548249, 'InitialPhase': 4.511823479832223, 'Polarization': 0.37243437483196656}, {'Amplitude': 1.5489932581229267e-23, 'EclipticLatitude': 0.08657927554360081, 'EclipticLongitude': 1.4212506785944994, 'Frequency': 0.003991707739793144, 'FrequencyDerivative': 2.5450811164809185e-19, 'Inclination': 0.690337600245047, 'InitialPhase': 4.776416946743948, 'Polarization': 2.7285199010849364}, {'Amplitude': 4.8915063162637816e-23, 'EclipticLatitude': 0.05850127222848512, 'EclipticLongitude': -1.5771336670989364, 'Frequency': 0.003991792400468529, 'FrequencyDerivative': 2.0393291309682316e-16, 'Inclination': 1.7946876616624199, 'InitialPhase': 0.21032908786235874, 'Polarization': 2.008309906700344}, {'Amplitude': 5.306481784161249e-23, 'EclipticLatitude': -0.06779363346484671, 'EclipticLongitude': -1.4982937212805427, 'Frequency': 0.0039976305216734365, 'FrequencyDerivative': 2.1745651981347919e-16, 'Inclination': 1.2246094355547388, 'InitialPhase': 3.635593311475872, 'Polarization': 2.94222383865636}, {'Amplitude': 3.170134031682973e-23, 'EclipticLatitude': -0.1638144315627871, 'EclipticLongitude': 2.2805313360736426, 'Frequency': 0.003991603900857941, 'FrequencyDerivative': 3.643939463427972e-16, 'Inclination': 1.333210234648431, 'InitialPhase': 5.651938249501982, 'Polarization': 0.49149396794277106}, {'Amplitude': 3.116936690610746e-23, 'EclipticLatitude': 1.1885433908903278, 'EclipticLongitude': 2.5170459939269474, 'Frequency': 0.003996828976830102, 'FrequencyDerivative': 1.5988168007227264e-20, 'Inclination': 1.0557899828457835, 'InitialPhase': 3.310627116151158, 'Polarization': 1.4377394054608283}, {'Amplitude': 2.3960477109986378e-23, 'EclipticLatitude': 0.1401272397490411, 'EclipticLongitude': 1.8337277998440449, 'Frequency': 0.003991510665799794, 'FrequencyDerivative': 1.5939672325245582e-16, 'Inclination': 1.2550432067088833, 'InitialPhase': 2.6421123815181846, 'Polarization': 1.6813994201476268}, {'Amplitude': 1.5436638529228088e-23, 'EclipticLatitude': -0.2994235536663984, 'EclipticLongitude': -1.777601458181164, 'Frequency': 0.003998405240231357, 'FrequencyDerivative': 2.2463289935659335e-16, 'Inclination': 2.5828881093784606, 'InitialPhase': 1.9768429966929055, 'Polarization': 2.909766221678198}, {'Amplitude': 2.3000388360870768e-23, 'EclipticLatitude': -0.3998074332863729, 'EclipticLongitude': -1.8791198152633801, 'Frequency': 0.003991779496126577, 'FrequencyDerivative': 1.1126971612921899e-18, 'Inclination': 2.0624847223908365, 'InitialPhase': 1.2557671610080516, 'Polarization': 1.3183175758247412}, {'Amplitude': 1.967779472994796e-23, 'EclipticLatitude': 0.60011337590119, 'EclipticLongitude': 1.5876658533766896, 'Frequency': 0.0039950335208913225, 'FrequencyDerivative': 2.8091633584119826e-17, 'Inclination': 2.0899583946535594, 'InitialPhase': 5.773970924058125, 'Polarization': 1.9873414967276564}]
# 2 searches found
# found_sources = [{'Amplitude': 4.109569184608354e-22, 'EclipticLatitude': 0.8719265085917549, 'EclipticLongitude': 0.48606457487968857, 'Frequency': 0.0039952211701314135, 'FrequencyDerivative': 1.0431129985759679e-16, 'Inclination': 1.029183450319783, 'InitialPhase': 5.471513112648915, 'Polarization': 1.0906241224884106}, {'Amplitude': 1.115343508540407e-22, 'EclipticLatitude': -1.0178575092392967, 'EclipticLongitude': -2.3529350317020983, 'Frequency': 0.0039966757231684, 'FrequencyDerivative': 1.3983248496408984e-16, 'Inclination': 0.2658436253432126, 'InitialPhase': 0.3329002340127984, 'Polarization': 1.691687559283213}]
# 2 searches reduced search space
# found_sources = [{'Amplitude': 4.109569184608354e-22, 'EclipticLatitude': 0.8719265085917549, 'EclipticLongitude': 0.48606457487968857, 'Frequency': 0.0039952211701314135, 'FrequencyDerivative': 1.0431129985759679e-16, 'Inclination': 1.029183450319783, 'InitialPhase': 5.471513112648915, 'Polarization': 1.0906241224884106}, {'Amplitude': 1.1665007374305327e-22, 'EclipticLatitude': -1.182687405357272, 'EclipticLongitude': -2.670730631099633, 'Frequency': 0.003994620326281378, 'FrequencyDerivative': 8.356951542815989e-17, 'Inclination': 1.9397835783608024, 'InitialPhase': 2.496234407368622, 'Polarization': 2.5120989956796698}]
# # found_sourcesrun = [{'Amplitude': 1.5364465535602838e-21, 'EclipticLatitude': 0.22842923568790388, 'EclipticLongitude': 3.9876628102916634, 'Frequency': 0.0034068399355733194, 'FrequencyDerivative': 1.7265154770131203e-16, 'Inclination': 1.140895523332838, 'InitialPhase': 1.3231144225092324, 'Polarization': 2.7082821014391674}, {'Amplitude': 1.8299785587003754e-21, 'EclipticLatitude': -0.5717437471770699, 'EclipticLongitude': 5.013423186770661, 'Frequency': 0.0021468057546732903, 'FrequencyDerivative': 2.5679108781904502e-17, 'Inclination': 2.381977466585251, 'InitialPhase': 1.7795097897215049, 'Polarization': 1.361815899661127}, {'Amplitude': 1.3740545095358953e-21, 'EclipticLatitude': 0.4586355029484734, 'EclipticLongitude': 2.326341135168954, 'Frequency': 0.002171903889305508, 'FrequencyDerivative': 2.656711011190014e-17, 'Inclination': 2.6231811067672464, 'InitialPhase': 2.246734702442237, 'Polarization': 0.4505374916175284}, {'Amplitude': 4.85984439130295e-22, 'EclipticLatitude': -0.6103079146464152, 'EclipticLongitude': 3.8165981027013838, 'Frequency': 0.007219358007407393, 'FrequencyDerivative': 3.215874944181044e-15, 'Inclination': 0.3611841295412159, 'InitialPhase': 2.8258044166152523, 'Polarization': 2.426471942484511}, {'Amplitude': 6.653710374755356e-22, 'EclipticLatitude': -0.359628257716563, 'EclipticLongitude': 4.9343500673177365, 'Frequency': 0.009211757029070448, 'FrequencyDerivative': 4.897974394062554e-15, 'Inclination': 0.8418300854577668, 'InitialPhase': 2.120331638811925, 'Polarization': 1.2267534109224667}, {'Amplitude': 3.312821037152804e-22, 'EclipticLatitude': 0.7326377959505177, 'EclipticLongitude': 6.056532678360872, 'Frequency': 0.004022512317404639, 'FrequencyDerivative': 7.928660484261939e-17, 'Inclination': 2.687294151927051, 'InitialPhase': 1.6080815997044122, 'Polarization': 1.9588214370682089}, {'Amplitude': 4.686184942845765e-22, 'EclipticLatitude': 0.20352748572849222, 'EclipticLongitude': 5.007749923410212, 'Frequency': 0.009211757978548653, 'FrequencyDerivative': 4.895550870611016e-15, 'Inclination': 1.1415963061394274, 'InitialPhase': 3.14159265358979, 'Polarization': 2.2644816658566564}, {'Amplitude': 5.18155313331795e-22, 'EclipticLatitude': 0.19478251263514385, 'EclipticLongitude': 5.219095549872743, 'Frequency': 0.009445485983463868, 'FrequencyDerivative': 7.023800557316083e-15, 'Inclination': 0.5777619297133901, 'InitialPhase': 3.141592653589793, 'Polarization': 2.0928411308295223}, {'Amplitude': 5.876055535475725e-22, 'EclipticLatitude': -0.46618198434204294, 'EclipticLongitude': 4.289702767312378, 'Frequency': 0.003077973937546665, 'FrequencyDerivative': 1.025303603104508e-16, 'Inclination': 0.7494287470138635, 'InitialPhase': 0.056843419329868985, 'Polarization': 2.8300541786136124}, {'Amplitude': 4.0839934840674097e-22, 'EclipticLatitude': 0.8720661480722399, 'EclipticLongitude': 0.4861465045465657, 'Frequency': 0.003995221091227311, 'FrequencyDerivative': 1.0688851160394643e-16, 'Inclination': 1.0242692798927258, 'InitialPhase': 2.319609919910819, 'Polarization': 2.6582674105848603}]
# # target_sourcesrun = [{'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}, {'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}, {'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}, {'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}, {'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}, {'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}, {'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}, {'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}, {'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}, {'Amplitude': 3.971727e-22, 'EclipticLatitude': 0.870896, 'EclipticLongitude': 0.486536, 'Frequency': 0.003995221, 'FrequencyDerivative': 1.106162e-16, 'Inclination': 1.009082, 'InitialPhase': 5.44632, 'Polarization': 4.229721}]
# # found_sources19 = [{'Amplitude': 7.158007512978195e-23, 'EclipticLatitude': -0.17879074789424193, 'EclipticLongitude': 4.603147784376865, 'Frequency': 0.00399097103940054, 'FrequencyDerivative': 2.29199470942113e-16, 'Inclination': 2.3395967762010357, 'InitialPhase': 1.6414925267236022, 'Polarization': 1.409625714095948}, {'Amplitude': 5.417316688438823e-23, 'EclipticLatitude': 0.09299982963031442, 'EclipticLongitude': 4.709278461472791, 'Frequency': 0.003991793230604934, 'FrequencyDerivative': 1.85302117460904e-16, 'Inclination': 1.8243956188626644, 'InitialPhase': 0.3119979003693603, 'Polarization': 1.9954273407108052}, {'Amplitude': 2.6586721928556783e-23, 'EclipticLatitude': -0.23347086835824327, 'EclipticLongitude': 4.667459964579811, 'Frequency': 0.003991718511582032, 'FrequencyDerivative': 1e-20, 'Inclination': 2.052394297268333, 'InitialPhase': 2.713778167855026, 'Polarization': 2.6371523097400362}, {'Amplitude': 1.4603651612868757e-23, 'EclipticLatitude': 0.5088538078457533, 'EclipticLongitude': 5.061705367435939, 'Frequency': 0.003991299832643548, 'FrequencyDerivative': 1.422149445394665e-16, 'Inclination': 3.1293675414452404, 'InitialPhase': 0.9692647620129494, 'Polarization': 1.6681641616303722}, {'Amplitude': 3.059691174442375e-23, 'EclipticLatitude': 0.4317147381501876, 'EclipticLongitude': 4.997451644492614, 'Frequency': 0.003989147910089541, 'FrequencyDerivative': 1.8726122581411511e-16, 'Inclination': 1.1753859885192686, 'InitialPhase': 1.6649392950095585e-10, 'Polarization': 1.7744204264185908}, {'Amplitude': 1.0744473573747681e-23, 'EclipticLatitude': -0.04537250792371965, 'EclipticLongitude': 4.75590442505031, 'Frequency': 0.003991905382027591, 'FrequencyDerivative': 1.1687840960558461e-16, 'Inclination': 3.1085189766581087, 'InitialPhase': 2.4798219832899684, 'Polarization': 0.23152794448077565}, {'Amplitude': 6.3758855191988285e-24, 'EclipticLatitude': -0.5711395779756441, 'EclipticLongitude': 4.719785340796603, 'Frequency': 0.003991184090209936, 'FrequencyDerivative': 5.197211730593386e-16, 'Inclination': 0.03112334053551802, 'InitialPhase': 0.1840427721876706, 'Polarization': 2.665465621636476}, {'Amplitude': 2.695963984061525e-23, 'EclipticLatitude': 0.698020789183839, 'EclipticLongitude': 4.544599641191273, 'Frequency': 0.00398921600077084, 'FrequencyDerivative': 2.748523299944693e-17, 'Inclination': 1.566152596953124, 'InitialPhase': 1.6976257328944115, 'Polarization': 0.17969203457850272}, {'Amplitude': 1.7140045975600926e-23, 'EclipticLatitude': 0.5353309026749867, 'EclipticLongitude': 4.661338785721523, 'Frequency': 0.003990068561124402, 'FrequencyDerivative': 1.448478884234823e-20, 'Inclination': 1.0325571150373603, 'InitialPhase': 2.8574124537257295, 'Polarization': 1.2248048144455297}, {'Amplitude': 8.600732555086791e-24, 'EclipticLatitude': -0.464957284959853, 'EclipticLongitude': 4.537557676317077, 'Frequency': 0.0039912439454365765, 'FrequencyDerivative': 1.811792639171847e-16, 'Inclination': 1.009201467440333, 'InitialPhase': 1.7030607900434433e-19, 'Polarization': 3.141592653589793}]
# maxpGBsearch = [[{'Amplitude': 5.857579179174051e-23, 'EclipticLatitude': -0.08784800634433576, 'EclipticLongitude': 2.102933791009969, 'Frequency': 0.006220280040884535, 'FrequencyDerivative': 7.503323987637861e-16, 'Inclination': 0.5104961907537644, 'InitialPhase': 3.8035280736734642, 'Polarization': 0.07816858489674128}]]
# maxpGBsearch = [[{'Amplitude': 2.1808846772946646e-22, 'EclipticLatitude': -0.5283847702958109, 'EclipticLongitude': -2.5389708751966173, 'Frequency': 0.0012531301631903925, 'FrequencyDerivative': 3.6807471100533655e-19, 'Inclination': 0.9737933771974954, 'InitialPhase': 3.95558129985279, 'Polarization': 2.8786856073448313}]]
# maxpGBsearch = [[{'Amplitude': 1.5454899498120303e-22, 'EclipticLatitude': 0.3239828534237335, 'EclipticLongitude': -2.7552505031608576, 'Frequency': 0.0013596198681856582, 'FrequencyDerivative': 6.7913371080273325e-18, 'Inclination': 0.9607561157105087, 'InitialPhase': 0.8835396431378469, 'Polarization': 2.498067711496641}]]
# maxpGBsearch = [[{'Amplitude': 5.79066183343146e-23, 'EclipticLatitude': -0.0878417628451765, 'EclipticLongitude': 2.102946967840941, 'Frequency': 0.006220280042028942, 'FrequencyDerivative': 7.503066211114356e-16, 'Inclination': 0.4881731730267816, 'InitialPhase': 0.7708964347916855, 'Polarization': 1.703361386089743}]]
# maxpGBsearch = [[{'Amplitude': 2.1769714187643506e-22, 'EclipticLatitude': -0.5274539110096076, 'EclipticLongitude': -2.538643237956386, 'Frequency': 0.001253130102785549, 'FrequencyDerivative': 1.540240386132294e-18, 'Inclination': 0.9727231980521276, 'InitialPhase': 0.8034948714833561, 'Polarization': 1.3059518704908937}]]
# maxpGBsearch = [[{'Amplitude': 3.4354064952566296e-22, 'EclipticLatitude': -0.5620249399525753, 'EclipticLongitude': -2.5410710550204922, 'Frequency': 0.0012531297336659981, 'FrequencyDerivative': 1.1347668907852654e-17, 'Inclination': 1.300864172138984, 'InitialPhase': 1.6801339406921025, 'Polarization': 1.7829668872856153}]]
maxpGBsearch = [[{'Amplitude': 2.1806906125113854e-22, 'EclipticLatitude': -0.5283216519943494, 'EclipticLongitude': -2.5390590647151567, 'Frequency': 0.0012531301743037816, 'FrequencyDerivative': 1.5344145457167598e-20, 'Inclination': 0.9737610507888282, 'InitialPhase': 3.963309240828965, 'Polarization': 2.8821924743499507}]]
# two signals
# maxpGBsearch = [[{'Amplitude': 5.650971275229933e-23, 'EclipticLatitude': -0.11623084966211612, 'EclipticLongitude': 2.1020434784406845, 'Frequency': 0.006220280231371811, 'FrequencyDerivative': 7.4488630026342e-16, 'Inclination': 0.3354188141532851, 'InitialPhase': 0.24642299320918085, 'Polarization': 1.3910468453167297}]]
# simultaneous
# maxpGBsearch = [[{'Amplitude': 2.062657723130844e-22, 'EclipticLatitude': -0.5301457680576338, 'EclipticLongitude': -2.5397130334982148, 'Frequency': 0.0012531301530795797, 'FrequencyDerivative': 5.322073428412324e-20, 'Inclination': 0.9230137170021916, 'InitialPhase': 4.001790784490439, 'Polarization': 2.9039001571193293}, {'Amplitude': 1.3315183227559572e-22, 'EclipticLatitude': -0.23145135178640172, 'EclipticLongitude': 1.3950802171550274, 'Frequency': 0.001253130094404677, 'FrequencyDerivative': 1.2483943968407192e-18, 'Inclination': 0.5031237729421293, 'InitialPhase': 1.8934905131551893, 'Polarization': 1.4483055437852574}]]
# sequential
# maxpGBsearch = [[{'Amplitude': 3.4340841609241628e-22, 'EclipticLatitude': -0.5642312761041255, 'EclipticLongitude': -2.540164501994462, 'Frequency': 0.0012531300751455023, 'FrequencyDerivative': 5.762233382334956e-20, 'Inclination': 1.299562538559146, 'InitialPhase': 1.69817397754988, 'Polarization': 1.781921395023023}, {'Amplitude': 1.1851659270590206e-22, 'EclipticLatitude': -0.25537481681007557, 'EclipticLongitude': 1.401133480829003, 'Frequency': 0.0012531299798796004, 'FrequencyDerivative': 3.3962492707384217e-20, 'Inclination': 0.6215859432165541, 'InitialPhase': 2.584930634802693, 'Polarization': 1.808144181137842}]]
# two signals search
# maxpGBsearch = [[{'Amplitude': 2.058659878434126e-22, 'EclipticLatitude': -0.5296350392160815, 'EclipticLongitude': -2.5396889596026955, 'Frequency': 0.0012531300651996446, 'FrequencyDerivative': 2.238804389295947e-18, 'Inclination': 0.9217091079486288, 'InitialPhase': 0.805347902851053, 'Polarization': 1.310212623598777}, {'Amplitude': 1.380834345334383e-22, 'EclipticLatitude': -0.23119253368043985, 'EclipticLongitude': 1.3948998681838116, 'Frequency': 0.001253130132149432, 'FrequencyDerivative': 9.868074590632422e-20, 'Inclination': 0.5626330817112827, 'InitialPhase': 5.1201456274601505, 'Polarization': 3.0607225778403317}]]
# maxpGBsearch = [[found_sources[1]]]
# found_sources = maxpGBsearch[0]
# low SNR 0.4
# maxpGBsearch = [[{'Amplitude': 8.130389665859798e-23, 'EclipticLatitude': -0.1612727865959, 'EclipticLongitude': 1.4432442279097897, 'Frequency': 0.0011045698481549845, 'FrequencyDerivative': 4.3663644378721594e-17, 'Inclination': 0.8760756907394487, 'InitialPhase': 3.014063485216334, 'Polarization': 1.8417354795091958}]]
#low SNR 0.3
# maxpGBsearch = [[{'Amplitude': 6.562776565327098e-23, 'EclipticLatitude': -0.14435385908719103, 'EclipticLongitude': 1.4515027376840077, 'Frequency': 0.0011045699201325618, 'FrequencyDerivative': 5.1245329565638194e-17, 'Inclination': 0.9047061450390467, 'InitialPhase': 0.006805230475244557, 'Polarization': 0.29069015470709964}]]
#low SNR 0.25
# maxpGBsearch = [[{'Amplitude': 6.037224261657494e-23, 'EclipticLatitude': -0.12922571652037737, 'EclipticLongitude': 1.4563047936833247, 'Frequency': 0.001104569849763138, 'FrequencyDerivative': 6.05714132012712e-17, 'Inclination': 0.9594247426082358, 'InitialPhase': 6.25009957316676, 'Polarization': 0.24055729773942267}]]
# found sources large window
# found_sources = [{'Amplitude': 4.109569184608354e-22, 'EclipticLatitude': 0.8719265085917549, 'EclipticLongitude': 0.48606457487968857, 'Frequency': 0.0039952211701314135, 'FrequencyDerivative': 1.0431129985759679e-16, 'Inclination': 1.029183450319783, 'InitialPhase': 5.471513112648915, 'Polarization': 1.0906241224884106}, {'Amplitude': 1.0803296578347787e-22, 'EclipticLatitude': -1.0179511288962544, 'EclipticLongitude': -2.3529400135381033, 'Frequency': 0.0039966759161445325, 'FrequencyDerivative': 1.342287594092441e-16, 'Inclination': 0.05714523526873107, 'InitialPhase': 0.19652634417343698, 'Polarization': 1.6157812641947626}, {'Amplitude': 7.392771321450061e-23, 'EclipticLatitude': 1.1921347847452342, 'EclipticLongitude': -2.866880402412061, 'Frequency': 0.003994619763755466, 'FrequencyDerivative': 9.652804956924732e-17, 'Inclination': 0.9960614898489382, 'InitialPhase': 2.0771763346360217, 'Polarization': 0.2602403484155356}, {'Amplitude': 6.052000597238673e-23, 'EclipticLatitude': -1.1631606640677326, 'EclipticLongitude': -2.4036068448739054, 'Frequency': 0.0039945897062114094, 'FrequencyDerivative': 6.337842485651739e-17, 'Inclination': 1.834189552461693, 'InitialPhase': 5.430323004640078, 'Polarization': 2.2250796760763114}, {'Amplitude': 6.111342145528687e-23, 'EclipticLatitude': -0.5048141424281883, 'EclipticLongitude': -1.8904187985228338, 'Frequency': 0.003993697471568838, 'FrequencyDerivative': 2.8560118450156315e-16, 'Inclination': 1.0428569619306647, 'InitialPhase': 3.9366078589266555, 'Polarization': 0.07131834736022112}, {'Amplitude': 3.1606218513833566e-23, 'EclipticLatitude': 0.2789703710219093, 'EclipticLongitude': 1.5604400217310923, 'Frequency': 0.0039950005098807, 'FrequencyDerivative': 6.427860454177372e-17, 'Inclination': 1.8172151552210758, 'InitialPhase': 2.6406393088881392, 'Polarization': 0.36853170848803696}, {'Amplitude': 2.025273370716966e-23, 'EclipticLatitude': -0.9591831345772711, 'EclipticLongitude': -3.053968784554868, 'Frequency': 0.00399672627261989, 'FrequencyDerivative': 3.006309504285046e-16, 'Inclination': 2.338038388457288, 'InitialPhase': 0.0, 'Polarization': 2.2700843125142045}]
# found_sources = [{'Amplitude': 4.121081985812462e-22, 'EclipticLatitude': 0.8720912997092143, 'EclipticLongitude': 0.486193256097311, 'Frequency': 0.003995221049705896, 'FrequencyDerivative': 1.0837989076641528e-16, 'Inclination': 1.0327323917486553, 'InitialPhase': 2.3202157687117384, 'Polarization': 2.659572961678935}, {'Amplitude': 1.5040603666335883e-22, 'EclipticLatitude': 1.0394646625641053, 'EclipticLongitude': -2.142446904035906, 'Frequency': 0.00399667671467272, 'FrequencyDerivative': 1.105081314114411e-16, 'Inclination': 1.9526161285228707, 'InitialPhase': 3.461589075762154, 'Polarization': 0.12582380766565493}, {'Amplitude': 1.264511080132232e-22, 'EclipticLatitude': -0.9859491231188499, 'EclipticLongitude': -2.529917267730215, 'Frequency': 0.003996712446900811, 'FrequencyDerivative': 9.338729169935121e-20, 'Inclination': 1.1403434998592155, 'InitialPhase': 0.31777310567149253, 'Polarization': 2.6904533118798613}, {'Amplitude': 1.1709284748235206e-22, 'EclipticLatitude': -1.1823578113186826, 'EclipticLongitude': -2.6717703722081003, 'Frequency': 0.003994620355463592, 'FrequencyDerivative': 8.401124965830958e-17, 'Inclination': 1.9389856617468404, 'InitialPhase': 5.65425649384106, 'Polarization': 0.9401872064188714}, {'Amplitude': 6.144572693838021e-23, 'EclipticLatitude': -0.5045599115646083, 'EclipticLongitude': -1.890197934558693, 'Frequency': 0.0039936972016079725, 'FrequencyDerivative': 2.9265422900918785e-16, 'Inclination': 1.0490386225781874, 'InitialPhase': 3.9023251300452775, 'Polarization': 0.06712552589177627}, {'Amplitude': 9.937091181958379e-23, 'EclipticLatitude': 1.1946307669091916, 'EclipticLongitude': -1.9792204752046842, 'Frequency': 0.0039966790273277325, 'FrequencyDerivative': 2.2054029511360456e-19, 'Inclination': 1.6673888659296807, 'InitialPhase': 3.7463458892212884, 'Polarization': 2.374074396897342}]
# small window
# 3.9945 - 3.9955
# found_sources = [{'Amplitude': 4.08091270139556e-22, 'EclipticLatitude': 0.8720294908527731, 'EclipticLongitude': 0.48611245457890506, 'Frequency': 0.003995221037609946, 'FrequencyDerivative': 1.0853165500503126e-16, 'Inclination': 1.0235918972353715, 'InitialPhase': 5.45676679259367, 'Polarization': 1.0874299088412087}, {'Amplitude': 4.241825789139581e-23, 'EclipticLatitude': -1.0829288148278788, 'EclipticLongitude': -2.8880314941141756, 'Frequency': 0.003994655528129478, 'FrequencyDerivative': 1e-20, 'Inclination': 3.138285557134122, 'InitialPhase': 6.271817715240895, 'Polarization': 0.009742582711271345}, {'Amplitude': 5.890739738384561e-23, 'EclipticLatitude': 1.0430343332229675, 'EclipticLongitude': -2.7043753127944035, 'Frequency': 0.003994652194752191, 'FrequencyDerivative': 8.254203907567565e-17, 'Inclination': 1.2413923074046642, 'InitialPhase': 3.976944820641031, 'Polarization': 1.648961773210886}, {'Amplitude': 4.052133117939597e-23, 'EclipticLatitude': -1.3031395106469301, 'EclipticLongitude': 0.7867711859149087, 'Frequency': 0.003994843705923819, 'FrequencyDerivative': 3.290260838486189e-19, 'Inclination': 1.5705158045184489, 'InitialPhase': 4.91001005392994, 'Polarization': 1.8095694071827446}]
# found_sources = [{'Amplitude': 4.08091270139556e-22, 'EclipticLatitude': 0.8720294908527731, 'EclipticLongitude': 0.48611245457890506, 'Frequency': 0.003995221037609946, 'FrequencyDerivative': 1.0853165500503126e-16, 'Inclination': 1.0235918972353715, 'InitialPhase': 5.45676679259367, 'Polarization': 1.0874299088412087}, {'Amplitude': 1.1702671226530929e-22, 'EclipticLatitude': -1.1827612284663045, 'EclipticLongitude': -2.670933609481371, 'Frequency': 0.003994620014958145, 'FrequencyDerivative': 9.373396650887945e-17, 'Inclination': 1.9401898236659252, 'InitialPhase': 5.615785193393405, 'Polarization': 0.9422765624096732}]
# found_sources = [{'Amplitude': 4.091252936013729e-22, 'EclipticLatitude': 0.8719925893796358, 'EclipticLongitude': 0.48611841268552825, 'Frequency': 0.0039952210896465855, 'FrequencyDerivative': 1.0703768603283968e-16, 'Inclination': 1.026056054251666, 'InitialPhase': 5.462174515586841, 'Polarization': 1.0880240462829833}, {'Amplitude': 1.1696512269972643e-22, 'EclipticLatitude': -1.1826720001456954, 'EclipticLongitude': -2.6706071863149625, 'Frequency': 0.003994619916795375, 'FrequencyDerivative': 9.639821852041419e-17, 'Inclination': 1.9405696619404715, 'InitialPhase': 5.607323710845002, 'Polarization': 0.9418562049325352}]
# found_sources = [{'Amplitude': 4.079729023951356e-22, 'EclipticLatitude': 0.8720028645326899, 'EclipticLongitude': 0.4861184116681905, 'Frequency': 0.0039952210800729485, 'FrequencyDerivative': 1.0729362427985777e-16, 'Inclination': 1.0234054245950561, 'InitialPhase': 5.460969734629043, 'Polarization': 1.0878440488246315}, {'Amplitude': 1.1705162052088404e-22, 'EclipticLatitude': -1.1826812683131485, 'EclipticLongitude': -2.6709159052066718, 'Frequency': 0.003994619913706621, 'FrequencyDerivative': 9.679235392481017e-17, 'Inclination': 1.9399892694512724, 'InitialPhase': 5.609087244892229, 'Polarization': 0.9419296285733922}, {'Amplitude': 1.2885439983397928e-23, 'EclipticLatitude': -0.2160888045559383, 'EclipticLongitude': -1.4636553007260287, 'Frequency': 0.0039953822920595515, 'FrequencyDerivative': 1.344649486483458e-19, 'Inclination': 0.2752472958012907, 'InitialPhase': 1.6595119937657414, 'Polarization': 1.762114328420676}, {'Amplitude': 1.8521127184966895e-23, 'EclipticLatitude': 0.3812186492244738, 'EclipticLongitude': 1.256068131079525, 'Frequency': 0.003994855309629712, 'FrequencyDerivative': 1.353386568754485e-15, 'Inclination': 1.8797565282129438, 'InitialPhase': 0.0, 'Polarization': 1.9493864228996722}, {'Amplitude': 3.0843284877290966e-22, 'EclipticLatitude': -0.2826645193353009, 'EclipticLongitude': 1.4808405521022285, 'Frequency': 0.003995983891137299, 'FrequencyDerivative': 3.998367893417067e-14, 'Inclination': 1.2773900117985248, 'InitialPhase': 4.427814115747999, 'Polarization': 0.9895902856309813}, {'Amplitude': 1.3941165323321472e-23, 'EclipticLatitude': 1.4578858842164413, 'EclipticLongitude': -1.4667050425253096, 'Frequency': 0.003994019592492943, 'FrequencyDerivative': 5.98283827336233e-17, 'Inclination': 2.895997278971493, 'InitialPhase': 4.199368021290813, 'Polarization': 2.407379947374118}]
# found_sources = [{'Amplitude': 4.079729023951356e-22, 'EclipticLatitude': 0.8720028645326899, 'EclipticLongitude': 0.4861184116681905, 'Frequency': 0.0039952210800729485, 'FrequencyDerivative': 1.0729362427985777e-16, 'Inclination': 1.0234054245950561, 'InitialPhase': 5.460969734629043, 'Polarization': 1.0878440488246315}, {'Amplitude': 1.1705162052088404e-22, 'EclipticLatitude': -1.1826812683131485, 'EclipticLongitude': -2.6709159052066718, 'Frequency': 0.003994619913706621, 'FrequencyDerivative': 9.679235392481017e-17, 'Inclination': 1.9399892694512724, 'InitialPhase': 5.609087244892229, 'Polarization': 0.9419296285733922}, {'Amplitude': 1.2885439983397928e-23, 'EclipticLatitude': -0.2160888045559383, 'EclipticLongitude': -1.4636553007260287, 'Frequency': 0.0039953822920595515, 'FrequencyDerivative': 1.344649486483458e-19, 'Inclination': 0.2752472958012907, 'InitialPhase': 1.6595119937657414, 'Polarization': 1.762114328420676}, {'Amplitude': 1.8521127184966895e-23, 'EclipticLatitude': 0.3812186492244738, 'EclipticLongitude': 1.256068131079525, 'Frequency': 0.003994855309629712, 'FrequencyDerivative': 1.353386568754485e-15, 'Inclination': 1.8797565282129438, 'InitialPhase': 0.0, 'Polarization': 1.9493864228996722}, {'Amplitude': 3.0843284877290966e-22, 'EclipticLatitude': -0.2826645193353009, 'EclipticLongitude': 1.4808405521022285, 'Frequency': 0.003995983891137299, 'FrequencyDerivative': 3.998367893417067e-14, 'Inclination': 1.2773900117985248, 'InitialPhase': 4.427814115747999, 'Polarization': 0.9895902856309813}, {'Amplitude': 1.3941165323321472e-23, 'EclipticLatitude': 1.4578858842164413, 'EclipticLongitude': -1.4667050425253096, 'Frequency': 0.003994019592492943, 'FrequencyDerivative': 5.98283827336233e-17, 'Inclination': 2.895997278971493, 'InitialPhase': 4.199368021290813, 'Polarization': 2.407379947374118}, {'Amplitude': 1.492854360123542e-23, 'EclipticLatitude': -0.20263042623361782, 'EclipticLongitude': 1.5502602981860276, 'Frequency': 0.003994999671458462, 'FrequencyDerivative': 8.254673943490034e-17, 'Inclination': 2.012993024314833, 'InitialPhase': 2.595715267486518, 'Polarization': 0.33177715919609146}, {'Amplitude': 1.2635482628044953e-23, 'EclipticLatitude': 0.4130973814425737, 'EclipticLongitude': -0.8738984913855745, 'Frequency': 0.003995399634350178, 'FrequencyDerivative': 4.1351087575970064e-16, 'Inclination': 1.8384584899716157, 'InitialPhase': 3.981352021823958, 'Polarization': 1.7277846594839554}]
# while
# found_sources = [{'Amplitude': 4.079729023951356e-22, 'EclipticLatitude': 0.8720028645326899, 'EclipticLongitude': 0.4861184116681905, 'Frequency': 0.0039952210800729485, 'FrequencyDerivative': 1.0729362427985777e-16, 'Inclination': 1.0234054245950561, 'InitialPhase': 5.460969734629043, 'Polarization': 1.0878440488246315}, {'Amplitude': 1.1705162052088404e-22, 'EclipticLatitude': -1.1826812683131485, 'EclipticLongitude': -2.6709159052066718, 'Frequency': 0.003994619913706621, 'FrequencyDerivative': 9.679235392481017e-17, 'Inclination': 1.9399892694512724, 'InitialPhase': 5.609087244892229, 'Polarization': 0.9419296285733922}, {'Amplitude': 2.3648307475712127e-23, 'EclipticLatitude': 0.24973495547255922, 'EclipticLongitude': 1.5520627725385836, 'Frequency': 0.003995000325578178, 'FrequencyDerivative': 8.210112456087282e-17, 'Inclination': 1.9857061043535384, 'InitialPhase': 2.5061784908179114, 'Polarization': 0.4668885557463254}, {'Amplitude': 2.078513556387228e-23, 'EclipticLatitude': 0.009302285249180517, 'EclipticLongitude': -1.4678510005512306, 'Frequency': 0.00399537934744601, 'FrequencyDerivative': 8.922723172067281e-17, 'Inclination': 0.9997108283183515, 'InitialPhase': 2.20684418221725, 'Polarization': 2.1510593971971494}, {'Amplitude': 2.3009137999076302e-23, 'EclipticLatitude': -1.4620742621532563, 'EclipticLongitude': -2.0177735980924876, 'Frequency': 0.003994020727203907, 'FrequencyDerivative': 2.9326534966910615e-17, 'Inclination': 1.01452098829612, 'InitialPhase': 5.842707916403813, 'Polarization': 1.8952814865950574}, {'Amplitude': 8.809448653203343e-23, 'EclipticLatitude': -0.2091459628466528, 'EclipticLongitude': 1.482052276495815, 'Frequency': 0.003995809866982243, 'FrequencyDerivative': 9.268455145460366e-14, 'Inclination': 0.5396658347048873, 'InitialPhase': 3.0449279758986716, 'Polarization': 0.38958482085232776}, {'Amplitude': 1.3846071189503718e-23, 'EclipticLatitude': -0.015062600723075503, 'EclipticLongitude': 1.9977522235127187, 'Frequency': 0.0039950661341102865, 'FrequencyDerivative': 3.201401195478979e-14, 'Inclination': 1.3755184410838133, 'InitialPhase': 2.297042501340329, 'Polarization': 0.9686693671594404}]
# optimized
# found_sources = [{'Amplitude': 3.9677348309561324e-22, 'EclipticLatitude': 0.8706663925878446, 'EclipticLongitude': 0.4864856376491269, 'Frequency': 0.003995221078648298, 'FrequencyDerivative': 1.074654099218549e-16, 'Inclination': 1.003291296310796, 'InitialPhase': 5.449444041401877, 'Polarization': 1.0868222545350241}, {'Amplitude': 1.2608384858561515e-22, 'EclipticLatitude': -1.185729083264578, 'EclipticLongitude': -2.6825275334464163, 'Frequency': 0.003994620257517766, 'FrequencyDerivative': 8.737281362904867e-17, 'Inclination': 1.902022201001059, 'InitialPhase': 5.621548359543076, 'Polarization': 0.953714407264287}, {'Amplitude': 2.843731580031459e-23, 'EclipticLatitude': 0.17063412744965353, 'EclipticLongitude': 1.5712041844800808, 'Frequency': 0.0039949987267497555, 'FrequencyDerivative': 1.1538608864748923e-16, 'Inclination': 1.9427405894898233, 'InitialPhase': 2.286430201434658, 'Polarization': 0.5026719162742546}, {'Amplitude': 2.222162420298637e-23, 'EclipticLatitude': -0.014661688533237226, 'EclipticLongitude': -1.4660897173675091, 'Frequency': 0.003995378909066518, 'FrequencyDerivative': 1.0254090587723359e-16, 'Inclination': 0.982476642481787, 'InitialPhase': 2.2278095805769316, 'Polarization': 2.182965959230866}, {'Amplitude': 2.443611563080014e-23, 'EclipticLatitude': -1.4618031858993956, 'EclipticLongitude': -1.9420454717268312, 'Frequency': 0.003994019505769568, 'FrequencyDerivative': 9.437101202803238e-17, 'Inclination': 1.0579755980776597, 'InitialPhase': 5.7968795135054805, 'Polarization': 1.8223057532979217}, {'Amplitude': 9.48027129302248e-23, 'EclipticLatitude': -0.22312586367228956, 'EclipticLongitude': 1.4464589913829462, 'Frequency': 0.003995800418265822, 'FrequencyDerivative': 9.297672219911292e-14, 'Inclination': 0.4949586216978129, 'InitialPhase': 3.3374448402826093, 'Polarization': 0.33877611887923553}, {'Amplitude': 1.3846071189503718e-23, 'EclipticLatitude': -0.015062600723075503, 'EclipticLongitude': 1.9977522235127187, 'Frequency': 0.0039950661341102865, 'FrequencyDerivative': 3.201401195478979e-14, 'Inclination': 1.3755184410838133, 'InitialPhase': 2.297042501340329, 'Polarization': 0.9686693671594404}]
# pop 8
#found_sources = [{'Amplitude': 3.9733841529745045e-22, 'EclipticLatitude': 0.8706193180380079, 'EclipticLongitude': 0.4868664263075614, 'Frequency': 0.003995221059654803, 'FrequencyDerivative': 1.0934359208926452e-16, 'Inclination': 1.0045035448484279, 'InitialPhase': 2.3069032092004123, 'Polarization': 2.6553854089368647}, {'Amplitude': 1.2518948740021577e-22, 'EclipticLatitude': -1.1858150967677368, 'EclipticLongitude': -2.682631307248918, 'Frequency': 0.00399461989272679, 'FrequencyDerivative': 9.802590130056444e-17, 'Inclination': 1.9071582823150672, 'InitialPhase': 2.443991210402824, 'Polarization': 2.524512740602364}, {'Amplitude': 2.2778660022182263e-23, 'EclipticLatitude': -0.020033846654579877, 'EclipticLongitude': -1.467228783866619, 'Frequency': 0.003995379182745883, 'FrequencyDerivative': 9.668607676384704e-17, 'Inclination': 0.9981284353717432, 'InitialPhase': 2.2620803970769106, 'Polarization': 2.179927018739246}, {'Amplitude': 2.5032393353832746e-23, 'EclipticLatitude': 0.14692326933236322, 'EclipticLongitude': 1.5732237256003172, 'Frequency': 0.0039949984234949444, 'FrequencyDerivative': 1.3304435386633334e-16, 'Inclination': 2.0387570164650017, 'InitialPhase': 2.2755856794210803, 'Polarization': 0.517684250634588}, {'Amplitude': 1.140003286421287e-23, 'EclipticLatitude': -0.35370013343872175, 'EclipticLongitude': -1.6544876148967642, 'Frequency': 0.003995085989736267, 'FrequencyDerivative': 1.657691449365034e-15, 'Inclination': 1.9947521879107803, 'InitialPhase': 2.1949373990730487, 'Polarization': 2.787426174384708}, {'Amplitude': 1.7560329908350063e-23, 'EclipticLatitude': 1.1673302328772148, 'EclipticLongitude': 0.054390645602732324, 'Frequency': 0.003994055682954232, 'FrequencyDerivative': 1.2452890391516023e-19, 'Inclination': 2.195965197657791, 'InitialPhase': 2.661747733201045, 'Polarization': 1.2917392485253638}, {'Amplitude': 9.195847842242905e-24, 'EclipticLatitude': -1.4588669948701634, 'EclipticLongitude': -2.1720044580476228, 'Frequency': 0.003994016023749778, 'FrequencyDerivative': 1.7848780180150961e-16, 'Inclination': 0.5776332095957758, 'InitialPhase': 1.5934082149022513, 'Polarization': 0.018023880549174257}, {'Amplitude': 1.3014877346399066e-23, 'EclipticLatitude': -1.1186543836182052, 'EclipticLongitude': -2.4887798439825097, 'Frequency': 0.003994174327120818, 'FrequencyDerivative': 3.426264409182625e-16, 'Inclination': 1.727693731912528, 'InitialPhase': 2.7644941573107173, 'Polarization': 0.9114844778043095}]
# 3.9965 - 3.9975 found sources
# found_sources = [{'Amplitude': 1.3193908464998936e-22, 'EclipticLatitude': -1.0180158762218303, 'EclipticLongitude': -2.3525546025237745, 'Frequency': 0.003996675600988893, 'FrequencyDerivative': 1.429993058592808e-16, 'Inclination': 0.6324098371137679, 'InitialPhase': 5.472624184261258, 'Polarization': 1.1264185956876631}, {'Amplitude': 5.373167674103956e-23, 'EclipticLatitude': -0.06058597427214923, 'EclipticLongitude': -1.4980227952309695, 'Frequency': 0.0039976302585943676, 'FrequencyDerivative': 2.2468475627445273e-16, 'Inclination': 1.234102701454988, 'InitialPhase': 3.6298354756627877, 'Polarization': 2.949206485578923}, {'Amplitude': 2.894245043610553e-23, 'EclipticLatitude': -0.4006057347422682, 'EclipticLongitude': 1.8564312823438724, 'Frequency': 0.003996567161691237, 'FrequencyDerivative': 4.051814603007406e-14, 'Inclination': 2.096582021970788, 'InitialPhase': 4.2407504417170125, 'Polarization': 0.5037630610377921}, {'Amplitude': 3.5865639021547016e-23, 'EclipticLatitude': -1.180315833499288, 'EclipticLongitude': 2.713964964253206, 'Frequency': 0.00399682396359513, 'FrequencyDerivative': 1.3942897460563206e-16, 'Inclination': 1.9215353969762918, 'InitialPhase': 2.9891619120375634, 'Polarization': 1.4518253339718061}, {'Amplitude': 7.394658008831879e-24, 'EclipticLatitude': -0.06840708862901682, 'EclipticLongitude': -1.510166881384511, 'Frequency': 0.003997284247341794, 'FrequencyDerivative': 6.553769317842602e-20, 'Inclination': 0.028501273528823836, 'InitialPhase': 3.285710858818603, 'Polarization': 1.0337908605128148}]
# found_sources = [{'Amplitude': 1.7510613538120334e-22, 'EclipticLatitude': -1.0192683275620062, 'EclipticLongitude': -2.350718675536724, 'Frequency': 0.003996675352787596, 'FrequencyDerivative': 1.5139757310327714e-16, 'Inclination': 0.9597111050835103, 'InitialPhase': 1.6718791011749683, 'Polarization': 2.370918294913893}, {'Amplitude': 4.9751208361442293e-23, 'EclipticLatitude': -1.177830138214736, 'EclipticLongitude': 2.722444456031965, 'Frequency': 0.00399682720006738, 'FrequencyDerivative': 8.771998230895396e-17, 'Inclination': 1.97302972008595, 'InitialPhase': 0.31883344851049383, 'Polarization': 3.0233806725362378}, {'Amplitude': 6.0502612169381655e-24, 'EclipticLatitude': -0.460672559731572, 'EclipticLongitude': 1.589505286128624, 'Frequency': 0.003997372196706388, 'FrequencyDerivative': 2.127505554707063e-16, 'Inclination': 1.735674843936081, 'InitialPhase': 4.849815528656579, 'Polarization': 2.4200560745703554}, {'Amplitude': 1.496896814581308e-23, 'EclipticLatitude': -0.019159730592215457, 'EclipticLongitude': -1.5184009286388875, 'Frequency': 0.003997280708194298, 'FrequencyDerivative': 1.1194290026277205e-16, 'Inclination': 0.9946102604453, 'InitialPhase': 4.695690253112564, 'Polarization': 1.8314674785542293}, {'Amplitude': 6.998698497967915e-24, 'EclipticLatitude': -0.04736401762331348, 'EclipticLongitude': -1.5815185569281762, 'Frequency': 0.003996754044603449, 'FrequencyDerivative': 1.040146797248039e-17, 'Inclination': 0.8264565120853646, 'InitialPhase': 4.482587162345392, 'Polarization': 2.8214607057360186}, {'Amplitude': 5.3762238773603e-23, 'EclipticLatitude': -0.060820759097666155, 'EclipticLongitude': -1.498009226190154, 'Frequency': 0.003997630265868493, 'FrequencyDerivative': 2.2446719053639957e-16, 'Inclination': 1.2344438327469707, 'InitialPhase': 3.6304313859043877, 'Polarization': 2.9492655452487377}, {'Amplitude': 5.7002100900321295e-24, 'EclipticLatitude': -0.7219885404615313, 'EclipticLongitude': -2.1430058237453133, 'Frequency': 0.003997711419694218, 'FrequencyDerivative': 2.585470835870089e-16, 'Inclination': 0.0623872574203187, 'InitialPhase': 3.190774806492449, 'Polarization': 0.04088407809367899}, {'Amplitude': 6.859636546160481e-24, 'EclipticLatitude': -0.4071166242504351, 'EclipticLongitude': 2.0559003093473507, 'Frequency': 0.003996011240694136, 'FrequencyDerivative': 1.5010185886498453e-19, 'Inclination': 1.0472696460734363, 'InitialPhase': 1.71961129969503, 'Polarization': 0.9595110678312323}, {'Amplitude': 5.175803893425514e-24, 'EclipticLatitude': -0.11098539794746572, 'EclipticLongitude': 3.019749530976217, 'Frequency': 0.003996035744244042, 'FrequencyDerivative': 1.8469359324270092e-17, 'Inclination': 2.3602365590548366, 'InitialPhase': 2.6489983988703703, 'Polarization': 0.23851625065642815}, {'Amplitude': 1.0843281132833893e-23, 'EclipticLatitude': 0.33068951070819524, 'EclipticLongitude': -2.1529916057446945, 'Frequency': 0.003996151760101549, 'FrequencyDerivative': 9.515548206546911e-17, 'Inclination': 1.6711838817576194, 'InitialPhase': 4.37199120019937, 'Polarization': 1.9778743681379183}]
# found_sources = [{'Amplitude': 1.7560237486055259e-22, 'EclipticLatitude': -1.0192448183379537, 'EclipticLongitude': -2.350656863293792, 'Frequency': 0.00399667533995668, 'FrequencyDerivative': 1.5164924284872915e-16, 'Inclination': 0.9621827919705668, 'InitialPhase': 1.6700174436274142, 'Polarization': 2.3707365365642756}, {'Amplitude': 4.989412601437577e-23, 'EclipticLatitude': -1.1774710636219279, 'EclipticLongitude': 2.722504272327886, 'Frequency': 0.00399682715510749, 'FrequencyDerivative': 8.956203203599157e-17, 'Inclination': 1.9734663937479489, 'InitialPhase': 3.454469414790749, 'Polarization': 1.4554272881129424}, {'Amplitude': 1.2537248867582883e-23, 'EclipticLatitude': -0.01512710397853698, 'EclipticLongitude': -1.5188657491543964, 'Frequency': 0.0039972805485313035, 'FrequencyDerivative': 1.1572222447467635e-16, 'Inclination': 0.7551265315420236, 'InitialPhase': 0.9856521199234776, 'Polarization': 3.141592653589793}, {'Amplitude': 7.554281389045426e-24, 'EclipticLatitude': 0.011908469219195427, 'EclipticLongitude': -1.5825321550283402, 'Frequency': 0.003996752046259143, 'FrequencyDerivative': 5.665222864836095e-17, 'Inclination': 0.89516476593619, 'InitialPhase': 0.9214846806606066, 'Polarization': 1.1517723555411419}, {'Amplitude': 5.375823613868071e-23, 'EclipticLatitude': -0.059520167092247916, 'EclipticLongitude': -1.4979281889590126, 'Frequency': 0.003997630171346563, 'FrequencyDerivative': 2.257537464830536e-16, 'Inclination': 1.2338569231632452, 'InitialPhase': 0.47757645789016445, 'Polarization': 1.379680964671422}, {'Amplitude': 8.568031437847458e-24, 'EclipticLatitude': -0.5489952437540906, 'EclipticLongitude': -2.046634003577373, 'Frequency': 0.003997687176509116, 'FrequencyDerivative': 8.223722922267948e-17, 'Inclination': 0.8466577401148874, 'InitialPhase': 3.7882434134671166, 'Polarization': 1.4598517632485832}, {'Amplitude': 6.669476779175032e-24, 'EclipticLatitude': -0.41203126110425947, 'EclipticLongitude': 2.0613971026234106, 'Frequency': 0.003996011095982604, 'FrequencyDerivative': 1.1059579624914732e-19, 'Inclination': 1.0252371608521584, 'InitialPhase': 1.6716288430304367, 'Polarization': 0.9609713737452854}, {'Amplitude': 5.156944876431513e-24, 'EclipticLatitude': -0.10745068215448333, 'EclipticLongitude': 3.021064491257313, 'Frequency': 0.003996035677820067, 'FrequencyDerivative': 2.1298455565525138e-17, 'Inclination': 2.360684148410426, 'InitialPhase': 2.682912288664525, 'Polarization': 0.2204737155565217}, {'Amplitude': 7.624822267375493e-24, 'EclipticLatitude': -0.19733236233887655, 'EclipticLongitude': 1.3658024921374485, 'Frequency': 0.0039976382379903015, 'FrequencyDerivative': 5.825892282310141e-19, 'Inclination': 1.7135452830719313, 'InitialPhase': 1.5986149099060554, 'Polarization': 2.792974723339009}, {'Amplitude': 1.0088748678874333e-23, 'EclipticLatitude': -0.3128966726603072, 'EclipticLongitude': -2.1112308481608304, 'Frequency': 0.00399614536188065, 'FrequencyDerivative': 3.358660415587361e-16, 'Inclination': 1.7163765470046464, 'InitialPhase': 3.6973690164967232, 'Polarization': 2.1475366892182257}]
# 3.9955 - 3.9965
# found_sources = [{'Amplitude': 4.1122160565447513e-22, 'EclipticLatitude': 0.8722627178988783, 'EclipticLongitude': 0.48654113291089507, 'Frequency': 0.0039952212328801095, 'FrequencyDerivative': 1.0325821075505724e-16, 'Inclination': 1.0337801752936657, 'InitialPhase': 5.471782551667702, 'Polarization': 1.085645029043506}, {'Amplitude': 1.3196087341219487e-22, 'EclipticLatitude': -1.01802415675539, 'EclipticLongitude': -2.352845358405966, 'Frequency': 0.003996675714076671, 'FrequencyDerivative': 1.4023967877038525e-16, 'Inclination': 0.6321538008842804, 'InitialPhase': 2.3413920047205843, 'Polarization': 2.6965969833231753}, {'Amplitude': 4.907896845237119e-23, 'EclipticLatitude': -0.9083115358404816, 'EclipticLongitude': 1.830946807621947, 'Frequency': 0.003996986130665041, 'FrequencyDerivative': 1.0020926813443886e-20, 'Inclination': 1.9049727112005492, 'InitialPhase': 4.4254524562509605, 'Polarization': 2.03493839501435}, {'Amplitude': 2.80765572753344e-23, 'EclipticLatitude': 0.44991418886259016, 'EclipticLongitude': 1.4372067700793787, 'Frequency': 0.003995, 'FrequencyDerivative': 7.821411239324986e-17, 'Inclination': 1.9770691405054526, 'InitialPhase': 0.008617516017601046, 'Polarization': 2.135208893996428}]
# new tol
# found_sources = [{'Amplitude': 4.078548524665942e-22, 'EclipticLatitude': 0.8719598285541009, 'EclipticLongitude': 0.48610736535582433, 'Frequency': 0.003995221114188715, 'FrequencyDerivative': 1.0624083369431835e-16, 'Inclination': 1.023230688156412, 'InitialPhase': 2.332809482837094, 'Polarization': 2.6643510855317127}, {'Amplitude': 8.599847426548067e-23, 'EclipticLatitude': -1.0083972289641516, 'EclipticLongitude': -3.0428977235295163, 'Frequency': 0.003994685927983612, 'FrequencyDerivative': 1.8332697826643778e-17, 'Inclination': 1.9139534641852995, 'InitialPhase': 5.114695209956535, 'Polarization': 2.605693433561302}, {'Amplitude': 4.6256730400469185e-23, 'EclipticLatitude': 1.1722276880402125, 'EclipticLongitude': -2.871446178228877, 'Frequency': 0.003994619105280041, 'FrequencyDerivative': 1.086012039379694e-16, 'Inclination': 0.9479641813421728, 'InitialPhase': 4.8371517203795085, 'Polarization': 1.7036543858894597}, {'Amplitude': 1.643553600616391e-23, 'EclipticLatitude': 0.16199423244320238, 'EclipticLongitude': 1.4165120905649182, 'Frequency': 0.003994969780666866, 'FrequencyDerivative': 2.5348637389615347e-17, 'Inclination': 2.4954310757426548, 'InitialPhase': 3.701982859408883, 'Polarization': 0.5289520237169353}]
# slice and dice without outsiders
# found_sources = [{'Amplitude': 3.954437185335676e-22, 'EclipticLatitude': 0.8705283227819178, 'EclipticLongitude': 0.48686210801377316, 'Frequency': 0.003995221067702464, 'FrequencyDerivative': 1.0889344296696101e-16, 'Inclination': 1.0044969922292915, 'InitialPhase': 2.2984551395638126, 'Polarization': 2.651756350361183}, {'Amplitude': 1.2539437693613833e-22, 'EclipticLatitude': -1.1856213219729361, 'EclipticLongitude': -2.683928479256104, 'Frequency': 0.00399462026922502, 'FrequencyDerivative': 8.769232289332856e-17, 'Inclination': 1.9068676016928692, 'InitialPhase': 2.4797972297671795, 'Polarization': 2.524582279457128}, {'Amplitude': 2.4065471042304014e-23, 'EclipticLatitude': 0.13558691114892432, 'EclipticLongitude': 1.5773316510537079, 'Frequency': 0.003995000156412037, 'FrequencyDerivative': 8.54693032650248e-17, 'Inclination': 2.0697606285415384, 'InitialPhase': 2.410026294141748, 'Polarization': 0.5353316544631529}, {'Amplitude': 2.1352291472674563e-23, 'EclipticLatitude': -0.02247515378134304, 'EclipticLongitude': -1.4648490465114308, 'Frequency': 0.003995378272058502, 'FrequencyDerivative': 1.1921575041379484e-16, 'Inclination': 0.9409699470043055, 'InitialPhase': 5.27492090265447, 'Polarization': 0.5916292927580608}, {'Amplitude': 1.3838203303205143e-23, 'EclipticLatitude': -0.02115513595420595, 'EclipticLongitude': -1.601928897252455, 'Frequency': 0.00399515473218528, 'FrequencyDerivative': 1.0156155118046183e-20, 'Inclination': 2.0216072557759834, 'InitialPhase': 1.1195504942202392, 'Polarization': 1.2139853659653514}, {'Amplitude': 1.9662728051288063e-23, 'EclipticLatitude': 0.32244724528551977, 'EclipticLongitude': -3.141592653589793, 'Frequency': 0.0039949685759942165, 'FrequencyDerivative': 9.570885840902786e-14, 'Inclination': 1.420762074316341, 'InitialPhase': 0.7280971343638819, 'Polarization': 2.900984054674992}, {'Amplitude': 2.450403228241954e-23, 'EclipticLatitude': 1.4600819824490758, 'EclipticLongitude': -1.4867428828144233, 'Frequency': 0.003994018890486554, 'FrequencyDerivative': 7.27066897106673e-17, 'Inclination': 2.068581292366034, 'InitialPhase': 5.335276341634204, 'Polarization': 1.7887059369037583}, {'Amplitude': 3.16655881881124e-22, 'EclipticLatitude': -0.36440785933466047, 'EclipticLongitude': 1.4778115680497415, 'Frequency': 0.0039959999999999996, 'FrequencyDerivative': 3.691848331425048e-14, 'Inclination': 1.2254897063352792, 'InitialPhase': 4.18392142231699, 'Polarization': 0.875022662275395}, {'Amplitude': 1.3062754628688594e-23, 'EclipticLatitude': -1.2473959799399659, 'EclipticLongitude': 2.7715207737627443, 'Frequency': 0.003994093500845711, 'FrequencyDerivative': 3.917743059398532e-20, 'Inclination': 1.3965202196689626, 'InitialPhase': 3.4356454300704917, 'Polarization': 2.674370951557413}, {'Amplitude': 1.5754587435594722e-23, 'EclipticLatitude': -0.6530606376730765, 'EclipticLongitude': -2.052637213940567, 'Frequency': 0.003994058956397211, 'FrequencyDerivative': 1.0097581829635257e-20, 'Inclination': 1.6418600927298175, 'InitialPhase': 3.263520624039774, 'Polarization': 1.1342667679686258}]
# found_sources = [{'Amplitude': 3.9689132736574507e-22, 'EclipticLatitude': 0.8707505533420358, 'EclipticLongitude': 0.48662742593757335, 'Frequency': 0.003995221168064815, 'FrequencyDerivative': 1.0495922668329137e-16, 'Inclination': 1.0030685433927617, 'InitialPhase': 2.312257571052444, 'Polarization': 2.655389100178587}, {'Amplitude': 1.2582248855092968e-22, 'EclipticLatitude': -1.1861909664486303, 'EclipticLongitude': -2.682255219986959, 'Frequency': 0.0039946198558335555, 'FrequencyDerivative': 9.899011309497499e-17, 'Inclination': 1.9034799708471297, 'InitialPhase': 2.438578058105272, 'Polarization': 2.526411018962479}, {'Amplitude': 2.890363929791713e-23, 'EclipticLatitude': 0.16361092519076423, 'EclipticLongitude': 1.5706937501136518, 'Frequency': 0.0039949980260170475, 'FrequencyDerivative': 1.3484819876946333e-16, 'Inclination': 1.9355321177427731, 'InitialPhase': 2.2076817294694266, 'Polarization': 0.505660865986898}, {'Amplitude': 1.2947038507777098e-23, 'EclipticLatitude': 1.1416157111708225, 'EclipticLongitude': -1.0097100690135445, 'Frequency': 0.003995130547135328, 'FrequencyDerivative': 3.6748793067569695e-14, 'Inclination': 1.5417163414944728, 'InitialPhase': 2.809507415256849, 'Polarization': 2.13508409583371}, {'Amplitude': 2.5045671136878724e-23, 'EclipticLatitude': -0.05210562118966992, 'EclipticLongitude': -1.459162174180093, 'Frequency': 0.003995380435174599, 'FrequencyDerivative': 7.297345996825393e-17, 'Inclination': 1.0847271487749357, 'InitialPhase': 5.654535579595753, 'Polarization': 0.6403848178881576}, {'Amplitude': 2.450403228241954e-23, 'EclipticLatitude': 1.4600819824490758, 'EclipticLongitude': -1.4867428828144233, 'Frequency': 0.003994018890486554, 'FrequencyDerivative': 7.27066897106673e-17, 'Inclination': 2.068581292366034, 'InitialPhase': 5.335276341634204, 'Polarization': 1.7887059369037583}, {'Amplitude': 8.432124140961827e-23, 'EclipticLatitude': -0.18573110882718108, 'EclipticLongitude': 1.060882704635124, 'Frequency': 0.003995717179439931, 'FrequencyDerivative': 9.61652656496001e-14, 'Inclination': 1.3580547994571575, 'InitialPhase': 6.14128922422039, 'Polarization': 2.905461713580476},{'Amplitude': 1.6069668047950477e-23, 'EclipticLatitude': 0.695553972366851, 'EclipticLongitude': -2.2813168281189378, 'Frequency': 0.003995209132034196, 'FrequencyDerivative': 1e-13, 'Inclination': 1.5806174383653766, 'InitialPhase': 5.691666472290874, 'Polarization': 0.03247000597996169}]
# found_sources = [{'Amplitude': 3.960763861925662e-22, 'EclipticLatitude': 0.87082501612068, 'EclipticLongitude': 0.4865306596743239, 'Frequency': 0.003995221147594301, 'FrequencyDerivative': 1.0550882754170559e-16, 'Inclination': 1.0014118629031084, 'InitialPhase': 2.3101731093466102, 'Polarization': 2.6552029377086903}, {'Amplitude': 1.2557721986337792e-22, 'EclipticLatitude': -1.186159654973634, 'EclipticLongitude': -2.6819160691010735, 'Frequency': 0.00399461981336803, 'FrequencyDerivative': 1.0037720328610126e-16, 'Inclination': 1.9047524347599567, 'InitialPhase': 2.434749861107631, 'Polarization': 2.526440577068304}, {'Amplitude': 2.3049589958174292e-23, 'EclipticLatitude': -0.11985084854915119, 'EclipticLongitude': -1.465346904354099, 'Frequency': 0.00399537915885506, 'FrequencyDerivative': 1.0682357470038055e-16, 'Inclination': 1.0107536282956504, 'InitialPhase': 2.2241909841352405, 'Polarization': 2.1471839347098385}, {'Amplitude': 2.7831892808451777e-23, 'EclipticLatitude': 0.16333318203734196, 'EclipticLongitude': 1.5720208015057582, 'Frequency': 0.003994998432218014, 'FrequencyDerivative': 1.2521477898223552e-16, 'Inclination': 1.9435027987144409, 'InitialPhase': 5.402174068645186, 'Polarization': 2.077720898308443}, {'Amplitude': 1.1025488670938274e-23, 'EclipticLatitude': 1.1673600818167074, 'EclipticLongitude': 0.057088882641640915, 'Frequency': 0.003994055648329161, 'FrequencyDerivative': 5.811149163003285e-19, 'Inclination': 3.141592653589793, 'InitialPhase': 4.011764949371346, 'Polarization': 0.6149910828045372}, {'Amplitude': 1.8296715009558222e-23, 'EclipticLatitude': 1.4763765365636556, 'EclipticLongitude': -1.4936674770977738, 'Frequency': 0.0039940143831784704, 'FrequencyDerivative': 2.1978625263064225e-16, 'Inclination': 1.829917463042957, 'InitialPhase': 4.894312836190836, 'Polarization': 1.891773422360853}, {'Amplitude': 1.0951489778321648e-23, 'EclipticLatitude': -1.1620567772915402, 'EclipticLongitude': 2.208658684486738, 'Frequency': 0.003994365782934417, 'FrequencyDerivative': 9.398702356575727e-14, 'Inclination': 0.8912620469775321, 'InitialPhase': 1.1248245768258496, 'Polarization': 0.33861953048977755}, {'Amplitude': 1.479345603740554e-23, 'EclipticLatitude': -0.8110628214488717, 'EclipticLongitude': -1.9835724175601763, 'Frequency': 0.00399405689262681, 'FrequencyDerivative': 1.135000009580916e-20, 'Inclination': 1.736437404437111, 'InitialPhase': 5.777394508225548, 'Polarization': 2.0847997431184604}]
# -20 to -14.5 = df/dt
# found_sources = [{'Amplitude': 3.9593411057982525e-22, 'EclipticLatitude': 0.8706035660730967, 'EclipticLongitude': 0.4867711389943583, 'Frequency': 0.003995221069402352, 'FrequencyDerivative': 1.0866829218214829e-16, 'Inclination': 1.0056149026804564, 'InitialPhase': 2.3014245768388166, 'Polarization': 2.653261596089537}, {'Amplitude': 1.24956315244586e-22, 'EclipticLatitude': -1.18622632093688, 'EclipticLongitude': -2.6828225533029064, 'Frequency': 0.003994619649028894, 'FrequencyDerivative': 1.0517947828832982e-16, 'Inclination': 1.9095721945209743, 'InitialPhase': 2.418942202423396, 'Polarization': 2.5244526060642483}, {'Amplitude': 1.3676610853313325e-23, 'EclipticLatitude': -0.018521276207872868, 'EclipticLongitude': -1.466764263927183, 'Frequency': 0.003995379125385088, 'FrequencyDerivative': 9.597030021659842e-17, 'Inclination': 0.010723818509195321, 'InitialPhase': 3.9165627042143085, 'Polarization': 3.0176937135707163}, {'Amplitude': 2.370847512918587e-23, 'EclipticLatitude': 0.13658490306299995, 'EclipticLongitude': 1.573869274356241, 'Frequency': 0.003994998419266302, 'FrequencyDerivative': 1.3198396363145704e-16, 'Inclination': 2.0931186332759637, 'InitialPhase': 2.234863721938927, 'Polarization': 0.5345964955705121}, {'Amplitude': 1.3279642690930703e-23, 'EclipticLatitude': -0.00171035385637249, 'EclipticLongitude': -1.6052081984746336, 'Frequency': 0.0039951550071572264, 'FrequencyDerivative': 1.2232904563154208e-20, 'Inclination': 2.0619567363134363, 'InitialPhase': 1.222346190648877, 'Polarization': 1.1977681031993352}, {'Amplitude': 1.4133785710858455e-23, 'EclipticLatitude': 1.2913444244053214, 'EclipticLongitude': 0.1445509942836649, 'Frequency': 0.003994082456377896, 'FrequencyDerivative': 1.4818904094633514e-16, 'Inclination': 2.31879974688655, 'InitialPhase': 2.615171521065686, 'Polarization': 0.3295668252765084}, {'Amplitude': 1.4234622912121888e-23, 'EclipticLatitude': 1.3829226364408154, 'EclipticLongitude': 0.9125530827150881, 'Frequency': 0.0039940814643270785, 'FrequencyDerivative': 1.988304873137983e-16, 'Inclination': 1.9700328063080959, 'InitialPhase': 1.876886550628801, 'Polarization': 2.3768487770660007}, {'Amplitude': 1.87191708625015e-23, 'EclipticLatitude': -1.122168161334067, 'EclipticLongitude': -1.4794955904688387, 'Frequency': 0.003994009296128951, 'FrequencyDerivative': 4.3798592478926975e-16, 'Inclination': 1.5513903560363402, 'InitialPhase': 2.859349925688812, 'Polarization': 1.310845535800841},{'Amplitude': 1.0610583492839024e-23, 'EclipticLatitude': -0.42836394662975397, 'EclipticLongitude': -0.6829651112254829, 'Frequency': 0.003994524925022837, 'FrequencyDerivative': 6.877624346349927e-16, 'Inclination': 1.6217388536201427, 'InitialPhase': 0.5344281759141147, 'Polarization': 0.25867501796914344}]
# slice and dice with outsiders
# found_sources = [{'Amplitude': 3.96289718306154e-22, 'EclipticLatitude': 0.8708928120472197, 'EclipticLongitude': 0.48665946453175035, 'Frequency': 0.003995221254538698, 'FrequencyDerivative': 1.0269904779836001e-16, 'Inclination': 1.001973215647365, 'InitialPhase': 2.319958680868902, 'Polarization': 2.654537339387363}, {'Amplitude': 1.2619710443797457e-22, 'EclipticLatitude': -1.1856184090524093, 'EclipticLongitude': -2.683104698487379, 'Frequency': 0.0039946203259605215, 'FrequencyDerivative': 8.553014677510986e-17, 'Inclination': 1.9015560309268713, 'InitialPhase': 2.483133469961604, 'Polarization': 2.5270559934483035}, {'Amplitude': 2.864376702953742e-23, 'EclipticLatitude': 0.1630042875267411, 'EclipticLongitude': 1.5724345857890278, 'Frequency': 0.003994998564933117, 'FrequencyDerivative': 1.1941881242605984e-16, 'Inclination': 1.9322708810269942, 'InitialPhase': 2.250537209576398, 'Polarization': 0.5092726324336495}, {'Amplitude': 2.3687747804856843e-23, 'EclipticLatitude': 1.459302793276365, 'EclipticLongitude': -1.5544211392843303, 'Frequency': 0.003994021200055487, 'FrequencyDerivative': 4.4173057035466856e-18, 'Inclination': 1.9986779156687422, 'InitialPhase': 5.8398053449499185, 'Polarization': 1.524542363540894}, {'Amplitude': 1.0997731423376216e-23, 'EclipticLatitude': 1.1094981126424124, 'EclipticLongitude': -0.9623293571222385, 'Frequency': 0.003995160830733279, 'FrequencyDerivative': 3.482976053683434e-14, 'Inclination': 1.697683986561341, 'InitialPhase': 3.8252099672398296, 'Polarization': 2.142227458770483}, {'Amplitude': 5.0625386251849795e-22, 'EclipticLatitude': -0.23503441908984832, 'EclipticLongitude': 1.6361307415835924, 'Frequency': 0.003995994570650909, 'FrequencyDerivative': 3.5099578525260406e-14, 'Inclination': 1.3398146930752173, 'InitialPhase': 1.8213670084595, 'Polarization': 0.47617840650086124}, {'Amplitude': 2.4539601480458732e-23, 'EclipticLatitude': -0.05780570833613541, 'EclipticLongitude': -1.4596419484178105, 'Frequency': 0.003995378959176098, 'FrequencyDerivative': 1.0846895889086299e-16, 'Inclination': 1.0768336533465692, 'InitialPhase': 2.3677586761599523, 'Polarization': 2.2248720169531526}, {'Amplitude': 1.6597338193344268e-23, 'EclipticLatitude': -1.2835634078448852, 'EclipticLongitude': -2.956056071255409, 'Frequency': 0.003994054152190724, 'FrequencyDerivative': 2.0361138761171463e-16, 'Inclination': 1.3622526771807204, 'InitialPhase': 5.649426320028763, 'Polarization': 2.204871630524406}, {'Amplitude': 4.864426062085493e-23, 'EclipticLatitude': -0.9620890981246751, 'EclipticLongitude': -2.4125645239308144, 'Frequency': 0.003995102747038606, 'FrequencyDerivative': 9.999999999826903e-14, 'Inclination': 1.794724968706271, 'InitialPhase': 0.6474549451054277, 'Polarization': 2.443162184273923}]
# found_sources = found_sources[:-1]

# plt.figure(figsize=fig_size)
# plt.plot(range(len(energies)),-energies, marker = 'o', markersize = 3, label='search')
# plt.hlines(search1.loglikelihood([search1.pGB]),0,len(energies),'k', label='true')
# plt.xlabel('Iteration')
# plt.ylabel('Log-likelihood')
# plt.legend()
# plt.show()


# tdi_fs_long_subtracted = deepcopy(tdi_fs_long)
# tdi_fs_subtracted = deepcopy(tdi_fs)
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
# for pGBadding in [pGBadded7]:#, pGBadded2, pGBadded3, pGBadded4]:
#     Xs_added, Ys_added, Zs_added = GB.get_fd_tdixyz(template=pGBadding, oversample=4, simulator="synthlisa")
#     source_added = dict({"X": Xs_added, "Y": Ys_added, "Z": Zs_added})
#     index_low = np.searchsorted(tdi_fs["X"].f, Xs_added.f[0])
#     index_high = index_low+len(Xs_added)
#     # tdi_fs['X'] = tdi_fs['X'] #+ Xs_added
#     for k in ["X", "Y", "Z"]:
#         tdi_fs[k].data[index_low:index_high] = tdi_fs[k].data[index_low:index_high] + source_added[k].data
# tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft(dt=dt)) for k, n in [["X", 1], ["Y", 2], ["Z", 3]]]))

# create two sets of found sources. found_sources_in with signals inside the boundary and founce_sources_out with outside sources
found_sources_in = []
found_sources_out = []
for i in range(len(found_sources)):
    if found_sources[i]['Frequency'] > lower_frequency and found_sources[i]['Frequency'] < upper_frequency:
        found_sources_in.append(found_sources[i])
    else:
        found_sources_out.append(found_sources[i])
for i in range(len(found_sources_out)):
    Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_out[i], oversample=4, simulator="synthlisa")
    source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
    index_low = np.searchsorted(tdi_fs["X"].f, Xs_subtracted.f[0])
    index_high = index_low+len(Xs_subtracted)
    for k in ["X", "Y", "Z"]:
        tdi_fs[k].data[index_low:index_high] = tdi_fs[k].data[index_low:index_high] - source_subtracted[k].data
    tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft()) for k in ["X", "Y", "Z"]]))

search1 = Search(tdi_fs,Tobs)
search1.plot()

total_boundaries = deepcopy(search1.boundaries)
amplitudes = []
for i in range(len(found_sources_in)):
    amplitudes.append(found_sources_in[i]['Amplitude'])
total_boundaries['Amplitude'] = [np.min(amplitudes),np.max(amplitudes)]
amplitudes_length = np.log10(total_boundaries['Amplitude'][1]) - np.log10(total_boundaries['Amplitude'][0])
total_boundaries['Amplitude'] = [np.log10(total_boundaries['Amplitude'][0]) - amplitudes_length/5, np.log10(total_boundaries['Amplitude'][1]) + amplitudes_length/5]


number_of_signals_optimize = len(found_sources_in)
start = time.time()
maxpGBsearch, pGB = search1.optimize([found_sources_in], boundaries= total_boundaries)
print(time.time()-start)
# start = time.time()
# maxpGBsearch, pGB = search1.differential_evolution_search(search1.boundaries['Frequency'], initial_guess=found_sources[:number_of_signals])
# print(time.time()-start)

# maxpGBsearch = [maxpGBsearch]
# found_sources = maxpGBsearch[0]

# subtraction = maxpGBsearch[0][1]
# subtraction = found_sources[0]
# Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=subtraction, oversample=4, simulator="synthlisa")
# source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
# index_low = np.searchsorted(tdi_fs["X"].f, Xs_subtracted.f[0])
# index_high = index_low+len(Xs_subtracted)
# for k in ["X", "Y", "Z"]:
#     tdi_fs[k].data[index_low:index_high] = tdi_fs[k].data[index_low:index_high] - source_subtracted[k].data
# tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft()) for k in ["X", "Y", "Z"]]))


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
index_low = np.searchsorted(p.get('Frequency')[indexes], lower_frequency-padding)
index_high = np.searchsorted(p.get('Frequency')[indexes], upper_frequency+padding)
range_index = np.logical_and(tdi_fs.f > lower_frequency-padding, tdi_fs.f < upper_frequency+padding)

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

fig = plt.figure(figsize=fig_size)
ax1 = plt.subplot(111)
for i in range(len(found_sources)):
    if i == 0:
        ax1.semilogy(found_sources[i]['Frequency']*10**3,found_sources[i]['Amplitude'], color = colors[i], marker = 'o', markersize = 5, label= 'initial guess')
    else:
        ax1.semilogy(found_sources[i]['Frequency']*10**3,found_sources[i]['Amplitude'], color = colors[i], marker = 'o', markersize = 5)
for i in range(len(maxpGBsearch)):
    if i == 0:
        ax1.semilogy(maxpGBsearch[i]['Frequency']*10**3,maxpGBsearch[i]['Amplitude'], color = colors[i], marker = '+', markersize = 8, label= 'global fit')
    else:
        ax1.semilogy(maxpGBsearch[i]['Frequency']*10**3,maxpGBsearch[i]['Amplitude'], color = colors[i], marker = '+', markersize = 8)
for i in range(len(pGB_injected)):   
    if i == 0: 
        ax1.semilogy(pGB_injected[i]['Frequency']*10**3,pGB_injected[i]['Amplitude'],color='grey', marker = 'o', markersize = 8, zorder=1, label= 'true')
    else:
        ax1.semilogy(pGB_injected[i]['Frequency']*10**3,pGB_injected[i]['Amplitude'],color='grey', marker = 'o', markersize = 8, zorder=1)
ax1.axvline(10**3*(search1.boundaries['Frequency'][0]+padding), color= 'red', label='Boundaries')
ax1.axvline(10**3*(search1.boundaries['Frequency'][1]-padding), color= 'red')
plt.xlim((lower_frequency-padding)*10**3, (upper_frequency+padding)*10**3)
ax1.xaxis.set_major_locator(plt.MaxNLocator(4))
plt.xlabel('f [mHz]')
plt.ylabel('Amplitude')
# plt.legend(loc='upper right')
fig.savefig('/home/stefan/LDC/LDC/pictures/global fit amplitude frequency 3,9945-3,9955 mHz2.png',dpi=300,bbox_inches='tight')
plt.show()

plt.figure(figsize=fig_size)
ax1 = plt.subplot(111)
for i in range(len(found_sources)):
    Xs, Ys, Zs = GB.get_fd_tdixyz(template= found_sources[i], oversample=4, simulator="synthlisa")
    ax1.plot(found_sources[i]['EclipticLatitude'],found_sources[i]['EclipticLongitude'], marker = 'o', markersize = 3)
for i in range(len(pGB_injected)):  
    ax1.plot(pGB_injected[i]['EclipticLatitude'],pGB_injected[i]['EclipticLongitude'],color='grey', marker = 'o', markersize = 5, zorder=1)
plt.show()

# start = time.time()
# maxpGB, pGB = search1.optimize([found_sources])
# print(time.time()-start)

# found_sources = maxpGB
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
search1 = Search(tdi_fs,Tobs)
pGB = search1.pGB
# maxpGB, pGB =  search1.optimize([[maxpGB]])
# maxpGB = maxpGB[0]
# search1.plot(maxpGB)
loglikelihood = search1.loglikelihood
print('log-lieklihood:', loglikelihood([pGB]))
SNR = search1.SNR
boundaries = search1.boundaries
best_value = loglikelihood([maxpGB])
boundaries_reduced = Reduce_boundaries(maxpGB, boundaries,ratio=0.1)

search1.plot(pGBadded=maxpGB, added_label='MLE')

fig = plt.figure(figsize=fig_size)
ax1 = plt.subplot(111)
ax1.semilogy(tdi_fs.f[range_index]*10**3,np.abs(tdi_fs['X'][range_index])**2,'k',zorder= 1, linewidth = 2, label = 'Data')
# ax1.semilogy(tdi_fs_long_subtracted.f[range_index],np.abs(tdi_fs_long_subtracted['X'][range_index])**2,'b',zorder= 5)
for i in range(len(pGB_injected)):
    Xs, Ys, Zs = GB.get_fd_tdixyz(template= pGB_injected[i], oversample=4, simulator="synthlisa")
    a,Xs = xr.align(search1.dataX, Xs, join='left',fill_value=0)
    ax1.semilogy(Xs.f*10**3,np.abs(Xs)**2,label= str(np.round(pGB_injected[i]['Frequency'],0)), color='grey', linewidth = 3)
for i in range(len(maxpGBsearch)):
    Xs, Ys, Zs = GB.get_fd_tdixyz(template= maxpGBsearch[i], oversample=4, simulator="synthlisa")
    ax1.semilogy(Xs.f*10**3,np.abs(Xs)**2,'--', color= colors[i], linewidth = 1.8)
# Xs, Ys, Zs = GB.get_fd_tdixyz(template= pGB, oversample=4, simulator="synthlisa")
# a,Xs = xr.align(search1.dataX, Xs, join='left',fill_value=0)
# ax1.semilogy(Xs.f,np.abs(tdi_fs['X'][range_index]-Xs)**2,label= 'residual')
plt.xlim((lower_frequency-padding)*10**3, (upper_frequency+padding)*10**3)
# ax1.axhline(y=0)
ax1.xaxis.set_major_locator(plt.MaxNLocator(4))
# plt.legend()
plt.xlabel('f [mHz]')
plt.ylabel('|X|')
# fig.savefig('/home/stefan/LDC/LDC/pictures/global fit strain 3,9945-3,9955 mHz2.png',dpi=300,bbox_inches='tight')
plt.show()

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

def confidence_ellipse(x,y, ax, n_std=3.0, facecolor='none', **kwargs):
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
    cov = np.cov(x, y)
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
# index_parameter1 = parameters.index('Inclination')
# index_parameter2 = parameters.index('Amplitude')
index_parameter1 = parameters.index('EclipticLatitude')
index_parameter2 = parameters.index('EclipticLongitude')
cov2d[0,0] = covariance_matrix[index_parameter1,index_parameter1]
cov2d[0,1] = covariance_matrix[index_parameter1,index_parameter2]
cov2d[1,0] = covariance_matrix[index_parameter2,index_parameter1]
cov2d[1,1] = covariance_matrix[index_parameter2,index_parameter2]
# cov2d = scipy.linalg.inv([[2,1],[1,1]])
x = maxpGB01[parameters[index_parameter1]]
y = maxpGB01[parameters[index_parameter2]]
samples = np.random.multivariate_normal([x,y],cov2d,10000)
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
confidence_ellipse(samples_x, samples_y, ax,n_std=5, edgecolor='red')
ax.scatter(samples_x,samples_y)
# ax.quiver(x,y,v[0,0], v[0,1],color='r',scale = 1/lambda_[0], scale_units='xy')
# ax.quiver(x,y,v[1,0], v[1,1],scale = 1/lambda_[1], scale_units='xy')
# ax.scatter(x+v[0,0]*lambda_[0]*5, y+v[1,0]*lambda_[0]*5)
# ax.scatter(x+v[0,1]*lambda_[1]*5, y+v[1,1]*lambda_[1]*5)
ax.scatter(x, y, color='g',label='Found parameters')
plt.xlabel(parameters[index_parameter1])
plt.ylabel(parameters[index_parameter2])
transf = v @ np.diag(lambda_)
scalematrix = np.max(np.abs(transf), axis=1)
scalematrix = np.sqrt(np.diag(cov2d))
xlim = scalematrix[0]
ylim = scalematrix[1]
# xlim = np.max(v[:,0]*lambda_[0])
# ylim = np.max(v[:,1]*lambda_[1])
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
scalematrix = np.sqrt(np.diag(covariance_matrix))

# for parameter in parameters:
#     scalematrix[parameters.index(parameter)] = np.sqrt(covariance_matrix)[parameters.index(parameter)][parameters.index(parameter)]

maxpGB01_low = deepcopy(maxpGB01)
maxpGB01_high = deepcopy(maxpGB01)
boundaries_reduced_fisher = {}
sigma_multiplyer = 2
for parameter in parameters:
    maxpGB01_low[parameter] = maxpGB01[parameter] - scalematrix[parameters.index(parameter)] * sigma_multiplyer 
    maxpGB01_high[parameter] = maxpGB01[parameter] + scalematrix[parameters.index(parameter)] * sigma_multiplyer 
    if parameter in [ 'InitialPhase', 'Polarization']:
        maxpGB01_low[parameter] = maxpGB01[parameter] - scalematrix[parameters.index(parameter)]  * 0.001
        maxpGB01_high[parameter] = maxpGB01[parameter] + scalematrix[parameters.index(parameter)] * 0.001
    if parameter in [ 'Frequency']:
        maxpGB01_low[parameter] = maxpGB01[parameter] - scalematrix[parameters.index(parameter)]  * 1
        maxpGB01_high[parameter] = maxpGB01[parameter] + scalematrix[parameters.index(parameter)] * 1
    if parameter == 'FrequencyDerivative':
        print('scale',scalematrix[parameters.index(parameter)])
        if scalematrix[parameters.index(parameter)] > 0.07:
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
if boundaries_reduced['FrequencyDerivative'][1] < split_fd+0.5:
    boundaries_reduced['FrequencyDerivative'][0] = -18
    boundaries_reduced['FrequencyDerivative'][1] = -16
elif boundaries_reduced['FrequencyDerivative'][0] < split_fd+0.5:
    boundaries_reduced['FrequencyDerivative'][0] = -17.5
# boundaries_reduced['FrequencyDerivative'][0] = -18
# if boundaries_reduced['FrequencyDerivative'][0] < split_fd and boundaries_reduced['FrequencyDerivative'][1] > split_fd:
#     boundaries_reduced1 = deepcopy(boundaries_reduced)
#     boundaries_reduced1['FrequencyDerivative'] = [boundaries_reduced['FrequencyDerivative'][0], split_fd]
#     boundaries_reduced2 = deepcopy(boundaries_reduced)
#     boundaries_reduced2['FrequencyDerivative'] = [split_fd, boundaries_reduced['FrequencyDerivative'][1]]
# else:
halfway = boundaries_reduced['FrequencyDerivative'][0]+(boundaries_reduced['FrequencyDerivative'][1]-boundaries_reduced['FrequencyDerivative'][0])*0.99
boundaries_reduced1 = deepcopy(boundaries_reduced)
boundaries_reduced1['FrequencyDerivative'] = [boundaries_reduced['FrequencyDerivative'][0], boundaries_reduced['FrequencyDerivative'][1]]
boundaries_reduced2 = deepcopy(boundaries_reduced)
boundaries_reduced2['FrequencyDerivative'] = [halfway, boundaries_reduced['FrequencyDerivative'][1]]

print(boundaries_reduced)
print('booundary1:', boundaries_reduced1)
print('booundary2:', boundaries_reduced2)
rmse = 2
train_size = 0
test_size = 500
added_trainig_size = 1000
j = 0
samples = np.random.rand(6000,8)
while rmse > 0.5 and j < 5:
    j += 1
    train_size += added_trainig_size
    # train_size = 3000
    if j == 1:
        resolution = train_size + test_size
    else:
        resolution = added_trainig_size
    # start = time.time()
    # samples = sampler(resolution,parameters,maxpGB, boundaries_reduced1, uniform=False, twoD=False)
    # print('sample time of', resolution, 'samples ',time.time() - start)
    start = time.time()
    samples_likelihood = np.zeros(resolution)
    samples_likelihood2 = np.zeros(resolution)
    for i in range(resolution):
        samples_p = scaletooriginal(samples[i+j*added_trainig_size], boundaries_reduced1)
        # samples_likelihood[i] = search1.loglikelihoodsdf([samples_p])
        samples_likelihood[i] = search1.loglikelihood([samples_p])
    print('sample time of', resolution, 'samples ',time.time() - start)

    samples_flat = np.zeros((resolution  , len(parameters)))
    i = 0
    boundary_ratio = 1# (boundaries_reduced1['FrequencyDerivative'][1]-boundaries_reduced1['FrequencyDerivative'][0])/(boundaries_reduced2['FrequencyDerivative'][1]-boundaries_reduced1['FrequencyDerivative'][0])
    for parameter in parameters:
        if parameter == 'FrequencyDerivative':
            samples_flat[:, i] = samples[j*added_trainig_size:j*added_trainig_size+resolution,parameters.index(parameter)]*boundary_ratio
        else:
            samples_flat[:, i] = samples[j*added_trainig_size:j*added_trainig_size+resolution,parameters.index(parameter)]
        i += 1
    if j == 1:
        train_y = samples_likelihood[test_size:]
        # train_y2 = samples_likelihood2[test_size:]
        test_y = samples_likelihood[:test_size]
        # test_y2 = samples_likelihood2[:test_size]
        train_x = samples_flat[test_size:]
        test_x = samples_flat[:test_size]
    else:
        train_y = np.append(train_y, samples_likelihood)
        # train_y2 = np.append(train_y2, samples_likelihood2)
        # test_y = np.append(test_y, samples_likelihood[train_size:])
        train_x = np.append(train_x, samples_flat,axis=0)
        # test_x = np.append(test_x, samples_flat[train_size:],axis=0)

    nu = np.mean(train_y)
    sigma = np.std(train_y)
    train_y_normalized = (train_y - nu) / sigma
    # nu2 = np.mean(train_y2)
    # sigma2 = np.std(train_y2)
    # train_y_normalized2 = (train_y2 - nu2) / sigma2
    kernel = RBF(length_scale=[1,2,5,1,1,1,1,1],length_scale_bounds=[(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,100),(0.1,20)])
    start = time.time()
    gpr = GaussianProcessRegressor(kernel=kernel, random_state=0).fit(train_x, train_y_normalized)
    # gpr2 = GaussianProcessRegressor(kernel=kernel, random_state=0).fit(train_x, train_y_normalized2)
    # gpr = traingpmodelsk(train_x, train_y, kernel)
    print('train',time.time() - start)
    start = time.time()
    observed_pred_sk = gpr.predict(test_x)
    # observed_pred_sk2 = gpr2.predict(test_x)
    print('eval time of ', test_size, 'samples: ',time.time() - start)
    observed_pred_sk_scaled = observed_pred_sk*sigma +nu
    # observed_pred_sk_scaled2 = observed_pred_sk2*sigma2 +nu2
    rmse = np.sqrt(mean_squared_error(test_y,observed_pred_sk_scaled))
    # rmse2 = np.sqrt(mean_squared_error(test_y2,observed_pred_sk_scaled2))
    print("RMSE ",rmse)
    print('training size', len(train_y))

# maxpGBmod = deepcopy(maxpGB)
# maxpGBmod['FrequencyDerivative'] = 10**(np.random.rand()*(boundaries_reduced2['FrequencyDerivative'][1]-boundaries_reduced2['FrequencyDerivative'][0])+boundaries_reduced2['FrequencyDerivative'][0])
# start = time.time()
# parameter = "Frequency"
# samples2 = sampler(resolution,parameters,maxpGBmod, boundaries_reduced2, uniform=False, twoD=False)
# print('sample time of', resolution, 'samples ',time.time() - start)
# train_x = np.zeros((resolution, len(parameters)))
# i = 0
# for parameter in parameters:
#     if parameter == 'FrequencyDerivative':
#         train_x[:, i] = samples2[parameter]*(1-boundary_ratio)+boundary_ratio
#     else:
#         train_x[:, i] = samples2[parameter]
#     i += 1
# train_y = samples2["Likelihood"][:train_size]
# test_y = samples2["Likelihood"][-500:]
# test_x = train_x[-500:]
# train_x = train_x[:train_size]
# nu2 = np.mean(train_y)
# sigma2 = np.std(train_y)
# train_y_scaled = (train_y - nu2) / sigma2
# test_y_scaled = (test_y - nu2) / sigma2
# kernel = RBF(length_scale=[1,2,5,1,1,1,1,1],length_scale_bounds=[(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,30),(0.1,30)])
# start = time.time()
# # gpr2 = traingpmodelsk(train_x, train_y_scaled, kernel)
# gpr2 = GaussianProcessRegressor(kernel=kernel, random_state=0).fit(train_x, train_y_scaled)
# print('train',time.time() - start)
# start = time.time()
# observed_pred_sk = gpr2.predict(test_x)
# print('eval',time.time() - start)
# observed_pred_sk_scaled2 = observed_pred_sk*sigma2 +nu2
# print("RMSE ",np.sqrt(mean_squared_error(test_y,observed_pred_sk_scaled2)))

# first sampling for calculating prior
resolution = 1*10 ** 6
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
numberinlowerbatch = resolution#next(test_x_m[:,4],boundary_ratio)
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
# def Evaluate2(i):
#     prediction = gpr2.predict(test_x_m2[(i)*partial_length:(i+1)*partial_length])
#     return prediction
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

# start = time.time()
# for n in range(int(numberinupperbatch/partial_length)):
#     observed_pred_sk2[n*partial_length:(n+1)*partial_length] = gpr2.predict(test_x_m2[(n)*partial_length:(n+1)*partial_length])
# try:
#     observed_pred_sk2[int(numberinupperbatch/partial_length)*partial_length:] = gpr2.predict(test_x_m2[int(numberinupperbatch/partial_length)*partial_length:])
# except:
#     pass
# observed_pred_sk2 = np.asarray(observed_pred_sk2)
# observed_pred_sk2 = observed_pred_sk2.reshape(numberinupperbatch)
# observed_pred_mean[-numberinupperbatch:] = observed_pred_sk2*sigma2 +nu2
# print('eval time', time.time()-start)

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
    resolution = 1*10**6
    numPoints = 2**6+1
    if round == 1:
        resolution = 1*10**6
        numPoints = 2**6+1
    if round == 2:
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
    if round != 3:
        start = time.time()
        # test_x_m = np.random.normal(loc= 0.5,scale=0.2,size=(resolution,len(parameters)))
        # std_scale = 0.4
        # test_x_m = scipy.stats.truncnorm.rvs((0-0.5)/std_scale,(1-0.5)/std_scale,loc=0.5,scale=std_scale,size=(resolution,len(parameters)))
        new_index = test_x_m[:,4].argsort()
        test_x_m = test_x_m[new_index]
        probability = probability[new_index]
        numberinlowerbatch = resolution# next(test_x_m[:,4],boundary_ratio)
        numberinupperbatch = len(probability)-numberinlowerbatch
        test_x_m1 = test_x_m #test_x_m[:numberinlowerbatch]
        numberinlowerbatch = len(test_x_m1)
        # test_x_m2 = test_x_m[numberinlowerbatch:]
        print('sample time', time.time()-start)
        partial_length = 1*10**3
        start = time.time()
        observed_pred_mean = np.zeros(len(probability))
        observed_pred_sk = np.zeros(numberinlowerbatch)
        # observed_pred_sk2 = np.zeros(numberinupperbatch)
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

        # !!!!!!! mistake in former versions of splitting the two datasets

        # start = time.time()
        # for n in range(int(numberinupperbatch/partial_length)):
        #     observed_pred_sk2[n*partial_length:(n+1)*partial_length] = gpr2.predict(test_x_m2[(n)*partial_length:(n+1)*partial_length])
        # try:
        #     observed_pred_sk2[int(numberinupperbatch/partial_length)*partial_length:] = gpr2.predict(test_x_m2[int(numberinupperbatch/partial_length)*partial_length:])
        # except:
        #     pass
        # observed_pred_sk2 = np.asarray(observed_pred_sk2)
        # observed_pred_sk2 = observed_pred_sk2.reshape(numberinupperbatch)
        # observed_pred_mean[-numberinupperbatch:] = observed_pred_sk2*sigma2 +nu2
        # print('eval time', time.time()-start)


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
        flatsamples_normalized = np.exp(flatsamplesparameters[:,0]-best_value)#/normalizer
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
# mcmc_samples = test_x_m
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
start = time.time()
df = pd.DataFrame(data=mcmc_samples_rescaled, columns=parameters)
df.to_csv('/home/stefan/Repositories/ldc1_evaluation_data/submission/ETH_LDC_1253/GW'+str(int(np.round(maxpGB['Frequency']*10**8)))+'seed42.csv',index=False)
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

dr = '/home/stefan/Repositories/ldc1_evaluation_data/submission/ETH_LDC_1253/GW125313seed40'
ETH_data3 = {}
dat = np.genfromtxt(dr+".csv", delimiter=',', names=True)
print (np.shape(dat))
datS2 = np.zeros((len(dat),len(parameters)))


datS2[:, 0] = dat['EclipticLongitude']
datS2[:, 1] = dat['EclipticLatitude']
datS2[:, 2] = dat['Frequency']
datS2[:, 3] = np.log10(np.abs(dat['FrequencyDerivative']))
# datS2[:, 4] = np.arccos(np.cos(dat['Inclination']))
datS2[:, 4] = dat['Inclination']
datS2[:, 5] = np.log10(dat['Amplitude'])
datS2[:, 6] = dat['Polarization']
datS2[:, 7] = dat['InitialPhase']

# rng = []
# for i in range(len(lbls)):
#     minrange = min(datS[:,i].min(), datS2[:,i].min())
#     maxrange = max(datS[:,i].max(), datS2[:,i].max())
#     oner = ( minrange, maxrange)
#     rng.append(oner)

fig =  corner.corner(datS[:,:6],  bins=40, hist_kwargs={'density':True, 'lw':3}, plot_datapoints=False, fill_contours=False,  show_titles=True, \
                        color='#348ABD', truths= tr_s[:6], truth_color='k', use_math_test=True, labels= lbls[:6],\
                        levels=[0.9,0.63], title_kwargs={"fontsize": 12})
             
corner.corner(datS2[:,:6],  bins=40, hist_kwargs={'density':True, 'lw':3},fig=fig, plot_datapoints=False, fill_contours=False,  show_titles=True, \
                        color='g', truths= tr_s[:6], truth_color='k', use_math_test=True, labels= lbls[:6],\
                        levels=[0.9,0.63], title_kwargs={"fontsize": 12}, range=rng[:6], label_size= 60)


# Extract the axes
ndim = len(parameters)-2
axes = np.array(fig.axes).reshape((ndim, ndim))
# # Loop over the diagonal
# for i in range(ndim):
#     ax = axes[i, i]
#     ax.axvline(maxvalues[i], color="g", label= 'Found parameters')
# # Loop over the histograms
# for yi in range(ndim):
#     for xi in range(yi):
#         ax = axes[yi, xi]
#         ax.axvline(maxvalues[xi], color="g")
#         ax.axhline(maxvalues[yi], color="g")
#         ax.scatter(maxvalues[xi], maxvalues[yi],s=130, color="g")

# corner.ove(fig, maxvalues[None], marker="s", color="C1")

# maxvalues_previous = np.zeros(len(parameters))
# maxpGBsearch_previous = {'Amplitude': 5.857579179174051e-23, 'EclipticLatitude': -0.08784800634433576, 'EclipticLongitude': 2.102933791009969, 'Frequency': 0.006220280040884535, 'FrequencyDerivative': 7.503323987637861e-16, 'Inclination': 0.5104961907537644, 'InitialPhase': 3.8035280736734642, 'Polarization': 0.07816858489674128}
# i = 0
# for parameter in ['EclipticLongitude','EclipticLatitude','Frequency','FrequencyDerivative','Inclination','Amplitude','InitialPhase','Polarization']:
#     if parameter in ['Amplitude','FrequencyDerivative']:
#         maxvalues_previous[i] = np.log10(maxpGBsearch_previous[parameter])
#     else:
#         maxvalues_previous[i] = maxpGBsearch_previous[parameter]
#     i += 1
# for i in range(ndim):
#     ax = axes[i, i]
#     ax.axvline(maxvalues_previous[i], color="r", label= 'Found parameters')
# # Loop over the histograms
# for yi in range(ndim):
#     for xi in range(yi):
#         ax = axes[yi, xi]
#         ax.axvline(maxvalues_previous[xi], color="r")
#         ax.axhline(maxvalues_previous[yi], color="r")
#         ax.plot(maxvalues_previous[xi], maxvalues_previous[yi],'sr')

legend_elements = [Line2D([0], [0], color='k', ls='-',lw=2, label='True'),
                   Line2D([0], [0], color='#348ABD', ls='-',lw=2, label='Signal Source'),
                   Line2D([0], [0], color='g', ls='-',lw=2, label='Other seed Signal Source')]
                #    Line2D([0], [0], color='r', ls='-',lw=2, label='Search2')]

axes[0,3].legend(handles=legend_elements, loc='upper left')
# plt.legend(loc='upper right')
for ax in fig.get_axes():
    ax.tick_params(axis='both', labelsize=14)
fig.subplots_adjust(right=1.5,top=1.5)
fig.set_size_inches(fig_size_squared[0],fig_size_squared[1], forward=True)
fig.savefig('/home/stefan/LDC/LDC/pictures/corner3.png',dpi=300,bbox_inches='tight')
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
