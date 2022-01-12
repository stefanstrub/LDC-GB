#%%
from re import escape
from gpytorch.likelihoods import likelihood
from matplotlib.markers import MarkerStyle
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rcParams
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch
from getdist import plots, MCSamples
import getdist
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
import seaborn as sns

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

import matplotlib.font_manager
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
    def __init__(self,tdi_fs,Tobs, lower_frequency, upper_frequency, recombination=0.75) -> None:
        self.tdi_fs = tdi_fs
        self.GB = fastGB.FastGB(delta_t=dt, T=Tobs)
        self.reduced_frequency_boundaries = None
        self.recombination = recombination

        tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft()) for k in ["X", "Y", "Z"]]))
        f, psdX =  scipy.signal.welch(tdi_ts["X"], fs=1.0/dt, window='hanning', nperseg=len(tdi_ts["X"])/1)
        f, psdY =  scipy.signal.welch(tdi_ts["Y"], fs=1.0/dt, window='hanning', nperseg=len(tdi_ts["X"])/1)
        f, psdZ =  scipy.signal.welch(tdi_ts["Z"], fs=1.0/dt, window='hanning', nperseg=len(tdi_ts["X"])/1)
        psd = psdX + psdY + psdZ
        indexes = np.logical_and(f>lower_frequency-padding, f<upper_frequency+padding)
        psd = psd[indexes]
        
        # plt.figure()
        # plt.plot(f[indexes], psdX[indexes])
        # plt.plot(f[indexes], psd)
        # plt.show()

        frequencyrange =  [lower_frequency-padding, upper_frequency+padding]
        indexes = np.logical_and(tdi_fs['X'].f > frequencyrange[0], tdi_fs['X'].f < frequencyrange[1]) 
        self.dataX = tdi_fs["X"][indexes]
        self.dataY = tdi_fs["Y"][indexes]
        self.dataZ = tdi_fs["Z"][indexes]

        self.DAf = (self.dataZ - self.dataX)/np.sqrt(2.0)
        self.DEf = (self.dataZ - 2.0*self.dataY + self.dataX)/np.sqrt(6.0)

        # plt.figure()
        # plt.plot(f[indexes], np.abs(self.DAf))
        # plt.plot(f[indexes], np.abs(self.DEf))
        # plt.plot(f[indexes], np.abs(self.DAf)+np.abs(self.DEf))
        # plt.show()

        fmin, fmax = float(self.dataX.f[0]), float(self.dataX.f[-1] + self.dataX.attrs["df"])
        freq = np.array(self.dataX.sel(f=slice(fmin, fmax)).f)
        Nmodel = get_noise_model(noise_model, freq)
        self.Sn = Nmodel.psd(freq=freq, option="X")
        self.SA = Nmodel.psd(freq=freq, option="A")
        self.SE = Nmodel.psd(freq=freq, option="E")

        f_0 = fmin
        f_transfer = 19.1*10**-3
        snr = 7
        amplitude_lower = 2*snr/(Tobs * np.sin(f_0/ f_transfer)/self.SA[0])**0.5
        snr = 2000
        amplitude_upper = 2*snr/(Tobs * np.sin(f_0/ f_transfer)/self.SA[0])**0.5
        amplitude = [amplitude_lower, amplitude_upper]
        print('lower frequency', lower_frequency)
        # print('amplitude boundaries', amplitude)
        # print('amplitude boundaries previous', np.sqrt(np.max(psd))/1000,np.sqrt(np.max(psd))/10)
        # print('amplitude boundaries previous', np.max(np.abs(self.DAf.data))/10**7,np.max(np.abs(self.DAf.data))/10**5)
        # amplitude = np.sqrt(np.max(psd))

        # indexes = np.argsort(p.get('Frequency'))
        # index_low = np.searchsorted(p.get('Frequency')[indexes], frequencyrange[0])
        # index_high = np.searchsorted(p.get('Frequency')[indexes], frequencyrange[1])

        # pGB = deepcopy(pGBadded)
        # self.pGB = {'Amplitude': 3.676495e-22, 'EclipticLatitude': 0.018181, 'EclipticLongitude': 1.268061, 'Frequency': 0.01158392, 'FrequencyDerivative': 8.009579e-15, 'Inclination': 0.686485, 'InitialPhase': 4.201455, 'Polarization': 2.288223}
        

        self.boundaries = {
            "Amplitude": [np.log10(amplitude[0]),np.log10(amplitude[1])],
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
        self.pGBs = {}
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

        # self.DAf = (2/3*self.dataX - self.dataY - self.dataZ)/3.0
        # self.DEf = (self.dataZ - self.dataY)/np.sqrt(3.0)

        # Xs, Ys, Zs = (
        #     Xs[lowerindex:higherindex],
        #     Ys[lowerindex:higherindex],
        #     Zs[lowerindex:higherindex],
        # )

        # diff = np.abs(self.dataX - Xs.values) ** 2 + np.abs(self.dataY - Ys.values) ** 2 + np.abs(self.dataZ - Zs.values) ** 2
        # p1 = float(np.sum(diff / (self.Sn + noise)) * Xs.df) / 2.0
        # p1 = -p1
        # diff = np.abs(self.dataX) ** 2 + np.abs(self.dataY) ** 2 + np.abs(self.dataZ) ** 2
        # null_pGBs = deepcopy(self.pGBs)
        # null_pGBs['Amplitude'] = 4*10**-25
        # print('pGB', self.pGB, self.loglikelihood([self.pGB]))

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

        # Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=self.pGB, oversample=4, simulator="synthlisa")
        # index_low = np.searchsorted(Xs.f, self.dataX.f[0])
        # Xs = Xs[index_low : index_low + len(self.dataX)]
        # Ys = Ys[index_low : index_low + len(self.dataY)]
        # Zs = Zs[index_low : index_low + len(self.dataZ)]

        # Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=self.pGB, oversample=8, simulator="synthlisa")
        # index_low = np.searchsorted(Xs.f, self.dataX.f[0])
        # Xs = Xs[index_low : index_low + len(self.dataX)]
        # Ys = Ys[index_low : index_low + len(self.dataY)]
        # Zs = Zs[index_low : index_low + len(self.dataZ)]
        # ax1.plot(Xs.f * 1000, Xs.values.real, label="VGB2", marker=".", zorder=5)
                    
        # Af = (Zs - Xs)/np.sqrt(2.0)

        ax1.semilogy(self.DAf.f* 1000, np.abs(self.DAf), label='Data')
        # ax1.semilogy(Af.f* 1000, np.abs(Af.data), marker='.', label='Injection')

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
            initial_guess01 = np.zeros((len(parameters)-1)*number_of_signals)
            for signal in range(number_of_signals):
                pGBstart01 = scaleto01(initial_guess[signal], self.boundaries_reduced)
                for count, parameter in enumerate(parameters-['Polarization']):
                    initial_guess01[count+(len(parameters)-1)*signal] = pGBstart01[parameter]
            start = time.time()
            res, energies = differential_evolution(self.function_evolution, bounds=bounds, disp=False, strategy='best1exp', popsize=10,tol= 1e-8 , maxiter=1000, recombination= self.recombination, mutation=(0.5,1), x0=initial_guess01)
            print('time',time.time()-start)
        else:
            start = time.time()
            res, energies = differential_evolution(self.function_evolution, bounds=bounds, disp=False, strategy='best1exp', popsize=8,tol= 1e-10 , maxiter=1000, recombination= self.recombination, mutation=(0.5,1))
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
        number_of_signals_optimize = len(pGBmodes[0])
        for signal in range(number_of_signals_optimize):
            bounds += ((0,1),(0,1),(0,1),(0,1),(0,1),(-0.25,1.25),(0,1),(0,1))
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
            # print('optimized loglikelihood', self.loglikelihood(maxpGB),self.loglikelihood([self.pGB]))
        maxpGB = current_maxpGB
        print('final optimized loglikelihood', self.loglikelihood(maxpGB),maxpGB[0]['Frequency'])
        return maxpGB

    def function(self, pGBs01, boundaries_reduced):
        pGBs = []
        for signal in range(int(len(pGBs01)/8)):
            pGBs.append({})
            i = 0
            for parameter in parameters:
                if parameter in ["EclipticLatitude"]:
                    pGBs[signal][parameter] = np.arcsin((pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
                elif parameter in ["Inclination"]:
                    shifted_inclination = pGBs01[signal*8:signal*8+8][i]
                    if pGBs01[signal*8:signal*8+8][i] < 0:
                        shifted_inclination = pGBs01[signal*8:signal*8+8][i] + 1
                    if pGBs01[signal*8:signal*8+8][i] > 1:
                        shifted_inclination = pGBs01[signal*8:signal*8+8][i] - 1
                    pGBs[signal][parameter] = np.arccos((shifted_inclination * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
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
                    step_size[parameter] = 1e-10
                    # if parameter == 'Frequency':
                    #     step_size[parameter] = 0.00001
                else:
                    step_size[parameter] = 0.001/np.sqrt(inner_product[parameter][parameter])
                if step_size[parameter] > 1e-10:
                    step_size[parameter] = 1e-10
                pGB_low = maxpGB01[parameter] #- step_size[parameter]/2
                pGB_high = maxpGB01[parameter] + step_size[parameter]
                print(parameter, step_size[parameter],i)
                # print(parameter, pGB_low, pGB_high)
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
            # print(step_size['Amplitude'],inner_product['Amplitude']['Amplitude'],step_size['Frequency'],inner_product['Frequency']['Frequency'])
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
pGBadded8 = {}
pGBadded8['Amplitude'] = 1.36368e-22
pGBadded8['EclipticLatitude'] = -0.2
pGBadded8['EclipticLongitude'] = 1.4
pGBadded8['Frequency'] = 0.00120457
pGBadded8['FrequencyDerivative'] = 1e-19
pGBadded8['Inclination'] = 0.5
pGBadded8['InitialPhase'] = 3
pGBadded8['Polarization'] = 2
pGBadded9 = {}
pGBadded9['Amplitude'] = 1.36368e-21
pGBadded9['EclipticLatitude'] = -0.2
pGBadded9['EclipticLongitude'] = 1.4
pGBadded9['Frequency'] = 0.00130457
pGBadded9['FrequencyDerivative'] = 1e-19
pGBadded9['Inclination'] = 0.3
pGBadded9['InitialPhase'] = 3
pGBadded9['Polarization'] = 2
pGBadded10 = {}
pGBadded10['Amplitude'] = 1.36368e-21*0.5
pGBadded10['EclipticLatitude'] = -0.2
pGBadded10['EclipticLongitude'] = 1.4
pGBadded10['Frequency'] = 0.00140457
pGBadded10['FrequencyDerivative'] = 1e-19
pGBadded10['Inclination'] = 0.4
pGBadded10['InitialPhase'] = 3
pGBadded10['Polarization'] = 2
pGBadded11 = {'Amplitude': 1.36368e-22, 'EclipticLatitude': -0.529009, 'EclipticLongitude': -2.51031, 'Frequency': 0.00125313, 'FrequencyDerivative': 9.159587298288947e-19, 'Inclination': 0.244346, 'InitialPhase': 2.64414439, 'Polarization': 2.229426357}
pGBadded11['Frequency'] -= 2*10**-6
# pGBadded = {'Amplitude': 6.486747e-24, 'EclipticLatitude': -0.318064, 'EclipticLongitude': 4.369758, 'Frequency': 0.01620542, 'FrequencyDerivative': 8.508619e-14, 'Inclination': 0.866814, 'InitialPhase': 6.28136, 'Polarization': 5.978979}

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

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
sangria_fn_noiseless = DATAPATH + "/LDC1-3_VGB_v2_FD_noiseless.hdf5"
FD5 = LISAhdf5(sangria_fn_noiseless)
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
tdi_ts_noiseless = xr.Dataset(dict([(k, TimeSeries(td[:int(len(td[:,1])/reduction), n], dt=dt)) for k, n in [["X", 1], ["Y", 2], ["Z", 3]]]))
# tdi_ts = xr.Dataset(dict([(k,TimeSeries(tdi_ts[k][:,1], dt=dt)) for k in ["X", "Y", "Z"]]))
tdi_fs_noiseless = xr.Dataset(dict([(k, tdi_ts_noiseless[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))


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

target_frequencies = p.get('Frequency')
lower_frequency = target_frequencies[0]-1e-6
upper_frequency = target_frequencies[0]+1e-6
range_index = np.logical_and(tdi_fs.f > lower_frequency, tdi_fs.f < upper_frequency)

fig = plt.figure()
tdi_ts_1 = xr.Dataset(dict([(k, TimeSeries(td[:int(len(td[:,1])), n], dt=dt)) for k, n in [["X", 1], ["Y", 2], ["Z", 3]]]))
tdi_fs_1 = xr.Dataset(dict([(k, tdi_ts_1[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
range_index_1 = np.logical_and(tdi_fs_1.f > lower_frequency, tdi_fs_1.f < upper_frequency)
plt.plot(tdi_fs_1.f[range_index_1],tdi_fs_1['X'][range_index_1])
plt.plot(tdi_fs.f[range_index],tdi_fs['X'][range_index])
#subtract noisless data
for k in ["X", "Y", "Z"]:
    tdi_ts[k].data = tdi_ts[k].data - tdi_ts_noiseless[k].data
tdi_fs = xr.Dataset(dict([(k, tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
plt.plot(tdi_fs.f[range_index],tdi_fs['X'][range_index])

# add signals using fastGB
pGB_injected = []
for i in range(len(p.get('Amplitude'))):
    pGBs = {}
    for parameter in parameters:
        pGBs[parameter] = p.get(parameter)[i]
    pGB_injected.append(pGBs)
for pGBadding in pGB_injected:#, pGBadded2, pGBadded3, pGBadded4]:
    Xs_added, Ys_added, Zs_added = GB.get_fd_tdixyz(template=pGBadding, oversample=4, simulator="synthlisa")
    source_added = dict({"X": Xs_added, "Y": Ys_added, "Z": Zs_added})
    index_low = np.searchsorted(tdi_fs["X"].f, Xs_added.f[0])
    index_high = index_low+len(Xs_added)
    # tdi_fs['X'] = tdi_fs['X'] #+ Xs_added
    for k in ["X", "Y", "Z"]:
        tdi_fs[k].data[index_low:index_high] = tdi_fs[k].data[index_low:index_high] + source_added[k].data
tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft(dt=dt)) for k, n in [["X", 1], ["Y", 2], ["Z", 3]]]))

plt.plot(tdi_fs.f[range_index],tdi_fs['X'][range_index])
plt.show()


noise_model = "MRDv1"
Nmodel = get_noise_model(noise_model, np.logspace(-5, -1, 100))
Npsd = Nmodel.psd()

# for pGBadding in [pGBadded7,pGBadded8,pGBadded9, pGBadded10, pGBadded11]:#, pGBadded2, pGBadded3, pGBadded4]:
#     for parameter in parameters:
#         values = p.get(parameter)
#         values = np.append(values, pGBadding[parameter])
#         unit = p.units[parameter]
#         p.addPar(parameter,values,unit)
#     Xs_added, Ys_added, Zs_added = GB.get_fd_tdixyz(template=pGBadding, oversample=4, simulator="synthlisa")
#     source_added = dict({"X": Xs_added, "Y": Ys_added, "Z": Zs_added})
#     index_low = np.searchsorted(tdi_fs["X"].f, Xs_added.f[0])
#     index_high = index_low+len(Xs_added)
#     # tdi_fs['X'] = tdi_fs['X'] #+ Xs_added
#     for k in ["X", "Y", "Z"]:
#         tdi_fs[k].data[index_low:index_high] = tdi_fs[k].data[index_low:index_high] + source_added[k].data
# tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft(dt=dt)) for k, n in [["X", 1], ["Y", 2], ["Z", 3]]]))



pGB = {}
ind = 0
found_sources = []
target_sources = []
first_start = time.time()
np.random.seed(40) #40
# for ind in range(1,len(p.get('Frequency'))):
number_of_signals = 1
signals_per_subtraction = 1

print(pGBadded7["FrequencyDerivative"] * Tobs)
print('smear f', 300*pGBadded7["Frequency"] * 10**3 / 10**9)

snr_peak = 10
def snr_prior(snr):
    return 3*snr/(4*snr_peak**2*(1+snr/(4*snr_peak))**5)
snr = np.linspace(0,200,1000)
# plt.figure()
# plt.plot(snr, snr_prior(snr))
# plt.show()

def transform_uniform_to_snr_prior(x, c):
    return -2 * np.sqrt(2) * np.sqrt((c**6*x - np.sqrt(c**12* x**3 + c**12 * x**2))**(1/3)/x - c**4/(c**6 *x - np.sqrt(c**12 *x**3 + c**12* x**2))**(1/3)) - 1/2* np.sqrt(-1*(32 *(c**6 *x - np.sqrt(c**12 *x**3 + c**12* x**2))**(1/3))/x + (32* c**4)/(c**6 *x - np.sqrt(c**12* x**3 + c**12 *x**2))**(1/3) - (2048* c**3 - (2048 *(c**3 *x + c**3))/x)/(16 *np.sqrt(2) * np.sqrt((c**6* x - np.sqrt(c**12* x**3 + c**12* x**2))**(1/3)/x - c**4/(c**6 *x - np.sqrt(c**12* x**3 + c**12 *x**2))**(1/3)))) - 4 * c


import scipy.stats as st

class my_pdf(st.rv_continuous):
    def _pdf(self,x):
        return 3*x/(4*snr_peak**2*(1+x/(4*snr_peak))**5)  # Normalized over its range, in this case [0,1]

my_cv = my_pdf(a=0, b=10000, name='my_pdf')
my_cv.pdf(1)
my_cv.rvs(size= 100)


class MLP_search():
    def __init__(self,tdi_fs, Tobs, signals_per_window) -> None:
        self.tdi_fs = tdi_fs
        self.Tobs = Tobs
        self.signals_per_window = signals_per_window

    def search(self, lower_frequency, upper_frequency):
        found_sources = []
        # indexes = np.argsort(p.get('Frequency'))
        # index_low = np.searchsorted(p.get('Frequency')[indexes], lower_frequency)
        # index_high = np.searchsorted(p.get('Frequency')[indexes], upper_frequency)
        # pGB_injected = []
        # for i in range(len(p.get('Amplitude')[indexes][index_low:index_high])):
        #     pGBs = {}
        #     for parameter in parameters:
        #         pGBs[parameter] = p.get(parameter)[indexes][index_low:index_high][i]
        #     pGB_injected.append(pGBs)
        tdi_fs_search = deepcopy(self.tdi_fs)
        # previous_found_sources = [{'Amplitude': 4.084935966774485e-22, 'EclipticLatitude': 0.8719934546490874, 'EclipticLongitude': 0.48611009683797857, 'Frequency': 0.003995221087430858, 'FrequencyDerivative': 1.0704703957490903e-16, 'Inclination': 1.0245091695238984, 'InitialPhase': 2.320136113624083, 'Polarization': 2.65883774239409}, {'Amplitude': 1.170377953453263e-22, 'EclipticLatitude': -1.1827019140449202, 'EclipticLongitude': -2.6708716710257203, 'Frequency': 0.003994619937260686, 'FrequencyDerivative': 9.604827167870394e-17, 'Inclination': 1.9399867466326164, 'InitialPhase': 2.468693959968005, 'Polarization': 2.5128702009090644}]

        current_SNR = 100
        ind = 0
        while current_SNR > 10 and ind < self.signals_per_window:
            ind += 1
            
            search1 = Search(tdi_fs_search,self.Tobs, lower_frequency, upper_frequency)

            start = time.time()

            # print('SNR ',np.round(search1.SNR([search1.pGB])))
            # print('SNR2', np.round(search1.loglikelihood([search1.pGB])))
            # print('SNR2', np.round(search1.loglikelihoodsdf([search1.pGB])))
            # print('SNRm', np.round(search1.SNRm([search1.pGB]),3))
            # print('SNRflat', np.round(search1.loglikelihoodflat([search1.pGB])))
            # search1.plot()#pGBadded=pGBadded5)
            # print(pGBadded7["FrequencyDerivative"] * self.Tobs)
            # print('smear f', 300*pGBadded7["Frequency"] * 10**3 / 10**9)
            # print(search1.reduced_frequency_boundaries)
            for i in range(1):
                # if i > 0:
                #     search1.recombination = 0.75
                #     maxpGBsearch_new, energies =  search1.differential_evolution_search(search1.boundaries['Frequency'], initial_guess=maxpGBsearch[0])
                # else:
                maxpGBsearch_new, energies =  search1.differential_evolution_search(search1.boundaries['Frequency'])
                # maxpGBsearch = [[previous_found_sources[ind-1]]]
                # print(search1.loglikelihood(maxpGBsearch[0]), ', ', ind, '. Source')
                print('SNRm of found signal', np.round(search1.SNRm(maxpGBsearch_new[0]),3))
                print('which signal per window', ind, i)
                print(maxpGBsearch_new[0][0]['Frequency'] > lower_frequency and maxpGBsearch_new[0][0]['Frequency'] < upper_frequency)
                # if maxpGBsearch[0][0]['Frequency'] > lower_frequency and maxpGBsearch[0][0]['Frequency'] < upper_frequency:
                new_SNR = search1.SNRm(maxpGBsearch_new[0])[0]
                if i == 0:
                    current_SNR = deepcopy(new_SNR)
                    maxpGBsearch = deepcopy(maxpGBsearch_new)
                if new_SNR > current_SNR:
                    current_SNR = deepcopy(new_SNR)
                    maxpGBsearch = deepcopy(maxpGBsearch_new)
                print('current SNR', current_SNR)
                if current_SNR < 9:
                    break

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
                    tdi_fs_subtracted = deepcopy(tdi_fs)
                    for i in range(len(found_sources_out)):
                        Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_out[i], oversample=4, simulator="synthlisa")
                        source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
                        index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
                        index_high = index_low+len(Xs_subtracted)
                        for k in ["X", "Y", "Z"]:
                            tdi_fs_subtracted[k].data[index_low:index_high] = tdi_fs_subtracted[k].data[index_low:index_high] - source_subtracted[k].data
                        tdi_ts_subtracted = xr.Dataset(dict([(k, tdi_fs_subtracted[k].ts.ifft()) for k in ["X", "Y", "Z"]]))

                    search_out_subtracted = Search(tdi_fs_subtracted,self.Tobs, lower_frequency, upper_frequency)

                    total_boundaries = deepcopy(search1.boundaries)
                    amplitudes = []
                    for i in range(len(found_sources_in)):
                        amplitudes.append(found_sources_in[i]['Amplitude'])
                    total_boundaries['Amplitude'] = [np.min(amplitudes),np.max(amplitudes)]
                    amplitudes_length = np.log10(total_boundaries['Amplitude'][1]) - np.log10(total_boundaries['Amplitude'][0])
                    total_boundaries['Amplitude'] = [np.log10(total_boundaries['Amplitude'][0]) - amplitudes_length/5, np.log10(total_boundaries['Amplitude'][1]) + amplitudes_length/5]

                    
                    start = time.time()
                    found_sources_in = search_out_subtracted.optimize([found_sources_in], boundaries= total_boundaries)
                    print(time.time()-start)

            found_sources = found_sources_in + found_sources_out

            #subtract the found sources from original
            tdi_fs_search = deepcopy(self.tdi_fs)
            for i in range(len(found_sources)):
                Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources[i], oversample=4, simulator="synthlisa")
                source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
                index_low = np.searchsorted(tdi_fs_search["X"].f, Xs_subtracted.f[0])
                index_high = index_low+len(Xs_subtracted)
                for k in ["X", "Y", "Z"]:
                    tdi_fs_search[k].data[index_low:index_high] = tdi_fs_search[k].data[index_low:index_high] - source_subtracted[k].data
        return found_sources

lower_frequency = 0.0039945
upper_frequency = 0.0039955
lower_frequency2 = 0.0039965
upper_frequency2 = 0.0039975

padding = 0.5e-6

save_name = 'fastGB'
# LDC1-3 ##########################################
target_frequencies = p.get('Frequency')
frequencies = []
window_length = 10**-6 # Hz
for i in range(len(target_frequencies)):
    window_shift = ((np.random.random(1)-0.5)*window_length*0.5)[0]
    frequencies.append([target_frequencies[i]-window_length/2+window_shift,target_frequencies[i]+window_length/2+window_shift])
frequencies = frequencies[:10]
number_of_windows = len(frequencies)
MLP = MLP_search(tdi_fs, Tobs, signals_per_window = 1)
# found_sources = MLP.search(frequencies[0][0], frequencies[0][1])
# start = time.time()
# pool = mp.Pool(mp.cpu_count())
# found_sources_mp= pool.starmap(MLP.search, frequencies)
# pool.close()
# pool.join()
# print('time to search ', number_of_windows, 'windows: ', time.time()-start)

found_sources_mp = np.load('/home/stefan/LDC/LDC/pictures/found_sources_ldc1-3'+save_name+'.npy', allow_pickle= True)
# found_sources_mp = found_sources_mp[:-1]
# np.save('/home/stefan/LDC/LDC/pictures/found_sources_ldc1-3'+save_name+'.npy', found_sources_mp)

# LDC1-4 #####################################
# frequencies = []
# frequencies_even = []
# frequencies_odd = []
# window_length = 1*10**-6 # Hz
# search_range = [0.00398, 0.0041]
# search_range = [0.0039885, 0.0040045]
# search_range = [0.0039935, 0.0039965]
# number_of_windows = int(round((search_range[1]-search_range[0])/window_length))
# for i in range(number_of_windows):
#     frequencies.append([search_range[0]+window_length*i,search_range[0]+window_length*(i+1)])
# frequencies_even = frequencies[::2]
# frequencies_odd = frequencies[1::2]


# MLP = MLP_search(tdi_fs, Tobs, signals_per_window = 3)
# start = time.time()
# pool = mp.Pool(mp.cpu_count())
# found_sources_mp_even = pool.starmap(MLP.search, frequencies_even)
# pool.close()
# pool.join()
# print('time to search ', number_of_windows, 'windows: ', time.time()-start)

# #np.save('/home/stefan/LDC/LDC/pictures/found_sources_ldc1-4 even 3.npy', found_sources_mp_even)

# found_sources_in = []
# found_sources_out = []
# for i in range(len(found_sources_mp_even)):
#     found_sources_in.append([])
#     found_sources_out.append([])
#     for j in range(len(found_sources_mp_even[i])):
#         if found_sources_mp_even[i][j]['Frequency'] > frequencies_even[i][0] and found_sources_mp_even[i][j]['Frequency'] < frequencies_even[i][1]:
#             found_sources_in[i].append(found_sources_mp_even[i][j])
#         else:
#             found_sources_out[i].append(found_sources_mp_even[i][j])

# #subtract the found sources from original
# tdi_fs_subtracted = deepcopy(tdi_fs)
# for i in range(len(found_sources_in)):
#     for j in range(len(found_sources_in[i])):
#         Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_in[i][j], oversample=4, simulator="synthlisa")
#         source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
#         index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
#         index_high = index_low+len(Xs_subtracted)
#         for k in ["X", "Y", "Z"]:
#             tdi_fs_subtracted[k].data[index_low:index_high] = tdi_fs_subtracted[k].data[index_low:index_high] - source_subtracted[k].data

# indexes = np.argsort(p.get('Frequency'))
# pGB_injected = []
# for j in range(len(found_sources_mp_even)):
#     index_low = np.searchsorted(p.get('Frequency')[indexes], frequencies_even[j][0])
#     index_high = np.searchsorted(p.get('Frequency')[indexes], frequencies_even[j][1])
#     pGB_injected_window = []
#     for i in range(len(p.get('Amplitude')[indexes][index_low:index_high])):
#         pGBs = {}
#         for parameter in parameters:
#             pGBs[parameter] = p.get(parameter)[indexes][index_low:index_high][i]
#         pGB_injected_window.append(pGBs)
#     pGB_injected.append(pGB_injected_window)

# # #plot strains
# for i in range(len(found_sources_mp_even)):
#     lower_frequency = frequencies_even[i][0]
#     upper_frequency = frequencies_even[i][1]
#     indexes = np.argsort(p.get('Frequency'))
#     range_index = np.logical_and(tdi_fs.f > lower_frequency-padding, tdi_fs.f < upper_frequency+padding)

#     indexes = np.logical_and(tdi_fs['X'].f > lower_frequency-padding, tdi_fs['X'].f < upper_frequency+padding) 
#     dataX = tdi_fs["X"][indexes]

#     fig = plt.figure(figsize=fig_size)
#     ax1 = plt.subplot(111)
#     ax1.semilogy(tdi_fs.f[range_index]*10**3,np.abs(tdi_fs['X'][range_index])**2,'k',zorder= 1, linewidth = 2, label = 'Data')
#     # ax1.semilogy(tdi_fs_long_subtracted.f[range_index],np.abs(tdi_fs_long_subtracted['X'][range_index])**2,'b',zorder= 5)


#     for j in range(len( pGB_injected[i])):
#         Xs, Ys, Zs = GB.get_fd_tdixyz(template= pGB_injected[i][j], oversample=4, simulator="synthlisa")
#         a,Xs = xr.align(dataX, Xs, join='left',fill_value=0)
#         ax1.semilogy(Xs.f*10**3,np.abs(Xs)**2,label= str(np.round(pGB_injected[i][j]['Frequency'],0)), color='grey', linewidth = 3)

#     for j in range(len(found_sources_in[i])):
#         Xs, Ys, Zs = GB.get_fd_tdixyz(template= found_sources_in[i][j], oversample=4, simulator="synthlisa")
#         ax1.semilogy(Xs.f*10**3,np.abs(Xs)**2,'--', color= colors[j], linewidth = 1.8)
#     # Xs, Ys, Zs = GB.get_fd_tdixyz(template= found_sources_mp2[i][0], oversample=4, simulator="synthlisa")
#     # ax1.semilogy(Xs.f*10**3,np.abs(Xs)**2,'-', color= colors[i], linewidth = 1.8)
#     # Xs, Ys, Zs = GB.get_fd_tdixyz(template= pGB, oversample=4, simulator="synthlisa")
#     # a,Xs = xr.align(search1.dataX, Xs, join='left',fill_value=0)
#     # ax1.semilogy(Xs.f,np.abs(tdi_fs['X'][range_index]-Xs)**2,label= 'residual')
#     plt.xlim((lower_frequency-padding)*10**3, (upper_frequency+padding)*10**3)
#     # ax1.axhline(y=0)
#     ax1.xaxis.set_major_locator(plt.MaxNLocator(4))
#     # plt.legend()
#     plt.xlabel('f [mHz]')
#     plt.ylabel('|X|')
#     fig.savefig('/home/stefan/LDC/LDC/pictures/strain first ldc1-4 '+ str(lower_frequency*1000) +'.png',dpi=300,bbox_inches='tight')
# # plt.show()  

# # #plot subtraction
# for i in range(len(frequencies_odd)):
#     lower_frequency = frequencies_odd[i][0]
#     upper_frequency = frequencies_odd[i][1]
#     indexes = np.argsort(p.get('Frequency'))
#     range_index = np.logical_and(tdi_fs.f > lower_frequency-padding, tdi_fs.f < upper_frequency+padding)

#     indexes = np.logical_and(tdi_fs['X'].f > lower_frequency-padding, tdi_fs['X'].f < upper_frequency+padding) 
#     dataX = tdi_fs["X"][indexes]

#     fig = plt.figure(figsize=fig_size)
#     ax1 = plt.subplot(111)
#     ax1.semilogy(tdi_fs.f[range_index]*10**3,np.abs(tdi_fs['X'][range_index])**2,'k',zorder= 1, linewidth = 2, label = 'Data')
#     ax1.semilogy(tdi_fs.f[range_index]*10**3,np.abs(tdi_fs_subtracted['X'][range_index])**2,zorder= 1, linewidth = 2, label = 'Subtracted data')
#     # ax1.semilogy(tdi_fs_long_subtracted.f[range_index],np.abs(tdi_fs_long_subtracted['X'][range_index])**2,'b',zorder= 5)



#     # Xs, Ys, Zs = GB.get_fd_tdixyz(template= found_sources_mp2[i][0], oversample=4, simulator="synthlisa")
#     # ax1.semilogy(Xs.f*10**3,np.abs(Xs)**2,'-', color= colors[i], linewidth = 1.8)
#     # Xs, Ys, Zs = GB.get_fd_tdixyz(template= pGB, oversample=4, simulator="synthlisa")
#     # a,Xs = xr.align(search1.dataX, Xs, join='left',fill_value=0)
#     # ax1.semilogy(Xs.f,np.abs(tdi_fs['X'][range_index]-Xs)**2,label= 'residual')
#     plt.xlim((lower_frequency-padding)*10**3, (upper_frequency+padding)*10**3)
#     # ax1.axhline(y=0)
#     ax1.xaxis.set_major_locator(plt.MaxNLocator(4))
#     # plt.legend()
#     plt.xlabel('f [mHz]')
#     plt.ylabel('|X|')
#     plt.legend()
#     fig.savefig('/home/stefan/LDC/LDC/pictures/subtraction ldc1-4 solved '+ str(lower_frequency*1000) +'.png',dpi=300,bbox_inches='tight')
# # plt.show()  

# MLP = MLP_search(tdi_fs_subtracted, Tobs, signals_per_window = 10)
# start = time.time()
# pool = mp.Pool(mp.cpu_count())
# found_sources_mp_odd = pool.starmap(MLP.search, frequencies_odd)
# pool.close()
# pool.join()
# print('time to search ', number_of_windows, 'windows: ', time.time()-start)

# found_sources_in = []
# found_sources_out = []
# for i in range(len(found_sources_mp_odd)):
#     found_sources_in.append([])
#     found_sources_out.append([])
#     for j in range(len(found_sources_mp_odd[i])):
#         if found_sources_mp_odd[i][j]['Frequency'] > frequencies_odd[i][0] and found_sources_mp_odd[i][j]['Frequency'] < frequencies_odd[i][1]:
#             found_sources_in[i].append(found_sources_mp_odd[i][j])
#         else:
#             found_sources_out[i].append(found_sources_mp_odd[i][j])

# #subtract the found sources from original
# tdi_fs_subtracted = deepcopy(tdi_fs)
# for i in range(len(found_sources_in)):
#     for j in range(len(found_sources_in[i])):
#         Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_in[i][j], oversample=4, simulator="synthlisa")
#         source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
#         index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
#         index_high = index_low+len(Xs_subtracted)
#         for k in ["X", "Y", "Z"]:
#             tdi_fs_subtracted[k].data[index_low:index_high] = tdi_fs_subtracted[k].data[index_low:index_high] - source_subtracted[k].data

# MLP = MLP_search(tdi_fs_subtracted, Tobs, signals_per_window = 10)
# start = time.time()
# pool = mp.Pool(mp.cpu_count())
# found_sources_mp_even = pool.starmap(MLP.search, frequencies_even)
# pool.close()
# pool.join()
# print('time to search ', number_of_windows, 'windows: ', time.time()-start)

# found_sources_mp =[]
# for i in range(number_of_windows):
#     ind = int(i/2)
#     if i % 2 == 0:
#         found_sources_mp.append(found_sources_mp_even[ind])
#     else:
#         found_sources_mp.append(found_sources_mp_odd[ind])

# np.save('/home/stefan/LDC/LDC/pictures/found_sources_ldc1-4.npy', found_sources_mp)
# found_sources_mp = np.load('/home/stefan/LDC/LDC/pictures/found_sources_ldc1-4.npy', allow_pickle= True)

# np.random.seed(42)
# start = time.time()
# pool = mp.Pool(mp.cpu_count())
# found_sources_mp2 = pool.starmap(MLP_search, frequencies)
# pool.close()
# pool.join()
# print('time to search ', len(target_frequencies), 'models: ', time.time()-start)

# pGB_injected = []
# for i in range(len(p.get('Amplitude'))):
#     pGBs = {}
#     for parameter in parameters:
#         pGBs[parameter] = p.get(parameter)[i]
#     pGB_injected.append(pGBs)
# found_sources_mp = [[{'Amplitude': 1.5472182628836777e-22, 'EclipticLatitude': 0.32392807398767537, 'EclipticLongitude': -2.7551030697541634, 'Frequency': 0.001359619859748239, 'FrequencyDerivative': 6.741316095619172e-18, 'Inclination': 0.9617478105498454, 'InitialPhase': 0.8814766948044847, 'Polarization': 2.497986422716522}],[{'Amplitude': 2.1806906125113854e-22, 'EclipticLatitude': -0.5283216519943494, 'EclipticLongitude': -2.5390590647151567, 'Frequency': 0.0012531301743037816, 'FrequencyDerivative': 1.5344145457167598e-20, 'Inclination': 0.9737610507888282, 'InitialPhase': 3.963309240828965, 'Polarization': 2.8821924743499507}]]

indexes = np.argsort(p.get('Frequency'))
pGB_injected = []
for j in range(len(found_sources_mp)):
    index_low = np.searchsorted(p.get('Frequency')[indexes], frequencies[j][0])
    index_high = np.searchsorted(p.get('Frequency')[indexes], frequencies[j][1])
    pGB_injected_window = []
    for i in range(len(p.get('Amplitude')[indexes][index_low:index_high])):
        pGBs = {}
        for parameter in parameters:
            pGBs[parameter] = p.get(parameter)[indexes][index_low:index_high][i]
        pGB_injected_window.append(pGBs)
    pGB_injected.append(pGB_injected_window)
    
found_sources_in = []
found_sources_out = []
for i in range(len(found_sources_mp)):
    found_sources_in.append([])
    found_sources_out.append([])
    for j in range(len(found_sources_mp[i])):
        if found_sources_mp[i][j]['Frequency'] > frequencies[i][0] and found_sources_mp[i][j]['Frequency'] < frequencies[i][1]:
            found_sources_in[i].append(found_sources_mp[i][j])
        else:
            found_sources_out[i].append(found_sources_mp[i][j])


#check SNR
for i in range(len(found_sources_mp)):
    search1 = Search(tdi_fs,Tobs, frequencies[i][0], frequencies[i][1])
    for j in range(len( pGB_injected[i])):
        print(search1.SNRm([pGB_injected[i][j]])[0])
    for j in range(len(found_sources_in[i])):
        print('found', search1.SNRm([found_sources_in[i][j]])[0])


#plot strains
for i in range(len(found_sources_mp)):
    lower_frequency = frequencies[i][0]
    upper_frequency = frequencies[i][1]
    indexes = np.argsort(p.get('Frequency'))
    range_index = np.logical_and(tdi_fs.f > lower_frequency-padding, tdi_fs.f < upper_frequency+padding)

    indexes = np.logical_and(tdi_fs['X'].f > lower_frequency-padding, tdi_fs['X'].f < upper_frequency+padding) 
    dataX = tdi_fs["X"][indexes]

    fig = plt.figure(figsize=fig_size)
    ax1 = plt.subplot(111)
    ax1.semilogy(tdi_fs.f[range_index]*10**3,np.abs(tdi_fs['X'][range_index])**2,'k',zorder= 1, linewidth = 2, label = 'Data')
    # ax1.semilogy(tdi_fs_long_subtracted.f[range_index],np.abs(tdi_fs_long_subtracted['X'][range_index])**2,'b',zorder= 5)


    for j in range(len( pGB_injected[i])):
        Xs, Ys, Zs = GB.get_fd_tdixyz(template= pGB_injected[i][j], oversample=4, simulator="synthlisa")
        a,Xs = xr.align(dataX, Xs, join='left',fill_value=0)
        ax1.semilogy(Xs.f*10**3,np.abs(Xs)**2,label= str(np.round(pGB_injected[i][j]['Frequency'],0)), color='grey', linewidth = 3)

    for j in range(len(found_sources_in[i])):
        Xs, Ys, Zs = GB.get_fd_tdixyz(template= found_sources_in[i][j], oversample=4, simulator="synthlisa")
        ax1.semilogy(Xs.f*10**3,np.abs(Xs)**2,'--', color= colors[j], linewidth = 1.8)
    # Xs, Ys, Zs = GB.get_fd_tdixyz(template= found_sources_mp2[i][0], oversample=4, simulator="synthlisa")
    # ax1.semilogy(Xs.f*10**3,np.abs(Xs)**2,'-', color= colors[i], linewidth = 1.8)
    # Xs, Ys, Zs = GB.get_fd_tdixyz(template= pGB, oversample=4, simulator="synthlisa")
    # a,Xs = xr.align(search1.dataX, Xs, join='left',fill_value=0)
    # ax1.semilogy(Xs.f,np.abs(tdi_fs['X'][range_index]-Xs)**2,label= 'residual')
    plt.xlim((lower_frequency-padding)*10**3, (upper_frequency+padding)*10**3)
    # ax1.axhline(y=0)
    ax1.xaxis.set_major_locator(plt.MaxNLocator(4))
    # plt.legend()
    plt.xlabel('f [mHz]')
    plt.ylabel('|X|')
    # fig.savefig('/home/stefan/LDC/LDC/pictures/strain added'+ str(int(np.round(lower_frequency*10**8))) +save_name+'.png',dpi=300,bbox_inches='tight')
plt.show()


prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

fig = plt.figure(figsize=fig_size)
ax1 = plt.subplot(111)
parameter1 = 'EclipticLatitude'
parameter2 = 'EclipticLongitude'
for i in range(len(found_sources_mp)):
    for j in range(len(found_sources_in[i])):
        if i == 0:
            ax1.plot(found_sources_in[i][j][parameter1],found_sources_in[i][j][parameter2], color = colors[i], marker = '+', markersize = 8, label= 'global fit')
        else:
            ax1.plot(found_sources_in[i][j][parameter1],found_sources_in[i][j][parameter2], color = colors[i], marker = '+', markersize = 8)
for i in range(len(pGB_injected)):   
    for j in range(len( pGB_injected[i])):
        if i == 0: 
            ax1.plot(pGB_injected[i][j][parameter1],pGB_injected[i][j][parameter2],color='grey', marker = 'o', markersize = 8, zorder=1, label= 'true')
        else:
            ax1.plot(pGB_injected[i][j][parameter1],pGB_injected[i][j][parameter2],color='grey', marker = 'o', markersize = 8, zorder=1)
# ax1.axvline(10**3*(search1.boundaries[parameter1][0]+padding), color= 'red', label='Boundaries')
# ax1.axvline(10**3*(search1.boundaries[parameter1][1]-padding), color= 'red')
# plt.xlim((lower_frequency-padding), (upper_frequency+padding))
ax1.xaxis.set_major_locator(plt.MaxNLocator(4))
plt.xlabel(parameter1)
plt.ylabel(parameter2)
plt.legend(loc='upper right')
fig.savefig('/home/stefan/LDC/LDC/pictures/global fit sky'+ str(int(np.round(lower_frequency*10**7))) +save_name+'.png',dpi=300,bbox_inches='tight')
plt.show()

# np.save('/home/stefan/LDC/LDC/pictures/found_sources_ldc1-4', found_sources_mp)
# found_sources_mp = [[{'Amplitude': 1.5472707659844358e-22, 'EclipticLatitude': 0.3239131241698293, 'EclipticLongitude': -2.7550743649206386, 'Frequency': 0.0013596198589500806, 'FrequencyDerivative': 6.760619501138522e-18, 'Inclination': 0.9617263677786225, 'InitialPhase': 0.8816156199962063, 'Polarization': 2.4981058468684956}]]

class Posterior_computer():
    def __init__(self, tdi_fs, Tobs, frequencies, maxpGB) -> None:
        self.tdi_fs = tdi_fs
        self.Tobs = Tobs
        self.frequencies = frequencies
        self.maxpGB = maxpGB
        self.search1 = Search(self.tdi_fs,self.Tobs, self.frequencies[0], self.frequencies[1])

    def reduce_boundaries(self):
        fisher_information_matrix = self.search1.fisher_information(self.maxpGB)
        FIM = np.zeros((len(parameters),len(parameters)))
        for i,parameter1 in enumerate(parameters):
            for j,parameter2 in enumerate(parameters):
                FIM[i,j] = fisher_information_matrix[parameter1][parameter2]
        covariance_matrix = scipy.linalg.inv(FIM)
        maxpGB01 = scaleto01(self.maxpGB, self.search1.boundaries)

        lambda_, v = scipy.linalg.eig(covariance_matrix)
        transf = np.zeros((len(lambda_),len(lambda_)))
        for i in range(len(lambda_)):
            transf[i] = v[i] * np.sqrt(lambda_[i])
        covariance2 = v @ np.diag(lambda_) @ scipy.linalg.inv(v)
        transf = v @ np.diag(np.sqrt(lambda_))
        scalematrix = np.max(np.abs(transf), axis=1)
        scalematrix = np.sqrt(np.diag(covariance_matrix))

        maxpGB01_low = deepcopy(maxpGB01)
        maxpGB01_high = deepcopy(maxpGB01)
        boundaries_reduced_fisher = {}
        sigma_multiplyer = 2
        for parameter in parameters:
            maxpGB01_low[parameter] = maxpGB01[parameter] - scalematrix[parameters.index(parameter)] * sigma_multiplyer 
            maxpGB01_high[parameter] = maxpGB01[parameter] + scalematrix[parameters.index(parameter)] * sigma_multiplyer 
            if parameter in [ 'Amplitude', 'Inclination']:
                maxpGB01_low[parameter] = maxpGB01[parameter] - scalematrix[parameters.index(parameter)] * 2.5 
                maxpGB01_high[parameter] = maxpGB01[parameter] + scalematrix[parameters.index(parameter)] * 2.5 
            if parameter in [ 'InitialPhase', 'Polarization']:
                maxpGB01_low[parameter] = maxpGB01[parameter] - 0.001
                maxpGB01_high[parameter] = maxpGB01[parameter] + 0.001
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
            maxpGB_fisher_low = (maxpGB01_low[parameter] * (self.search1.boundaries[parameter][1] - self.search1.boundaries[parameter][0])) + self.search1.boundaries[parameter][0]
            maxpGB_fisher_high = (maxpGB01_high[parameter] * (self.search1.boundaries[parameter][1] - self.search1.boundaries[parameter][0])) + self.search1.boundaries[parameter][0]
            boundaries_reduced_fisher[parameter] = [maxpGB_fisher_low, maxpGB_fisher_high]

        self.boundaries_reduced = deepcopy(boundaries_reduced_fisher)
        split_fd = -17
        if self.boundaries_reduced['FrequencyDerivative'][1] < split_fd+0.5:
            self.boundaries_reduced['FrequencyDerivative'][0] = -18
            self.boundaries_reduced['FrequencyDerivative'][1] = -16
        elif self.boundaries_reduced['FrequencyDerivative'][0] < split_fd+0.5:
            self.boundaries_reduced['FrequencyDerivative'][0] = -18
        print(self.boundaries_reduced)
    
    def train_model(self):
        rmse = 2
        train_size = 0
        test_size = 500
        added_trainig_size = 1000
        j = 0
        samples = np.random.rand(6000,8)
        while rmse > 0.5 and j < 5:
            j += 1
            train_size += added_trainig_size
            if j == 1:
                resolution = train_size + test_size
            else:
                resolution = added_trainig_size
            start = time.time()
            samples_likelihood = np.zeros(resolution)
            samples_likelihood2 = np.zeros(resolution)
            for i in range(resolution):
                samples_p = scaletooriginal(samples[i+j*added_trainig_size], self.boundaries_reduced)
                samples_likelihood[i] = self.search1.loglikelihood([samples_p])
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
                test_y = samples_likelihood[:test_size]
                train_x = samples_flat[test_size:]
                test_x = samples_flat[:test_size]
            else:
                train_y = np.append(train_y, samples_likelihood)
                train_x = np.append(train_x, samples_flat,axis=0)

            self.mu = np.mean(train_y)
            self.sigma = np.std(train_y)
            train_y_normalized = (train_y - self.mu) / self.sigma
            kernel = RBF(length_scale=[1,2,5,1,1,1,1,1],length_scale_bounds=[(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,100),(0.1,20)])
            start = time.time()
            self.gpr = GaussianProcessRegressor(kernel=kernel, random_state=0).fit(train_x, train_y_normalized)
            print('train',time.time() - start)
            start = time.time()
            observed_pred_sk = self.gpr.predict(test_x)
            print('eval time of ', test_size, 'samples: ',time.time() - start)
            observed_pred_sk_scaled = observed_pred_sk*self.sigma + self.mu
            rmse = np.sqrt(mean_squared_error(test_y,observed_pred_sk_scaled))
            print("RMSE ",rmse)
            print('training size', len(train_y))

    def evaluate(self, x):
        partial_length = 1*10**3
        start = time.time()
        observed_pred_mean = np.zeros(len(x))
        observed_pred_sk = np.zeros(len(x))
        for n in range(int(len(x)/partial_length)):
            observed_pred_sk[n*partial_length:(n+1)*partial_length] = self.gpr.predict(x[(n)*partial_length:(n+1)*partial_length])
        try:
            observed_pred_sk[int(len(x)/partial_length)*partial_length:] = self.gpr.predict(x[int(len(x)/partial_length)*partial_length:])
        except:
            pass
        observed_pred_sk = np.asarray(observed_pred_sk)
        observed_pred_sk = observed_pred_sk.reshape(len(x))
        observed_pred_mean[:len(x)] = observed_pred_sk[:len(x)]*self.sigma + self.mu
        print('eval time', time.time()-start)
        return observed_pred_mean

    def calculate_posterior(self,resolution = 1*10**6, proposal= None, temperature = 1):
        self.resolution = resolution
        numPoints = 2**6+1
        start = time.time()
        if proposal is None:
            test_x_m = np.random.uniform(size=(self.resolution,len(parameters)))
            probability = np.ones(self.resolution)
        else:
            test_x_m = np.zeros((self.resolution,len(parameters)))
            ax = np.linspace(-0.15,1.15,numPoints)
            mypdf,axes = fastKDE.pdf(proposal[:,1],proposal[:,2], axes=[ax,ax])
            dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
            data, pdfs = dist(self.resolution)
            test_x_m[:,1] = data[1]
            test_x_m[:,2] = data[0]
            probability = pdfs
            mypdf,axes = fastKDE.pdf(proposal[:,0],proposal[:,5], axes=[ax,ax])
            dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
            data, pdfs = dist(self.resolution)
            test_x_m[:,0] = data[1]
            test_x_m[:,5] = data[0]
            # plt.figure()
            # plt.scatter(np.log10(test_x_m[:10000,0]),np.arccos(test_x_m[:10000,5]* (boundaries_reduced['Inclination'][1] - boundaries_reduced['Inclination'][0]) + boundaries_reduced['Inclination'][0]), c=pdfs[:10000])
            probability *= pdfs
            mypdf,axes = fastKDE.pdf(proposal[:,3],proposal[:,4], axes=[ax,ax])
            dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
            data, pdfs = dist(self.resolution)
            test_x_m[:,3] = data[1]
            test_x_m[:,4] = data[0]
            probability *= pdfs

            test_x_m[:,6] = np.random.uniform(size=self.resolution)
            test_x_m[:,7] = np.random.uniform(size=self.resolution)

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

        observed_pred_mean = self.evaluate(test_x_m)
        flatsamples = np.zeros(len(test_x_m))
        flatsamplesparameters = np.zeros((len(test_x_m),len(parameters)+1))
        i = 0
        flatsamples[:] = observed_pred_mean
        flatsamplesparameters[:,1:] = test_x_m
        flatsamplesparameters[:,0] = observed_pred_mean

        maxindx = np.unravel_index(flatsamplesparameters[:,0].argmax(), flatsamplesparameters[:,0].shape)
        max_parameters = flatsamplesparameters[maxindx[0],1:]
        max_loglike = flatsamplesparameters[:,0].max()
        maxpGBpredicted = scaletooriginal(max_parameters, self.boundaries_reduced)
        if self.search1.loglikelihood([maxpGBpredicted]) > self.search1.loglikelihood([self.maxpGB]):
            maxpGB = maxpGBpredicted
        best_value = self.search1.loglikelihood([self.maxpGB])
        print("pred", max_loglike, "true", self.search1.loglikelihood([scaletooriginal(max_parameters, self.boundaries_reduced)]), "max", self.search1.loglikelihood([self.maxpGB]), self.maxpGB)

        # np.random.shuffle(flatsamplesparameters)
        start = time.time()
        normalizer = np.sum(np.exp(flatsamplesparameters[:,0]-best_value))
        flatsamples_normalized = np.exp(flatsamplesparameters[:,0]-best_value)
        # flatsamples_normalized = flatsamplesparameters[:,0]
        mcmc_samples = []
        mcmc_samples.append(flatsamplesparameters[0,1:])
        previous_p = flatsamples_normalized[0]
        if previous_p == 0:
            previous_p == 1e-300
        current_parameters = flatsamplesparameters[0,1:]
        # probability = scipy.stats.multivariate_normal.pdf(flatsamplesparameters[:,1:], mean=[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5], cov=[std_scale,std_scale,std_scale,std_scale,std_scale,std_scale,std_scale,std_scale])
        
        previous_probability = probability[0]
        accepted = 0
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

        return mcmc_samples

    def plot_corner(self, mcmc_samples, pGB, save_bool = False, save_chain = False):
        start = time.time()
        mcmc_samples_rescaled = np.zeros(np.shape(mcmc_samples))
        i = 0
        for parameter in parameters:
            if parameter in ["EclipticLatitude"]:
                mcmc_samples_rescaled[:,i] = np.arcsin((mcmc_samples[:,parameters.index(parameter)] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
            elif parameter in ["Inclination"]:
                mcmc_samples_rescaled[:,i] = np.arccos((mcmc_samples[:,parameters.index(parameter)] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
            elif parameter in ['Amplitude',"FrequencyDerivative"]:
                mcmc_samples_rescaled[:,i] = 10**((mcmc_samples[:,parameters.index(parameter)] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
            else:
                mcmc_samples_rescaled[:,i] = (mcmc_samples[:,parameters.index(parameter)] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0]
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
        datS[:,2] = mcmc_samples_rescaled[:,3]*10**3
        datS[:,3] = np.log10(mcmc_samples_rescaled[:,4])
        datS[:,4] = mcmc_samples_rescaled[:,5]
        datS[:,5] = np.log10(mcmc_samples_rescaled[:,0])
        datS[:,6] = mcmc_samples_rescaled[:,6]
        datS[:,7] = mcmc_samples_rescaled[:,7]

        lbls = [r'\lambda', r'\beta', 'f(mHz)', r'\dot{f}', r'\iota', r'A', r'\phi', r'\Phi']

        tr_s = np.zeros(len(parameters))
        maxvalues = np.zeros(len(parameters))
        i = 0
        # pGB = search1.pGB
        for parameter in ['EclipticLongitude','EclipticLatitude','Frequency','FrequencyDerivative','Inclination','Amplitude','InitialPhase','Polarization']:
            if parameter in ['Amplitude','FrequencyDerivative']:
                tr_s[i] = np.log10(pGB[parameter])
                maxvalues[i] = np.log10(self.maxpGB[parameter])
            elif parameter in ['Frequency']:
                tr_s[i] = pGB[parameter]*10**3
                maxvalues[i] = self.maxpGB[parameter]*10**3
            else:
                tr_s[i] = pGB[parameter]
                maxvalues[i] = self.maxpGB[parameter]
            i += 1
        if tr_s[0] > np.pi:
            tr_s[0] -= 2*np.pi


        rng = []
        for i in range(len(lbls)):
            minrange = min(datS[:,i].min(), tr_s[i])
            maxrange = max(datS[:,i].max(), tr_s[i])
            range_width = np.abs(maxrange - minrange)
            oner = ( minrange- range_width/10, maxrange + range_width/10)
            rng.append(oner)
        # Get the getdist MCSamples objects for the samples, specifying same parameter
        # names and labels; if not specified weights are assumed to all be unity
        ndim = 6
        names = ['EclipticLongitude','EclipticLatitude','Frequency','FrequencyDerivative','Inclination','Amplitude']
        labels =  lbls[:6]
        samples = MCSamples(samples=datS[:,:6],names = names, labels = labels)

        g = plots.get_subplot_plotter()
        samples.updateSettings({'contours': [0.68, 0.95]})
        g.settings.num_plot_contours = 3
        g.triangle_plot([samples], shaded=True,title_limit=1)
        #markers vertical
        for i in range(ndim):
            for ax in g.subplots[i:,i]:
                ax.axvline(tr_s[i], color='orange', ls='--')
                ax.axvline(maxvalues[i], color='green', ls='--')
            i += 1
        #markers horizontal
        for i in range(ndim):
            for ax in g.subplots[i,:i]:
                ax.axhline(tr_s[i], color='orange', ls='--')
                ax.axhline(maxvalues[i], color='green', ls='--')
            i += 1
        save_frequency = self.frequencies[0]
        save_frequency = pGB['Frequency']
        if save_bool:
            g.export('/home/stefan/LDC/LDC/pictures/corner_frequency'+ str(int(np.round(save_frequency*10**8))) +'added'+save_name+'.png')
        if save_chain:
            df = pd.DataFrame(data=mcmc_samples_rescaled, columns=parameters)
            df.to_csv('/home/stefan/Repositories/ldc1_evaluation_data/submission/ETH_2/GW'+str(int(np.round(save_frequency*10**8)))+save_name+'.csv',index=False)
        plt.show()


def compute_posterior(tdi_fs, Tobs, frequencies, maxpGB, pGB_true):
    start = time.time()
    posterior1 = Posterior_computer(tdi_fs, Tobs, frequencies, maxpGB)
    posterior1.reduce_boundaries()
    posterior1.train_model()
    mcmc_samples = posterior1.calculate_posterior(resolution = 1*10**5, temperature= 10)
    # posterior1.plot_corner(mcmc_samples, pGB_injected[1][0])
    mcmc_samples = posterior1.calculate_posterior(resolution = 1*10**5, proposal= mcmc_samples, temperature= 2)
    # posterior1.plot_corner(mcmc_samples, pGB_injected[1][0])
    mcmc_samples = posterior1.calculate_posterior(resolution = 1*10**5, proposal= mcmc_samples, temperature= 1)
    mcmc_samples = posterior1.calculate_posterior(resolution = 5*10**5, proposal= mcmc_samples, temperature= 1)
    print('time to compute posterior: ', time.time()-start)
    posterior1.plot_corner(mcmc_samples, pGB_true, save_bool= True, save_chain= True)
    return mcmc_samples

# LDC1-3 ######################
start = time.time()
number_of_total_signals = 0
for i in range(len(found_sources_in)):
    for j in range(len(found_sources_in[i])):
        mcmc_samples = compute_posterior(tdi_fs, Tobs, frequencies[i], found_sources_in[i][j], pGB_injected[i][j])
        number_of_total_signals += 1
print('time to calculate posterior for ', number_of_total_signals, 'signals: ', time.time()-start)

# LDC1-4 ####################
# for i in range(len(found_sources_in)):
#     for j in range(len(found_sources_in[i])):
#         #subtract the found sources of neighbours and own window from original except the signal itself
#         tdi_fs_subtracted = deepcopy(tdi_fs)
#         for m in range(3):
#             if i-1+m < 0:
#                 pass
#             elif i-1+m > len(found_sources_in)-1:
#                 pass
#             else:
#                 for n in range(len(found_sources_in[i-1+m])):
#                     if j != n or m != 1:
#                         print(i,j,m,n)
#                         Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_in[i-1+m][n], oversample=4, simulator="synthlisa")
#                         source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
#                         index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
#                         index_high = index_low+len(Xs_subtracted)
#                         for k in ["X", "Y", "Z"]:
#                             tdi_fs_subtracted[k].data[index_low:index_high] = tdi_fs_subtracted[k].data[index_low:index_high] - source_subtracted[k].data

#         mcmc_samples = compute_posterior(tdi_fs_subtracted, Tobs, frequencies[i], found_sources_in[i][j], pGB_injected[i][j])

cmap = plt.get_cmap("Dark2")
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']


labels = {'EclipticLongitude': r'$\lambda$','EclipticLatitude': r'$\beta$','Frequency': '$f$ mHz','FrequencyDerivative': r'$\dot{f}$ Hz$^2$','Inclination': r'cos $\iota$','Amplitude': r'$A$','Polarization': r'$\psi$','InitialPhase': r'$\phi_0$'}

fig = plt.figure(figsize=fig_size_squared)
for n in range(3):
    ax1 = plt.subplot(3,1,n+1)
    if n == 0:
        parameter1 = 'EclipticLongitude'
        parameter2 = 'EclipticLatitude'
    if n == 1:
        parameter1 = 'Inclination'
        parameter2 = 'Amplitude'
    if n == 2:
        parameter1 = 'Frequency'
        parameter2 = 'FrequencyDerivative'
    for i in range(len(pGB_injected)):   
        for j in range(len( pGB_injected[i])):
            if parameter1 == 'Inclination':
                ax1.plot(np.cos(pGB_injected[i][j][parameter1]),pGB_injected[i][j][parameter2],color='lightgrey', marker = 'o', markersize = 8, zorder=1, label= 'true')
            else:
                ax1.plot(pGB_injected[i][j][parameter1],pGB_injected[i][j][parameter2],color='lightgrey', marker = 'o', markersize = 8, zorder=1, label= 'true')
    for i in range(len(found_sources_mp)):
        for j in range(len(found_sources_in[i])):
            save_frequency = pGB_injected[i][j]['Frequency']
            df = pd.read_csv('/home/stefan/Repositories/ldc1_evaluation_data/submission/ETH_2/GW'+str(int(np.round(save_frequency*10**8)))+save_name+'.csv')
            if parameter1 == 'Inclination':
                data1 = np.cos(df[parameter1][:10000])
            else:
                data1 = df[parameter1][:10000]
            # ax1.scatter(data1[:100000],df[parameter2][:100000], color = cmap(i), alpha = 0.03, marker='.', s = 0.05)
            ax1 = sns.kdeplot(x=data1[:1000], y=df[parameter2][:1000],fill=False, levels=[0.1,0.5])


            # ax = np.linspace(-0.15,1.15,numPoints)
            # mypdf,axes = fastKDE.pdf(proposal[:,1],proposal[:,2], axes=[ax,ax])
            # dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
            # data, pdfs = dist(self.resolution)
    
    if n > 0:
        ax1.set_yscale('log')
    ax1.xaxis.set_major_locator(plt.MaxNLocator(6))
    ax1.set_xlabel(labels[parameter1])
    ax1.set_ylabel(labels[parameter2])

    # for j in range(len(found_sources_in[i])):
    #     if i == 0:
    #         ax1.plot(found_sources_in[i][j][parameter1],found_sources_in[i][j][parameter2], color = colors[i], marker = '+', markersize = 8, label= 'global fit')
    #     else:
    #         ax1.plot(found_sources_in[i][j][parameter1],found_sources_in[i][j][parameter2], color = colors[i], marker = '+', markersize = 8)

plt.tight_layout()
fig.savefig('/home/stefan/LDC/LDC/pictures/global fit all '+save_name+'.png',dpi=300,bbox_inches='tight')
plt.show()




import matplotlib.font_manager

injected_frequencies = []
m = 0
for i in range(len(pGB_injected)):   
    for j in range(len( pGB_injected[i])):
        injected_frequencies.append(pGB_injected[i][j]['Frequency'])
        m += 1
pGB_injected_reshaped = np.reshape(pGB_injected, m)
def closest(list, Number):
    aux = []
    for valor in list:
        aux.append(abs(Number-valor))

    return aux.index(min(aux))

lbls = [ r'$\log A$', r'$\sin \beta$',r'$\lambda$', '$f - f_{True}$ (nHz)', '$\log \dot{f}$ (Hz/s)', r'$\cos \iota$', r'$\phi$', r'$\Phi$']
g = plots.get_subplot_plotter(subplot_size_ratio=9/16*0.7, subplot_size=8)
g.settings.line_styles = 'tab10'
g.settings.solid_colors='tab10'
boundaries = {
    "EclipticLatitude": [-1.0, 1.0],
    "EclipticLongitude": [-np.pi, np.pi],
    "Inclination": [-1.0, 1.0],
    "InitialPhase": [0.0, 2.0 * np.pi],
    "Polarization": [0.0, 1.0 * np.pi],
}
names = parameters
parameter_pairs = [['EclipticLongitude', 'EclipticLatitude'],['Inclination', 'Amplitude'],['Frequency', 'FrequencyDerivative']]
samples = []
pGB_injected_sorted_index = []
m = 0
for i in range(len(found_sources_mp)):
    for j in range(len(found_sources_in[i])):
        save_frequency = pGB_injected[i][j]['Frequency']
        df = pd.read_csv('/home/stefan/Repositories/ldc1_evaluation_data/submission/ETH_2/GW'+str(int(np.round(save_frequency*10**8)))+'.csv')
        df['Inclination'] = np.cos(df['Inclination'].values)
        df['EclipticLatitude'] = np.sin(df['EclipticLatitude'].values)
        df['FrequencyDerivative'] = np.log10(df['FrequencyDerivative'].values)
        df['Amplitude'] = np.log10(df['Amplitude'].values)
        pGB_injected_sorted_index.append(closest(injected_frequencies, df['Frequency'][0]))
        df['Frequency'] = (df['Frequency'] - pGB_injected_reshaped[pGB_injected_sorted_index[-1]]['Frequency'] + m*2e-9) * 1e9
        samples.append(MCSamples(samples=df.to_numpy(), names = names, labels = lbls))
        samples[-1].updateSettings({'contours': [0.68, 0.95]})
        m += 1
pGB_injected_sorted = []
for i in range(len(found_sources_mp)):
    pGB_injected_sorted.append(pGB_injected_reshaped[pGB_injected_sorted_index[i]])

g.settings.num_plot_contours = 2
# 3D (scatter) triangle plot
# you can adjust the scaling factor if font sizes are too small when
# making many subplots in a fixed size (default=2 would give smaller fonts)
g.settings.scaling_factor = 1
g.plots_2d(samples, param_pairs=parameter_pairs,legend_labels=[],lws=1.5)
for n, ax in enumerate(g.subplots[:,0]):
    parameter1, parameter2 = parameter_pairs[n]
    m = 0
    for i in range(len(pGB_injected_sorted)):   
        pGB_injected_scaled = deepcopy(pGB_injected_sorted[i])
        pGB_injected_scaled['Inclination'] = np.cos(pGB_injected_scaled['Inclination'])
        pGB_injected_scaled['EclipticLatitude'] = np.sin(pGB_injected_scaled['EclipticLatitude'])
        pGB_injected_scaled['FrequencyDerivative'] = np.log10(pGB_injected_scaled['FrequencyDerivative'])
        pGB_injected_scaled['Amplitude'] = np.log10(pGB_injected_scaled['Amplitude'])
        pGB_injected_scaled['Frequency'] = m * 2
        m += 1
        ax.plot(pGB_injected_scaled[parameter1],pGB_injected_scaled[parameter2],color='black', marker = '+',zorder=1, markersize = 10, label= 'true')
        ax.plot(pGB_injected_scaled[parameter1],pGB_injected_scaled[parameter2], marker = '+',zorder=1.1, markersize = 15,alpha = 0.5, label= 'true', linewidth = 4)
    try:
        ax.set_xlim(boundaries[parameter1])
    except:
        xlim = ax.get_xlim()
        x_length = xlim[1]-xlim[0]
        ax.set_xlim([xlim[0]-x_length*0.02, xlim[1]+x_length*0.02])
    try:
        ax.set_ylim(boundaries[parameter2])
    except:
        ylim = ax.get_ylim()
        y_length = ylim[1]-ylim[0]
        ax.set_ylim([ylim[0]-y_length*0.02, ylim[1]+y_length*0.02])
    if parameter2 in ['FrequencyDerivative']:
        ax.axhline(y=np.log10(1/Tobs**2/100), color='grey', linestyle = '--', zorder = 0.5)
    # if parameter2 in ['Amplitude', 'FrequencyDerivative']:
    #     ax.set_yscale('log')
g.export('/home/stefan/LDC/LDC/pictures/global fit all '+save_name+'log frequency.png')


stock = lambda A, amp, angle, phase: A * angle + amp * np.sin(angle + phase)
theta = np.linspace(0., 2 * np.pi, 250) # x-axis
np.random.seed(100)
noise = 0.2 * np.random.random(250)
y = stock(.1, .2, theta, 1.2) + noise # y-axis
fig = plt.figure(figsize=fig_size)
plt.subplots_adjust(bottom = 0., left = 0, top = 1., right = 1,wspace= 0.1)
# Create first axes, the top-left plot with green plot
sub1 = fig.add_subplot(2,2,1) # two rows, two columns, fist cell
sub1.plot(theta, y, color = 'green')
sub1.set_xlim(1, 2)
sub1.set_ylim(0, 1)
sub1.set_ylabel('y', labelpad = 15)
# Create second axes, the top-left plot with orange plot
sub2 = fig.add_subplot(2,2,2) # two rows, two columns, second cell
sub2.plot(theta, y, color = 'orange')
sub2.set_xlim(5, 6)
sub2.set_ylim(0, 1)
sub2.yaxis.set_ticks([])
# Create third axes, a combination of third and fourth cell
sub3 = fig.add_subplot(2,2,(3,4)) # two rows, two colums, combined third and fourth cell
sub3.plot(theta, y, color = 'darkorchid', alpha = .7)
sub3.set_xlim(0, 6.5)
sub3.set_ylim(0, 1)
sub3.set_xlabel(r'$\theta$ (rad)', labelpad = 15)
sub3.set_ylabel('y', labelpad = 15)

sub1.spines['right'].set_visible(False)
sub2.spines['left'].set_visible(False)

d = .015 # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=sub1.transAxes, color='k', clip_on=False)
sub1.plot((1-d,1+d), (-d,+d), **kwargs)
sub1.plot((1-d,1+d),(1-d,1+d), **kwargs)
kwargs.update(transform=sub2.transAxes)  # switch to the bottom axes
sub2.plot((-d,+d), (1-d,1+d), **kwargs)
sub2.plot((-d,+d), (-d,+d), **kwargs)
# Create blocked area in third axes
sub3.fill_between((1,2), 0, 1, facecolor='green', alpha=0.2) # blocked area for first axes
sub3.fill_between((5,6), 0, 1, facecolor='orange', alpha=0.2) # blocked area for second axes
# Create left side of Connection patch for first axes
con1 = ConnectionPatch(xyA=(1, .2), coordsA=sub1.transData, xyB=(1, .3), coordsB=sub3.transData, color = 'green')
fig.add_artist(con1)
con2 = ConnectionPatch(xyA=(2, .2), coordsA=sub1.transData, xyB=(2, .3), coordsB=sub3.transData, color = 'green')
fig.add_artist(con2)
con3 = ConnectionPatch(xyA=(5, .4), coordsA=sub2.transData, xyB=(5, .5), coordsB=sub3.transData, color = 'orange')
fig.add_artist(con3)
con4 = ConnectionPatch(xyA=(6, .4), coordsA=sub2.transData, xyB=(6, .9), coordsB=sub3.transData, color = 'orange')
fig.add_artist(con4)
plt.savefig('/home/stefan/LDC/LDC/pictures/zoom_effect_2.png', dpi = 300, bbox_inches = 'tight', pad_inches = .1)

# compute_posterior_input = [[tdi_fs, Tobs, frequencies[1], found_sources_in[1][0]]]
# start = time.time()
# pool = mp.Pool(mp.cpu_count())
# found_sources_mp= pool.starmap(compute_posterior, compute_posterior_input)
# pool.close()
# pool.join()
# print('time to search ', number_of_windows, 'windows: ', time.time()-start)
# %%
