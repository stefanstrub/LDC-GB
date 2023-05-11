#%%
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Ellipse
from getdist import plots, MCSamples
import scipy
from scipy.optimize import differential_evolution
import numpy as np
import xarray as xr
# from getdist import plots, MCSamples
import time
from copy import deepcopy
import multiprocessing as mp
import pandas as pd
import os
import h5py
from KDEpy import FFTKDE
import sys
sys.path.append('/cluster/home/sstrub/Repositories/LDC/lib/lib64/python3.8/site-packages/ldc-0.1-py3.8-linux-x86_64.egg')

from ldc.lisa.noise import get_noise_model
from ldc.lisa.noise import AnalyticNoise
from ldc.common.series import TimeSeries
import ldc.waveform.fastGB as fastGB
from ldc.common.tools import window
# from ldc.waveform.fastGB import fastGB
# import ldc.waveform.fastGB as fastGB

# from ldc.common.tools import compute_tdi_snr

from fastkde import fastKDE
# from sklearn.metrics import mean_squared_error
# from sklearn.gaussian_process import GaussianProcessRegressor
# from sklearn.gaussian_process.kernels import RBF
# from chainconsumer import ChainConsumer


try:
    import cupy as xp

    gpu_available = True

except (ImportError, ModuleNotFoundError) as e:
    import numpy as xp

    gpu_available = False

from gbgpu.gbgpu import GBGPU

from gbgpu.utils.constants import *

from sources import *

import subprocess as sp
import os

def get_gpu_memory():
    command = "nvidia-smi --query-gpu=memory.free --format=csv"
    memory_free_info = sp.check_output(command.split()).decode('ascii').split('\n')[:-1][1:]
    memory_free_values = [int(x.split()[0]) for i, x in enumerate(memory_free_info)]
    return memory_free_values

gb_gpu = GBGPU(use_gpu=gpu_available)

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
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

# def window(tm, offs=1000.0):
#     xl = offs
#     ind_r = np.argwhere(tm[-1] - tm <= offs)[0][0]
#     xr = tm[ind_r]
#     kap = 0.005
#     winl = 0.5 * (1.0 + np.tanh(kap * (tm - xl)))
#     winr = 0.5 * (1.0 - np.tanh(kap * (tm - xr)))
#     return winl * winr
def pGB_dict_to_gpu_input(pGB):
    params = []
    for i in range(len(pGB)):
        params.append(
            [pGB[i]['Amplitude'],
            pGB[i]['Frequency'],
            pGB[i]['FrequencyDerivative'],
            0,
            pGB[i]['InitialPhase'],
            pGB[i]['Inclination'],
            pGB[i]['Polarization'],
            pGB[i]['EclipticLongitude'],
            pGB[i]['EclipticLatitude']]
        )
    params = np.asanyarray(params)
    params = np.swapaxes(params, 0,1)
    return params


# def get_noise_from_frequency_domain(tdi_fs, number_of_windows=100):
#     tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft()) for k in ["X", "Y", "Z"]]))
#     tdi_ts["E"] = deepcopy(tdi_ts["X"])
#     tdi_ts["A"] = (tdi_ts["Z"] - tdi_ts["X"])/np.sqrt(2.0)
#     tdi_ts["E"].values = (tdi_ts["Z"] - 2.0*tdi_ts["Y"] + tdi_ts["X"])/np.sqrt(6.0)
#     tdi_ts["T"] = (tdi_ts["Z"] + tdi_ts["Y"] + tdi_ts["X"])/np.sqrt(3.0)
#     f, psdA =  scipy.signal.welch(tdi_ts["A"], fs=1.0/dt, nperseg=len(tdi_ts["X"])/number_of_windows, average='mean', window= 'boxcar')
#     f, psdE =  scipy.signal.welch(tdi_ts["E"], fs=1.0/dt, nperseg=len(tdi_ts["X"])/number_of_windows)
#     # f2, psdE2 =  scipy.signal.welch(tdi_ts["E"], fs=1.0/dt, nperseg=len(tdi_ts["X"]), scaling='spectrum')
#     f, psdT =  scipy.signal.welch(tdi_ts["T"], fs=1.0/dt, nperseg=len(tdi_ts["X"])/number_of_windows)
#     return f, psdA, psdE, psdT

# def median_windows(y, window_size):
#     medians = deepcopy(y)
#     for i in range(int(len(y)/window_size*2)-1):
#         start_index = int(i/2*window_size)
#         end_index = int((i/2+1)*window_size)
#         median = np.median(y[start_index:end_index])
#         outliers = np.abs(medians[start_index:end_index]) > median
#         medians[start_index:end_index][outliers] = median
#     return medians

# def smooth_psd(psd, f):
#     smoothed = median_windows(psd, 20)
#     smoothed[:40] = psd[:40]
#     index_cut = np.searchsorted(f, 0.0008)
#     index_cut_lower = np.searchsorted(f, 3*10**-4)
#     psd_fit = np.ones_like(smoothed)
#     psd_fit_low = scipy.signal.savgol_filter(smoothed, 10, 1)
#     psd_fit_high = scipy.signal.savgol_filter(smoothed, 70, 1)
#     psd_fit[:index_cut] = psd_fit_low[:index_cut] 
#     psd_fit[index_cut:] = psd_fit_high[index_cut:] 
#     # psd_fit[:index_cut+100] = scipy.signal.savgol_filter(smoothed[:index_cut+100], 10, 1)
#     # psd_fit[index_cut-2:] = scipy.signal.savgol_filter(smoothed[index_cut-2:], 70, 1)
#     psd_fit[:index_cut_lower] = smoothed[:index_cut_lower]
#     return psd_fit

# def objective(trial,changeableparameters, maxpGB2, signal):
#     for parameter in changeableparameters:
#         parametervalue = trial.suggest_uniform(parameter, boundaries[parameter][0], boundaries[parameter][1])
#         maxpGB2[signal][parameter] = parametervalue
#         if parameter in ["EclipticLatitude"]:
#             maxpGB2[signal][parameter] = np.arcsin(parametervalue)
#         elif parameter in ["Inclination"]:
#             maxpGB2[signal][parameter] = np.arccos(parametervalue)
#         elif parameter in parameters_log_uniform:
#             maxpGB2[signal][parameter] = 10**(parametervalue)
#         # elif parameter in ['Amplitude']:
#         #     maxpGB2[signal][parameter] = 10**(parametervalue)/100
#     p = loglikelihood(maxpGB2)
#     return p

# def scaletooriginal_array(previous_pGB, boundaries):
#     start = time.time()
#     no_array = False
#     if len(np.shape(previous_pGB)) == 1:
#         no_array = True
#         previous_pGB = np.array([previous_pGB])
#     original_pGB = np.zeros(np.shape(previous_pGB))
#     i = 0
#     for parameter in parameters:
#         if parameter in ["EclipticLatitude"]:
#             original_pGB[:,i] = np.arcsin((previous_pGB[:,parameters.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
#         elif parameter in ["Inclination"]:
#             original_pGB[:,i] = np.arccos((previous_pGB[:,parameters.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
#         elif parameter in parameters_log_uniform:
#             original_pGB[:,i] = 10**((previous_pGB[:,parameters.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
#         else:
#             original_pGB[:,i] = (previous_pGB[:,parameters.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]
#         i += 1
#     if no_array:
#         original_pGB = original_pGB[0]
#     # print('time rescale', time.time()-start)
#     return original_pGB
    
# def scaletooriginal(previous_pGB, boundaries):
#     original_pGB = {}
#     for parameter in parameters:
#         if parameter in ["EclipticLatitude"]:
#             original_pGB[parameter] = np.arcsin((previous_pGB[parameters.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
#         elif parameter in ["Inclination"]:
#             original_pGB[parameter] = np.arccos((previous_pGB[parameters.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
#         elif parameter in parameters_log_uniform:
#             original_pGB[parameter] = 10**((previous_pGB[parameters.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
#         else:
#             original_pGB[parameter] = (previous_pGB[parameters.index(parameter)] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]
#     return original_pGB

# def scaletooriginalparameter(previous_max, boundaries):
#     maxpGB = {}
#     for parameter in parameters:
#         if parameter in ["EclipticLatitude"]:
#             maxpGB[parameter] = np.arcsin((previous_max[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
#         elif parameter in ["Inclination"]:
#             maxpGB[parameter] = np.arccos((previous_max[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
#         elif parameter in parameters_log_uniform:
#             maxpGB[parameter] = 10**((previous_max[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0])
#         else:
#             maxpGB[parameter] = (previous_max[parameter] * (boundaries[parameter][1] - boundaries[parameter][0])) + boundaries[parameter][0]
#     return maxpGB

# def scaleto01(previous_max, boundaries):
#     maxpGB = {}
#     for parameter in parameters:
#         if parameter in ["EclipticLatitude"]:
#             maxpGB[parameter] = (np.sin(previous_max[parameter]) - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0])
#         elif parameter in ["Inclination"]:
#             maxpGB[parameter] = (np.cos(previous_max[parameter]) - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0])
#         elif parameter in parameters_log_uniform:
#             maxpGB[parameter] = (np.log10(previous_max[parameter]) - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0])
#         else:
#             maxpGB[parameter] = (previous_max[parameter] - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0])
#     return maxpGB

# def reduce_boundaries(maxpGB, boundaries, ratio=0.1):
#     boundaries_reduced = deepcopy(boundaries)
#     for parameter in parameters:
#         length = boundaries[parameter][1] - boundaries[parameter][0]
#         if parameter == "EclipticLatitude":
#             boundaries_reduced[parameter] = [
#                 np.sin(maxpGB[parameter]) - length * ratio / 2,
#                 np.sin(maxpGB[parameter]) + length * ratio / 2,
#             ]
#         elif parameter == "Inclination":
#             boundaries_reduced[parameter] = [
#                 np.cos(maxpGB[parameter]) - length * ratio / 2,
#                 np.cos(maxpGB[parameter]) + length * ratio / 2,
#             ]
#         elif parameter == "Frequency":
#             boundaries_reduced[parameter] = [
#                 maxpGB[parameter] - length * ratio / 2,
#                 maxpGB[parameter] + length * ratio / 2,
#             ]
#         elif parameter == "FrequencyDerivative":
#             boundaries_reduced[parameter] = [
#                 maxpGB[parameter] - length * ratio / 2,
#                 maxpGB[parameter] + length * ratio / 2,
#             ]
#         elif parameter == "Amplitude":
#             boundaries_reduced[parameter] = [
#                 np.log10(maxpGB[parameter]) - length * ratio / 2,
#                 np.log10(maxpGB[parameter]) + length * ratio / 2,
#             ]
#         elif parameter in ["InitialPhase",'Polarization']:
#             boundaries_reduced[parameter] = [
#                 maxpGB[parameter] - length * ratio / 2,
#                 maxpGB[parameter] + length * ratio / 2,
#             ]
#         elif parameter == "EclipticLongitude":
#             boundaries_reduced[parameter] = [
#                 maxpGB[parameter] - length * ratio / 2,
#                 maxpGB[parameter] + length * ratio / 2,
#             ]
#         if boundaries_reduced[parameter][0] < boundaries[parameter][0]:
#             boundaries_reduced[parameter][0] = boundaries[parameter][0]
#         if boundaries_reduced[parameter][1] > boundaries[parameter][1]:
#             boundaries_reduced[parameter][1] = boundaries[parameter][1]
#     return boundaries_reduced

# def CoordinateMC(n, pGBs, boundaries, parameters_recorded, loglikelihood, n_trials=50):
#     maxpGB = []
#     parameters_recorded1 = []
#     no_improvement_counter = 0
#     for i in range(number_of_signals):
#         maxpGB.append(deepcopy(pGBs))
#         parameters_recorded1.append({})
#         for parameter in parameters_no_amplitude:
#             parametervalue = np.random.uniform(boundaries[parameter][0], boundaries[parameter][1])
#             maxpGB[i][parameter] = parametervalue
#             if parameter in ["EclipticLatitude"]:
#                 maxpGB[i][parameter] = np.arcsin(parametervalue)
#             elif parameter in ["Inclination"]:
#                 maxpGB[i][parameter] = np.arccos(parametervalue)
#             elif parameter in ['Amplitude',"FrequencyDerivative"]:
#                 maxpGB[i][parameter] = 10**(parametervalue)
#             parameters_recorded1[i][parameter] = []
#             parameters_recorded1[i][parameter].append(maxpGB[i][parameter])
#         parameters_recorded1[i]["Loglikelihood"] = []
#     # if n == 0:
#     #     maxpGB = deepcopy(pGBs)
#     try:
#         best_value
#     except:
#         best_value = loglikelihood(maxpGB)
#     previous_best = deepcopy(best_value)
#     maxpGB2 = deepcopy(maxpGB)
#     for i in range(100):
#         # if i > 48:
#         #     n_trials = 50
#         parameter1 = parameters_no_amplitude[i % 7]
#         parameter2 = parameters_no_amplitude[np.random.randint(0, 6)]
#         parameter3 = parameters_no_amplitude[np.random.randint(0, 6)]
#         # parameter2 = 'InitialPhase'
#         # parameter1 = 'Inclination'
#         while parameter2 == parameter1:
#             parameter2 = parameters_no_amplitude[np.random.randint(0, 6)]
#         while parameter3 == parameter1 or parameter3 == parameter2:
#             parameter3 = parameters_no_amplitude[np.random.randint(0, 6)]
#         # if parameter1 == 'Frequency':
#         #     parameter2 = 'Polarization'
#         changeableparameters = [parameter1, parameter2]
#         signal = i % number_of_signals
#         for j in range(n_trials):
#             k = 0
#             for parameter in changeableparameters:
#                 parametervalue = np.random.uniform(boundaries[parameter][0], boundaries[parameter][1],1)[0]
#                 # if k == 0:
#                 #     parametervalue = np.random.uniform(j%number_per_row * 1/number_per_row, (j%number_per_row+1)*1/number_per_row)
#                 # else:
#                 #     parametervalue = np.random.uniform(int(j/number_per_row) * 1/number_per_row, (int(j/number_per_row)+1)*1/number_per_row)
#                 # parametervalue = parametervalue * (boundaries[parameter][1] - boundaries[parameter][0]) + boundaries[parameter][0]
#                 maxpGB2[signal][parameter] = parametervalue
#                 if parameter in ["EclipticLatitude"]:
#                     maxpGB2[signal][parameter] = np.arcsin(parametervalue)
#                 elif parameter in ["Inclination"]:
#                     maxpGB2[signal][parameter] = np.arccos(parametervalue)
#                 elif parameter in ['Amplitude',"FrequencyDerivative"]:
#                     maxpGB2[signal][parameter] = 10**(parametervalue)
#                 k += 1
#             suggestion = loglikelihood(maxpGB2)
#             if suggestion > previous_best:
#                 no_improvement_counter = 0
#                 previous_best = suggestion
#                 for parameter in changeableparameters:
#                     maxpGB[signal][parameter] = maxpGB2[signal][parameter]
#         if no_improvement_counter != 0:
#             no_improvement_counter += 1
#         if no_improvement_counter > 16:
#             past_mean = 0
#             sum_count = 0
#             for l in range(n):
#                 try:
#                     past_mean += parameters_recorded[l]["Loglikelihood"][i]
#                     sum_count += 1
#                 except:
#                     pass
#             try:
#                 past_mean = past_mean / sum_count
#                 if previous_best > past_mean:
#                     no_improvement_counter = 0
#                 else:
#                     # print("no improvement")
#                     break
#             except:
#                 # no_improvement_counter = 0
#                 break

#         if i in [30,40,50,60,70]:
#             past_mean = 0
#             sum_count = 0
#             for l in range(n):
#                 try:
#                     past_mean += parameters_recorded[l][0]["Loglikelihood"][i]
#                     sum_count += 1
#                 except:
#                     pass
#             try:
#                 past_mean = past_mean / sum_count
#                 if previous_best > past_mean:
#                     pass
#                 else:
#                     break
#             except:
#                 pass

#         # start = time.time()
#         # print(n, i, previous_best, loglikelihood(maxpGB), maxpGB)
#         parameters_recorded1[0]["Loglikelihood"].append(loglikelihood(maxpGB))
#         for i in range(number_of_signals):
#             for parameter in parameters_no_amplitude:
#                 parameters_recorded1[i][parameter].append(maxpGB[i][parameter])

#         maxpGB2 = deepcopy(maxpGB)
#     if previous_best > best_value:
#         best_value = previous_best
#     parameters_recorded[n] = parameters_recorded1
#     return parameters_recorded1

# class Distribution(object):
#     """
#     draws samples from a one dimensional probability distribution,
#     by means of inversion of a discrete inverstion of a cumulative density function

#     the pdf can be sorted first to prevent numerical error in the cumulative sum
#     this is set as default; for big density functions with high contrast,
#     it is absolutely necessary, and for small density functions,
#     the overhead is minimal

#     a call to this distibution object returns indices into density array
#     """
#     def __init__(self, pdf, sort = True, interpolation = True, transform = lambda x: x):
#         self.shape          = pdf.shape
#         self.pdf            = pdf.ravel()
#         self.sort           = sort
#         self.interpolation  = interpolation
#         self.transform      = transform

#         #a pdf can not be negative
#         assert(np.all(pdf>=0))

#         #sort the pdf by magnitude
#         if self.sort:
#             self.sortindex = np.argsort(self.pdf, axis=None)
#             self.pdf = self.pdf[self.sortindex]
#         #construct the cumulative distribution function
#         self.cdf = np.cumsum(self.pdf)
#     @property
#     def ndim(self):
#         return len(self.shape)
#     @property
#     def sum(self):
#         """cached sum of all pdf values; the pdf need not sum to one, and is imlpicitly normalized"""
#         return self.cdf[-1]
#     def __call__(self, N):
#         """draw """
#         #pick numbers which are uniformly random over the cumulative distribution function
#         choice = np.random.uniform(high = self.sum, size = N)
#         #find the indices corresponding to this point on the CDF
#         index = np.searchsorted(self.cdf, choice)
#         pdfs = self.pdf[index]
#         #if necessary, map the indices back to their original ordering
#         if self.sort:
#             index = self.sortindex[index]
#         #map back to multi-dimensional indexing
#         index = np.unravel_index(index, self.shape)
#         index = np.vstack(index)
#         #is this a discrete or piecewise continuous distribution?
#         if self.interpolation:
#             index = index + np.random.uniform(size=index.shape)
#         return self.transform(index), pdfs

# def moving_average(a, n=3) :
#     ret = np.cumsum(a, dtype=float)
#     ret[n:] = ret[n:] - ret[:-n]
#     return ret[n - 1:] / n

# class Search():
#     def __init__(self,tdi_fs,Tobs, lower_frequency, upper_frequency, recombination=0.75, noise=None):
#         self.tdi_fs = tdi_fs
#         self.GB = fastGB.FastGB(delta_t=dt, T=Tobs)
#         self.reduced_frequency_boundaries = None
#         self.recombination = recombination
#         self.lower_frequency = lower_frequency
#         self.upper_frequency = upper_frequency
#         self.padding = (upper_frequency - lower_frequency)/2
#         # self.tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft()) for k in ["X", "Y", "Z"]]))
#         # psd = psdX + psdY + psdZ
#         # psd = psd[indexes]
#         self.frequency_T_threshold = 19.1*10**-3/2
        
#         self.use_T_component = False
#         if self.upper_frequency + self.padding > self.frequency_T_threshold:
#             self.use_T_component = True
#         # plt.figure()
#         # plt.plot(f[indexes], psdX[indexes])
#         # plt.plot(f[indexes], psd)
#         # plt.show()

#         frequencyrange =  [self.lower_frequency-self.padding, self.upper_frequency+self.padding]
#         indexes = np.logical_and(tdi_fs['X'].f > frequencyrange[0], tdi_fs['X'].f < frequencyrange[1]) 
#         self.dataX = tdi_fs["X"][indexes]
#         self.dataY = tdi_fs["Y"][indexes]
#         self.dataZ = tdi_fs["Z"][indexes]
#         self.dataX_full_f = tdi_fs["X"]
#         self.dataY_full_f = tdi_fs["Y"]
#         self.dataZ_full_f = tdi_fs["Z"]

#         self.DAf_full_f = (self.dataZ_full_f - self.dataX_full_f)/np.sqrt(2.0)
#         self.DEf_full_f = (self.dataZ_full_f - 2.0*self.dataY_full_f + self.dataX_full_f)/np.sqrt(6.0)
#         self.DTf_full_f = (self.dataZ_full_f + self.dataY_full_f + self.dataX_full_f)/np.sqrt(3.0)
#         self.DAf = self.DAf_full_f[indexes]
#         self.DEf = self.DEf_full_f[indexes]
#         self.DTf = self.DTf_full_f[indexes]
#         self.df = self.DAf.f[1].values-self.DAf.f[0].values

#         # self.tdi_ts["E"] = deepcopy(self.tdi_ts["X"])
#         # self.tdi_ts["A"] = (self.tdi_ts["Z"] - self.tdi_ts["X"])/np.sqrt(2.0)
#         # self.tdi_ts["E"].values = (self.tdi_ts["Z"] - 2.0*self.tdi_ts["Y"] + self.tdi_ts["X"])/np.sqrt(6.0)
#         # self.tdi_ts["T"] = (self.tdi_ts["Z"] + self.tdi_ts["Y"] + self.tdi_ts["X"])/np.sqrt(3.0)
#         # f, self.psdA =  scipy.signal.welch(self.tdi_ts["A"], fs=1.0/dt, nperseg=len(self.tdi_ts["X"]))
#         # f, self.psdE =  scipy.signal.welch(self.tdi_ts["E"], fs=1.0/dt, nperseg=len(self.tdi_ts["X"]))
#         # f2, self.psdE2 =  scipy.signal.welch(self.tdi_ts["E"], fs=1.0/dt, nperseg=len(self.tdi_ts["X"]), scaling='spectrum')
#         # f, self.psdT =  scipy.signal.welch(self.tdi_ts["T"], fs=1.0/dt, nperseg=len(self.tdi_ts["X"]))
#         # indexes = np.logical_and(f > frequencyrange[0], f < frequencyrange[1]) 
#         # indexes2 = np.logical_and(f2 > frequencyrange[0], f2 < frequencyrange[1]) 
#         # self.f_noise = f[indexes]
#         # self.f_noise2 = f2[indexes2]
#         # self.psdA = self.psdA[indexes]
#         # self.psdE = self.psdE[indexes]
#         # self.psdE2 = self.psdE2[indexes2]
#         # self.psdT = self.psdT[indexes]
#         # self.tdifE = self.tdi_ts["E"].ts.fft(win=window)
#         # # self.tdifE =  scipy.fft.fft(self.tdi_ts["E"].values)
#         # # f =  scipy.fft.fftshift(scipy.fft.fftfreq(self.tdi_ts["E"].t.values.shape[-1], d=1/self.tdi_ts["X"].t[1].values))
#         # indexes = np.logical_and(self.tdifE.f > frequencyrange[0], self.tdifE.f < frequencyrange[1]) 
#         # self.tdifE = self.tdifE[indexes]

#         # plt.figure()
#         # plt.plot(f[indexes], np.abs(self.DAf))
#         # plt.plot(f[indexes], np.abs(self.DEf))
#         # plt.plot(f[indexes], np.abs(self.DAf)+np.abs(self.DEf))
#         # plt.show()
#         # print('frequencyrange',frequencyrange)
#         if noise is None:
#             fmin, fmax = float(self.dataX.f[0]), float(self.dataX.f[-1] + self.dataX.attrs["df"])
#             freq = np.array(self.dataX.sel(f=slice(fmin, fmax)).f)
#             Nmodel = get_noise_model(noise_model, freq)
#             # self.Sn_full_f = Nmodel.psd(freq=self.DAf_full_f.f, option="X")
#             self.SA_full_f = Nmodel.psd(freq=self.DAf_full_f.f, option="A")
#             self.SE_full_f = Nmodel.psd(freq=self.DEf_full_f.f, option="E")
#             # self.ST_full_f = Nmodel.psd(freq=self.DAf_full_f.f, option="T")
#             self.Sn = Nmodel.psd(freq=freq, option="X")
#             self.SA = Nmodel.psd(freq=freq, option="A")
#             # self.SE = Nmodel.psd(freq=freq, option="E")
#             self.ST = Nmodel.psd(freq=freq, option="T")
#         else:
#             self.SA_full_f = noise['A']
#             self.SE_full_f = noise['E']
#             self.ST_full_f = noise['T']
#             fmin, fmax = float(self.dataX.f[0]), float(self.dataX.f[-1] + self.dataX.attrs["df"])
#             freq = np.array(self.dataX.sel(f=slice(fmin, fmax)).f)
#             self.SA = self.SA_full_f[indexes]
#             self.SE = self.SE_full_f[indexes]
#             self.ST = self.ST_full_f[indexes]

#         self.data_GPU = [xp.array(self.DAf_full_f),
#                 xp.array(self.DEf_full_f),
#         ]
        
#         self.PSD_GPU =  [xp.array(self.SA_full_f),
#                 xp.array(self.SA_full_f),
#         ]

#         # self.dd =  4.0*self.dataX.df*np.sum((np.absolute(self.DAf.data)**2 + np.absolute(self.DEf.data)**2) /self.SA)

#         f_0 = fmin
#         f_transfer = 19.1*10**-3
#         snr = 4
#         amplitude_lower = 2*snr/(Tobs * np.sin(f_0/ f_transfer)**2/self.SA[0])**0.5
#         snr = 1100
#         amplitude_upper = 2*snr/(Tobs * np.sin(f_0/ f_transfer)**2/self.SA[0])**0.5
#         amplitude = [amplitude_lower, amplitude_upper]
#         # print('lower frequency', lower_frequency)
#         # print('amplitude boundaries', amplitude)
#         # print('amplitude boundaries previous', np.sqrt(np.max(psd))/1000,np.sqrt(np.max(psd))/10)
#         # print('amplitude boundaries previous', np.max(np.abs(self.DAf.data))/10**7,np.max(np.abs(self.DAf.data))/10**5)
#         # amplitude = np.sqrt(np.max(psd))

#         # indexes = np.argsort(p.get('Frequency'))
#         # index_low = np.searchsorted(p.get('Frequency')[indexes], frequencyrange[0])
#         # index_high = np.searchsorted(p.get('Frequency')[indexes], frequencyrange[1])

#         # pGB = deepcopy(pGBadded)
#         # self.pGB = {'Amplitude': 3.676495e-22, 'EclipticLatitude': 0.018181, 'EclipticLongitude': 1.268061, 'Frequency': 0.01158392, 'FrequencyDerivative': 8.009579e-15, 'Inclination': 0.686485, 'InitialPhase': 4.201455, 'Polarization': 2.288223}
        
#         # fd_range = [np.log10(frequency_derivative(lower_frequency,0.1)),np.log10(frequency_derivative(lower_frequency,M_chirp_upper_boundary))]
#         # fd_range = [frequency_derivative(lower_frequency,0.1),frequency_derivative(lower_frequency,M_chirp_upper_boundary)]
#         # fd_range = [frequency_derivative_tyson_lower(lower_frequency),frequency_derivative_tyson(lower_frequency)]
#         fd_range = [frequency_derivative_tyson_lower(lower_frequency),frequency_derivative(upper_frequency,M_chirp_upper_boundary)]
#         self.boundaries = {
#             "Amplitude": [np.log10(amplitude[0]),np.log10(amplitude[1])],
#             # "Amplitude": [-23.5,-21],
#             # "Amplitude": [np.log10(self.pGB['Amplitude'])-2,np.log10(self.pGB['Amplitude'])+1],
#             "EclipticLatitude": [-1.0, 1.0],
#             "EclipticLongitude": [-np.pi, np.pi],
#             # "Frequency": [self.pGB["Frequency"] * 0.99995, self.pGB["Frequency"] * 1.00015],
#             # "Frequency": [self.pGB["Frequency"] - 3e-7, self.pGB["Frequency"] + 3e-7],
#             "Frequency": frequencyrange,
#             "FrequencyDerivative": fd_range,
#             # "FrequencyDerivative": [np.log10(5e-6*self.pGB['Frequency']**(13/3)),np.log10(8e-8*self.pGB['Frequency']**(11/3))],
#             "Inclination": [-1.0, 1.0],
#             "InitialPhase": [0.0, 2.0 * np.pi],
#             "Polarization": [0.0, 1.0 * np.pi],
#         }
#         # print('fd_rage',self.boundaries['FrequencyDerivative'], 'at f=', lower_frequency)

#         # print('fd * Tobs', 10**self.boundaries['FrequencyDerivative'][1] * Tobs)
#         # print('smear f', 300*lower_frequency * 10**3 / 10**9)

#         if self.boundaries['FrequencyDerivative'][0] > self.boundaries['FrequencyDerivative'][1]:
#             c = self.boundaries['FrequencyDerivative'][0]
#             self.boundaries['FrequencyDerivative'][0] = self.boundaries['FrequencyDerivative'][1]
#             self.boundaries['FrequencyDerivative'][1] = c
        
#         previous_max = np.random.rand(8)
#         # previous_max[0] = np.random.rand(1)*0.1 +0.5
#         # previous_max[3] = np.random.rand(1)*0.1 +0.5
#         i = 0
#         self.pGBs = {}
#         for parameter in parameters:
#             # if parameter in ["FrequencyDerivative"]:
#             #     i -= 1
#             if parameter in ["EclipticLatitude"]:
#                 self.pGBs[parameter] = np.arcsin((previous_max[i] * (self.boundaries[parameter][1] - self.boundaries[parameter][0])) + self.boundaries[parameter][0])
#             elif parameter in ["Inclination"]:
#                 self.pGBs[parameter] = np.arccos((previous_max[i] * (self.boundaries[parameter][1] - self.boundaries[parameter][0])) + self.boundaries[parameter][0])
#             elif parameter in parameters_log_uniform:
#                 self.pGBs[parameter] = 10**((previous_max[i] * (self.boundaries[parameter][1] - self.boundaries[parameter][0])) + self.boundaries[parameter][0])
#             else:
#                 self.pGBs[parameter] = (previous_max[i] * (self.boundaries[parameter][1] - self.boundaries[parameter][0])) + self.boundaries[parameter][0]
#             i += 1
#         # self.pGBs = {'Amplitude': 4.0900673126042746e-22, 'EclipticLatitude': 0.8718477251317046, 'EclipticLongitude': 0.48599945403230693, 'Frequency': 0.003995220986111426, 'FrequencyDerivative': 1.0985841703423861e-16, 'Inclination': 1.0262955111380103, 'InitialPhase': 5.453865686076588, 'Polarization': 1.089057196561609}

#         # cutoff_ratio = 1000
#         # Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=self.pGBs, oversample=4)
#         # psd_signal = np.abs(Xs.values) ** 2 + np.abs(Ys.values) ** 2 + np.abs(Zs.values) ** 2
#         # highSNR = psd_signal > np.max(psd_signal) / cutoff_ratio
#         # lowerindex = np.where(highSNR)[0][0] - 30
#         # higherindex = np.where(highSNR)[0][-1] + 30
#         # self.dataX = tdi_fs["X"].isel(f=slice(Xs.kmin, Xs.kmin + len(Xs)))[lowerindex:higherindex]
#         # self.dataY = tdi_fs["Y"].isel(f=slice(Ys.kmin, Ys.kmin + len(Ys)))[lowerindex:higherindex]
#         # self.dataZ = tdi_fs["Z"].isel(f=slice(Zs.kmin, Zs.kmin + len(Zs)))[lowerindex:higherindex]

#         # self.DAf = (2/3*self.dataX - self.dataY - self.dataZ)/3.0
#         # self.DEf = (self.dataZ - self.dataY)/np.sqrt(3.0)

#         # Xs, Ys, Zs = (
#         #     Xs[lowerindex:higherindex],
#         #     Ys[lowerindex:higherindex],
#         #     Zs[lowerindex:higherindex],
#         # )

#         # diff = np.abs(self.dataX - Xs.values) ** 2 + np.abs(self.dataY - Ys.values) ** 2 + np.abs(self.dataZ - Zs.values) ** 2
#         # p1 = float(np.sum(diff / (self.Sn + noise)) * Xs.df) / 2.0
#         # p1 = -p1
#         # diff = np.abs(self.dataX) ** 2 + np.abs(self.dataY) ** 2 + np.abs(self.dataZ) ** 2
#         # null_pGBs = deepcopy(self.pGBs)
#         # null_pGBs['Amplitude'] = 4*10**-25
#         # print('pGB', self.pGB, self.loglikelihood([self.pGB]))

#     def update_noise(self, pGB=None):
#         if pGB != None:
#             Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB, oversample=4)
#             Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#             Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#             Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]

#             Af = (Zs_total - Xs_total)/np.sqrt(2.0)
#             # Ef = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)
#             Tf = (Zs_total + Ys_total + Xs_total)/np.sqrt(3.0)
#         else:
#             Af = 0
#             Tf = 0
#         self.SA = np.ones_like(self.SA)*(np.mean(np.abs(self.DAf-Af)).data)**2/(dt*len(self.tdi_fs['X']))*2
#         # self.SE = Nmodel.psd(freq=freq, option="E")
#         self.ST = np.ones_like(self.ST)*(np.mean(np.abs(self.DTf-Tf)).data)**2/(dt*len(self.tdi_fs['X']))*2

#         #GPU
#         self.SA_full_f = np.ones_like(self.SA_full_f)*(np.mean(np.abs(self.DAf_full_f-Af)).data)**2/(dt*len(self.tdi_fs['X']))*2
#         # self.SA_full_f = np.ones_like(self.SA_full_f)*(np.mean(np.abs(DAf_original-Af)).data)**2/(dt*len(self.tdi_fs['X']))*2
#         # self.SA_full_f = np.ones_like(self.SA_full_f)*(np.mean(np.abs(tdi_fs_residual['A']-Af)).data)**2/(dt*len(self.tdi_fs['X']))*2
        
        
#         # self.ST_full_f = np.ones_like(self.ST_full_f)*(np.mean(np.abs(self.DTf_full_f-Tf)).data)**2
#         self.PSD_GPU =  [xp.array(self.SA_full_f),
#                 xp.array(self.SA_full_f),
#         ]

#     def f_statistic(self, N_frequency, N_sky):
#         # We construct a global proposal density using the single
#         # source F statistic to compute the individual likelihoods
#         F_stat = []
#         frequency = []
#         eclipticlatitude = []
#         eclipticlongitude = []
#         pGBf = {}
#         for parameter in parameters:
#             pGBf[parameter] = 0
#         pGBf['Amplitude'] = 1e-24
#         pGBf['FrequencyDerivative'] = 0
#         frequency_boundaries = [self.lower_frequency,self.upper_frequency]
#         for n in range(N_sky):
#             eclipticlatitude.append(self.boundaries['EclipticLatitude'][0]+(self.boundaries['EclipticLatitude'][1]-self.boundaries['EclipticLatitude'][0])*n/N_sky)
#             eclipticlongitude.append(self.boundaries['EclipticLongitude'][0]+(self.boundaries['EclipticLongitude'][1]-self.boundaries['EclipticLongitude'][0])*n/N_sky)
#         for k in range(N_frequency):
#             F_stat.append([])
#             frequency.append(frequency_boundaries[0] + (frequency_boundaries[1]-frequency_boundaries[0])*k/(N_frequency-1))
#             for l in range(N_sky):
#                 F_stat[-1].append([])
#                 for m in range(N_sky):
#                     F_stat[-1][-1].append(self.F_fd0(frequency[-1],eclipticlatitude[l],eclipticlongitude[m],pGBf))
#         F_stat = np.asarray(F_stat)
#         return F_stat, frequency, eclipticlatitude, eclipticlongitude

#     def F_fd0(self, f0, theta, phi, pGBs):
#         g = []
#         pGBs['Frequency'] = f0
#         pGBs['EclipticLongitude'] = theta
#         pGBs['EclipticLatitude'] = phi
#         pGBs['InitialPhase'] = 0
#         pGBs['Inclination'] = np.pi/2
#         pGBs['Polarization'] = 0
#         g.append(self.GB.get_fd_tdixyz(template=pGBs, oversample=4))
#         pGBs['InitialPhase'] = np.pi
#         pGBs['Inclination'] = np.pi/2
#         pGBs['Polarization'] = np.pi/4
#         g.append(self.GB.get_fd_tdixyz(template=pGBs, oversample=4))
#         pGBs['InitialPhase'] = 3*np.pi/2
#         pGBs['Inclination'] = np.pi/2
#         pGBs['Polarization'] = 0
#         g.append(self.GB.get_fd_tdixyz(template=pGBs, oversample=4))
#         pGBs['InitialPhase'] = np.pi/2
#         pGBs['Inclination'] = np.pi/2
#         pGBs['Polarization'] = np.pi/4
#         g.append(self.GB.get_fd_tdixyz(template=pGBs, oversample=4))
#         g2 = []
#         for i in range(4):
#             g2.append([])
#             for j in range(3):
#                 g2[i].append(xr.align(self.dataX, g[i][j], join='left',fill_value=0)[1])
#         g = g2
#         data = [self.dataX,self.dataY,self.dataZ]
#         f = 0
#         for i in range(4):
#             for j in range(4):
#                 if i != j:
#                     f += 1/2* self.scalarproduct(g[i],g[j])**(-1)*self.scalarproduct(data,g[j])*self.scalarproduct(data,g[j])
#         return f

#     def F(self, intrinsic_parameter_values):
#         g = []
#         pGBs = {}
#         pGBs['Amplitude'] = 1e-24
#         for parameter in intrinsic_parameters:
#             pGBs[parameter] = intrinsic_parameter_values[parameter]
#         pGBs['InitialPhase'] = 0
#         pGBs['Inclination'] = np.pi/2
#         pGBs['Polarization'] = 0
#         g.append(self.GB.get_fd_tdixyz(template=pGBs, oversample=4))
#         pGBs['InitialPhase'] = np.pi
#         pGBs['Inclination'] = np.pi/2
#         pGBs['Polarization'] = np.pi/4
#         g.append(self.GB.get_fd_tdixyz(template=pGBs, oversample=4))
#         pGBs['InitialPhase'] = 3*np.pi/2
#         pGBs['Inclination'] = np.pi/2
#         pGBs['Polarization'] = 0
#         g.append(self.GB.get_fd_tdixyz(template=pGBs, oversample=4))
#         pGBs['InitialPhase'] = np.pi/2
#         pGBs['Inclination'] = np.pi/2
#         pGBs['Polarization'] = np.pi/4
#         g.append(self.GB.get_fd_tdixyz(template=pGBs, oversample=4))
#         g2 = []
#         for i in range(4):
#             g2.append([])
#             for j in range(3):
#                 g2[i].append(xr.align(self.dataX, g[i][j], join='left',fill_value=0)[1])
#         g = g2
#         data = [self.dataX,self.dataY,self.dataZ]
#         f = 0
#         for i in range(4):
#             for j in range(4):
#                 if i != j:
#                     f += 1/2* self.scalarproduct(g[i],g[j])**(-1)*self.scalarproduct(data,g[j])*self.scalarproduct(data,g[j])
#         return f

#     def scalarproduct(self, a, b):
#         diff = np.real(a[0].values * np.conjugate(b[0].values)) ** 2 + np.real(a[1].values * np.conjugate(b[1].values)) ** 2 + np.real(a[2].values * np.conjugate(b[2].values)) ** 2
#         res = 4*float(np.sum(diff / self.Sn) * self.dataX.df)
#         return res

#     def plot(self, maxpGBs=None, pGBadded=None, second_data = None,  found_sources_in= [], found_sources_not_matched= [], pGB_injected = [], pGB_injected_matched = [], added_label='Injection2', saving_label =None, vertical_lines = []):
#         plt.figure(figsize=fig_size)
#         fig, [ax1, ax2] = plt.subplots(2, 1, sharex=True, figsize=fig_size)
#         # plt.plot(dataX_training.f*1000,dataX_training.values, label='data')
#         # ax1.plot(self.dataX.f * 1000, self.dataX.values.real, label="data", marker="o", zorder=5)

#         # Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=self.pGB, oversample=4)
#         # index_low = np.searchsorted(Xs.f, self.dataX.f[0])
#         # Xs = Xs[index_low : index_low + len(self.dataX)]
#         # Ys = Ys[index_low : index_low + len(self.dataY)]
#         # Zs = Zs[index_low : index_low + len(self.dataZ)]

#         # Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=self.pGB, oversample=8)
#         # index_low = np.searchsorted(Xs.f, self.dataX.f[0])
#         # Xs = Xs[index_low : index_low + len(self.dataX)]
#         # Ys = Ys[index_low : index_low + len(self.dataY)]
#         # Zs = Zs[index_low : index_low + len(self.dataZ)]
#         # ax1.plot(Xs.f * 1000, Xs.values.real, label="VGB2", marker=".", zorder=5)
                    
#         # Af = (Zs - Xs)/np.sqrt(2.0)
#         ax1.plot(self.DAf.f*10**3,self.DAf,'k',zorder= 1, linewidth = 2, label = 'Data')
#         ax2.plot(self.DEf.f*10**3,np.abs(self.DEf),'k',zorder= 1, linewidth = 2, label = 'Data')
#         ax1.plot(self.DAf_full_f.f*10**3,self.DAf_full_f,'b.',zorder= 1, linewidth = 2, label = 'Data')
#         # ax2.plot(self.DEf_full_f.f*10**3,np.abs(self.DEf_full_f),'b.',zorder= 1, linewidth = 2, label = 'Data')
#         # ax2.plot(self.f_noise*10**3,np.sqrt(np.abs(self.SA)),'r',zorder= 1, linewidth = 2, label = 'Noise')
#         # ax2.plot(self.DAf_full_f.f*10**3,np.sqrt(np.abs(self.SA_full_f)),'b',zorder= 1, linewidth = 2, label = 'Noise')
        
#         # ax2.plot(self.f_noise*10**3,np.sqrt(self.psdE),'g',zorder= 1, linewidth = 2, label = 'Noise')
#         # ax2.plot(self.f_noise*10**3,np.sqrt(self.psdE*self.df),'r.',zorder= 1, linewidth = 2, label = 'Noise')
#         # ax2.plot(self.f_noise2*10**3,np.sqrt(self.psdE2),'g.',zorder= 1, linewidth = 2, label = 'Noise')
#         # ax2.plot(self.tdifE.f*10**3,np.abs(self.tdifE),'b',zorder= 1, linewidth = 2, label = 'Noise')
#         # ax1.plot(tdi_fs_long_subtracted.f[range_index],np.abs(tdi_fs_long_subtracted['X'][range_index])**2,'b',zorder= 5)


#         if second_data != None:
#             a,Xs = xr.align(self.dataX, second_data['X'], join='left',fill_value=0)
#             a,Ys = xr.align(self.dataY, second_data['Y'], join='left',fill_value=0)
#             a,Zs = xr.align(self.dataZ, second_data['Z'], join='left',fill_value=0)
#             Af = (Zs - Xs)/np.sqrt(2.0)
#             Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
#             ax1.plot(Af.f*10**3,Af,'k--',zorder= 1, linewidth = 2, label = 'Data subtracted')
#             ax2.plot(Ef.f*10**3,np.abs(Ef),'k--',zorder= 1, linewidth = 2, label = 'Data subtracted')

#         for j in range(len( pGB_injected)):
#             Xs, Ys, Zs = GB.get_fd_tdixyz(template= pGB_injected[j], oversample=4)
#             a,Xs = xr.align(self.dataX, Xs, join='left',fill_value=0)
#             a,Ys = xr.align(self.dataY, Ys, join='left',fill_value=0)
#             a,Zs = xr.align(self.dataZ, Zs, join='left',fill_value=0)
#             Af = (Zs - Xs)/np.sqrt(2.0)
#             Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
#             ax1.plot(Af.f*10**3,Af.data, color='grey', linewidth = 5, alpha = 0.5)
#             ax2.plot(Ef.f*10**3,np.abs(Ef.data), color='grey', linewidth = 5, alpha = 0.5)

#         for j in range(len(pGB_injected_matched)):
#             Xs, Ys, Zs = GB.get_fd_tdixyz(template= pGB_injected_matched[j], oversample=4)
#             a,Xs = xr.align(self.dataX, Xs, join='left',fill_value=0)
#             a,Ys = xr.align(self.dataY, Ys, join='left',fill_value=0)
#             a,Zs = xr.align(self.dataZ, Zs, join='left',fill_value=0)
#             Af = (Zs - Xs)/np.sqrt(2.0)
#             Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
#             ax1.plot(Af.f*10**3,Af.data, color=colors[j%10], linewidth = 5, alpha = 0.5)
#             ax2.plot(Ef.f*10**3,np.abs(Ef.data), color=colors[j%10], linewidth = 5, alpha = 0.5)


#         if pGBadded != None:
#             Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=pGBadded, oversample=4)
#             index_low = np.searchsorted(Xs.f, self.dataX.f[0])
#             Xs = Xs[index_low : index_low + len(self.dataX)]
#             Ys = Ys[index_low : index_low + len(self.dataY)]
#             Zs = Zs[index_low : index_low + len(self.dataZ)]
#             Af = (Zs - Xs)/np.sqrt(2.0)
#             Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
#             ax1.plot(Af.f* 1000, Af.data, marker='.', label=added_label)
#             ax2.plot(Ef.f* 1000, np.abs(Ef.data), marker='.', label=added_label)

#         for j in range(len(found_sources_in)):
#             Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=found_sources_in[j], oversample=4)
#             index_low = np.searchsorted(Xs.f, self.dataX.f[0])
#             Xs = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#             Zs = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
#             Af = (Zs - Xs)/np.sqrt(2.0)
#             Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
#             ax1.plot(Af.f* 1000, Af.data,'--', color= colors[j%10], linewidth = 1.6)
#             ax1.plot(Ef.f* 1000, Ef.data,'--', color= colors[j%10], linewidth = 1.6)
#             ax2.plot(Ef.f* 1000, np.abs(Ef.data),'--', color= colors[j%10], linewidth = 1.6)

#             params = pGB_dict_to_gpu_input([found_sources_in[j]])
#             gb.run_wave(*params, N=128, dt=dt, T=Tobs, oversample=4)
#             A = gb.A[0]
#             E = gb.E[0]
#             ax1.plot(gb.freqs[0].get()* 1000, A.get(),'.', color= colors[j%10], linewidth = 1.6)
#             ax1.plot(gb.freqs[0].get()* 1000, E.get(),'.', color= colors[j%10], linewidth = 1.6)
#             ax2.plot(gb.freqs[0].get()* 1000, np.abs(E.get()), color= colors[j%10], linewidth = 1.6)
#             # tdi_fs = {'X': Xs, 'Y': Ys, 'Z': Zs}
#             # tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft()) for k in ["X", "Y", "Z"]]))
#             # f, psdX =  scipy.signal.welch(tdi_ts["X"], fs=1.0/dt, nperseg=len(tdi_ts["X"]))
#             # ax2.plot(f* 1000, np.sqrt(psdX),'-.', color= 'b', linewidth = 1.6)
            

#         for j in range(len(found_sources_not_matched)):
#             Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=found_sources_not_matched[j], oversample=4)
#             index_low = np.searchsorted(Xs.f, self.dataX.f[0])
#             Xs = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#             Zs = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
#             Af = (Zs - Xs)/np.sqrt(2.0)
#             Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
#             ax1.plot(Af.f* 1000,Af.data,'.', color= colors[j%10], linewidth = 1.6)
#             ax2.plot(Ef.f* 1000, np.abs(Ef.data),'.', color= colors[j%10], linewidth = 1.6)

#         # ax1.plot(Xs_added2.f * 1000, Xs_added2.values.real, label="VGB2", marker=".", zorder=5)
#         ax1.axvline(self.lower_frequency* 1000, color= 'red', label='Boundaries')
#         ax1.axvline(self.upper_frequency* 1000, color= 'red')
#         ax2.axvline(self.lower_frequency* 1000, color= 'red')
#         ax2.axvline(self.upper_frequency* 1000, color= 'red')
#         for j in range(len(vertical_lines)):
#             ax1.axvline(vertical_lines[j]* 1000, color= 'red')
#             ax2.axvline(vertical_lines[j]* 1000, color= 'red')
#         # ax2.axvline(self.lower_frequency* 1000- 4*32*10**-6, color= 'green')
#         # ax2.axvline(self.upper_frequency* 1000+ 4*32*10**-6, color= 'green')
#         # if self.reduced_frequency_boundaries != None:
#         #     ax1.axvline(self.reduced_frequency_boundaries[0]* 1000, color= 'green', label='Reduced Boundaries')
#         #     ax1.axvline(self.reduced_frequency_boundaries[1]* 1000, color= 'green')

#         # ax1.plot(Xs.f * 1000, dataX.values.real - Xs.values.real, label="residual", alpha=0.8, color="red", marker=".")
#         plt.xlabel('f (mHz)')
#         ax1.set_ylabel('real A')    
#         ax2.set_ylabel('|E|') 
#         # ax1.set_yscale('log')  
#         ax2.set_yscale('log')   
#         ax1.set_xlim((self.lower_frequency-self.padding)*10**3, (self.upper_frequency+self.padding)*10**3)
#         ax2.set_xlim((self.lower_frequency-self.padding)*10**3, (self.upper_frequency+self.padding)*10**3)
#         # ax1.set_xlim((self.lower_frequency)*10**3, (self.upper_frequency)*10**3)
#         # ax2.set_xlim((self.lower_frequency)*10**3, (self.upper_frequency)*10**3)
#         ax1.xaxis.set_major_locator(plt.MaxNLocator(4))
#         ax2.xaxis.set_major_locator(plt.MaxNLocator(4))
#         # plt.legend()
#         # plt.pause(1)
#         if saving_label != None:
#             plt.savefig(saving_label,dpi=300,bbox_inches='tight')
#         # plt.pause(1)
#         plt.show()
#         # print("p true", self.loglikelihood([pGB]), "null hypothesis", self.loglikelihood([null_pGBs]))


#     def SNR_split(self, pGBs):
#         for i in range(len(pGBs)):
#             Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs[i], oversample=4)
#             if i == 0:
#                 Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
#             else:
#                 Xs_total += xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total += xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total += xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            
#         Af = (Zs_total - Xs_total)/np.sqrt(2.0)
#         Ef = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)
#         if self.upper_frequency > self.frequency_T_threshold:
#             Tf = (Zs_total + Ys_total + Xs_total)/np.sqrt(3.0)
#             hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2)/self.SA + np.absolute(Tf.data)**2 /self.ST)
#             SNR2 = np.sum( np.real(self.DAf * np.conjugate(Af.data) + self.DEf * np.conjugate(Ef.data))/self.SA + np.real(self.DTf * np.conjugate(Tf.data))/self.ST )
#         else:
#             hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /self.SA)
#             SNR2 = np.sum( np.real(self.DAf * np.conjugate(Af.data) + self.DEf * np.conjugate(Ef.data))/self.SA )
#         SNR = 4.0*Xs.df* hh
#         SNR2 = 4.0*Xs.df* SNR2
#         SNR3 = SNR2 / np.sqrt(SNR)
#         return SNR3.values

#     def SNR(self, pGBs):
#         for i in range(len(pGBs)):
#             Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs[i], oversample=4)
#             if i == 0:
#                 Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
#             else:
#                 Xs_total += xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total += xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total += xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            
#         Af = (Zs_total - Xs_total)/np.sqrt(2.0)
#         Ef = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)
#         if self.use_T_component:
#             Tf = (Zs_total + Ys_total + Xs_total)/np.sqrt(3.0)
#             hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2)/self.SA + np.absolute(Tf.data)**2 /self.ST)
#             SNR2 = np.sum( np.real(self.DAf * np.conjugate(Af.data) + self.DEf * np.conjugate(Ef.data))/self.SA + np.real(self.DTf * np.conjugate(Tf.data))/self.ST )
#         else:
#             SNR2 = np.sum( np.real(self.DAf * np.conjugate(Af.data) + self.DEf * np.conjugate(Ef.data))/self.SA )
#             hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /self.SA)
#         SNR = 4.0*Xs.df* hh
#         SNR2 = 4.0*Xs.df* SNR2
#         SNR3 = SNR2 / np.sqrt(SNR)
#         # plotIt = False
#         # if plotIt:
#         #     fig, ax = plt.subplots(nrows=3, sharex=True) 
#         #     ax[0].plot(Af.f, np.abs(self.DAf))
#         #     ax[0].plot(Af.f, np.abs(Af.data))
            
#         #     ax[1].plot(Af.f, np.abs(self.DEf))
#         #     ax[1].plot(Af.f, np.abs(Ef.data))
#         #     ax[2].plot(Af.f, np.abs(self.DTf))
#         #     ax[2].plot(Af.f, np.abs(Tf.data))
#         #     plt.show()
#         return SNR3.values

#     def SNR_AET_compute(self, pGBs):
#         for i in range(len(pGBs)):
#             Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs[i], oversample=4)
#             if i == 0:
#                 Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
#             else:
#                 Xs_total += xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total += xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total += xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            
#         Af = (Zs_total - Xs_total)/np.sqrt(2.0)
#         Ef = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)
#         Tf = (Zs_total + Ys_total + Xs_total)/np.sqrt(3.0)

#         tdi = dict({"A":Af, "E":Ef, "T":Tf, "X":Xs})
#         tdi_data = dict({"A":self.DAf, "E":self.DEf, "T":self.DTf})
#         hh = compute_tdi_snr(tdi, Nmodel, AET=True, fmin=Af.f[0], fmax=Af.f[-1])["tot2"]
#         SNR2 = compute_tdi_snr(tdi, Nmodel, data= tdi_data, AET=True, fmin=Af.f[0], fmax=Af.f[-1])["tot2"]
#         SNR3 = SNR2 / np.sqrt(hh)
#         return SNR3

#     def SNR_XYZ(self, pGBs):
#         for i in range(len(pGBs)):
#             Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs[i], oversample=4)
#             if i == 0:
#                 Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
#             else:
#                 Xs_total += xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total += xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total += xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            
#         hh = np.sum((np.absolute(Xs_total.data)**2 + np.absolute(Ys_total.data)**2 + np.absolute(Zs_total.data)**2) /self.Sn)
#         SNR2 = np.sum( np.real(self.dataX * np.conjugate(Xs_total.data) + self.dataY * np.conjugate(Ys_total.data)+ self.dataZ * np.conjugate(Zs_total.data))/self.Sn)
#         SNR = 4.0*Xs.df* hh
#         SNR2 = 4.0*Xs.df* SNR2
#         SNR3 = SNR2 / np.sqrt(SNR)
#         return SNR3.values

#     def SNR_XYZ_Sa(self, pGBs):
#         for i in range(len(pGBs)):
#             Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs[i], oversample=4)
#             if i == 0:
#                 Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
#             else:
#                 Xs_total += xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total += xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total += xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            
#         hh = np.sum((np.absolute(Xs_total.data)**2 + np.absolute(Ys_total.data)**2 + np.absolute(Zs_total.data)**2) /self.SA)
#         SNR2 = np.sum( np.real(self.dataX * np.conjugate(Xs_total.data) + self.dataY * np.conjugate(Ys_total.data)+ self.dataZ * np.conjugate(Zs_total.data))/self.SA)
#         SNR = 4.0*Xs.df* hh
#         SNR2 = 4.0*Xs.df* SNR2
#         SNR3 = SNR2 / np.sqrt(SNR)
#         return SNR3.values

#     def SNR_noise_matrix(self, pGBs):
#         for i in range(len(pGBs)):
#             Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs[i], oversample=4)
#             if i == 0:
#                 Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
#             else:
#                 Xs_total += xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total += xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total += xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
#         noise_model = "SciRDv1"
#         Nmodel = get_noise_model(noise_model, np.logspace(-5, -1, 100))

#         tdi = dict({"X":Xs_total, "Y":Ys_total, "Z":Zs_total})
#         tdi_data = dict({"X":self.dataX, "Y":self.dataY, "Z":self.dataZ})
#         hh = compute_tdi_snr(tdi, Nmodel)["tot2"]
#         SNR2 = compute_tdi_snr(tdi, Nmodel, data= tdi_data)["tot2"]
#         SNR3 = SNR2 / np.sqrt(hh)
#         return SNR3

#     def SNR_AE(self, pGBs):
#         for i in range(len(pGBs)):
#             Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs[i], oversample=4)
#             if i == 0:
#                 Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
#             else:
#                 Xs_total += xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total += xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total += xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            
#         Af = (Zs_total - Xs_total)/np.sqrt(2.0)
#         Ef = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)
#         SNR2 = np.sum( np.real(self.DAf * np.conjugate(Af.data) + self.DEf * np.conjugate(Ef.data))/self.SA )
#         hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /self.SA)
#         SNR = 4.0*Xs.df* hh
#         SNR2 = 4.0*Xs.df* SNR2
#         SNR3 = SNR2 / np.sqrt(SNR)
#         return SNR3.values

#     def SNR2(self, pGBs):
#         for i in range(len(pGBs)):
#             Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=pGBs[i], oversample=4)
#             index_low = np.searchsorted(Xs.f, self.dataX.f[0])
#             if i == 0:
#                 Xs_total = Xs[index_low : index_low + len(self.dataX)]
#                 Ys_total = Ys[index_low : index_low + len(self.dataY)]
#                 Zs_total = Zs[index_low : index_low + len(self.dataZ)]
#             else:
#                 Xs_total += Xs[index_low : index_low + len(self.dataX)]
#                 Ys_total += Ys[index_low : index_low + len(self.dataY)]
#                 Zs_total += Zs[index_low : index_low + len(self.dataZ)]
#             if len(Xs_total) < len(self.dataX):
#                 a,Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)
#                 a,Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)
#                 a,Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)
#         Af = (Zs_total - Xs_total)/np.sqrt(2.0)
#         Ef = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)

#         diff = np.abs(Af.values) ** 2 + np.abs(Ef.values) ** 2
#         SNR = -float(np.sum(diff / self.SA) * self.dataX.df) /2

#         return SNR#/10000

#     def loglikelihoodXYZ(self, pGBs):
#         for i in range(len(pGBs)):
#             Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=pGBs[i], oversample=4)
#             index_low = np.searchsorted(Xs.f, self.dataX.f[0])
#             if i == 0:
#                 Xs_total = Xs[index_low : index_low + len(self.dataX)]
#                 Ys_total = Ys[index_low : index_low + len(self.dataY)]
#                 Zs_total = Zs[index_low : index_low + len(self.dataZ)]
#             else:
#                 Xs_total += Xs[index_low : index_low + len(self.dataX)]
#                 Ys_total += Ys[index_low : index_low + len(self.dataY)]
#                 Zs_total += Zs[index_low : index_low + len(self.dataZ)]
#             if len(Xs_total) < len(self.dataX):
#                 a,Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)
#                 a,Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)
#                 a,Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)

#         # Af = (Zs_total - Xs_total)/np.sqrt(2.0)
#         # Ef = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)
#         # diff = np.abs(self.DAf - Af.values) ** 2 + np.abs(self.DEf - Ef.values) ** 2
#         diff = np.abs(self.dataX - Xs_total.values) ** 2 + np.abs(self.dataY - Ys_total.values) ** 2 + np.abs(self.dataZ - Zs_total.values) ** 2
#         # p1 = -float(np.sum(diff / Sn)*Xs.attrs['df'])/2.0
#         p1 = float(np.sum(diff / self.Sn) * Xs_total.df) / 2.0
#         # p1 = np.exp(p1)
#         return -p1#/10000

#     def loglikelihood(self, pGBs):
#         for i in range(len(pGBs)):
#             Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs[i], oversample=4)
#             if i == 0:
#                 Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
#             else:
#                 Xs_total += xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total += xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total += xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            
#         Af = (Zs_total - Xs_total)/np.sqrt(2.0)
#         Ef = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)
#         Tf = (Zs_total + Ys_total + Xs_total)/np.sqrt(3.0)
#         if self.use_T_component:
#             Tf = (Zs_total + Ys_total + Xs_total)/np.sqrt(3.0)
#             hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2)/self.SA + np.absolute(Tf.data)**2 /self.ST)
#             SNR2 = np.sum( np.real(self.DAf * np.conjugate(Af.data) + self.DEf * np.conjugate(Ef.data))/self.SA + np.real(self.DTf * np.conjugate(Tf.data))/self.ST )
#         else:
#             SNR2 = np.sum( np.real(self.DAf * np.conjugate(Af.data) + self.DEf * np.conjugate(Ef.data))/self.SA )
#             hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /self.SA)
#         # dd = np.sum((np.absolute(self.DAf.data)**2 + np.absolute(self.DEf.data)**2) /self.SA)
#         plotIt = False
#         if plotIt:
#             fig, ax = plt.subplots(nrows=2, sharex=True) 
#             ax[0].plot(Af.f, np.abs(self.DAf))
#             ax[0].plot(Af.f, np.abs(Af.data))
            
#             ax[1].plot(Af.f, np.abs(self.DEf))
#             ax[1].plot(Af.f, np.abs(Ef.data))
#             plt.show()
#         logliks = 4.0*Xs.df*( SNR2 - 0.5 * hh )
#         return logliks.values
    
#     def loglikelihood_gpu_signal(self, params):
#         A, E = gb.inject_signal(*params, dt=dt, T=Tobs, oversample=4)
#         A = A[:-1]
#         E = E[:-1]
        

#         SNR2 = np.sum( np.real(self.DAf_full_f * np.conjugate(A) + self.DEf_full_f * np.conjugate(E))/self.SA_full_f )
#         hh = np.sum((np.absolute(A)**2 + np.absolute(E)**2) /self.SA_full_f)

#         logliks = 4.0*self.df*( SNR2 - 0.5 * hh )
#         return logliks.values
        
#     def loglikelihood_gpu(self, parameters, start_freq_ind=0):
#         # N_index = np.searchsorted(self.N_values,int(len(self.dataX)))
#         N = 256
#         gb.d_d = 0
#         # parameters[4] *= -1
#         partial_length = 1*10**4
#         full_length = parameters.shape[-1]
#         like = np.zeros(parameters.shape[-1])
#         if len(parameters.shape) == 2:
#             parameters = np.array([parameters])
#         for n in range(int(full_length/partial_length)):
#             like[n*partial_length:(n+1)*partial_length] = gb.get_ll(parameters[:,:,(n)*partial_length:(n+1)*partial_length], self.data_GPU, self.PSD_GPU, N=N, oversample=4, dt=dt, T=Tobs, start_freq_ind=start_freq_ind)
#         try:
#             like[int(full_length/partial_length)*partial_length:] = gb.get_ll(parameters[:,:,int(full_length/partial_length)*partial_length:], self.data_GPU, self.PSD_GPU, N=N, oversample=4, dt=dt, T=Tobs, start_freq_ind=start_freq_ind)
#         except:
#             pass
#         # like = gb.get_ll(parameters, self.data_GPU, self.PSD_GPU, N=N, dt=dt, T=Tobs)

#         return like

#     def loglikelihood_dd(self, pGBs):
#         for i in range(len(pGBs)):
#             Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs[i], oversample=4)
#             if i == 0:
#                 Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
#             else:
#                 Xs_total += xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total += xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total += xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            
#         Af = (Zs_total - Xs_total)/np.sqrt(2.0)
#         Ef = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)
#         Tf = (Zs_total + Ys_total + Xs_total)/np.sqrt(3.0)
#         if self.use_T_component:
#             Tf = (Zs_total + Ys_total + Xs_total)/np.sqrt(3.0)
#             hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2)/self.SA + np.absolute(Tf.data)**2 /self.ST)
#             SNR2 = np.sum( np.real(self.DAf * np.conjugate(Af.data) + self.DEf * np.conjugate(Ef.data))/self.SA + np.real(self.DTf * np.conjugate(Tf.data))/self.ST )
#             dd = np.sum((np.absolute(self.DAf.data)**2 + np.absolute(self.DEf.data)**2) /self.SA + np.absolute(self.DTf.data)**2 /self.ST)
#         else:
#             SNR2 = np.sum( np.real(self.DAf * np.conjugate(Af.data) + self.DEf * np.conjugate(Ef.data))/self.SA )
#             hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /self.SA)
#             dd = np.sum((np.absolute(self.DAf.data)**2 + np.absolute(self.DEf.data)**2) /self.SA)
#         plotIt = False
#         if plotIt:
#             fig, ax = plt.subplots(nrows=2, sharex=True) 
#             ax[0].plot(Af.f, np.abs(self.DAf))
#             ax[0].plot(Af.f, np.abs(Af.data))
            
#             ax[1].plot(Af.f, np.abs(self.DEf))
#             ax[1].plot(Af.f, np.abs(Ef.data))
#             plt.show()
#         logliks = 4.0*Xs.df*( SNR2 - 0.5 * hh - 0.5 * dd)
#         return logliks.values

#     def loglikelihood_noise_matrix(self, pGBs):
#         for i in range(len(pGBs)):
#             Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs[i], oversample=4)
#             if i == 0:
#                 Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
#             else:
#                 Xs_total += xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total += xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total += xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            

#         tdi_signal = dict({"X":Xs_total, "Y":Ys_total, "Z":Zs_total})
#         tdi_data = dict({"X":self.dataX, "Y":self.dataY, "Z":self.dataZ})
#         hh = compute_tdi_snr(tdi_signal, Nmodel)["tot2"]
#         SNR2 = compute_tdi_snr(tdi_signal, Nmodel, data= tdi_data)["tot2"]
#         logliks = SNR2 - 0.5 * hh
#         return logliks

#     def intrinsic_SNR(self, pGBs):
#         for i in range(len(pGBs)):
#             Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs[i], oversample=4)
#             if i == 0:
#                 Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
#             else:
#                 Xs_total += xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total += xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total += xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            
#         Af = (Zs_total - Xs_total)/np.sqrt(2.0)
#         Ef = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)
#         if self.use_T_component:
#             Tf = (Zs_total + Ys_total + Xs_total)/np.sqrt(3.0)
#             hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2)/self.SA + np.absolute(Tf.data)**2 /self.ST)
#         else:
#             hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /self.SA)
#         SNR = 4.0*Xs.df* hh
#         return np.sqrt(SNR)

#     def intrinsic_SNR_T(self, pGBs):
#         for i in range(len(pGBs)):
#             Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs[i], oversample=4)
#             if i == 0:
#                 Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
#             else:
#                 Xs_total += xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total += xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total += xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            
#         Af = (Zs_total - Xs_total)/np.sqrt(2.0)
#         Ef = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)
#         Tf = (Zs_total + Ys_total + Xs_total)/np.sqrt(3.0)
#         hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /np.absolute(Tf.data)**2)
#         SNR = 4.0*Xs.df* hh
#         return np.sqrt(SNR)

#     def intrinsic_SNR_old(self, pGBs):
#         for i in range(len(pGBs)):
#             Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs[i], oversample=4)
#             if i == 0:
#                 Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
#             else:
#                 Xs_total += xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total += xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total += xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            
#         Af = (Zs_total - Xs_total)/np.sqrt(2.0)
#         Ef = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)
#         hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /self.SA)
#         SNR = 4.0*Xs.df* hh
#         return np.sqrt(SNR)

#     def SNR_with_rolling_mean(self, pGBs):
#         for i in range(len(pGBs)):
#             Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs[i], oversample=4)
#             if i == 0:
#                 Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
#             else:
#                 Xs_total += xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total += xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total += xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            
#         Af = (Zs_total - Xs_total)/np.sqrt(2.0)
#         Ef = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)
#         residual_A = self.DAf - Af.data
#         residual_E = self.DEf - Ef.data
#         average_length = 7
#         noise_rolling_mean_A = moving_average(np.abs(residual_A.values)**2, n=average_length)

#         hh = np.sum((np.absolute(Af.data[int((average_length-1)/2):len(Af.data)-int((average_length-1)/2)])**2 + np.absolute(Ef.data[int((average_length-1)/2):len(Af.data)-int((average_length-1)/2)])**2) /noise_rolling_mean_A)
#         SNR = 4.0*Xs.df* hh
#         return np.sqrt(SNR)

#     def differential_evolution_search(self, frequency_boundaries, initial_guess = None):
#         bounds = []
#         for signal in range(number_of_signals):
#             for i in range(7):
#                 bounds.append((0,1))

#         maxpGB = []
#         self.boundaries_reduced = deepcopy(self.boundaries)
#         self.boundaries_reduced['Frequency'] = frequency_boundaries
#         if initial_guess != None:
#             initial_guess01 = np.zeros((len(parameters)-1)*number_of_signals)
#             for signal in range(number_of_signals):
#                 pGBstart01 = scaleto01(initial_guess[signal], self.boundaries_reduced)

#                 for count, parameter in enumerate(parameters_no_amplitude):
#                     if pGBstart01[parameter] < 0:
#                         pGBstart01[parameter] = 0
#                     if pGBstart01[parameter] > 1:
#                         pGBstart01[parameter] = 1
#                     initial_guess01[count+(len(parameters_no_amplitude))*signal] = pGBstart01[parameter]
#             start = time.time()
#             res = differential_evolution(self.function_evolution, bounds=bounds, disp=False, strategy='best1exp', popsize=8,tol= 1e-8 , maxiter=1000, recombination= self.recombination, mutation=(0.5,1), x0=initial_guess01)
#             print('time',time.time()-start)
#         else:
#             start = time.time()
#             res = differential_evolution(self.function_evolution, bounds=bounds, disp=False, strategy='best1exp', popsize=8, tol= 1e-8 , maxiter=1000, recombination= self.recombination, mutation=(0.5,1))
#             print('time',time.time()-start)
#         for signal in range(number_of_signals):
#             pGB01 = [0.5] + res.x[signal*7:signal*7+7].tolist()
#             maxpGB.append(scaletooriginal(pGB01,self.boundaries_reduced))
#         print(res)
#         print(maxpGB)
#         print('log-likelihood',self.loglikelihood(maxpGB))
#         # print(pGB)
#         return [maxpGB], res.nfev

#     def differential_evolution_search_F(self, frequency_boundaries):
#         bounds = []
#         for signal in range(number_of_signals):
#             for i in range(4):
#                 bounds.append((0,1))

#         maxpGB = []
#         self.boundaries_reduced = deepcopy(self.boundaries)
#         self.boundaries_reduced['Frequency'] = frequency_boundaries
#         start = time.time()
#         res, energies = differential_evolution(self.function_evolution_F, bounds=bounds, disp=True, strategy='best1exp', popsize=5,tol= 1e-6 , maxiter=300, recombination= self.recombination, mutation=(0.5,1))
#         print('time',time.time()-start)
#         for signal in range(number_of_signals):
#             pGB01 = [0.5] + res.x[signal*4:signal*4+4].tolist() + [0.5,0.5,0.5] 
#             maxpGB.append(scaletooriginal(pGB01,self.boundaries_reduced))
#         print(res)
#         print(maxpGB)
#         print(self.loglikelihood(maxpGB))
#         # print(pGB)
#         return [maxpGB], energies

#     def searchCD(self):
#         # np.random.seed(42)

#         parameters_recorded = [None] * 10

#         n_trials = 50
#         number_of_evaluations = 0
#         for n in range(len(parameters_recorded)):
#             start = time.time()
#             parameters_recorded[n] = CoordinateMC(n, self.pGBs, self.boundaries, parameters_recorded, self.SNR, n_trials=n_trials)

#             print('n',n+1,'time', int(time.time()-start), np.round(parameters_recorded[n][0]['Loglikelihood'][-1],2),len(parameters_recorded[n][0]['Loglikelihood']))
#             number_of_evaluations += len(parameters_recorded[n][0]['Loglikelihood']) * n_trials
#         # pbar = tqdm(total=len(parameters_recorded))
#         # pool = mp.Pool(mp.cpu_count())
#         # parameters_recorded = pool.map(CoordinateMC, [n for n in range(len(parameters_recorded))])
#         # pool.close()
#         # pool.join()
#         # pbar.close()

#         best_run = 0
#         loglikelihoodofruns = np.zeros(len(parameters_recorded))
#         for i in range(len(parameters_recorded)):
#             loglikelihoodofruns[i] = parameters_recorded[i][0]['Loglikelihood'][-1]
#         best_value = np.max(loglikelihoodofruns)
#         best_run = np.argmax(loglikelihoodofruns)
#         good_runs = loglikelihoodofruns > 0
#         indices = (-loglikelihoodofruns).argsort()[:len(good_runs)]
#         pGBmodes = []
#         for i in range(len(good_runs)):
#             if good_runs[i]:
#                 pGBmodes.append({})
#         indices = (-loglikelihoodofruns).argsort()[:len(pGBmodes)]
#         pGBmodes = []
#         for n in indices:
#             pGBmodes.append([])
#             for i in range(number_of_signals):
#                 pGBmodes[-1].append({})
#                 for parameter in parameters_no_amplitude + ['Loglikelihood']:
#                     if parameter == 'Loglikelihood':
#                         pGBmodes[-1][i][parameter] = parameters_recorded[n][0][parameter][-1]
#                     else:
#                         pGBmodes[-1][i][parameter] = parameters_recorded[n][i][parameter][-1]
#         if len(pGBmodes) > 5:
#             for signal in range(number_of_signals):
#                 pGBmodes = pGBmodes[:5]
#         pGBmodes_opt = self.optimize_without_amplitude(pGBmodes)
#         return [pGBmodes_opt], number_of_evaluations

#     def optimize(self, pGBmodes, boundaries = None):
#         if boundaries == None:
#             boundaries = self.boundaries
#         bounds = ()
#         number_of_signals_optimize = len(pGBmodes[0])
#         for signal in range(number_of_signals_optimize):
#             bounds += ((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1))
#         for i in range(len(pGBmodes)):
#             maxpGB = []
#             boundaries_reduced = []
#             pGBs01 = []

#             # print('initial loglikelihood noise matrix', self.loglikelihood_noise_matrix(pGBmodes[0]),pGBmodes[0][0]['Frequency'])
#             # print('initial loglikelihood', self.loglikelihood(pGBmodes[0]),pGBmodes[0][0]['Frequency'])
#             for j in range(2):
#                 x = []
#                 for signal in range(number_of_signals_optimize):
#                     if j == 0:
#                         maxpGB.append({})
#                         boundaries_reduced.append({})
#                         for parameter in parameters:
#                             maxpGB[signal][parameter] = pGBmodes[i][signal][parameter]
#                     # print(maxpGB)
#                     # boundaries_reduced[signal] = reduce_boundaries(maxpGB[signal], boundaries,ratio=0.1)
#                     if j > 0:
#                         boundaries_reduced[signal] = reduce_boundaries(maxpGB[signal], boundaries,ratio=0.4)
#                     if j in [0]:
#                         boundaries_reduced[signal] = deepcopy(boundaries)
#                     pGBs01.append({})
#                     for parameter in parameters:
#                         if parameter in ["EclipticLatitude"]:
#                             pGBs01[signal][parameter] = (np.sin(maxpGB[signal][parameter]) - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
#                         elif parameter in ["Inclination"]:
#                             pGBs01[signal][parameter] = (np.cos(maxpGB[signal][parameter]) - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
#                         elif parameter in parameters_log_uniform:
#                             pGBs01[signal][parameter] = (np.log10(maxpGB[signal][parameter]) - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
#                         else:
#                             pGBs01[signal][parameter] = (maxpGB[signal][parameter] - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
#                     for parameter in parameters:
#                         x.append(pGBs01[signal][parameter])
#                 # print(loglikelihood(maxpGB))
#                 res = scipy.optimize.minimize(self.function, x, args=boundaries_reduced, method='SLSQP', bounds=bounds, tol=1e-5)
#                 # res = scipy.optimize.minimize(self.function, x, args=boundaries_reduced, method='Nelder-Mead', tol=1e-10)
#                 # res = scipy.optimize.least_squares(self.function, x, args=boundaries_reduced, bounds=bounds)
#                 for signal in range(number_of_signals_optimize):
#                     maxpGB[signal] = scaletooriginal(res.x[signal*8:signal*8+8],boundaries_reduced[signal])
#                 # print('optimized loglikelihood', loglikelihood(maxpGB),maxpGB)
#                 # print('boundaries reduced', boundaries_reduced)

#             best_value = self.loglikelihood(maxpGB)
#             if i == 0:
#                 current_best_value = best_value
#                 current_maxpGB = maxpGB
#             try:
#                 if current_best_value < best_value:
#                     current_best_value = best_value
#                     current_maxpGB = maxpGB
#             except:
#                 pass
#             # print('optimized loglikelihood', self.loglikelihood(maxpGB),self.loglikelihood([self.pGB]))
#         maxpGB = current_maxpGB
#         # print('final optimized loglikelihood noise matrix', self.loglikelihood_noise_matrix(maxpGB),maxpGB[0]['Frequency'])
#         # print('final optimized loglikelihood', self.loglikelihood(maxpGB),maxpGB[0]['Frequency'])
#         return maxpGB

#     def optimize_without_amplitude(self, pGBmodes, boundaries = None):
#         if boundaries == None:
#             boundaries = self.boundaries
#         bounds = ()
#         number_of_signals_optimize = len(pGBmodes[0])
#         for signal in range(number_of_signals_optimize):
#             bounds += ((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1))
#         for i in range(len(pGBmodes)):
#             maxpGB = []
#             pGBs01 = []

#             for j in range(2):
#                 x = []
#                 for signal in range(number_of_signals_optimize):
#                     if j == 0:
#                         maxpGB.append({})
#                         self.boundaries_reduced = {}
#                         for parameter in parameters_no_amplitude:
#                             maxpGB[signal][parameter] = pGBmodes[i][signal][parameter]
#                     # print(maxpGB)
#                     self.boundaries_reduced = reduce_boundaries(maxpGB[signal], boundaries,ratio=0.1)
#                     if j == 2:
#                         self.boundaries_reduced = reduce_boundaries(maxpGB[signal], boundaries,ratio=0.1)
#                     if j in [0,1]:
#                         self.boundaries_reduced = deepcopy(boundaries)
#                     pGBs01.append({})
#                     for parameter in parameters_no_amplitude:
#                         if parameter in ["EclipticLatitude"]:
#                             pGBs01[signal][parameter] = (np.sin(maxpGB[signal][parameter]) - self.boundaries_reduced[parameter][0]) / (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])
#                         elif parameter in ["Inclination"]:
#                             pGBs01[signal][parameter] = (np.cos(maxpGB[signal][parameter]) - self.boundaries_reduced[parameter][0]) / (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])
#                         elif parameter in parameters_log_uniform:
#                             pGBs01[signal][parameter] = (np.log10(maxpGB[signal][parameter]) - self.boundaries_reduced[parameter][0]) / (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])
#                         else:
#                             pGBs01[signal][parameter] = (maxpGB[signal][parameter] - self.boundaries_reduced[parameter][0]) / (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])
#                     for parameter in parameters_no_amplitude:
#                         x.append(pGBs01[signal][parameter])
#                 # print(loglikelihood(maxpGB))
#                 res = scipy.optimize.minimize(self.function_evolution, x, method='SLSQP', bounds=bounds, tol=1e-10)
#                 # res = scipy.optimize.minimize(self.function, x, args=self.boundaries_reduced, method='Nelder-Mead', tol=1e-10)
#                 # res = scipy.optimize.least_squares(self.function, x, args=self.boundaries_reduced, bounds=bounds)
#                 for signal in range(number_of_signals_optimize):
#                     pGB01 = [0.5] + res.x[signal*7:signal*7+7].tolist()
#                     maxpGB[signal] = scaletooriginal(pGB01,self.boundaries_reduced)
#                     # maxpGB[signal] = scaletooriginal(res.x[signal*7:signal*7+7],self.boundaries_reduced)
#                 # print('optimized loglikelihood', loglikelihood(maxpGB),maxpGB)
#                 # print('boundaries reduced', self.boundaries_reduced)

#             best_value = self.loglikelihood(maxpGB)
#             if i == 0:
#                 current_best_value = best_value
#                 current_maxpGB = maxpGB
#             try:
#                 if current_best_value < best_value:
#                     current_best_value = best_value
#                     current_maxpGB = maxpGB
#             except:
#                 pass
#             # print('optimized loglikelihood', self.loglikelihood(maxpGB),self.loglikelihood([self.pGB]))
#         maxpGB = current_maxpGB
#         print('final optimized loglikelihood', self.loglikelihood(maxpGB),maxpGB[0]['Frequency'])
#         return maxpGB

#     def calculate_Amplitude(self, pGBs):
#         for i in range(len(pGBs)):
#             Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs[i], oversample=4)
#             if i == 0:
#                 Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
#             else:
#                 Xs_total += xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total += xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total += xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            
#         Af = (Zs_total - Xs_total)/np.sqrt(2.0)
#         Ef = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)
#         SNR2 = np.sum( np.real(self.DAf * np.conjugate(Af.data) + self.DEf * np.conjugate(Ef.data))/self.SA )
#         hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /self.SA)
#         logliks = 4.0*Xs.df*( SNR2 - 0.5 * hh )
#         scalar_product_hh = 4.0*Xs.df* hh
#         scalar_product_dh = 4.0*Xs.df* SNR2
#         A = scalar_product_dh / scalar_product_hh
#         return A


#     def optimizeA(self, pGBmodes, boundaries = None):
#         if boundaries == None:
#             boundaries = self.boundaries
#         bounds = ()
#         number_of_signals_optimize = len(pGBmodes[0])
#         for signal in range(number_of_signals_optimize):
#             bounds += ((0,1))
#         for i in range(len(pGBmodes)):
#             maxpGB = []
#             boundaries_reduced = []
#             pGBs01 = []

#             for j in range(2):
#                 x = []
#                 pGBx = []
#                 for signal in range(number_of_signals_optimize):
#                     if j == 0:
#                         maxpGB.append({})
#                         boundaries_reduced.append({})
#                         for parameter in parameters:
#                             maxpGB[signal][parameter] = pGBmodes[i][signal][parameter]
#                     # print(maxpGB)
#                     boundaries_reduced[signal] = reduce_boundaries(maxpGB[signal], boundaries,ratio=0.1)
#                     if j == 2:
#                         boundaries_reduced[signal] = reduce_boundaries(maxpGB[signal], boundaries,ratio=0.1)
#                     if j in [0,1]:
#                         boundaries_reduced[signal] = deepcopy(boundaries)
#                     pGBs01.append({})
#                     for parameter in parameters:
#                         if parameter in ["EclipticLatitude"]:
#                             pGBs01[signal][parameter] = (np.sin(maxpGB[signal][parameter]) - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
#                         elif parameter in ["Inclination"]:
#                             pGBs01[signal][parameter] = (np.cos(maxpGB[signal][parameter]) - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
#                         elif parameter in parameters_log_uniform:
#                             pGBs01[signal][parameter] = (np.log10(maxpGB[signal][parameter]) - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
#                         else:
#                             pGBs01[signal][parameter] = (maxpGB[signal][parameter] - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
#                     for parameter in ['Amplitude']:
#                         x.append(pGBs01[signal][parameter])
#                     for parameter in parameters:
#                         pGBx.append(pGBs01[signal][parameter])
#                 self.pGBx = pGBx
#                 res = scipy.optimize.minimize(self.functiona, x, args=(pGBx, boundaries_reduced), method='trust-constr', bounds=[bounds], tol=1e-1)
#                 # res = scipy.optimize.minimize(self.function, x, args=boundaries_reduced, method='Nelder-Mead', tol=1e-10)
#                 # res = scipy.optimize.least_squares(self.function, x, args=boundaries_reduced, bounds=bounds)
#                 for signal in range(number_of_signals_optimize):
#                     pGB01 = deepcopy(pGBx)
#                     pGB01[signal*8] = res.x[signal]
#                     maxpGB[signal] = scaletooriginal(pGB01,boundaries_reduced[signal])
#                 # print('optimized loglikelihood', loglikelihood(maxpGB),maxpGB)
#                 # print('boundaries reduced', boundaries_reduced)

#             best_value = self.loglikelihood(maxpGB)
#             if i == 0:
#                 current_best_value = best_value
#                 current_maxpGB = maxpGB
#             try:
#                 if current_best_value < best_value:
#                     current_best_value = best_value
#                     current_maxpGB = maxpGB
#             except:
#                 pass
#             # print('optimized loglikelihood', self.loglikelihood(maxpGB),self.loglikelihood([self.pGB]))
#         maxpGB = current_maxpGB
#         print('final optimized loglikelihood', self.loglikelihood(maxpGB),maxpGB[0]['Frequency'])
#         return maxpGB

#     def functionA(self, a):
#         pGBs01 = self.pGBx
#         boundaries_reduced = deepcopy(self.boundaries)
#         pGBs = []
#         for signal in range(int(len(pGBs01)/8)):
#             pGBs.append({})
#             i = 0
#             for parameter in parameters:
#                 if parameter in ["EclipticLatitude"]:
#                     pGBs[signal][parameter] = np.arcsin((pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
#                 elif parameter in ["Inclination"]:
#                     shifted_inclination = pGBs01[signal*8:signal*8+8][i]
#                     if pGBs01[signal*8:signal*8+8][i] < 0:
#                         shifted_inclination = pGBs01[signal*8:signal*8+8][i] + 1
#                     if pGBs01[signal*8:signal*8+8][i] > 1:
#                         shifted_inclination = pGBs01[signal*8:signal*8+8][i] - 1
#                     pGBs[signal][parameter] = np.arccos((shifted_inclination * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
#                 elif parameter in ["FrequencyDerivative"]:
#                     pGBs[signal][parameter] = 10**((pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
#                 elif parameter in ['Amplitude']:
#                     pGBs[signal][parameter] = 10**((a[0] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
#                 else:
#                     pGBs[signal][parameter] = (pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0]
#                 i += 1
        
#         print(pGBs)
#         p = -self.loglikelihood(pGBs)
#         return p#/10**4

#     def functiona(self,a, pGBs01, boundaries_reduced):
#         pGBs = []
#         for signal in range(int(len(pGBs01)/8)):
#             pGBs.append({})
#             i = 0
#             for parameter in parameters:
#                 if parameter in ["EclipticLatitude"]:
#                     pGBs[signal][parameter] = np.arcsin((pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
#                 elif parameter in ["Inclination"]:
#                     shifted_inclination = pGBs01[signal*8:signal*8+8][i]
#                     if pGBs01[signal*8:signal*8+8][i] < 0:
#                         shifted_inclination = pGBs01[signal*8:signal*8+8][i] + 1
#                     if pGBs01[signal*8:signal*8+8][i] > 1:
#                         shifted_inclination = pGBs01[signal*8:signal*8+8][i] - 1
#                     pGBs[signal][parameter] = np.arccos((shifted_inclination * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
#                 elif parameter in ["FrequencyDerivative"]:
#                     pGBs[signal][parameter] = 10**((pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
#                 elif parameter in ['Amplitude']:
#                     pGBs[signal][parameter] = 10**((a[0] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
#                 else:
#                     pGBs[signal][parameter] = (pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0]
#                 i += 1
#         p = -self.loglikelihood(pGBs)
#         return p#/10**4

#     def function(self, pGBs01, boundaries_reduced):
#         pGBs = []
#         for signal in range(int(len(pGBs01)/8)):
#             pGBs.append({})
#             i = 0
#             for parameter in parameters:
#                 if parameter in ["EclipticLatitude"]:
#                     pGBs[signal][parameter] = np.arcsin((pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
#                 elif parameter in ["Inclination"]:
#                     shifted_inclination = pGBs01[signal*8:signal*8+8][i]
#                     if pGBs01[signal*8:signal*8+8][i] < 0:
#                         shifted_inclination = pGBs01[signal*8:signal*8+8][i] + 1
#                     if pGBs01[signal*8:signal*8+8][i] > 1:
#                         shifted_inclination = pGBs01[signal*8:signal*8+8][i] - 1
#                     pGBs[signal][parameter] = np.arccos((shifted_inclination * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
#                 elif parameter in parameters_log_uniform:
#                     pGBs[signal][parameter] = 10**((pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
#                 else:
#                     pGBs[signal][parameter] = (pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0]
#                 i += 1
#         p = -self.loglikelihood(pGBs)
#         return p#/10**4

#     def function_evolution(self, pGBs01):
#         pGBs = []
#         for signal in range(number_of_signals):
#             pGBs.append({})
#             i = 0
#             for parameter in parameters:
#                 if parameter in ["EclipticLatitude"]:
#                     pGBs[signal][parameter] = np.arcsin((pGBs01[signal*7:signal*7+7][i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
#                 elif parameter in ["Inclination"]:
#                     pGBs[signal][parameter] = np.arccos((pGBs01[signal*7:signal*7+7][i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
#                 elif parameter in ['Amplitude']:
#                     i -= 1
#                     pGBs[signal][parameter] = 10**((0.1 * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
#                 # elif parameter in ["FrequencyDerivative"]:
#                 #     pGBs[signal][parameter] = 10**((pGBs01[signal*7:signal*7+7][i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
#                 else:
#                     pGBs[signal][parameter] = (pGBs01[signal*7:signal*7+7][i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0]
#                 i += 1
#         p = -self.SNR(pGBs)
#         return p

#     def function_evolution_F(self, pGBs01):
#         pGBs =  {}
#         i = 0
#         for parameter in intrinsic_parameters:
#             if parameter in ["EclipticLatitude"]:
#                 pGBs[parameter] = np.arcsin((pGBs01[i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
#             elif parameter in ["Inclination"]:
#                 pGBs[parameter] = np.arccos((pGBs01[i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
#             elif parameter in ['Amplitude']:
#                 i -= 1
#                 pGBs[parameter] = 10**((0.1 * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
#             elif parameter in ["FrequencyDerivative"]:
#                 pGBs[parameter] = 10**((pGBs01[i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
#             else:
#                 pGBs[parameter] = (pGBs01[i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0]
#             i += 1
#         intrinsic_parameter_values = {}
#         for parameter in intrinsic_parameters:
#             intrinsic_parameter_values[parameter] = pGBs[parameter]
#         p = -self.F(intrinsic_parameter_values)
#         return p

#     def function_evolution_8(self, pGBs01):
#         pGBs = []
#         for signal in range(number_of_signals):
#             pGBs.append({})
#             i = 0
#             for parameter in parameters:
#                 if parameter in ["EclipticLatitude"]:
#                     pGBs[signal][parameter] = np.arcsin((pGBs01[signal*8:signal*8+8][i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
#                 elif parameter in ["Inclination"]:
#                     pGBs[signal][parameter] = np.arccos((pGBs01[signal*8:signal*8+8][i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
#                 elif parameter in parameters_log_uniform:
#                     pGBs[signal][parameter] = 10**((pGBs01[signal*8:signal*8+8][i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
#                 else:
#                     pGBs[signal][parameter] = (pGBs01[signal*8:signal*8+8][i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0]
#                 i += 1
#         p = -self.loglikelihood(pGBs)
#         return p#/10**4

#     def fisher_information(self, maxpGB):
#         maxpGB_changed = deepcopy(maxpGB)
#         maxpGB01 = scaleto01(maxpGB, self.boundaries)
#         maxpGB01_changed = deepcopy(maxpGB01)
#         step_size = {}
#         pGB_low = {}
#         pGB_high = {}
#         derivativeAf = {}
#         derivativeEf = {}
#         inner_product = {}
#         for i in range(1):
#             for parameter in parameters:
#                 if i == 0:
#                     step_size[parameter] = 1e-9
#                     # if parameter == 'Frequency':
#                     #     step_size[parameter] = 0.00001
#                 else:
#                     step_size[parameter] = 0.001/np.sqrt(inner_product[parameter][parameter])
#                 # if step_size[parameter] > 1e-9:
#                 #     step_size[parameter] = 1e-9
#                 pGB_low = maxpGB01[parameter] - step_size[parameter]/2
#                 pGB_high = maxpGB01[parameter] + step_size[parameter]/2
#                 # print(parameter, step_size[parameter],i)
#                 # print(parameter, pGB_low, pGB_high)
#                 if pGB_low < 0:
#                     pGB_low = 0
#                 if pGB_high > 1:
#                     pGB_high = 1
#                 maxpGB01_changed[parameter] = pGB_low
#                 maxpGB_changed = scaletooriginalparameter(maxpGB01_changed,self.boundaries)
#                 # print(maxpGB_changed)
#                 Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=maxpGB_changed, oversample=4)
#                 index_low = np.searchsorted(Xs.f, self.dataX.f[0])
#                 Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
#                 Af_low = (Zs_total - Xs_total)/np.sqrt(2.0)
#                 Ef_low = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)

#                 maxpGB01_changed[parameter] = pGB_high
#                 maxpGB_changed = scaletooriginalparameter(maxpGB01_changed,self.boundaries)
#                 Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=maxpGB_changed, oversample=4)
#                 index_low = np.searchsorted(Xs.f, self.dataX.f[0])
#                 Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
#                 Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
#                 Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
#                 Af_high = (Zs_total - Xs_total)/np.sqrt(2.0)
#                 Ef_high = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)

#                 derivativeAf[parameter] = (Af_high - Af_low)/step_size[parameter]
#                 derivativeEf[parameter] = (Ef_high - Ef_low)/step_size[parameter]

#                 maxpGB01_changed[parameter] = maxpGB01[parameter]

#             for parameter1 in parameters:
#                 inner_product[parameter1] = {}
#                 for parameter2 in parameters:
#                     AE = derivativeAf[parameter1]*np.conjugate(derivativeAf[parameter2]) + derivativeAf[parameter1]*np.conjugate(derivativeAf[parameter2])
#                     inner_product[parameter1][parameter2] = 4*float(np.real(np.sum(AE / self.SA) * self.dataX.df))
#             print(step_size['Amplitude'],inner_product['Amplitude']['Amplitude'],step_size['Frequency'],inner_product['Frequency']['Frequency'])
#         return inner_product

# def objective(n,tdi_fs,Tobs):
#     print('index',n)
#     search = Search(n,tdi_fs,Tobs)
#     # search.plot(search.pGBs)
#     pGBmodes =  search.search()
#     maxpGB, pGB =  search.optimize(pGBmodes)
#     return maxpGB, pGB

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
# parameters_log_uniform = ['Amplitude','FrequencyDerivative']
parameters_log_uniform = ['Amplitude']
parameters_no_amplitude = parameters[1:]
intrinsic_parameters = ['EclipticLatitude','EclipticLongitude','Frequency', 'FrequencyDerivative']

# get current directory
path = os.getcwd()
 
# parent directory
parent = os.path.dirname(path)
# grandparent directory
grandparent = os.path.dirname(parent)
Radler = False
if Radler:
    DATAPATH = grandparent+"/LDC/Radler/data"
    SAVEPATH = grandparent+"/LDC/pictures/LDC1-4/"
else:
    DATAPATH = grandparent+"/LDC/Sangria/data"
    SAVEPATH = grandparent+"/LDC/pictures/Sangria/"

if Radler:
    sangria_fn = DATAPATH + "/LDC1-4_GB_v2.hdf5"
    # sangria_fn = DATAPATH + "/LDC1-3_VGB_v2_FD_noiseless.hdf5"
    # sangria_fn = DATAPATH + "/LDC1-3_VGB_v2.hdf5"
else:
    sangria_fn = DATAPATH + "/LDC2_sangria_training_v2.h5"
fid = h5py.File(sangria_fn)

reduction = 1

# get TDI 
if Radler:
    td = np.array(fid["H5LISA/PreProcess/TDIdata"])
    td = np.rec.fromarrays(list(td.T), names=["t", "X", "Y", "Z"])
    dt = float(np.array(fid['H5LISA/GWSources/GalBinaries']['Cadence']))
    Tobs = float(int(np.array(fid['H5LISA/GWSources/GalBinaries']['ObservationDuration']))/reduction)
else:
    td = fid["obs/tdi"][()]
    td = np.rec.fromarrays(list(td.T), names=["t", "X", "Y", "Z"])
    td = td['t']
    dt = td["t"][1]-td["t"][0]
    
    td_mbhb = fid["sky/mbhb/tdi"][()]
    # cat_mbhb = fid["sky/mbhb/cat"]
    td_mbhb  = np.rec.fromarrays(list(td_mbhb .T), names=["t", "X", "Y", "Z"])
    td_mbhb  = td_mbhb ['t']
    # tdi_ts_mbhb = dict([(k, TimeSeries(td_mbhb[k][:int(len(td_mbhb[k][:])/reduction)], dt=dt)) for k in ["X", "Y", "Z"]])
    # tdi_fs_mbhb = xr.Dataset(dict([(k, tdi_ts_mbhb[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))

    # td_dgb = fid["sky/dgb/tdi"][()]
    # cat_dgb = fid["sky/dgb/cat"]
    # td_dgb  = np.rec.fromarrays(list(td_dgb .T), names=["t", "X", "Y", "Z"])
    # td_dgb  = td_dgb ['t']
    # tdi_ts_dgb = dict([(k, TimeSeries(td_dgb[k][:int(len(td_dgb[k][:])/reduction)], dt=dt)) for k in ["X", "Y", "Z"]])
    # tdi_fs_dgb = xr.Dataset(dict([(k, tdi_ts_dgb[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))

    # td_igb = fid["sky/igb/tdi"][()]
    # cat_igb = fid["sky/igb/cat"]
    # td_igb  = np.rec.fromarrays(list(td_igb .T), names=["t", "X", "Y", "Z"])
    # td_igb  = td_igb ['t']
    # tdi_ts_igb = dict([(k, TimeSeries(td_igb[k][:int(len(td_igb[k][:])/reduction)], dt=dt)) for k in ["X", "Y", "Z"]])
    # tdi_fs_igb = xr.Dataset(dict([(k, tdi_ts_igb[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))

    # td_vgb = fid["sky/vgb/tdi"][()]
    # cat_vgb = fid["sky/vgb/cat"]
    # td_vgb  = np.rec.fromarrays(list(td_vgb .T), names=["t", "X", "Y", "Z"])
    # td_vgb  = td_vgb ['t']
    # tdi_ts_vgb = dict([(k, TimeSeries(td_vgb[k][:int(len(td_vgb[k][:])/reduction)], dt=dt)) for k in ["X", "Y", "Z"]])
    # tdi_fs_vgb = xr.Dataset(dict([(k, tdi_ts_vgb[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))

    Tobs = float(int(np.array(fid['obs/config/t_max']))/reduction)
    for k in ["X", "Y", "Z"]:
        td[k] = td[k] - td_mbhb[k]


# Build timeseries and frequencyseries object for X,Y,Z
tdi_ts = dict([(k, TimeSeries(td[k][:int(len(td[k][:])/reduction)], dt=dt, t0=td.t[0])) for k in ["X", "Y", "Z"]])
tdi_fs = xr.Dataset(dict([(k, tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
GB = fastGB.FastGB(delta_t=dt, T=Tobs)  # in seconds


pGBadded20 = {}
pGBadded20['Amplitude'] = 1e-21
pGBadded20['EclipticLatitude'] = 0.2
pGBadded20['EclipticLongitude'] = 1.5
pGBadded20['Frequency'] = 0.0003
pGBadded20['FrequencyDerivative'] = 5*1e-17
pGBadded20['Inclination'] = 1.2
pGBadded20['InitialPhase'] = 3
pGBadded20['Polarization'] = 2


# add signal
add_signal = False
if add_signal:
    for pGBadding in [pGBadded20]: # faint
        Xs_added, Ys_added, Zs_added = GB.get_fd_tdixyz(template=pGBadding, oversample=4)
        source_added = dict({"X": Xs_added, "Y": Ys_added, "Z": Zs_added})
        index_low = np.searchsorted(tdi_fs["X"].f, Xs_added.f[0])
        index_high = index_low+len(Xs_added)
        # tdi_fs['X'] = tdi_fs['X'] #+ Xs_added
        for k in ["X", "Y", "Z"]:
            tdi_fs[k].data[index_low:index_high] = tdi_fs[k].data[index_low:index_high] + source_added[k].data
    tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft(dt=dt)) for k, n in [["X", 1], ["Y", 2], ["Z", 3]]]))

noise_model = "SciRDv1"
Nmodel = get_noise_model(noise_model, np.logspace(-5, -1, 100))

pGB = {}
ind = 0
found_sources = []
target_sources = []
first_start = time.time()
np.random.seed(42) #40
number_of_signals = 1
signals_per_subtraction = 1

f = 0.0115248
def frequency_derivative(f, Mc):
    G = 6.674*10**(-11)
    c = 3*10**8
    Mc = Mc * 2*10**30
    Mc_s = Mc*G/c**3
    return 96/(5*np.pi*Mc_s**2)*(np.pi*Mc_s*f)**(11/3)
def frequency_derivative_tyson(f):
    return 8*10**-7*f**(11/3)
def frequency_derivative_tyson_lower(f):
    return -5*10**-6*f**(13/3)
def frequency_derivative2(f, Mc_s):
    return 96/5*np.pi**(8/3)*Mc_s**(5/3)*(f)**(11/3)
print('frequency derivative', frequency_derivative(f,0.1),frequency_derivative(f,2),' at f=', f)
chandrasekhar_limit = 1.4
M_chirp_upper_boundary = (chandrasekhar_limit**2)**(3/5)/(2*chandrasekhar_limit)**(1/5)


def tdi_subtraction(tdi_fs,found_sources_mp_subtract, frequencies_search, ratio = 1):
    #subtract the found sources from original
    tdi_fs_subtracted2 = deepcopy(tdi_fs)
    for i in range(len(found_sources_mp_subtract)):
        # for j in range(len(found_sources_to_subtract[i])):
            Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_mp_subtract[i], oversample=4)
            source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
            index_low = np.searchsorted(tdi_fs_subtracted2["X"].f, Xs_subtracted.f[0])
            index_high = index_low+len(Xs_subtracted)
            for k in ["X", "Y", "Z"]:
                tdi_fs_subtracted2[k].data[index_low:index_high] -= source_subtracted[k].data * ratio
    return tdi_fs_subtracted2

load_category = False
if load_category:
    try:
        cat = np.load(SAVEPATH+'cat_sorted.npy', allow_pickle = True)
        print('cat sorted loaded')
    except:
        # get the source parameters
        # Radler
        if Radler:
            names = np.array(fid['H5LISA/GWSources/GalBinaries']) # radler
            params = [fid['H5LISA/GWSources/GalBinaries'][k] for k in names]
            reduced_names = []
            i = 0
            for p in params:
                i += 1
                if p.shape:
                    reduced_names.append(names[i-1])
            params = [np.array(p) for p in params if p.shape]
            names = reduced_names
        # Sangria
        else:
            names_dgb = fid["sky/dgb/cat"].dtype.names # Sangria
            params_dgb = [np.array(fid["sky/dgb/cat"][k]).squeeze() for k in names_dgb]
            names_igb = fid["sky/igb/cat"].dtype.names # Sangria
            params_igb = [np.array(fid["sky/igb/cat"][k]).squeeze() for k in names_igb]
            names_vgb = fid["sky/vgb/cat"].dtype.names # Sangria
            params_vgb = [np.array(fid["sky/vgb/cat"][k]).squeeze() for k in names_vgb]

        cat = np.rec.fromarrays(params_dgb, names=list(names_dgb))
        indexes = np.argsort(cat['Frequency'])
        cat = cat[indexes]
        np.save(SAVEPATH+'cat_sorted.npy',cat)

# LDC1-4 #####################################
frequencies = []
frequencies_even = []
frequencies_odd = []
# search_range = [0.00398, 0.0041]
# search_range = [0.0039885, 0.0040205]
# search_range = [0.0039935, 0.0039965]
f_Nyquist = 1/dt/2
search_range = [0.0003, f_Nyquist]
search_range = [0.000299, f_Nyquist]
# search_range = [0.0003, 0.03]
# search_range = [0.0019935, 0.0020135]
# search_range = [0.0029935, 0.0030135]
# window_length = 1*10**-7 # Hz

frequencies = create_frequency_windows(search_range, Tobs)


frequencies_half_shifted = []
for i in range(len(frequencies)-1):
    frequencies_half_shifted.append([(frequencies[i][1]-frequencies[i][0])/2 +frequencies[i][0],(frequencies[i+1][1]-frequencies[i+1][0])/2 +frequencies[i+1][0]])
# frequencies = frequencies_half_shifted
frequencies_even = frequencies[::2]
frequencies_odd = frequencies[1::2]

##### plot number of signals per frequency window
# frequencies_search = frequencies[::10]
# pGB_injected = []
# for j in range(len(frequencies_search)):
#     padding = (frequencies_search[j][1] - frequencies_search[j][0])/2 *0
#     index_low = np.searchsorted(cat_sorted['Frequency'], frequencies_search[j][0]-padding)
#     index_high = np.searchsorted(cat_sorted['Frequency'], frequencies_search[j][1]+padding)
#     if cat_sorted['Frequency'][index_high] < frequencies_search[j][1]:
#         index_high -= 1
#     indexesA = np.argsort(-cat_sorted[index_low:index_high]['Amplitude'])
#     pGB_injected_window = []
#     pGB_stacked = {}
#     for parameter in parameters:
#         pGB_stacked[parameter] = cat_sorted[parameter][index_low:index_high][indexesA]
#     for i in range(len(cat_sorted['Amplitude'][index_low:index_high])):
#         pGBs = {}
#         for parameter in parameters:
#             pGBs[parameter] = pGB_stacked[parameter][i]
#         pGB_injected_window.append(pGBs)
#     pGB_injected.append(pGB_injected_window)

# counts = np.zeros(len(pGB_injected))
# for i in range(len(pGB_injected)):
#     counts[i] = len(pGB_injected[i])

# frequencies_search = np.asarray(frequencies_search)
# figure = plt.figure()
# plt.loglog(frequencies_search[:,1],counts, '.')
# plt.xlabel('Frequency [Hz]')
# plt.ylabel('Number of signals')
# plt.show()
# figure = plt.figure()
# plt.loglog(frequencies_search[:,1],frequencies_search[:,1]-frequencies_search[:,0],  linewidth= 4, label= 'Frequency window width')
# plt.loglog(frequencies_search[:,1],np.ones(len(frequencies_search[:,1]))*4*32*10**-9*2, label= 'LISA rotation')
# plt.loglog(frequencies_search[:,1],frequencies_search[:,1]*3* 10**-4, label= 'Doppler modulation')
# plt.loglog(frequencies_search[:,1],frequency_derivative(frequencies_search[:,1],2)*Tobs, label= '$\dot{f}_{max} \cdot T_{obs}$')
# plt.xlabel('Frequency [Hz]')
# plt.ylabel('Frequency window witdh [Hz]')
# plt.ylim(bottom=(frequencies_search[0,1]-frequencies_search[0,0])/10**1)
# plt.legend()
# plt.show()


# save_name = 'not_anticorrelatedLDC1-4_2_optimized_second'
save_name = 'Sangria_12m'
# save_name = 'Sangria_1_dynamic_noise'
# save_name = 'not_anticorrelatedRadler_half_dynamic_noise'
# for i in range(65):
frequencies_search = frequencies
frequencies_search_full = deepcopy(frequencies_search)
# batch_index = int(sys.argv[1])
batch_index = 2
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.003977)
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.00399)
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.00404)
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.00264612)
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.007977)
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.016308)-1
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.00545)-1
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.001373)-1
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.02355)-1
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.01488)-1
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], cat[-1]['Frequency'])-5
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.0004)-1
batch_size = 1
start_index = batch_size*batch_index
# print('batch',batch_index, start_index)
# frequencies_search = frequencies_search[start_index:start_index+batch_size]

# print(i, frequencies_search[0])
### highest + padding has to be less than f Nyqist
while frequencies_search[-1][1] + (frequencies_search[-1][1] - frequencies_search[-1][0])/2 > f_Nyquist:
    frequencies_search = frequencies_search[:-1]
# frequencies_search = frequencies_search[70:80]
# frequencies_search = frequencies_search[25:]


do_print = False
if do_print:
    found_sources_mp = np.load(SAVEPATH+'found_sources_' +save_name+'.npy', allow_pickle = True)
    # found_sources_mp = np.load(SAVEPATH+'found_sources' +save_name+'optimized.npy', allow_pickle = True)

    found_sources_in_flat = []
    for i in range(len(found_sources_mp)):
        for j in range(len(found_sources_mp[i])):
            found_sources_in_flat.append(found_sources_mp[i][j])
    found_sources_in_flat = np.asarray(found_sources_in_flat)
    found_sources_in_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_in_flat]) for attribute in found_sources_in_flat[0].keys()}
    found_sources_in_flat_df = pd.DataFrame(found_sources_in_flat_array)
    found_sources_in_flat_df = found_sources_in_flat_df.sort_values('Frequency')
    found_sources_in = []
    for i in range(len(frequencies_search)):
        found_sources_in.append(found_sources_in_flat_df[(found_sources_in_flat_df['Frequency'] > frequencies_search[i][0]) & (found_sources_in_flat_df['Frequency'] < frequencies_search[i][1])].to_dict(orient='records'))


    # found_sources_in_flat = []
    # found_sources_in_flat_frequency = []
    # number_of_found_flat = 0
    # for i in range(len(found_sources_mp)):
    #     for j in range(len(found_sources_mp[i][3])):
    #         found_sources_in_flat.append(found_sources_mp[i][3][j])
    #         found_sources_in_flat_frequency.append(found_sources_in_flat[-1]['Frequency'])
    #         number_of_found_flat += 1
    # found_sources_in_flat_frequency = np.asarray(found_sources_in_flat_frequency)
    # found_sources_in_flat = np.asarray(found_sources_in_flat)
    # indexes_in = np.argsort(found_sources_in_flat_frequency)
    # found_sources_in_flat_frequency = found_sources_in_flat_frequency[indexes_in]
    # found_sources_in_flat = found_sources_in_flat[indexes_in]

    # found_sources_in = []
    # for i in range(len(frequencies_search)):
    #     lower_index = np.searchsorted(found_sources_in_flat_frequency,frequencies_search[i][0])
    #     higher_index = np.searchsorted(found_sources_in_flat_frequency,frequencies_search[i][1])
    #     if i == 102:
    #         print(lower_index, higher_index)
    #     found_sources_in.append(found_sources_in_flat[lower_index:higher_index])

    pGB_injected = []
    for j in range(len(frequencies_search)):
        padding = (frequencies_search[j][1] - frequencies_search[j][0])/2 *0
        index_low = np.searchsorted(cat['Frequency'], frequencies_search[j][0]-padding)
        index_high = np.searchsorted(cat['Frequency'], frequencies_search[j][1]+padding)
        try:
            if cat['Frequency'][index_high] < frequencies_search[j][1]:
                index_high -= 1
        except:
            pass
        indexesA = np.argsort(-cat[index_low:index_high]['Amplitude'])
        pGB_injected_window = []
        pGB_stacked = {}
        for parameter in parameters:
            pGB_stacked[parameter] = cat[parameter][index_low:index_high][indexesA]
        for i in range(len(cat['Amplitude'][index_low:index_high])):
            pGBs = {}
            for parameter in parameters:
                pGBs[parameter] = pGB_stacked[parameter][i]
            pGB_injected_window.append(pGBs)
        pGB_injected.append(pGB_injected_window)


    # #check SNR
    # for i in range(len(found_sources_in)):
    #     if i != 1:
    #         continue
    #     print('frequency range', frequencies_search[i][0],frequencies_search[i][1])
    #     search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
    #     for j in range(len( pGB_injected[i][:10])):
    #         #subtract the found sources from original
    #         tdi_fs_subtracted = deepcopy(tdi_fs)
    #         for n in range(len( pGB_injected[i][:j])):
    #             Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=pGB_injected[i][n], oversample=4)
    #             source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
    #             index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
    #             index_high = index_low+len(Xs_subtracted)
    #             for k in ["X", "Y", "Z"]:
    #                 tdi_fs_subtracted[k].data[index_low:index_high] = tdi_fs_subtracted[k].data[index_low:index_high] - source_subtracted[k].data
    #         search_subtracted = Search(tdi_fs_subtracted,Tobs, frequencies_search[i][0], frequencies_search[i][1])
    #         print('true subtracted',search_subtracted.SNR([pGB_injected[i][j]]), 'original data', search1.SNR([pGB_injected[i][j]]))
    #     for j in range(len(found_sources_in[i])):
    #         #subtract the found sources from original
    #         tdi_fs_subtracted = deepcopy(tdi_fs)
    #         for n in range(len( found_sources_in[i][:j])):
    #             Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_in[i][n], oversample=4)
    #             source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
    #             index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
    #             index_high = index_low+len(Xs_subtracted)
    #             for k in ["X", "Y", "Z"]:
    #                 tdi_fs_subtracted[k].data[index_low:index_high] = tdi_fs_subtracted[k].data[index_low:index_high] - source_subtracted[k].data
    #         search_subtracted = Search(tdi_fs_subtracted,Tobs, frequencies_search[i][0], frequencies_search[i][1])
    #         print('found subtracted',search_subtracted.SNR([found_sources_in[i][j]]), 'original data', search1.SNR([found_sources_in[i][j]]))
    #         # print('found', search1.SNR([found_sources_mp_even_all[i][j]]))

    # number_of_found_signals = 0
    # for i in range(len(found_sources_in)):
    #     number_of_found_signals += len(found_sources_in[i])

def hamiltonian_monte_carlo(n_samples, negative_log_prob, grad_log_prob, initial_position, path_len=0.1, step_size=0.5):
    """Run Hamiltonian Monte Carlo sampling.

    Parameters
    ----------
    n_samples : int
        Number of samples to return
    negative_log_prob : callable
        The negative log probability to sample from
    initial_position : np.array
        A place to start sampling from.
    path_len : float
        How long each integration path is. Smaller is faster and more correlated.
    step_size : float
        How long each integration step is. Smaller is slower and more accurate.

    Returns
    -------
    np.array
        Array of length `n_samples`.
    """
    # autograd magic
    dVdq = grad_log_prob

    # collect all our samples in a list
    samples = [initial_position]

    # Keep a single object for momentum resampling
    momentum = scipy.stats.norm(0, 1)

    # If initial_position is a 10d vector and n_samples is 100, we want
    # 100 x 10 momentum draws. We can do this in one call to momentum.rvs, and
    # iterate over rows
    size = (n_samples,) + initial_position.shape[:1]
    for p0 in momentum.rvs(size=size):
        # Integrate over our path to get a new position and momentum
        q_new, p_new = leapfrog(
            samples[-1],
            p0,
            dVdq,
            path_len=path_len,
            step_size=step_size,
        )

        # Check Metropolis acceptance criterion
        start_log_p = negative_log_prob(samples[-1]) - np.sum(momentum.logpdf(p0))
        new_log_p = negative_log_prob(q_new) - np.sum(momentum.logpdf(p_new))
        if np.log(np.random.rand()) < start_log_p - new_log_p:
            samples.append(q_new)
        else:
            samples.append(np.copy(samples[-1]))

    return np.array(samples[1:])

def leapfrog(q, p, dVdq, path_len, step_size):
    """Leapfrog integrator for Hamiltonian Monte Carlo.

    Parameters
    ----------
    q : np.floatX
        Initial position
    p : np.floatX
        Initial momentum
    dVdq : callable
        Gradient of the velocity
    path_len : float
        How long to integrate for
    step_size : float
        How long each integration step should be

    Returns
    -------
    q, p : np.floatX, np.floatX
        New position and momentum
    """
    q, p = np.copy(q), np.copy(p)

    p -= step_size * dVdq(q) / 2  # half step
    for _ in range(int(path_len / step_size) - 1):
        q += step_size * p  # whole step
        p -= step_size * dVdq(q)  # whole step
    q += step_size * p  # whole step
    p -= step_size * dVdq(q) / 2  # half step

    # momentum flip at end
    return q, -p

class Posterior_computer():
    def __init__(self, tdi_fs, Tobs, frequencies, maxpGB, noise=None) -> None:
        self.tdi_fs = tdi_fs
        self.Tobs = Tobs
        self.frequencies = frequencies
        self.maxpGB = maxpGB
        self.search1 = Search(self.tdi_fs,self.Tobs, self.frequencies[0], self.frequencies[1], noise=noise, gb_gpu=gb_gpu)

    def reduce_boundaries(self, plot_confidance = False):
        fisher_information_matrix = self.search1.fisher_information(self.maxpGB)
        FIM = np.zeros((len(parameters),len(parameters)))
        for i,parameter1 in enumerate(parameters):
            for j,parameter2 in enumerate(parameters):
                FIM[i,j] = fisher_information_matrix[parameter1][parameter2]
        covariance_matrix = scipy.linalg.inv(FIM)
        maxpGB01 = scaleto01(self.maxpGB, self.search1.boundaries, self.search1.parameters, self.search1.parameters_log_uniform)

        # lambda_, v = scipy.linalg.eig(covariance_matrix)
        # transf = np.zeros((len(lambda_),len(lambda_)))
        # for i in range(len(lambda_)):
            # transf[i] = v[i] * np.sqrt(lambda_[i])
        # covariance2 = v @ np.diag(lambda_) @ scipy.linalg.inv(v)
        # transf = v @ np.diag(np.sqrt(lambda_))
        # scalematrix = np.max(np.abs(transf), axis=1)
        scalematrix = np.sqrt(np.diag(covariance_matrix))
        # print('scalematrix', scalematrix)
        maxpGB01_low = {}
        maxpGB01_high = {}
        self.boundaries_reduced = {}
        sigma_multiplier = 4
        for parameter in parameters:
            maxpGB01_low[parameter] = maxpGB01[parameter] - scalematrix[parameters.index(parameter)] * sigma_multiplier 
            maxpGB01_high[parameter] = maxpGB01[parameter] + scalematrix[parameters.index(parameter)] * sigma_multiplier 
            if parameter in [ 'Amplitude', 'Inclination']:
                maxpGB01_low[parameter] = maxpGB01[parameter] - scalematrix[parameters.index(parameter)] * sigma_multiplier
                maxpGB01_high[parameter] = maxpGB01[parameter] + scalematrix[parameters.index(parameter)] * sigma_multiplier
            if parameter in [ 'InitialPhase', 'Polarization']:
                maxpGB01_low[parameter] = maxpGB01[parameter] - 0.001
                maxpGB01_high[parameter] = maxpGB01[parameter] + 0.001
            if parameter in [ 'Frequency']:
                mulitplier = 1
                if np.abs(self.search1.boundaries['Frequency'][1] - self.search1.boundaries['Frequency'][0])  > 1e-4:
                    mulitplier = 0.2
                maxpGB01_low[parameter] = maxpGB01[parameter] - scalematrix[parameters.index(parameter)]  * mulitplier
                maxpGB01_high[parameter] = maxpGB01[parameter] + scalematrix[parameters.index(parameter)] * mulitplier
            if parameter == 'FrequencyDerivative':
                # print('scale',scalematrix[parameters.index(parameter)])
                # if scalematrix[parameters.index(parameter)] > 0.07:
                print('scale', scalematrix[parameters.index(parameter)], 'fd')
                range_fd = 1
                if self.search1.boundaries['FrequencyDerivative'][1] > 1e-14:
                    range_fd = 0.02
                # if self.search1.boundaries['FrequencyDerivative'][1] > 1e-17:
                #     range_fd = 0.4
                maxpGB01_low[parameter] = maxpGB01[parameter] - range_fd
                maxpGB01_high[parameter] = maxpGB01[parameter] + range_fd
            if maxpGB01_low[parameter] > maxpGB01_high[parameter]:
                placeholder = deepcopy(maxpGB01_low[parameter])
                maxpGB01_low[parameter] = deepcopy(maxpGB01_high[parameter])
                maxpGB01_high[parameter] = deepcopy(placeholder)
            if maxpGB01_low[parameter] < 0:
                maxpGB01_low[parameter] = 0
            if maxpGB01_high[parameter] > 1:
                maxpGB01_high[parameter] = 1
            maxpGB_fisher_low = (maxpGB01_low[parameter] * (self.search1.boundaries[parameter][1] - self.search1.boundaries[parameter][0])) + self.search1.boundaries[parameter][0]
            maxpGB_fisher_high = (maxpGB01_high[parameter] * (self.search1.boundaries[parameter][1] - self.search1.boundaries[parameter][0])) + self.search1.boundaries[parameter][0]
            self.boundaries_reduced[parameter] = [maxpGB_fisher_low, maxpGB_fisher_high]
        # self.boundaries_reduced = deepcopy(boundaries_reduced_fisher)

        # correct Frequency Derivative
        # split_fd = -17
        # if self.boundaries_reduced['FrequencyDerivative'][1] < split_fd+1:
        #     self.boundaries_reduced['FrequencyDerivative'][0] = -18.5
        #     self.boundaries_reduced['FrequencyDerivative'][1] = -16
        # elif self.boundaries_reduced['FrequencyDerivative'][0] < split_fd+0.5:
        #     self.boundaries_reduced['FrequencyDerivative'][0] = -18.5
        # correct Inclination and Amplitude
        # split_inclination = 0.95
        # if self.boundaries_reduced['Inclination'][1] > split_inclination or self.boundaries_reduced['Inclination'][0] < -split_inclination:
        #     if self.boundaries_reduced['Inclination'][1] > split_inclination and self.boundaries_reduced['Inclination'][0] < -split_inclination:
        #         if maxpGB01['Inclination'] > 0.5:
        #             self.boundaries_reduced['Inclination'][0] = 0.5
        #             self.boundaries_reduced['Inclination'][1] = 1
        #         else:
        #             self.boundaries_reduced['Inclination'][0] = -1
        #             self.boundaries_reduced['Inclination'][1] = -0.5
        #     if self.boundaries_reduced['Inclination'][1] > split_inclination:
        #         self.boundaries_reduced['Inclination'][0] = 0.5
        #         self.boundaries_reduced['Inclination'][1] = 1
        #     if self.boundaries_reduced['Inclination'][0] < -split_inclination:
        #         self.boundaries_reduced['Inclination'][0] = -1
        #         self.boundaries_reduced['Inclination'][1] = -0.5
        #     parameter = 'Amplitude'
        #     maxpGB01_low[parameter]  = maxpGB01[parameter] - 0.08
        #     maxpGB01_high[parameter]  = maxpGB01[parameter] + 0.08
        #     if maxpGB01_low[parameter] > maxpGB01_high[parameter]:
        #         placeholder = deepcopy(maxpGB01_low[parameter])
        #         maxpGB01_low[parameter] = deepcopy(maxpGB01_high[parameter])
        #         maxpGB01_high[parameter] = deepcopy(placeholder)
        #     if maxpGB01_low[parameter] < 0:
        #         maxpGB01_low[parameter] = 0
        #     if maxpGB01_high[parameter] > 1:
        #         maxpGB01_high[parameter] = 1
        #     maxpGB_fisher_low = (maxpGB01_low[parameter] * (self.search1.boundaries[parameter][1] - self.search1.boundaries[parameter][0])) + self.search1.boundaries[parameter][0]
        #     maxpGB_fisher_high = (maxpGB01_high[parameter] * (self.search1.boundaries[parameter][1] - self.search1.boundaries[parameter][0])) + self.search1.boundaries[parameter][0]
        #     self.boundaries_reduced[parameter] = [maxpGB_fisher_low, maxpGB_fisher_high]
        # print('boundaries reduced', self.boundaries_reduced)
    
    def reduce_boundaries_from_kde(self, resolution=10**6, proposal=None):
        test_x_m, probability = self.get_samples(resolution=10**6, proposal=proposal)
        min_parameters = np.min(test_x_m, axis=0)
        max_parameters = np.max(test_x_m, axis=0)
        self.boundaries_reduced_kde = deepcopy(self.boundaries_reduced)
        for i, parameter in enumerate(parameters):
            if parameter in ['InitialPhase', 'Polarization']:
                continue
            length = self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0] 
            self.boundaries_reduced_kde[parameter][0] = min_parameters[i]*length + self.boundaries_reduced[parameter][0]
            self.boundaries_reduced_kde[parameter][1] = max_parameters[i]*length + self.boundaries_reduced[parameter][0]
        # self.boundaries_reduced = self.boundaries_reduced_kde
        return self.boundaries_reduced_kde
    
    def reduce_boundaries_from_samples(self, samples):
        min_parameters = np.min(samples, axis=0)
        max_parameters = np.max(samples, axis=0)
        self.boundaries_reduced_samples = deepcopy(self.boundaries_reduced)
        for i, parameter in enumerate(parameters):
            if parameter in ['InitialPhase', 'Polarization']:
                continue
            length = self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0] 
            self.boundaries_reduced_samples[parameter][0] = min_parameters[i]*length + self.boundaries_reduced[parameter][0]
            self.boundaries_reduced_samples[parameter][1] = max_parameters[i]*length + self.boundaries_reduced[parameter][0]
        return self.boundaries_reduced_samples
    
    def train_model(self):
        rmse = 2
        train_size = 0
        test_size = 500
        added_trainig_size = 1000
        j = 0
        samples = np.random.rand(6000,8)
        samples[0] = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]
        samples[test_size] = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]
        while rmse > 0.6 and j < 5:
            j += 1
            train_size += added_trainig_size
            if j == 1:
                resolution = train_size + test_size
            else:
                resolution = added_trainig_size
            start = time.time()
            samples_likelihood = np.zeros(resolution)
            samples_likelihood2 = np.zeros(resolution)
            incl = np.zeros(resolution)
            for i in range(resolution):
                samples_p = scaletooriginal(samples[i+j*added_trainig_size], self.boundaries_reduced, self.search1.parameters, self.search1.parameters_log_uniform)
                incl[i] = samples_p['Inclination']
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
            # samples_flat = samples_flat*2-1
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
            kernel = RBF(length_scale=[1,2,5,1,1,1,1,1],length_scale_bounds=[(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,100),(0.1,100)])
            start = time.time()
            self.gpr = GaussianProcessRegressor(kernel=kernel, random_state=0).fit(train_x, train_y_normalized)
            print('train',time.time() - start)
            start = time.time()
            observed_pred_sk = self.gpr.predict(test_x)
            print('eval time of ', test_size, 'samples: ',time.time() - start)
            observed_pred_sk_scaled = observed_pred_sk*self.sigma + self.mu
            rmse = np.sqrt(mean_squared_error(test_y,observed_pred_sk_scaled))
            print("RMSE ",rmse,'with training size', len(train_y))
            if rmse > 30:
                print('high RMSE')

            # fig = plt.figure(figsize=(15,6))
            # plt.scatter(train_x[:,0], train_x[:,5], c=train_y, cmap='gray')
            # plt.show()
            if rmse > 5 and j == 5:
                j = 0
                samples = np.random.rand(6000,8)
                samples[0] = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]
                samples[test_size] = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]
                train_size = 0

    def evaluate(self, x):
        partial_length = 1*10**3
        # start = time.time()
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
        # print('eval time', time.time()-start)
        return observed_pred_mean
    
    def get_loglikelihood_gpu(self, x):
        partial_length = 1*10**3
        # start = time.time()
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
        # print('eval time', time.time()-start)
        return observed_pred_mean

    def sampler(self, resolution=1000, path_len=0.01, step_size=0.01):
        x0 = np.asarray([0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5])
        samples = hamiltonian_monte_carlo(n_samples=resolution, negative_log_prob=self.logp_func, grad_log_prob=self.dlogp_func, initial_position= x0, path_len=0.1, step_size=0.01)
        return samples

    def sample_dist(self, data0, data1, numPoints, resolution):
        data = np.array([data0, data1]).T
        grid, mypdf = FFTKDE(bw=0.01, kernel='gaussian').fit(data).evaluate(numPoints)
        axes = np.unique(grid[:, 0]), np.unique(grid[:, 1])
        mypdf = mypdf.reshape(numPoints, numPoints).T
        mypdf = np.clip(mypdf, a_min=0, a_max=None)

        dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
        data, pdfs = dist(resolution)
        return data, pdfs
    
    def sample_dist_fastkde(self, data0, data1, numPoints, resolution):
        ax = np.linspace(-0.15,1.15, numPoints)
        mypdf,axes = fastKDE.pdf(data0, data1, axes=[ax,ax])
        while np.min(mypdf) < 0:
            ax_expand = (np.random.rand(1)*0.2+0.1)[0]
            ax = np.linspace(-ax_expand,1+ax_expand,numPoints)
            mypdf,axes = fastKDE.pdf(data0, data1, axes=[ax,ax])

        dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
        data, pdfs = dist(resolution)
        return data, pdfs

    def get_samples(self, resolution=1000, proposal= None):
        numPoints = 2**5+1
        if proposal is None:
            test_x_m = np.random.uniform(size=(resolution,len(parameters)))
            probability = np.ones(resolution)
        else:
            probability = np.ones(resolution)
            test_x_m = np.zeros((resolution,len(parameters)))
            # [test_x_m[:,2],test_x_m[:,1]], pdfs = self.sample_dist(proposal[:,1],proposal[:,2], numPoints, resolution)
            # probability *= pdfs

            # [test_x_m[:,5],test_x_m[:,0]], pdfs = self.sample_dist(proposal[:,0],proposal[:,5], numPoints, resolution)
            # probability *= pdfs

            # [test_x_m[:,4],test_x_m[:,3]], pdfs = self.sample_dist(proposal[:,3],proposal[:,4], numPoints, resolution)
            # probability *= pdfs

            # plt.figure()
            # plt.hist2d(test_x_m[:,5],test_x_m[:,0])
            # plt.show()

            [test_x_m[:,2],test_x_m[:,1]], pdfs = self.sample_dist_fastkde(proposal[:,1],proposal[:,2], numPoints, resolution)
            probability *= pdfs

            [test_x_m[:,5],test_x_m[:,0]], pdfs = self.sample_dist_fastkde(proposal[:,0],proposal[:,5], numPoints, resolution)
            probability *= pdfs

            [test_x_m[:,4],test_x_m[:,3]], pdfs = self.sample_dist_fastkde(proposal[:,3],proposal[:,4], numPoints, resolution)
            probability *= pdfs

            # data = np.array([proposal[:,0],proposal[:,5]]).T
            # grid, mypdf = FFTKDE(bw=0.1, kernel='gaussian').fit(data).evaluate(numPoints)
            # axes = np.unique(grid[:, 0]), np.unique(grid[:, 1])
            # mypdf = mypdf.reshape(numPoints, numPoints).T
            # mypdf = np.clip(mypdf, a_min=0, a_max=None)
            # dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
            # data, pdfs = dist(resolution)
            # test_x_m[:,0] = data[1]
            # test_x_m[:,5] = data[0]
            # plt.figure()
            # plt.scatter(np.log10(test_x_m[:10000,0]),np.arccos(test_x_m[:10000,5]* (boundaries_reduced['Inclination'][1] - boundaries_reduced['Inclination'][0]) + boundaries_reduced['Inclination'][0]), c=pdfs[:10000])
            # probability *= pdfs

            ### fastkde
            # ax = np.linspace(-0.15,1.15,numPoints)
            # mypdf,axes = fastKDE.pdf(proposal[:,1],proposal[:,2], axes=[ax,ax])
            # #a pdf can not be negative
            # while np.min(mypdf) < 0:
            #     ax_expand = (np.random.rand(1)*0.2+0.1)[0]
            #     ax = np.linspace(-ax_expand,1+ax_expand,numPoints)
            #     mypdf,axes = fastKDE.pdf(data[0], data[1], axes=[ax,ax])

            # dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
            # data, pdfs = dist(resolution)
            # test_x_m[:,1] = data[1]
            # test_x_m[:,2] = data[0]
            # probability *= pdfs

            # ax = np.linspace(-0.15,1.15,numPoints)
            # mypdf,axes = fastKDE.pdf(proposal[:,0],proposal[:,5], axes=[ax,ax])
            # #a pdf can not be negative
            # while np.min(mypdf) < 0:
            #     ax_expand = (np.random.rand(1)*0.2+0.1)[0]
            #     ax = np.linspace(-ax_expand,1+ax_expand,numPoints)
            #     mypdf,axes = fastKDE.pdf(data[0], data[1], axes=[ax,ax])
            # dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
            # data, pdfs = dist(resolution)
            # test_x_m[:,0] = data[1]
            # test_x_m[:,5] = data[0]
            # probability *= pdfs

            # ax = np.linspace(-0.15,1.15,numPoints)
            # mypdf,axes = fastKDE.pdf(proposal[:,3],proposal[:,4], axes=[ax,ax])
            # #a pdf can not be negative
            # while np.min(mypdf) < 0:
            #     ax_expand = (np.random.rand(1)*0.2+0.1)[0]
            #     ax = np.linspace(-ax_expand,1+ax_expand,numPoints)
            #     mypdf,axes = fastKDE.pdf(data[0], data[1], axes=[ax,ax])
            # dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
            # data, pdfs = dist(resolution)
            # test_x_m[:,3] = data[1]
            # test_x_m[:,4] = data[0]
            # probability *= pdfs




            # mypdf,axes = fastKDE.pdf(proposal[:,3],proposal[:,4], axes=[ax,ax])
            # while not(np.all(mypdf>=0)):
            #     proposal = self.calculate_posterior(resolution = self.previous_resolution, proposal= self.previous_mcmc_samples, temperature= self.previous_T)
            #     mypdf,axes = fastKDE.pdf(proposal[:,3],proposal[:,4], axes=[ax,ax])
            # data = np.array([proposal[:,3],proposal[:,4]]).T
            # grid, mypdf = FFTKDE(bw=0.1, kernel='gaussian').fit(data).evaluate(numPoints)
            # axes = np.unique(grid[:, 0]), np.unique(grid[:, 1])
            # mypdf = mypdf.reshape(numPoints, numPoints).T
            # mypdf = np.clip(mypdf, a_min=0, a_max=None)
            # dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
            # data, pdfs = dist(resolution)
            # test_x_m[:,3] = data[1]
            # test_x_m[:,4] = data[0]
            # probability *= pdfs

            test_x_m[:,6] = np.random.uniform(size=resolution)
            test_x_m[:,7] = np.random.uniform(size=resolution)

            for n in range(8):
                index = np.where(test_x_m[:,n] > 0)
                test_x_m = test_x_m[index]
                probability = probability[index]
                index = np.where(test_x_m[:,n] < 1)
                test_x_m = test_x_m[index]
                probability = probability[index]
        return test_x_m, probability

    def calculate_posterior(self,resolution = 1*10**6, proposal= None, temperature = 1, pGB_true=None):
        test_x_m, probability = self.get_samples(resolution, proposal)
        start = time.time()

        # fig =  corner.corner(test_x_m,  bins=40, hist_kwargs={'density':True, 'lw':3}, plot_datapoints=False, fill_contours=False,  show_titles=True, \
        #                 color='#348ABD',  truth_color='k', use_math_test=True, \
        #                  levels=[0.9], title_kwargs={"fontsize": 12})
        # plt.show()

        use_gpu = True
        start = time.time()
        if use_gpu:
            # test_x_m[0] = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]
            test_x_m_scaled = scaletooriginal_array(test_x_m, self.boundaries_reduced, self.search1.parameters, self.search1.parameters_log_uniform)
            test_x_m_scaled = np.swapaxes(test_x_m_scaled, 0,1)

            params = np.array(
                [test_x_m_scaled[parameters.index('Amplitude')],
                test_x_m_scaled[parameters.index('Frequency')],
                test_x_m_scaled[parameters.index('FrequencyDerivative')],
                np.zeros_like(test_x_m_scaled[parameters.index('FrequencyDerivative')]),
                test_x_m_scaled[parameters.index('InitialPhase')],
                test_x_m_scaled[parameters.index('Inclination')],
                test_x_m_scaled[parameters.index('Polarization')],
                test_x_m_scaled[parameters.index('EclipticLongitude')],
                test_x_m_scaled[parameters.index('EclipticLatitude')]]
            )
            # for i, parameter in enumerate(['Amplitude', 'Frequency', 'FrequencyDerivative', 'InitialPhase', 'Inclination', 'Polarization', 'EclipticLongitude', 'EclipticLatitude']):
            #     j = 0
            #     if i > 2:
            #         j = 1
            #     params[i+j,0] = self.search1.pGBs[parameter]
            # params = np.array(
            #     [amp_in, f0_in, fdot_in, fddot_in, phi0_in, iota_in, psi_in, lam_in, beta_sky_in,]
            # )

            # params[1,0] = 0.00443235
            observed_pred_mean = self.search1.loglikelihood_gpu(params)
        else:
            observed_pred_mean = self.evaluate(test_x_m)

        # print('time loglikelihood for ', resolution, ' signals: ', time.time()-start)
        flatsamples = np.zeros(len(test_x_m))
        flatsamplesparameters = np.zeros((len(test_x_m),len(parameters)+1))
        i = 0
        flatsamples[:] = observed_pred_mean
        flatsamplesparameters[:,1:] = test_x_m
        flatsamplesparameters[:,0] = observed_pred_mean

        maxindx = np.unravel_index(flatsamplesparameters[:,0].argmax(), flatsamplesparameters[:,0].shape)
        max_parameters = flatsamplesparameters[maxindx[0],1:]
        # max_loglike = flatsamplesparameters[:,0].max()
        maxpGBpredicted = scaletooriginal(max_parameters, self.boundaries_reduced, self.search1.parameters, self.search1.parameters_log_uniform)
        # print(self.search1.loglikelihood_gpu_signal(pGB_dict_to_gpu_input([maxpGBpredicted])), self.search1.loglikelihood_gpu_signal(pGB_dict_to_gpu_input([self.maxpGB])))
        # print(self.search1.loglikelihood_gpu(pGB_dict_to_gpu_input([maxpGBpredicted])), self.search1.loglikelihood_gpu(pGB_dict_to_gpu_input([self.maxpGB])))
        # print(self.search1.loglikelihood([maxpGBpredicted]), self.search1.loglikelihood([self.maxpGB]))
        if self.search1.loglikelihood([maxpGBpredicted]) > self.search1.loglikelihood([self.maxpGB]):
            maxpGB = maxpGBpredicted
            # print('better maxpGB is found',maxpGB,self.search1.loglikelihood([maxpGBpredicted]), self.search1.loglikelihood([self.maxpGB]), self.maxpGB)
        
        best_value = np.max(observed_pred_mean)
        # print("pred", max_loglike, "true", self.search1.loglikelihood([scaletooriginal(max_parameters, self.boundaries_reduced)]), "max", self.search1.loglikelihood([self.maxpGB]), self.maxpGB)

        # np.random.shuffle(flatsamplesparameters)
        start = time.time()
        # normalizer = np.sum(np.exp(flatsamplesparameters[:,0]-best_value))
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
        # print('time MHMC', time.time()-start)
        print('acceptance rate',np.round(accepted/len(probability)*100),'%')

        return mcmc_samples
    

    def predict(self,x,k=0):
        #x of shape (m)
        
        #returns the gp predictions where f is the true function and
        #df, ddf, If, IIf are its first and second derivate respectively antiderivates
        #the outputs are the predictions f_p,df_p,ddf_p,If_p,IIf_p where
        #f(x) = f_p(x), df(x) = df_p(x), ddf(x) = ddf_p(x), If(x) = If_p(x) + C1, 
        #IIf(x) = IIf_p(x) + C1*x + C2 with some constants C1,C2
        #set k = 0 for the normal prediction, K = 1,2 for the first or second derivates
        #and k = -1,-2 for the first or second antiderivates
    
        # x = x.reshape(-1,1)
    
        X = x - self.gpr.X_train_
        l = self.gpr.kernel_.length_scale
        A = self.gpr.alpha_

        K_trans = self.gpr.kernel_(x, self.gpr.X_train_)
        y_mean = K_trans @ self.gpr.alpha_

        f = np.prod(np.exp(-(X)**2 / (2*l**2)), axis=1)
        df = f * (-X / l ** 2).T
        
        if k == 0: 
            return f @ A
        elif k == 1: 
            return df @ A
        else:
            raise Exception('Unknown parameter k: {}'.format(k))

    def gradient(self,x):
        step_length = 0.01
        grad = np.ones_like(x)
        for i in range(len(x)):
            x_predict = np.copy(x)
            x_predict[i] = x[i]-step_length/2
            if x_predict[i] < 0:
                x_predict[i] = 0
            low = self.gpr.predict([x_predict])
            x_predict[i] = x[i]+step_length/2
            if x_predict[i] > 1:
                x_predict[i] = 1
            high = self.gpr.predict([x_predict])
            grad[i] = (high-low)/step_length
        return grad

    def predict2(self, x):
        # gets 'l' used in denominator of expected value of gradient for RBF kernel 
        k2_l = self.gpr.kernel_.length_scale

        # not necessary to do predict, but now y_pred has correct shape
        y_pred, sigma = self.gpr.predict(np.asanyarray(x).reshape(1,-1),  return_std=True)

        # allocate array to store gradient
        y_pred_grad = 0.0*y_pred

        # set of points where gradient is to be queried
        # x = np.atleast_2d(np.linspace(-5, 0.8, 1000)).T
        
        X = self.gpr.X_train_
        x_star = x

        # eval_gradient can't be true when eval site doesn't match X
        # this gives standard RBF kernel evaluations
        k_val= self.gpr.kernel_(X, np.atleast_2d(x_star), eval_gradient=False).ravel()

        # x_i - x_star / l^2
        x_diff_over_l_sq = ((X-x_star)/np.power(k2_l,2)).ravel()

        # pair-wise multiply
        intermediate_result = np.multiply(k_val, x_diff_over_l_sq)

        # dot product intermediate_result with the alphas
        final_result = np.dot(intermediate_result, self.gpr.alpha_)

        # store gradient at this point
        y_pred_grad[key] = final_result
        return y_pred_grad

    def logp_func(self,x):
        return -(self.evaluate([x])[0] * self.sigma + self.mu)

    def dlogp_func(self,x, loc=0, scale=1):
        return -(self.predict(x,k=1) * self.sigma)

    def plot_corner(self, mcmc_samples, pGB = {}, save_figure = False, save_chain = False, number_of_signal = 0, parameter_titles = False, rescaled = False):
        start = time.time()
        if not(rescaled):
            mcmc_samples_rescaled = scaletooriginal_array(mcmc_samples, self.boundaries_reduced, self.search1.parameters, self.search1.parameters_log_uniform)
            # mcmc_samples_rescaled = np.zeros(np.shape(mcmc_samples))
            # i = 0
            # for parameter in parameters:
            #     if parameter in ["EclipticLatitude"]:
            #         mcmc_samples_rescaled[:,i] = np.arcsin((mcmc_samples[:,parameters.index(parameter)] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
            #     elif parameter in ["Inclination"]:
            #         mcmc_samples_rescaled[:,i] = np.arccos((mcmc_samples[:,parameters.index(parameter)] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
            #     elif parameter in parameters_log_uniform:
            #         mcmc_samples_rescaled[:,i] = 10**((mcmc_samples[:,parameters.index(parameter)] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
            #     else:
            #         mcmc_samples_rescaled[:,i] = (mcmc_samples[:,parameters.index(parameter)] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0]
            #     i += 1
            # print('time rescale', time.time()-start)
        else:
            mcmc_samples_rescaled = mcmc_samples
        
        mcmc_samples_rescaled[:,parameters.index('Frequency')] *= 10**3 
        save_frequency = self.maxpGB['Frequency']

        if save_chain:
            df = pd.DataFrame(data=mcmc_samples_rescaled, columns=parameters)
            for parameter in ['Amplitude','FrequencyDerivative','EclipticLatitude','EclipticLongitude','Inclination','InitialPhase','Polarization']:
                df[parameter] = df[parameter].astype('float32')
            df.to_csv(SAVEPATH+'Chains_gpu/frequency'+str(int(np.round(save_frequency*10**12)))+'pHz'+save_name+'.csv',index=False)
        mcmc_samples_rescaled[:,parameters.index('Frequency')] /= 10**3 
        # start = time.time()
        # df = pd.DataFrame(data=mcmc_samples_rescaled, columns=parameters)
        # df.to_csv('/home/stefan/Repositories/ldc1_evaluation_data/submission/Stefan_LDC14/GW'+str(int(np.round(maxpGB['Frequency']*10**8)))+'.csv',index=False)
        # print('saving time', time.time()-start)

        if save_figure:
            # print('full time', time.time()-first_start)
            datS = np.zeros(np.shape(mcmc_samples))
            datS[:,0] = mcmc_samples_rescaled[:,2]
            datS[:,1] = np.sin(mcmc_samples_rescaled[:,1])
            datS[:,2] = mcmc_samples_rescaled[:,3]*10**3
            datS[:,3] = mcmc_samples_rescaled[:,4]
            datS[:,4] = np.cos(mcmc_samples_rescaled[:,5])
            datS[:,5] = np.log10(mcmc_samples_rescaled[:,0])
            datS[:,6] = mcmc_samples_rescaled[:,6]
            datS[:,7] = mcmc_samples_rescaled[:,7]
            lbls = [r'\lambda', r'\sin \beta', 'f$ $($mHz$)', r'\dot{f}$ $ ($Hz/s$)', r'\cos \iota', r'A', r'\phi', r'\Phi']

            if pGB:
                tr_s = np.zeros(len(parameters))
                i = 0
                for parameter in ['EclipticLongitude','EclipticLatitude','Frequency','FrequencyDerivative','Inclination','Amplitude','InitialPhase','Polarization']:
                    if parameter in parameters_log_uniform:
                        tr_s[i] = np.log10(pGB[parameter])
                    elif parameter in ['Frequency']:
                        tr_s[i] = pGB[parameter]*10**3
                    elif parameter in ['Inclination']:
                        tr_s[i] = np.cos(pGB[parameter])
                    elif parameter in ['EclipticLatitude']:
                        tr_s[i] = np.sin(pGB[parameter])
                    else:
                        tr_s[i] = pGB[parameter]
                    i += 1
                if tr_s[0] > np.pi:
                    tr_s[0] -= 2*np.pi
            maxvalues = np.zeros(len(parameters))
            i = 0
            for parameter in ['EclipticLongitude','EclipticLatitude','Frequency','FrequencyDerivative','Inclination','Amplitude','InitialPhase','Polarization']:
                if parameter in parameters_log_uniform:
                    maxvalues[i] = np.log10(self.maxpGB[parameter])
                elif parameter in ['Frequency']:
                    maxvalues[i] = self.maxpGB[parameter]*10**3
                elif parameter in ['Inclination']:
                    maxvalues[i] = np.cos(self.maxpGB[parameter])
                elif parameter in ['EclipticLatitude']:
                    maxvalues[i] = np.sin(self.maxpGB[parameter])
                else:
                    maxvalues[i] = self.maxpGB[parameter]
                i += 1

            # rng = []
            # for i in range(len(lbls)):
            #     minrange = min(datS[:,i].min(), tr_s[i])
            #     maxrange = max(datS[:,i].max(), tr_s[i])
            #     range_width = np.abs(maxrange - minrange)
            #     oner = ( minrange- range_width/10, maxrange + range_width/10)
            #     rng.append(oner)
            # Get the getdist MCSamples objects for the samples, specifying same parameter
            # names and labels; if not specified weights are assumed to all be unity
            ndim = 6
            names = ['EclipticLongitude','EclipticLatitude','Frequency','FrequencyDerivative','Inclination','Amplitude']
            labels =  lbls[:6]
            samples = MCSamples(samples=datS[:,:6],names = names, labels = labels)

            g = plots.get_subplot_plotter(subplot_size=0.9)
            samples.updateSettings({'contours': [0.68, 0.95]})
            g.settings.num_plot_contours = 2
            if parameter_titles:
                g.triangle_plot([samples], shaded=True, title_limit=2)
            else:
                g.triangle_plot([samples], shaded=True)
            
            #markers vertical
            for i in range(ndim):
                for ax in g.subplots[i:,i]:
                    if pGB:
                        xlim = ax.get_xlim()
                        ax.set_xlim(np.min([xlim[0], tr_s[i]]),np.max([xlim[1], tr_s[i]]))
                        xlim = ax.get_xlim()
                        if xlim[0] + (xlim[1]-xlim[0])*0.1 > tr_s[i]:
                            ax.set_xlim(xlim[0] - (xlim[1]-xlim[0])*0.1, xlim[1])
                        if xlim[1] - (xlim[1]-xlim[0])*0.1 < tr_s[i]:
                            ax.set_xlim(xlim[0], xlim[1] + (xlim[1]-xlim[0])*0.1)
                        ax.axvline(tr_s[i], color='red', lw = 1)
                    ax.axvline(maxvalues[i], color='green', ls='--', lw = 1)
                # i += 1
            #markers horizontal
            for i in range(ndim):
                for ax in g.subplots[i,:i]:
                    if pGB:
                        ylim = ax.get_ylim()
                        ax.set_ylim(np.min([ylim[0], tr_s[i]]),np.max([ylim[1], tr_s[i]]))
                        ylim = ax.get_ylim()
                        if ylim[0] + (ylim[1]-ylim[0])*0.1 > tr_s[i]:
                            ax.set_ylim(ylim[0] - (ylim[1]-ylim[0])*0.1, ylim[1])
                        if ylim[1] - (ylim[1]-ylim[0])*0.1 < tr_s[i]:
                            ax.set_ylim(ylim[0], ylim[1] + (ylim[1]-ylim[0])*0.1)
                        ax.axhline(tr_s[i], color='red', lw = 1)
                    ax.axhline(maxvalues[i], color='green', ls='--', lw = 1)
                # i += 1
            g.export(SAVEPATH+'Posteriors_gpu/frequency'+ str(int(np.round(save_frequency*10**9)))+save_name+str(parameter_titles)+'.png')
            # plt.show()

def compute_posterior(tdi_fs, Tobs, frequencies, maxpGB, pGB_true = [], number_of_signal = 0, noise=None):
    start = time.time()
    posterior1 = Posterior_computer(tdi_fs, Tobs, frequencies, maxpGB, noise=noise)
    print('time to initialize: ', time.time()-start)
    print(posterior1.search1.loglikelihood([maxpGB]))
    # posterior1.search1.update_noise(pGB=maxpGB)
    # print(posterior1.search1.loglikelihood([maxpGB]))
    # print('SNR found',posterior1.search1.SNR([maxpGB]))
    # posterior1.search1.update_noise(pGB=maxpGB)
    # print('SNR found full noise',posterior1.search1.SNR([maxpGB]))
    # if pGB_true:
    #     print('SNR injected',posterior1.search1.SNR([pGB_true]))
    # posterior1.search1.plot(found_sources_in=[maxpGB])
    start = time.time()
    posterior1.reduce_boundaries()
    print('time to reduce boundaries: ', time.time()-start)
    start = time.time()
    # posterior1.train_model()
    # print('GPU memory free', get_gpu_memory())
    # mcmc_samples = posterior1.calculate_posterior(resolution = 1*10**5, temperature= 10)
    temperature = [15, 10, 5, 3, 2, 1]
    resolution = [10**4, 10**4, 10**4, 10**4, 10**4, 10**4]
    posterior_round = 0
    mcmc_samples = posterior1.calculate_posterior(resolution = resolution[posterior_round], temperature= temperature[posterior_round])
    posterior_round += 1
    # posterior1.boundaries_reduced = posterior1.reduce_boundaries_from_samples(mcmc_samples)

    mcmc_samples = posterior1.calculate_posterior(resolution = resolution[posterior_round], proposal= mcmc_samples, temperature= temperature[posterior_round])
    posterior_round += 1
    
    mcmc_samples = posterior1.calculate_posterior(resolution = resolution[posterior_round], proposal= mcmc_samples, temperature= temperature[posterior_round])
    posterior_round += 1

    mcmc_samples = posterior1.calculate_posterior(resolution = resolution[posterior_round], proposal= mcmc_samples, temperature= temperature[posterior_round])
    posterior_round += 1

    mcmc_samples = posterior1.calculate_posterior(resolution = resolution[posterior_round], proposal= mcmc_samples, temperature= temperature[posterior_round])
    posterior_round += 1

    mcmc_samples_final = posterior1.calculate_posterior(resolution = resolution[posterior_round], proposal= mcmc_samples, temperature= temperature[posterior_round])
    # posterior1.plot_corner(mcmc_samples_final, pGB_true, save_figure= True)
    print('time to compute posterior: ', time.time()-start)
    start = time.time()
    posterior1.plot_corner(mcmc_samples_final, pGB_true, save_figure= False, save_chain= True, number_of_signal = 0, parameter_titles = False)
    # kde_samples, probability = posterior1.get_samples(resolution = 1*10**4, proposal= mcmc_samples)

    # plt.figure()
    # plt.scatter(kde_samples[:100000,0],np.arccos(kde_samples[:100000,5]* (posterior1.boundaries_reduced['Inclination'][1] - posterior1.boundaries_reduced['Inclination'][0]) + posterior1.boundaries_reduced['Inclination'][0]), c=probability[:100000])
    # plt.show()

    # posterior1.plot_corner(kde_samples, pGB_true, save_figure= True, save_chain= True, number_of_signal = 0, parameter_titles = False)
    # print('time to save posterior: ', time.time()-start)
    return mcmc_samples


print('GPU memory free', get_gpu_memory())

SAVEPATH_sangria = grandparent+"/LDC/pictures/Sangria/"
end_string = '_SNR_scaled_03_injected_snr5'
# found_sources = np.load(SAVEPATH+'/found_sources_matched' +save_name+end_string+'.npy', allow_pickle=True)
# found_sources_not_matched = np.load(SAVEPATH+'/found_sources_not_matched' +save_name+'.npy', allow_pickle=True)
# pGB_injected_not_matched = np.load(SAVEPATH+'/injected_not_matched_windows' +save_name+'.npy', allow_pickle=True)
# pGB_injected_matched = np.load(SAVEPATH+'/injected_matched_windows' +save_name+'.npy', allow_pickle=True)

# found_sources_mp = np.load(SAVEPATH+'/found_sourcesSangria_1_full.npy', allow_pickle=True)
found_sources_in_flat = np.load(SAVEPATH+'found_sources_' +save_name+'_flat.npy', allow_pickle = True)


amp = 1e-22  # amplitude
f0 = 2e-3  # f0
fdot = 1e-16  # fdot
fddot = 0.0
phi0 = 0.1  # initial phase
iota = 0.2  # inclination
psi = 0.3  # polarization angle
lam = 0.4  # ecliptic longitude
beta_sky = 0.5  # ecliptic latitude
window_length = 1e-5

pGB = {'Amplitude': amp, 'Frequency': f0, 'FrequencyDerivative': fdot, 'InitialPhase': phi0, 'Inclination': iota, 'Polarization': psi, 'EclipticLongitude': lam, 'EclipticLatitude': beta_sky}
# start = time.time()
# search1 = Search(tdi_fs,Tobs, f0-window_length, f0+window_length)
# print('search init time:', time.time() - start)

tdi_fs_added = deepcopy(tdi_fs)
Xs_added, Ys_added, Zs_added = GB.get_fd_tdixyz(template=pGB, oversample=4)
source_added = dict({"X": Xs_added, "Y": Ys_added, "Z": Zs_added})
index_low = np.searchsorted(tdi_fs_added["X"].f, Xs_added.f[0])
index_high = index_low+len(Xs_added)
for k in ["X", "Y", "Z"]:
    tdi_fs_added[k].data[index_low:index_high] = tdi_fs_added[k].data[index_low:index_high] + 2* source_added[k].data

# start = time.time()
# search1 = Search(tdi_fs_added,Tobs, f0-window_length, f0+window_length)
# print('search init time:', time.time() - start)


# print(search1.loglikelihood([pGB]))
# print(search1.loglikelihood_gpu(np.array([pGB_dict_to_gpu_input([pGB]),pGB_dict_to_gpu_input([pGB])])))
# print(search1.loglikelihood_gpu_signal(pGB_dict_to_gpu_input([pGB])))

# A_inj, E_inj = gb.inject_signal(*pGB_dict_to_gpu_input([pGB]), dt=dt,T=Tobs)
Af = (Zs_added - Xs_added)/np.sqrt(2.0)
Ef = (Zs_added - 2.0*Ys_added + Xs_added)/np.sqrt(6.0)

tdi_fs_added['A'] = (tdi_fs_added['Z'] - tdi_fs_added['X'])/np.sqrt(2.0)
indexes = np.logical_and(tdi_fs_added.f >= Af.f.values[0], tdi_fs_added.f <= Af.f.values[-1]) 


# plt.figure()
# plt.semilogy(Af.f,np.abs(A_inj[:-1][indexes]))
# plt.semilogy(tdi_fs_added.f,np.abs(tdi_fs_added['A']))
# plt.semilogy(Af.f,np.abs(Af))
# plt.semilogy(Af.f,np.abs(A_inj[:-1][indexes]-Af))
# # plt.semilogy(Af.f,np.abs(data[0][indexes].get()))
# plt.xlim(f0 - window_length, f0 + window_length)
# plt.show()


# found_sources_in_flat = []
# for i in range(len(found_sources)):
#     for j in range(len(found_sources[i])):
#         found_sources_in_flat.append(found_sources[i][j])
# for i in range(len(found_sources)):
#     for j in range(len(found_sources[i][3])):
#         found_sources_in_flat.append(found_sources[i][3][j])
found_sources_in_flat = np.asarray(found_sources_in_flat)
found_sources_in_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_in_flat]) for attribute in found_sources_in_flat[0].keys()}
found_sources_in_flat_df = pd.DataFrame(found_sources_in_flat_array)
found_sources_in_flat_df = found_sources_in_flat_df.sort_values('Frequency')
found_sources_in = []
for i in range(len(frequencies_search)):
    found_sources_in.append(found_sources_in_flat_df[(found_sources_in_flat_df['Frequency'] > frequencies_search[i][0]) & (found_sources_in_flat_df['Frequency'] < frequencies_search[i][1])].to_dict(orient='records'))
# found_sources_matched =found_sources_in

# pGB_injected = pGB_injected_matched

# number_of_signals = 0
# for i in range(len(found_sources_in)):
#     number_of_signals += len(found_sources_in[i])
# number_of_signals2 = 0
# for i in range(len(found_sources_not_matched)):
#     number_of_signals2 += len(found_sources_not_matched[i])

# parameter = 'EclipticLongitude'
# for i in range(len(pGB_injected)):
#     for j in range(len(pGB_injected[i])):
#         if pGB_injected[i][j][parameter] > np.pi:
#             pGB_injected[i][j][parameter] -= 2*np.pi

better_pGB = {'Amplitude': 3.7080776510756e-23, 'EclipticLatitude': -0.0864329194471405, 'EclipticLongitude': -1.5608489415225566, 'Frequency': 0.011097063538503463, 'FrequencyDerivative': 3.795997584356877e-15, 'Inclination': 1.3544536642993756, 'InitialPhase': 3.802341846303522, 'Polarization': 3.0450807858161113}
start = time.time()
number_of_signals = 0
frequencies_found = []
# LDC1-4 ####################
posterior_calculation_input = []
# batch = int(sys.argv[1])

# lower_window = 7210

# plt.figure()
# plt.loglog(np.array(frequencies)[:,1],frequency_derivative(np.array(frequencies)[:,1], M_chirp_upper_boundary))
# plt.show()

do_subtract = True
if do_subtract:
    start = time.time()
    # found_sources_flat = found_sources_in_flat

    # found_sources_flat = []
    # for i in range(len(found_sources_mp_subtract)):
    #     for j in range(len(found_sources_mp_subtract[i])):
    #         found_sources_flat.append(found_sources_mp_subtract[i][j])
    found_sources_flat = np.asarray(found_sources_in_flat)
    found_sources_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_flat]) for attribute in found_sources_flat[0].keys()}
    found_sources_flat_df = pd.DataFrame(found_sources_flat_array)
    found_sources_flat_df = found_sources_flat_df.sort_values('Frequency')
    # for i in range(len(frequencies_search_full)):
    #     found_sources_flat_df = found_sources_flat_df[(found_sources_flat_df['Frequency']< frequencies_search_full[i][0]) | (found_sources_flat_df['Frequency']> frequencies_search_full[i][1])]
    # found_sources_flat_df = found_sources_flat_df.sort_values('Frequency')
    found_sources_out_flat = found_sources_flat_df.to_dict(orient='records')
    tdi_fs_residual = tdi_subtraction(tdi_fs,found_sources_out_flat, frequencies_search_full)
    tdi_fs_partial_residual = tdi_subtraction(tdi_fs,found_sources_out_flat, frequencies_search_full, ratio = 0.7)

    print('subtraction time', time.time()-start)
    plot_subtraction = False
    if plot_subtraction:
        i = 4000
        lower_frequency = frequencies_search[i][0]
        upper_frequency = frequencies_search[i][1]
        search1 = Search(tdi_fs,Tobs, lower_frequency, upper_frequency)
        search1.plot(second_data= tdi_fs_residual)#, found_sources_in=found_sources_out_flat)
        # search1.plot(second_data= tdi_fs_residual, found_sources_in=found_sources_mp_o[start_index][0])
        
    # tdi_fs = deepcopy(tdi_fs_residual)

tdi_fs['A'] = (tdi_fs["Z"] - tdi_fs["X"])/np.sqrt(2.0)

tdi_fs_residual['A'] = (tdi_fs_residual["Z"] - tdi_fs_residual["X"])/np.sqrt(2.0)



def get_noise_from_frequency_domain(tdi_fs, number_of_windows=100):
    tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft()) for k in ["X", "Y", "Z"]]))
    tdi_ts["E"] = deepcopy(tdi_ts["X"])
    tdi_ts["A"] = (tdi_ts["Z"] - tdi_ts["X"])/np.sqrt(2.0)
    tdi_ts["E"].values = (tdi_ts["Z"] - 2.0*tdi_ts["Y"] + tdi_ts["X"])/np.sqrt(6.0)
    tdi_ts["T"] = (tdi_ts["Z"] + tdi_ts["Y"] + tdi_ts["X"])/np.sqrt(3.0)
    f, psdA =  scipy.signal.welch(tdi_ts["A"], fs=1.0/dt, nperseg=len(tdi_ts["X"])/number_of_windows, average='mean', window= 'boxcar')
    f, psdE =  scipy.signal.welch(tdi_ts["E"], fs=1.0/dt, nperseg=len(tdi_ts["X"])/number_of_windows)
    # f2, psdE2 =  scipy.signal.welch(tdi_ts["E"], fs=1.0/dt, nperseg=len(tdi_ts["X"]), scaling='spectrum')
    f, psdT =  scipy.signal.welch(tdi_ts["T"], fs=1.0/dt, nperseg=len(tdi_ts["X"])/number_of_windows)
    return f, psdA, psdE, psdT

number_of_windows = 500
start = time.time()
f, psdA, psdE, psdT = get_noise_from_frequency_domain(tdi_fs, number_of_windows=number_of_windows)
print(time.time()-start)
# psdA = np.interp(tdi_fs.f, f, psdA)
# psdE = np.interp(tdi_fs.f, f, psdE)
# psdT = np.interp(tdi_fs.f, f, psdT)

f_partial_residual, psdA_partial, psdE_partial, psdT_partial = get_noise_from_frequency_domain(tdi_fs_partial_residual, number_of_windows=number_of_windows)
# psdA_partial_residual = np.interp(tdi_fs.f, f_partial_residual, psdA_partial_residual)
# psdE_partial_residual = np.interp(tdi_fs.f, f_partial_residual, psdE_partial_residual)
# psdT_partial_residual = np.interp(tdi_fs.f, f_partial_residual, psdT_partial_residual)


f_res, psdA_residual, psdE_residual, psdT_residual = get_noise_from_frequency_domain(tdi_fs_residual, number_of_windows=number_of_windows)
# psdA_residual = np.interp(tdi_fs.f, f_res, psdA_residual)
# psdE_residual = np.interp(tdi_fs.f, f_res, psdE_residual)
# psdT_residual = np.interp(tdi_fs.f, f_res, psdT_residual)

i = 4000
lower_frequency = frequencies_search[i][0]
upper_frequency = frequencies_search[i][1]
search1 = Search(tdi_fs_residual,Tobs, lower_frequency, upper_frequency)

### mean in each frequency window
# noise_means = deepcopy(search1.SA_full_f)
# noise_means2 = []
# frequencies_means = []
# for i in range(len(frequencies)):
#     padding = (frequencies[i][1] - frequencies[i][0])/2
#     lower_index = np.searchsorted(search1.DAf_full_f.f, frequencies[i][0]-padding)
#     upper_index = np.searchsorted(search1.DAf_full_f.f, frequencies[i][1]+padding)
#     lower_index_inner = np.searchsorted(search1.DAf_full_f.f, frequencies[i][0])
#     upper_index_inner = np.searchsorted(search1.DAf_full_f.f, frequencies[i][1])
#     # indexes = np.logical_and(tdi_fs['X'].f >= frequencies[i][0]-padding, tdi_fs['X'].f < frequencies[i][1]+padding) 
#     DAf = search1.DAf_full_f[lower_index:upper_index]
#     noise_means2.append(np.median((np.abs(DAf).data)**2)/(Tobs)*2)
#     frequencies_means.append(search1.DAf_full_f.f[lower_index+int((upper_index-lower_index)/2)].data)
#     noise_means[lower_index_inner:upper_index_inner] = np.ones_like(search1.DAf_full_f[lower_index_inner:upper_index_inner])*np.median((np.abs(DAf).data)**2)/(Tobs)*2
    # noise_means[lower_index:upper_index] = ((np.abs(DAf).data))**2/(Tobs)*2




# psdA_fit = median_windows(psdA,20)
# psdE_fit = median_windows(psdE,20)
# psdT_fit = median_windows(psdT,20)
psdA_fit, smoothedA = smooth_psd(psdA_partial, f_res)
psdE_fit, smoothedE = smooth_psd(psdE_partial, f_res)
psdT_fit, smoothedT = smooth_psd(psdT_partial, f_res)

# psdA_fit = smooth_psd(noise_means, frequencies_means)

# psdA_fit = np.interp(tdi_fs.f, f_res, psdA_fit)
# psdE_fit = np.interp(tdi_fs.f, f_res, psdE_fit)
# psdT_fit = np.interp(tdi_fs.f, f_res, psdT_fit)

psdA_fit = scipy.interpolate.InterpolatedUnivariateSpline(f_res, psdA_fit)(tdi_fs.f)
psdE_fit = scipy.interpolate.InterpolatedUnivariateSpline(f_res, psdE_fit)(tdi_fs.f)
psdT_fit = scipy.interpolate.InterpolatedUnivariateSpline(f_res, psdT_fit)(tdi_fs.f)

# t = np.linspace(f_res[8], f_res[-2], 50)
# # t = np.append( f_res[1:10], t)
# spline = scipy.interpolate.LSQUnivariateSpline(f_res, smoothed_A, t)
# spline.set_smoothing_factor(10**5)



SA_full_f = Nmodel.psd(freq=tdi_fs.f, option="A")
df = tdi_fs['A'].f[1] - tdi_fs['A'].f[0]
fs = 1/dt

N = len(tdi_ts['X'])
tdi_fsX = np.fft.rfft(tdi_ts['X'].data)[0:N//2+1]/fs

noise_model = 'sangria'
Nmodel = get_noise_model(noise_model, f_res[1:])
Sn = Nmodel.psd(freq=f_res[1:], option="X")
SA = Nmodel.psd(freq=f_res[1:], option="A")
SE = Nmodel.psd(freq=f_res[1:], option="E")
ST = Nmodel.psd(freq=f_res[1:], option="T")

ldc_noise = AnalyticNoise(f_res[1:], model="sangria", wd=1)
SAa = ldc_noise.psd(f_res[1:], option='A')

plt.figure()
# plt.plot(xnew, np.sin(xnew), '-.', label='sin(x)')
# plt.plot(xnew, BSpline(*tck)(xnew), '-', label='s=0')
# plt.loglog(f_res, psdA)   
plt.loglog(f_res, psdA_residual, label='welch')  
# plt.loglog(f_res, psdA_partial)   
# plt.loglog(f_res, smoothed_A)     
# plt.loglog(tdi_fs.f, spline(tdi_fs.f))
# plt.loglog(frequencies_means, psd_fit, zorder=5)    
# plt.loglog(f_res[1:], SA, zorder=5)   
plt.loglog(f_res, smoothedA, zorder=4, label='median window smoothed')   
plt.loglog(tdi_fs.f, psdA_fit, zorder=5, label='estimate')   
plt.loglog(f_res[1:], SAa, 'k--', zorder=5, label='instrument')   
# plt.loglog(tdi_fs.f, noise_means, zorder=3)  
# plt.loglog(tdi_fs.f, np.abs(tdi_fs_residual['A']/dt)**2 /(fs*len(tdi_ts['X']))*2 ,'.', zorder =1)     
# plt.loglog(tdi_fs.f, tdi_fs['X'],'.', zorder =1)   
# plt.loglog(tdi_fs.f, tdi_fsX,'.', zorder =1)    
# plt.loglog(tdi_fs.f, (psdA_residual-np.abs(tdi_fs_residual['A'])**2 /(Tobs)*2)/psdA_residual,'.')    
# plt.loglog(f_res[peaks], psdA_residual[peaks],'x', color='red')     
# plt.loglog(f_res[lower_index_res:], y_pred)     
plt.legend()
plt.xlim(0.0003,0.1)
plt.ylim(10**-43,10**-37)
plt.show()


# psdA_residual = (np.abs(search1.DAf_full_f)).data**2*2/(dt*len(search1.tdi_fs['X']))
# psdE_residual = (np.abs(search1.DEf_full_f)).data**2/(dt*len(search1.tdi_fs['X']))*2
# psdT_residual = (np.abs(search1.DTf_full_f)).data**2/(dt*len(search1.tdi_fs['X']))*2


noise_data = {'A': psdA, 'E': psdE, 'T': psdT}
noise_fit = {'A': psdA_fit, 'E': psdE_fit, 'T': psdT_fit}
noise_residual = {'A': psdA_residual, 'E': psdE_residual, 'T': psdT_residual}
# noise_partial_residual = {'A': psdA_partial_residual, 'E': psdE_partial_residual, 'T': psdT_partial_residual}
noise = noise_fit

# noise['f'] = tdi_fs.f
# noise_df = pd.DataFrame(noise)
# noise_df.to_csv(SAVEPATH+'ETH_sangria_noise.csv')

search1 = Search(tdi_fs_residual,Tobs, lower_frequency, upper_frequency, noise=noise)

whitened_data = np.abs(tdi_fs_residual['A']/dt)**2 /(fs*len(tdi_ts['X']))*2 / noise['A']/dt
whitened_data = tdi_fs_residual['A']*2/dt / np.sqrt(noise['A']*fs*len(tdi_ts['X']))

psdA_residual_interpol = scipy.interpolate.InterpolatedUnivariateSpline(f_res, psdA_residual)(tdi_fs.f)
whitened_data = tdi_fs_residual['A']*2/dt / np.sqrt(psdA_residual_interpol*fs*len(tdi_ts['X']))
# whitened_data = search1.DAf_full_f / np.real(search1.DAf_full_f)

lower_index = np.searchsorted(search1.tdi_fs['X'].f, 0.0003)
upper_index = np.searchsorted(search1.tdi_fs['X'].f, 0.1)
print(np.mean(np.real(whitened_data[lower_index:upper_index])))
print(np.std(np.real(whitened_data[lower_index:upper_index])))
print(np.std(np.imag(whitened_data[lower_index:upper_index])))

plt.figure()
plt.hist(np.real(whitened_data[lower_index:upper_index]), bins=2000, density=True)
plt.xlim(-4,4)
plt.show()



i = 3000
# search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1], noise=noise)
plt.figure()
# plt.loglog(tdi_fs_residual.f, search1.SA_full_f, 'k')
# plt.loglog(tdi_fs_residual.f, psdA)
# plt.loglog(tdi_fs_residual.f, psdA_residual*(dt*len(search1.tdi_fs['X'])), zorder = 5)
# plt.loglog(tdi_fs_residual.f,  np.sqrt(psdA_residual), zorder = 5)
# plt.loglog(tdi_fs_residual.f,  np.sqrt(psdE_residual), zorder = 5)
# plt.loglog(tdi_fs_residual.f,  np.sqrt(noise_means), zorder = 5)
# plt.loglog(tdi_fs_residual.f, psdA_partial_residual)
plt.semilogx(tdi_fs_residual.f[100:], np.real(whitened_data[100:]))
# plt.loglog(tdi_fs.f, np.sqrt(psdA_interp),'-x')
# plt.loglog(f2, np.sqrt(psdA_residual2))
# plt.loglog(tdi_fs.f, np.abs(tdi_fs['A']))
# plt.semilogx(tdi_fs.f, np.abs(search1.DAf_full_f),'k')
# search1.update_noise()
# plt.loglog(tdi_fs.f, np.sqrt(search1.SA_full_f),'r')
# plt.loglog(tdi_fs_residual.f, np.abs(tdi_fs_residual['A']))
plt.xlim(0.0003, 0.1)
plt.show()


tdi_fs_original = deepcopy(tdi_fs)
dataX_original = tdi_fs["X"]
dataY_original = tdi_fs["Y"]
dataZ_original = tdi_fs["Z"]
DAf_original = (dataZ_original - dataX_original)/np.sqrt(2.0)

start_all = time.time()
for i in range(len(found_sources_in)):
    # if i < 2898:
    #     continue
    # if i%10 != 0:
    #     continue
    # if i in [0,len(found_sources_in)-1]:
    #     continue
    for j in range(len(found_sources_in[i])):
        start = time.time()
        # if j < 3 and i == 7213:
        #     continue

        print(i)
        #subtract the found sources of neighbours and own window from original except the signal itself
        tdi_fs_added = deepcopy(tdi_fs_residual)
        # tdi_fs_added_part = deepcopy(tdi_fs)

        Xs_added, Ys_added, Zs_added = GB.get_fd_tdixyz(template=found_sources_in[i][j], oversample=4)
        source_added = dict({"X": Xs_added, "Y": Ys_added, "Z": Zs_added})
        index_low = np.searchsorted(tdi_fs_added["X"].f, Xs_added.f[0])
        index_high = index_low+len(Xs_added)
        for k in ["X", "Y", "Z"]:
            tdi_fs_added[k].data[index_low:index_high] = tdi_fs_added[k].data[index_low:index_high] + source_added[k].data
            # tdi_fs_added_part[k].data[index_low:index_high] = tdi_fs_added_part[k].data[index_low:index_high] - source_added[k].data * 1


        plot_addition = False
        if plot_addition:
            # pGB_true = []
            # for k in range(len(pGB_injected[i])):
            #     pGB_true.append(pGB_injected[i][k].to_dict())
            # pGB_true_not_matched = []
            # for k in range(len(pGB_injected_not_matched[i])):
            #     pGB_true_not_matched.append(pGB_injected_not_matched[i].iloc[k].to_dict())
            search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1], noise=noise)
            search1.plot(second_data= tdi_fs_added, found_sources_in=found_sources_in[i])
            # print(search1.loglikelihood([found_sources_in[i][j]]))
            # search1.update_noise(pGB=found_sources_in[i][j])
            # print(search1.loglikelihood([found_sources_in[i][j]]))
            # search1.plot(second_data= tdi_fs_added, found_sources_in=found_sources_in[i])

        # print(search_subtracted.loglikelihood([found_sources_in[i][j]]),search_subtracted.loglikelihood_dd([found_sources_in[i][j]]))
        # print(search_subtracted.SNR([found_sources_in[i][j]]))
        # print(search_subtracted.intrinsic_SNR([found_sources_in[i][j]]))
        # j =1
        # print(search_subtracted.loglikelihood([found_sources_in[i][j]]),search_subtracted.loglikelihood_dd([found_sources_in[i][j]]))
        # print(search_subtracted.SNR([found_sources_in[i][j]]))
        # print(search_subtracted.intrinsic_SNR([found_sources_in[i][j]]))
        # search_subtracted.plot(found_sources_in=found_sources_in[i])
        # j=0
        # search_subtracted.update_noise(pGB=found_sources_in[i][j])
        # search_subtracted.plot(found_sources_in=found_sources_in[i])
        # print(search_subtracted.loglikelihood([found_sources_in[i][j]]),search_subtracted.loglikelihood_dd([found_sources_in[i][j]]))
        # print(search_subtracted.SNR([found_sources_in[i][j]]))
        # print(search_subtracted.intrinsic_SNR([found_sources_in[i][j]]))
        # j =1
        # print(search_subtracted.loglikelihood([found_sources_in[i][j]]),search_subtracted.loglikelihood_dd([found_sources_in[i][j]]))
        # print(search_subtracted.SNR([found_sources_in[i][j]]))
        # print(search_subtracted.intrinsic_SNR([found_sources_in[i][j]]))
        # search_subtracted.plot(second_data= tdi_fs_subtracted, found_sources_in=found_sources_in[i])
        # posterior_calculation_input.append((tdi_fs_subtracted, Tobs, frequencies_search[i], found_sources_in[i][j], pGB_injected[i][j]))

        # print('compute posterior of the signal',i,j, found_sources_in[i][j])
        # print(search_subtracted.loglikelihood([found_sources_in[i][j]]))
        # print(search_subtracted.loglikelihood_gpu(pGB_dict_to_gpu_input([found_sources_in[i][j]]),start_freq_ind=0))
        # print(search_subtracted.loglikelihood_gpu_signal(pGB_dict_to_gpu_input([found_sources_in[i][j]])))
        # tdi_fs_subtracted_part['A'] = (tdi_fs_subtracted_part['Z'] - tdi_fs_subtracted_part['X'])/np.sqrt(2.0)

        number_of_signals += 1
        frequencies_found.append(found_sources_in[i][j]['Frequency'])
        compute_posterior(tdi_fs_added, Tobs, frequencies_search[i], found_sources_in[i][j], pGBadded20, noise=noise)#, pGB_injected_matched[i][j].to_dict())
        print('time for one signal posterior', time.time() - start)
print('time to search ', number_of_signals, 'signals: ', time.time()-start_all, (time.time()-start_all)/number_of_signals, 's per signal')
start = time.time()

#{'Amplitude': 2.031616281777476e-23, 'EclipticLatitude': -0.15997150265171872, 'EclipticLongitude': -1.6348821584004698,
# 'Frequency': 0.003979533465546961, 'FrequencyDerivative': 1.0153651382637698e-16, 'Inclination': 2.455284370386729,
#  'InitialPhase': 2.01142405115654, 'Polarization': 0.7772759959676965, 'IntrinsicSNR': 27.674888111907023}
# for i in range(len(posterior_calculation_input)):
#     for j in range(len(posterior_calculation_input)):
#         if int(np.round(posterior_calculation_input[i][3]['Frequency']*10**12)) == int(np.round(posterior_calculation_input[j][3]['Frequency']*10**12)):
#             if i != j:
#                 print(i,j,int(np.round(posterior_calculation_input[i][3]['Frequency']*10**12)))
# pool = mp.Pool(mp.cpu_count())
# number_of_threads = 16 
# pool = mp.Pool(number_of_threads)
# batch = 0
# print(len(posterior_calculation_input[batch*number_of_threads:(batch+1)*number_of_threads]))
# mcmc_samples = pool.starmap(compute_posterior,posterior_calculation_input[batch*number_of_threads:(batch+1)*number_of_threads])
# pool.close()
# pool.join()
print('time to search ', len(posterior_calculation_input), 'signals: ', time.time()-start)
