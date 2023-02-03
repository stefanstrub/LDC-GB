#%%
from re import A
# from matplotlib.lines import _LineStyle
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scipy
from scipy.optimize import differential_evolution
import numpy as np
import xarray as xr
import time
from copy import deepcopy
import multiprocessing as mp
import pandas as pd
import os
import h5py
import sys
sys.path.append('/cluster/home/sstrub/Repositories/LDC/lib/lib64/python3.8/site-packages/ldc-0.1-py3.8-linux-x86_64.egg')

from astropy import units as u
import astropy.coordinates as coord

from ldc.lisa.noise import get_noise_model
from ldc.common.series import TimeSeries, window
import ldc.waveform.fastGB as fastGB
from ldc.common.tools import compute_tdi_snr


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
        elif parameter in parameters_log_uniform:
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
        elif parameter in parameters_log_uniform:
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
        elif parameter in parameters_log_uniform:
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
        elif parameter in parameters_log_uniform:
            maxpGB[parameter] = (np.log10(previous_max[parameter]) - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0])
        else:
            maxpGB[parameter] = (previous_max[parameter] - boundaries[parameter][0]) / (boundaries[parameter][1] - boundaries[parameter][0])
    return maxpGB

def reduce_boundaries(maxpGB, boundaries, ratio=0.1):
    boundaries_reduced = deepcopy(boundaries)
    for parameter in parameters:
        length = boundaries[parameter][1] - boundaries[parameter][0]
        if parameter == "EclipticLatitude":
            boundaries_reduced[parameter] = [
                np.sin(maxpGB[parameter]) - length * ratio / 2,
                np.sin(maxpGB[parameter]) + length * ratio / 2,
            ]
        elif parameter == "Inclination":
            boundaries_reduced[parameter] = [
                np.cos(maxpGB[parameter]) - length * ratio / 2,
                np.cos(maxpGB[parameter]) + length * ratio / 2,
            ]
        elif parameter == "Frequency":
            boundaries_reduced[parameter] = [
                maxpGB[parameter] - length * ratio / 2,
                maxpGB[parameter] + length * ratio / 2,
            ]
        elif parameter == "FrequencyDerivative":
            boundaries_reduced[parameter] = [
                maxpGB[parameter] - length * ratio / 2,
                maxpGB[parameter] + length * ratio / 2,
            ]
        elif parameter == "Amplitude":
            boundaries_reduced[parameter] = [
                np.log10(maxpGB[parameter]) - length * ratio / 2,
                np.log10(maxpGB[parameter]) + length * ratio / 2,
            ]
        elif parameter in ["InitialPhase",'Polarization']:
            boundaries_reduced[parameter] = [
                maxpGB[parameter] - length * ratio / 2,
                maxpGB[parameter] + length * ratio / 2,
            ]
        elif parameter == "EclipticLongitude":
            boundaries_reduced[parameter] = [
                maxpGB[parameter] - length * ratio / 2,
                maxpGB[parameter] + length * ratio / 2,
            ]
        if boundaries_reduced[parameter][0] < boundaries[parameter][0]:
            boundaries_reduced[parameter][0] = boundaries[parameter][0]
        if boundaries_reduced[parameter][1] > boundaries[parameter][1]:
            boundaries_reduced[parameter][1] = boundaries[parameter][1]
    return boundaries_reduced

def CoordinateMC(n, pGBs, boundaries, parameters_recorded, loglikelihood, n_trials=50):
    maxpGB = []
    parameters_recorded1 = []
    no_improvement_counter = 0
    for i in range(number_of_signals):
        maxpGB.append(deepcopy(pGBs))
        parameters_recorded1.append({})
        for parameter in parameters_no_amplitude:
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
    for i in range(100):
        # if i > 48:
        #     n_trials = 50
        parameter1 = parameters_no_amplitude[i % 7]
        parameter2 = parameters_no_amplitude[np.random.randint(0, 6)]
        parameter3 = parameters_no_amplitude[np.random.randint(0, 6)]
        # parameter2 = 'InitialPhase'
        # parameter1 = 'Inclination'
        while parameter2 == parameter1:
            parameter2 = parameters_no_amplitude[np.random.randint(0, 6)]
        while parameter3 == parameter1 or parameter3 == parameter2:
            parameter3 = parameters_no_amplitude[np.random.randint(0, 6)]
        # if parameter1 == 'Frequency':
        #     parameter2 = 'Polarization'
        changeableparameters = [parameter1, parameter2]
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
            for parameter in parameters_no_amplitude:
                parameters_recorded1[i][parameter].append(maxpGB[i][parameter])

        maxpGB2 = deepcopy(maxpGB)
    if previous_best > best_value:
        best_value = previous_best
    parameters_recorded[n] = parameters_recorded1
    return parameters_recorded1

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

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

class Search():
    def __init__(self,tdi_fs,Tobs, lower_frequency, upper_frequency, recombination=0.75):
        self.tdi_fs = tdi_fs
        self.GB = fastGB.FastGB(delta_t=dt, T=Tobs)
        self.reduced_frequency_boundaries = None
        self.recombination = recombination
        self.lower_frequency = lower_frequency
        self.upper_frequency = upper_frequency
        self.padding = (upper_frequency - lower_frequency)/2
        tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft()) for k in ["X", "Y", "Z"]]))
        # f, psdX =  scipy.signal.welch(tdi_ts["X"], fs=1.0/dt, window='hanning', nperseg=len(tdi_ts["X"])/1)
        # f, psdY =  scipy.signal.welch(tdi_ts["Y"], fs=1.0/dt, window='hanning', nperseg=len(tdi_ts["X"])/1)
        # f, psdZ =  scipy.signal.welch(tdi_ts["Z"], fs=1.0/dt, window='hanning', nperseg=len(tdi_ts["X"])/1)
        # psd = psdX + psdY + psdZ
        # indexes = np.logical_and(f>lower_frequency-self.padding, f<upper_frequency+self.padding)
        # psd = psd[indexes]
        self.frequency_T_threshold = 19.1*10**-3/2
        
        self.use_T_component = False
        if self.upper_frequency + self.padding > self.frequency_T_threshold:
            self.use_T_component = True
        # plt.figure()
        # plt.plot(f[indexes], psdX[indexes])
        # plt.plot(f[indexes], psd)
        # plt.show()

        frequencyrange =  [lower_frequency-self.padding, upper_frequency+self.padding]
        indexes = np.logical_and(tdi_fs['X'].f > frequencyrange[0], tdi_fs['X'].f < frequencyrange[1]) 
        self.dataX = tdi_fs["X"][indexes]
        self.dataY = tdi_fs["Y"][indexes]
        self.dataZ = tdi_fs["Z"][indexes]

        self.DAf = (self.dataZ - self.dataX)/np.sqrt(2.0)
        self.DEf = (self.dataZ - 2.0*self.dataY + self.dataX)/np.sqrt(6.0)
        self.DTf = (self.dataZ + self.dataY + self.dataX)/np.sqrt(3.0)

        # plt.figure()
        # plt.plot(f[indexes], np.abs(self.DAf))
        # plt.plot(f[indexes], np.abs(self.DEf))
        # plt.plot(f[indexes], np.abs(self.DAf)+np.abs(self.DEf))
        # plt.show()
        # print('frequencyrange',frequencyrange)
        fmin, fmax = float(self.dataX.f[0]), float(self.dataX.f[-1] + self.dataX.attrs["df"])
        freq = np.array(self.dataX.sel(f=slice(fmin, fmax)).f)
        Nmodel = get_noise_model(noise_model, freq)
        self.Sn = Nmodel.psd(freq=freq, option="X")
        self.SA = Nmodel.psd(freq=freq, option="A")
        self.SE = Nmodel.psd(freq=freq, option="E")
        self.ST = Nmodel.psd(freq=freq, option="T")

        f_0 = fmin
        f_transfer = 19.1*10**-3
        snr = 7
        amplitude_lower = 2*snr/(Tobs * np.sin(f_0/ f_transfer)**2/self.SA[0])**0.5
        snr = 1000
        amplitude_upper = 2*snr/(Tobs * np.sin(f_0/ f_transfer)**2/self.SA[0])**0.5
        amplitude = [amplitude_lower, amplitude_upper]
        # print('lower frequency', lower_frequency)
        # print('amplitude boundaries', amplitude)
        # print('amplitude boundaries previous', np.sqrt(np.max(psd))/1000,np.sqrt(np.max(psd))/10)
        # print('amplitude boundaries previous', np.max(np.abs(self.DAf.data))/10**7,np.max(np.abs(self.DAf.data))/10**5)
        # amplitude = np.sqrt(np.max(psd))

        # indexes = np.argsort(p.get('Frequency'))
        # index_low = np.searchsorted(p.get('Frequency')[indexes], frequencyrange[0])
        # index_high = np.searchsorted(p.get('Frequency')[indexes], frequencyrange[1])

        # pGB = deepcopy(pGBadded)
        # self.pGB = {'Amplitude': 3.676495e-22, 'EclipticLatitude': 0.018181, 'EclipticLongitude': 1.268061, 'Frequency': 0.01158392, 'FrequencyDerivative': 8.009579e-15, 'Inclination': 0.686485, 'InitialPhase': 4.201455, 'Polarization': 2.288223}
        
        fd_range = [np.log10(frequency_derivative(lower_frequency,0.1)),np.log10(frequency_derivative(lower_frequency,M_chirp_upper_boundary))]
        # fd_range = [frequency_derivative(lower_frequency,0.1),frequency_derivative(lower_frequency,M_chirp_upper_boundary)]
        # fd_range = [frequency_derivative_tyson_lower(lower_frequency),frequency_derivative_tyson(lower_frequency)]
        fd_range = [frequency_derivative_tyson_lower(lower_frequency),frequency_derivative(upper_frequency,M_chirp_upper_boundary)]
        self.boundaries = {
            "Amplitude": [np.log10(amplitude[0]),np.log10(amplitude[1])],
            # "Amplitude": [-23.5,-21],
            # "Amplitude": [np.log10(self.pGB['Amplitude'])-2,np.log10(self.pGB['Amplitude'])+1],
            "EclipticLatitude": [-1.0, 1.0],
            "EclipticLongitude": [-np.pi, np.pi],
            # "Frequency": [self.pGB["Frequency"] * 0.99995, self.pGB["Frequency"] * 1.00015],
            # "Frequency": [self.pGB["Frequency"] - 3e-7, self.pGB["Frequency"] + 3e-7],
            "Frequency": frequencyrange,
            "FrequencyDerivative": fd_range,
            # "FrequencyDerivative": [np.log10(5e-6*self.pGB['Frequency']**(13/3)),np.log10(8e-8*self.pGB['Frequency']**(11/3))],
            "Inclination": [-1.0, 1.0],
            "InitialPhase": [0.0, 2.0 * np.pi],
            "Polarization": [0.0, 1.0 * np.pi],
        }
        # print('fd_rage',self.boundaries['FrequencyDerivative'], 'at f=', lower_frequency)

        # print('fd * Tobs', 10**self.boundaries['FrequencyDerivative'][1] * Tobs)
        # print('smear f', 300*lower_frequency * 10**3 / 10**9)

        if self.boundaries['FrequencyDerivative'][0] > self.boundaries['FrequencyDerivative'][1]:
            c = self.boundaries['FrequencyDerivative'][0]
            self.boundaries['FrequencyDerivative'][0] = self.boundaries['FrequencyDerivative'][1]
            self.boundaries['FrequencyDerivative'][1] = c
        
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
            elif parameter in parameters_log_uniform:
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

    def f_statistic(self, N_frequency, N_sky):
        # We construct a global proposal density using the single
        # source F statistic to compute the individual likelihoods
        F_stat = []
        frequency = []
        eclipticlatitude = []
        eclipticlongitude = []
        pGBf = {}
        for parameter in parameters:
            pGBf[parameter] = 0
        pGBf['Amplitude'] = 1e-24
        pGBf['FrequencyDerivative'] = 0
        frequency_boundaries = [self.lower_frequency,self.upper_frequency]
        for n in range(N_sky):
            eclipticlatitude.append(self.boundaries['EclipticLatitude'][0]+(self.boundaries['EclipticLatitude'][1]-self.boundaries['EclipticLatitude'][0])*n/N_sky)
            eclipticlongitude.append(self.boundaries['EclipticLongitude'][0]+(self.boundaries['EclipticLongitude'][1]-self.boundaries['EclipticLongitude'][0])*n/N_sky)
        for k in range(N_frequency):
            F_stat.append([])
            frequency.append(frequency_boundaries[0] + (frequency_boundaries[1]-frequency_boundaries[0])*k/(N_frequency-1))
            for l in range(N_sky):
                F_stat[-1].append([])
                for m in range(N_sky):
                    F_stat[-1][-1].append(self.F_fd0(frequency[-1],eclipticlatitude[l],eclipticlongitude[m],pGBf))
        F_stat = np.asarray(F_stat)
        return F_stat, frequency, eclipticlatitude, eclipticlongitude

    def F_fd0(self, f0, theta, phi, pGBs):
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

    def F(self, intrinsic_parameter_values):
        g = []
        pGBs = {}
        pGBs['Amplitude'] = 1e-24
        for parameter in intrinsic_parameters:
            pGBs[parameter] = intrinsic_parameter_values[parameter]
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

    def plot(self, maxpGBs=None, pGBadded=None, second_data = None,  found_sources_in= [], found_sources_not_matched= [], pGB_injected = [], pGB_injected_matched = [], added_label='Injection2', saving_label =None, vertical_lines = []):
        plt.figure(figsize=fig_size)
        fig, [ax1, ax2] = plt.subplots(2, 1, sharex=True, figsize=fig_size)
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
        ax1.plot(self.DAf.f*10**3,self.DAf,'k',zorder= 1, linewidth = 2, label = 'Data')
        ax2.plot(self.DEf.f*10**3,np.abs(self.DEf),'k',zorder= 1, linewidth = 2, label = 'Data')
        # ax1.plot(tdi_fs_long_subtracted.f[range_index],np.abs(tdi_fs_long_subtracted['X'][range_index])**2,'b',zorder= 5)

        if second_data != None:
            a,Xs = xr.align(self.dataX, second_data['X'], join='left',fill_value=0)
            a,Ys = xr.align(self.dataY, second_data['Y'], join='left',fill_value=0)
            a,Zs = xr.align(self.dataZ, second_data['Z'], join='left',fill_value=0)
            Af = (Zs - Xs)/np.sqrt(2.0)
            Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
            ax1.plot(Af.f*10**3,Af,'r--',zorder= 1, linewidth = 2, label = 'Data subtracted')
            ax2.plot(Ef.f*10**3,np.abs(Ef),'r--',zorder= 1, linewidth = 2, label = 'Data subtracted')


        for j in range(len( pGB_injected)):
            Xs, Ys, Zs = GB.get_fd_tdixyz(template= pGB_injected[j], oversample=4, simulator="synthlisa")
            a,Xs = xr.align(self.dataX, Xs, join='left',fill_value=0)
            a,Ys = xr.align(self.dataY, Ys, join='left',fill_value=0)
            a,Zs = xr.align(self.dataZ, Zs, join='left',fill_value=0)
            Af = (Zs - Xs)/np.sqrt(2.0)
            Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
            ax1.plot(Af.f*10**3,Af.data, color='grey', linewidth = 5, alpha = 0.5)
            ax2.plot(Ef.f*10**3,np.abs(Ef.data), color='grey', linewidth = 5, alpha = 0.5)

        for j in range(len(pGB_injected_matched)):
            Xs, Ys, Zs = GB.get_fd_tdixyz(template= pGB_injected_matched[j], oversample=4, simulator="synthlisa")
            a,Xs = xr.align(self.dataX, Xs, join='left',fill_value=0)
            a,Ys = xr.align(self.dataY, Ys, join='left',fill_value=0)
            a,Zs = xr.align(self.dataZ, Zs, join='left',fill_value=0)
            Af = (Zs - Xs)/np.sqrt(2.0)
            Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
            ax1.plot(Af.f*10**3,Af.data, color=colors[j%10], linewidth = 5, alpha = 0.5)
            ax2.plot(Ef.f*10**3,np.abs(Ef.data), color=colors[j%10], linewidth = 5, alpha = 0.5)


        if pGBadded != None:
            Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=pGBadded, oversample=4, simulator="synthlisa")
            index_low = np.searchsorted(Xs.f, self.dataX.f[0])
            Xs = Xs[index_low : index_low + len(self.dataX)]
            Ys = Ys[index_low : index_low + len(self.dataY)]
            Zs = Zs[index_low : index_low + len(self.dataZ)]
            Af = (Zs - Xs)/np.sqrt(2.0)
            Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
            ax1.plot(Af.f* 1000, Af.data, marker='.', label=added_label)
            ax2.plot(Ef.f* 1000, np.abs(Ef.data), marker='.', label=added_label)

        for j in range(len(found_sources_in)):
            Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=found_sources_in[j], oversample=4, simulator="synthlisa")
            index_low = np.searchsorted(Xs.f, self.dataX.f[0])
            Xs = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
            Ys = xr.align(self.dataX, Ys, join='left',fill_value=0)[1]
            Zs = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            Af = (Zs - Xs)/np.sqrt(2.0)
            Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
            ax1.plot(Af.f* 1000, Af.data,'--', color= colors[j%10], linewidth = 1.6)
            ax2.plot(Ef.f* 1000, np.abs(Ef.data),'--', color= colors[j%10], linewidth = 1.6)

        for j in range(len(found_sources_not_matched)):
            Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=found_sources_not_matched[j], oversample=4, simulator="synthlisa")
            index_low = np.searchsorted(Xs.f, self.dataX.f[0])
            Xs = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
            Zs = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            Af = (Zs - Xs)/np.sqrt(2.0)
            Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
            ax1.plot(Af.f* 1000,Af.data,'.', color= colors[j%10], linewidth = 1.6)
            ax2.plot(Ef.f* 1000, np.abs(Ef.data),'.', color= colors[j%10], linewidth = 1.6)

        # ax1.plot(Xs_added2.f * 1000, Xs_added2.values.real, label="VGB2", marker=".", zorder=5)
        ax1.axvline(self.lower_frequency* 1000, color= 'red', label='Boundaries')
        ax1.axvline(self.upper_frequency* 1000, color= 'red')
        ax2.axvline(self.lower_frequency* 1000, color= 'red')
        ax2.axvline(self.upper_frequency* 1000, color= 'red')
        for j in range(len(vertical_lines)):
            ax1.axvline(vertical_lines[j]* 1000, color= 'red')
            ax2.axvline(vertical_lines[j]* 1000, color= 'red')
        # ax2.axvline(self.lower_frequency* 1000- 4*32*10**-6, color= 'green')
        # ax2.axvline(self.upper_frequency* 1000+ 4*32*10**-6, color= 'green')
        # if self.reduced_frequency_boundaries != None:
        #     ax1.axvline(self.reduced_frequency_boundaries[0]* 1000, color= 'green', label='Reduced Boundaries')
        #     ax1.axvline(self.reduced_frequency_boundaries[1]* 1000, color= 'green')

        # ax1.plot(Xs.f * 1000, dataX.values.real - Xs.values.real, label="residual", alpha=0.8, color="red", marker=".")
        plt.xlabel('f (mHz)')
        ax1.set_ylabel('real A')    
        ax2.set_ylabel('|E|') 
        # ax1.set_yscale('log')  
        ax2.set_yscale('log')   
        ax1.set_xlim((self.lower_frequency-self.padding)*10**3, (self.upper_frequency+self.padding)*10**3)
        ax2.set_xlim((self.lower_frequency-self.padding)*10**3, (self.upper_frequency+self.padding)*10**3)
        # ax1.set_xlim((self.lower_frequency)*10**3, (self.upper_frequency)*10**3)
        # ax2.set_xlim((self.lower_frequency)*10**3, (self.upper_frequency)*10**3)
        ax1.xaxis.set_major_locator(plt.MaxNLocator(4))
        ax2.xaxis.set_major_locator(plt.MaxNLocator(4))
        # plt.legend()
        plt.pause(1)
        if saving_label != None:
            plt.savefig(saving_label,dpi=300,bbox_inches='tight')
        plt.pause(1)
        plt.show()
        # print("p true", self.loglikelihood([pGB]), "null hypothesis", self.loglikelihood([null_pGBs]))

    def plotA(self, maxpGBs=None, pGBadded=None, found_sources_in= [], found_sources_not_matched= [], pGB_injected = [], pGB_injected_matched = [], added_label='Injection2', saving_label =None, vertical_lines = [], second_data= None):
        plt.figure(figsize=fig_size)
        fig, [ax1, ax2] = plt.subplots(2, 1, sharex=False, figsize=fig_size)
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
        parameter_x = 'Frequency'
        parameter_y = 'Amplitude'
                    
        # Af = (Zs - Xs)/np.sqrt(2.0)
        ax1.plot(self.DAf.f*10**3,np.abs(self.DAf),'k',zorder= 1, linewidth = 2, label = 'Data')
        # ax1.plot(tdi_fs_long_subtracted.f[range_index],np.abs(tdi_fs_long_subtracted['X'][range_index])**2,'b',zorder= 5)

        if second_data != None:
            a,Xs = xr.align(self.dataX, second_data['X'], join='left',fill_value=0)
            a,Ys = xr.align(self.dataY, second_data['Y'], join='left',fill_value=0)
            a,Zs = xr.align(self.dataZ, second_data['Z'], join='left',fill_value=0)
            Af = (Zs - Xs)/np.sqrt(2.0)
            Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
            ax1.plot(Af.f*10**3,Af,'c--',zorder= 1, linewidth = 2, label = 'Data subtracted')
            ax2.plot(Ef.f*10**3,np.abs(Ef),'c--',zorder= 1, linewidth = 2, label = 'Data subtracted')

        for j in range(len( pGB_injected)):
            Xs, Ys, Zs = GB.get_fd_tdixyz(template= pGB_injected[j], oversample=4, simulator="synthlisa")
            a,Xs = xr.align(self.dataX, Xs, join='left',fill_value=0)
            a,Ys = xr.align(self.dataY, Ys, join='left',fill_value=0)
            a,Zs = xr.align(self.dataZ, Zs, join='left',fill_value=0)
            Af = (Zs - Xs)/np.sqrt(2.0)
            Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
            ax1.plot(Af.f*10**3,np.abs(Af.data), color=colors[j%10], linewidth = 5, alpha = 0.5)
            ax2.plot(pGB_injected[j][parameter_x]*10**3,pGB_injected[j][parameter_y],'o', color=colors[j%10], markersize=7, alpha = 0.5)
            # ax3.plot(pGB_injected[j]['EclipticLongitude']*10**3,pGB_injected[j]['EclipticLatitude'],'o', color=colors[j%10], markersize=7, alpha = 0.5)
            # ax4.plot(pGB_injected[j]['Inclination']*10**3,pGB_injected[j]['FrequencyDerivative'],'o', color=colors[j%10], markersize=7, alpha = 0.5)
            # ax5.plot(pGB_injected[j]['InitialPhase']*10**3,pGB_injected[j]['Polarization'],'o', color=colors[j%10], markersize=7, alpha = 0.5)

        for j in range(len(pGB_injected_matched)):
            Xs, Ys, Zs = GB.get_fd_tdixyz(template= pGB_injected_matched[j], oversample=4, simulator="synthlisa")
            index_low = np.searchsorted(Xs.f, self.dataX.f[0])
            try:
                if np.abs(self.dataX.f[0] - Xs.f[index_low-1]) < np.abs(self.dataX.f[0] - Xs.f[index_low]):
                    index_low = index_low-1
            except:
                pass
            Xs = Xs[index_low : index_low + len(self.dataX)]
            Ys = Ys[index_low : index_low + len(self.dataY)]
            Zs = Zs[index_low : index_low + len(self.dataZ)]
            Af = (Zs - Xs)/np.sqrt(2.0)
            Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
            ax1.plot(Af.f*10**3,np.abs(Af.data), color=colors[j%10], linewidth = 5, alpha = 0.5)
            ax2.plot(pGB_injected_matched[j][parameter_x]*10**3,pGB_injected_matched[j][parameter_y],'o', color=colors[j%10], markersize=7, alpha = 0.5)
            # ax3.plot(pGB_injected_matched[j]['EclipticLongitude']*10**3,pGB_injected_matched[j]['EclipticLatitude'],'o', color=colors[j%10], markersize=7, alpha = 0.5)


        if pGBadded != None:
            Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=pGBadded, oversample=4, simulator="synthlisa")
            index_low = np.searchsorted(Xs.f, self.dataX.f[0])
            Xs = Xs[index_low : index_low + len(self.dataX)]
            Ys = Ys[index_low : index_low + len(self.dataY)]
            Zs = Zs[index_low : index_low + len(self.dataZ)]
            Af = (Zs - Xs)/np.sqrt(2.0)
            Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
            ax1.plot(Af.f* 1000, np.abs(Af.data), marker='.', label=added_label)
            # ax2.plot(Ef.f* 1000, np.abs(Ef.data), marker='.', label=added_label)

        # for j in range(len(found_sources_in)):
        #     Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=found_sources_in[j], oversample=4, simulator="synthlisa")
        #     index_low = np.searchsorted(Xs.f, self.dataX.f[0])
        #     Xs = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
        #     Zs = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
        #     Af = (Zs - Xs)/np.sqrt(2.0)
        #     Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
        #     ax1.plot(Af.f* 1000, np.abs(Af.data),'--', color= colors[j%10], linewidth = 1.6)
        #     ax2.plot(found_sources_in[j][parameter_x]*10**3,found_sources_in[j][parameter_y],'.', color=colors[j%10], markersize=7)
        #     # ax2.plot(Ef.f* 1000, np.abs(Ef.data),'--', color= colors[j%10], linewidth = 1.6)

        # for j in range(len(found_sources_in)):
        #     Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=found_sources_in[j], oversample=4, simulator="synthlisa")
        #     index_low = np.searchsorted(Xs.f, self.dataX.f[0])
        #     if j == 0:
        #         Xs = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
        #         Ys = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
        #         Zs = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
        #     # ax2.plot(found_sources_in[j][parameter_x]*10**3,found_sources_in[j][parameter_y],'o', color='grey', markersize=7)
        #     Af = (Zs - Xs)/np.sqrt(2.0)
        #     Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
        #     ax1.plot(Af.f* 1000, np.abs(Af.data),linewidth = 5, alpha = 0.5, color= 'grey')
        #     # ax2.plot(Ef.f* 1000, np.abs(Ef.data),'--', color= colors[j%10], linewidth = 1.6)

        for j in range(len(found_sources_in)):
            Xs, Ys, Zs = GB.get_fd_tdixyz(template=found_sources_in[j], oversample=4, simulator="synthlisa")
            index_low = np.searchsorted(Xs.f, self.dataX.f[0])
            try:
                if np.abs(self.dataX.f[0] - Xs.f[index_low-1]) < np.abs(self.dataX.f[0] - Xs.f[index_low]):
                    index_low = index_low-1
            except:
                pass
            Xs = Xs[index_low : index_low + len(self.dataX)]
            Ys = Ys[index_low : index_low + len(self.dataY)]
            Zs = Zs[index_low : index_low + len(self.dataZ)]
            # Xs = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
            # Ys = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
            # Zs = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            Af = (Zs - Xs)/np.sqrt(2.0)
            Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
            ax1.plot(Af.f* 1000, np.abs(Af.data),'--', color= colors[j%10], linewidth = 1.6)
            ax2.plot(found_sources_in[j][parameter_x]*10**3,found_sources_in[j][parameter_y],'+', color=colors[j%10], markersize=7)

        for j in range(len(found_sources_not_matched)):
            Xs, Ys, Zs = GB.get_fd_tdixyz(template=found_sources_not_matched[j], oversample=4, simulator="synthlisa")
            index_low = np.searchsorted(Xs.f, self.dataX.f[0])
            Xs = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
            Ys = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
            Zs = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            Af = (Zs - Xs)/np.sqrt(2.0)
            Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
            ax1.plot(Af.f* 1000, np.abs(Af.data),'--', color= colors[j%10], linewidth = 1.6)
            ax2.plot(found_sources_not_matched[j][parameter_x]*10**3,found_sources_not_matched[j][parameter_y],'+', color=colors[j%10], markersize=7)
            # ax3.plot(found_sources_not_matched[j]['EclipticLongitude']*10**3,found_sources_not_matched[j]['EclipticLatitude'],'+', color=colors[j%10], markersize=7)
            # ax2.plot(Ef.f* 1000, np.abs(Ef.data),'.', color= colors[j%10], linewidth = 1.6)

        # ax1.plot(Xs_added2.f * 1000, Xs_added2.values.real, label="VGB2", marker=".", zorder=5)

        # ax1.axvline(self.lower_frequency* 1000, color= 'red', label='Boundaries')
        # ax1.axvline(self.upper_frequency* 1000, color= 'red')
        # ax2.axvline(self.lower_frequency* 1000, color= 'red')
        # ax2.axvline(self.upper_frequency* 1000, color= 'red')
        # for j in range(len(vertical_lines)):
        #     ax1.axvline(vertical_lines[j]* 1000, color= 'red')
        #     ax2.axvline(vertical_lines[j]* 1000, color= 'red')

        # ax2.axvline(self.lower_frequency* 1000- 4*32*10**-6, color= 'green')
        # ax2.axvline(self.upper_frequency* 1000+ 4*32*10**-6, color= 'green')
        # if self.reduced_frequency_boundaries != None:
        #     ax1.axvline(self.reduced_frequency_boundaries[0]* 1000, color= 'green', label='Reduced Boundaries')
        #     ax1.axvline(self.reduced_frequency_boundaries[1]* 1000, color= 'green')

        # ax1.plot(Xs.f * 1000, dataX.values.real - Xs.values.real, label="residual", alpha=0.8, color="red", marker=".")
        plt.xlabel(parameter_x)
        if parameter_x == 'Frequency':
            plt.xlabel(parameter_x+' (mHz)')
        ax1.set_ylabel('|A|')    
        ax2.set_ylabel(parameter_y) 
        ax1.set_yscale('log')  
        ax2.set_yscale('log')   
        ax1.set_xlim((self.lower_frequency-self.padding)*10**3, (self.upper_frequency+self.padding)*10**3)
        ax2.set_xlim((self.lower_frequency-self.padding)*10**3, (self.upper_frequency+self.padding)*10**3)
        # ax2.set_ylim(10**-23,4*10**-23)
        # ax1.set_xlim((self.lower_frequency)*10**3, (self.upper_frequency)*10**3)
        # ax2.set_xlim((self.lower_frequency)*10**3, (self.upper_frequency)*10**3)
        ax1.xaxis.set_major_locator(plt.MaxNLocator(4))
        # ax2.xaxis.set_major_locator(plt.MaxNLocator(4))
        # plt.legend()
        plt.tight_layout()
        plt.pause(1)
        if saving_label != None:
            plt.savefig(saving_label,dpi=300,bbox_inches='tight')
        plt.pause(1)
        plt.show()
        # print("p true", self.loglikelihood([pGB]), "null hypothesis", self.loglikelihood([null_pGBs]))

    def SNR_split(self, pGBs):
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
        if self.upper_frequency > self.frequency_T_threshold:
            Tf = (Zs_total + Ys_total + Xs_total)/np.sqrt(3.0)
            hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2)/self.SA + np.absolute(Tf.data)**2 /self.ST)
            SNR2 = np.sum( np.real(self.DAf * np.conjugate(Af.data) + self.DEf * np.conjugate(Ef.data))/self.SA + np.real(self.DTf * np.conjugate(Tf.data))/self.ST )
        else:
            hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /self.SA)
            SNR2 = np.sum( np.real(self.DAf * np.conjugate(Af.data) + self.DEf * np.conjugate(Ef.data))/self.SA )
        SNR = 4.0*Xs.df* hh
        SNR2 = 4.0*Xs.df* SNR2
        SNR3 = SNR2 / np.sqrt(SNR)
        return SNR3.values

    def SNR(self, pGBs):
        for i in range(len(pGBs)):
            Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs[i], oversample=4, simulator="synthlisa")
            index_low = np.searchsorted(Xs.f, self.dataX.f[0])
            try:
                if np.abs(self.dataX.f[0] - Xs.f[index_low-1]) < np.abs(self.dataX.f[0] - Xs.f[index_low]):
                    index_low = index_low-1
            except:
                pass
            Xs_total = Xs[index_low : index_low + len(self.dataX)]
            Ys_total = Ys[index_low : index_low + len(self.dataY)]
            Zs_total = Zs[index_low : index_low + len(self.dataZ)]

            # if i == 0:
            #     Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
            #     Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
            #     Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            # else:
            #     Xs_total += xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
            #     Ys_total += xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
            #     Zs_total += xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            
        Af = (Zs_total - Xs_total)/np.sqrt(2.0)
        Ef = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)
        if self.use_T_component:
            Tf = (Zs_total + Ys_total + Xs_total)/np.sqrt(3.0)
            hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2)/self.SA + np.absolute(Tf.data)**2 /self.ST)
            SNR2 = np.sum( np.real(self.DAf * np.conjugate(Af.data) + self.DEf * np.conjugate(Ef.data))/self.SA + np.real(self.DTf * np.conjugate(Tf.data))/self.ST )
        else:
            SNR2 = np.sum( np.real(self.DAf * np.conjugate(Af.data) + self.DEf * np.conjugate(Ef.data))/self.SA )
            hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /self.SA)
        SNR = 4.0*Xs.df* hh
        SNR2 = 4.0*Xs.df* SNR2
        SNR3 = SNR2 / np.sqrt(SNR)
        # plotIt = False
        # if plotIt:
        #     fig, ax = plt.subplots(nrows=3, sharex=True) 
        #     ax[0].plot(Af.f, np.abs(self.DAf))
        #     ax[0].plot(Af.f, np.abs(Af.data))
            
        #     ax[1].plot(Af.f, np.abs(self.DEf))
        #     ax[1].plot(Af.f, np.abs(Ef.data))
        #     ax[2].plot(Af.f, np.abs(self.DTf))
        #     ax[2].plot(Af.f, np.abs(Tf.data))
        #     plt.show()
        return SNR3.values

    def SNR_AET_compute(self, pGBs):
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
        Tf = (Zs_total + Ys_total + Xs_total)/np.sqrt(3.0)

        tdi = dict({"A":Af, "E":Ef, "T":Tf, "X":Xs})
        tdi_data = dict({"A":self.DAf, "E":self.DEf, "T":self.DTf})
        hh = compute_tdi_snr(tdi, Nmodel, AET=True, fmin=Af.f[0], fmax=Af.f[-1])["tot2"]
        SNR2 = compute_tdi_snr(tdi, Nmodel, data= tdi_data, AET=True, fmin=Af.f[0], fmax=Af.f[-1])["tot2"]
        SNR3 = SNR2 / np.sqrt(hh)
        return SNR3

    def SNR_XYZ(self, pGBs):
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
            
        hh = np.sum((np.absolute(Xs_total.data)**2 + np.absolute(Ys_total.data)**2 + np.absolute(Zs_total.data)**2) /self.Sn)
        SNR2 = np.sum( np.real(self.dataX * np.conjugate(Xs_total.data) + self.dataY * np.conjugate(Ys_total.data)+ self.dataZ * np.conjugate(Zs_total.data))/self.Sn)
        SNR = 4.0*Xs.df* hh
        SNR2 = 4.0*Xs.df* SNR2
        SNR3 = SNR2 / np.sqrt(SNR)
        return SNR3.values

    def SNR_XYZ_Sa(self, pGBs):
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
            
        hh = np.sum((np.absolute(Xs_total.data)**2 + np.absolute(Ys_total.data)**2 + np.absolute(Zs_total.data)**2) /self.SA)
        SNR2 = np.sum( np.real(self.dataX * np.conjugate(Xs_total.data) + self.dataY * np.conjugate(Ys_total.data)+ self.dataZ * np.conjugate(Zs_total.data))/self.SA)
        SNR = 4.0*Xs.df* hh
        SNR2 = 4.0*Xs.df* SNR2
        SNR3 = SNR2 / np.sqrt(SNR)
        return SNR3.values

    def SNR_noise_matrix(self, pGBs):
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
        noise_model = "SciRDv1"
        Nmodel = get_noise_model(noise_model, np.logspace(-5, -1, 100))

        tdi = dict({"X":Xs_total, "Y":Ys_total, "Z":Zs_total})
        tdi_data = dict({"X":self.dataX, "Y":self.dataY, "Z":self.dataZ})
        hh = compute_tdi_snr(tdi, Nmodel)["tot2"]
        SNR2 = compute_tdi_snr(tdi, Nmodel, data= tdi_data)["tot2"]
        SNR3 = SNR2 / np.sqrt(hh)
        return SNR3

    def SNR_AE(self, pGBs):
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
        SNR = 4.0*Xs.df* hh
        SNR2 = 4.0*Xs.df* SNR2
        SNR3 = SNR2 / np.sqrt(SNR)
        return SNR3.values

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

    def loglikelihood(self, pGBs):
        # start = time.time()
        for i in range(len(pGBs)):
            Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs[i], oversample=4, simulator="synthlisa")
            # print('fastGB', time.time()-start)
            # start = time.time()
            if i == 0:
                Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
                Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
                Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            else:
                Xs_total += xr.align(self.dataX, Xs, join='left',fill_value=0)[1]
                Ys_total += xr.align(self.dataY, Ys, join='left',fill_value=0)[1]
                Zs_total += xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
        # print('align', time.time()-start)
        # start = time.time()
        Af = (Zs_total - Xs_total)/np.sqrt(2.0)
        Ef = (Zs_total - 2.0*Ys_total + Xs_total)/np.sqrt(6.0)
        Tf = (Zs_total + Ys_total + Xs_total)/np.sqrt(3.0)
        # print('AET', time.time()-start)
        # start = time.time()
        if self.use_T_component:
            Tf = (Zs_total + Ys_total + Xs_total)/np.sqrt(3.0)
            hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2)/self.SA + np.absolute(Tf.data)**2 /self.ST)
            SNR2 = np.sum( np.real(self.DAf * np.conjugate(Af.data) + self.DEf * np.conjugate(Ef.data))/self.SA + np.real(self.DTf * np.conjugate(Tf.data))/self.ST )
        else:
            SNR2 = np.sum( np.real(self.DAf * np.conjugate(Af.data) + self.DEf * np.conjugate(Ef.data))/self.SA )
            hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /self.SA)
        
        # print('SNR', time.time()-start)
        # start = time.time()
        logliks = 4.0*Xs.df*( SNR2 - 0.5 * hh )
        # print('log', time.time()-start)
        return logliks.values

    def loglikelihood_noise_matrix(self, pGBs):
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
            

        tdi_signal = dict({"X":Xs_total, "Y":Ys_total, "Z":Zs_total})
        tdi_data = dict({"X":self.dataX, "Y":self.dataY, "Z":self.dataZ})
        hh = compute_tdi_snr(tdi_signal, Nmodel)["tot2"]
        SNR2 = compute_tdi_snr(tdi_signal, Nmodel, data= tdi_data)["tot2"]
        logliks = SNR2 - 0.5 * hh
        return logliks

    def intrinsic_SNR(self, pGBs):
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
        if self.use_T_component:
            Tf = (Zs_total + Ys_total + Xs_total)/np.sqrt(3.0)
            hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2)/self.SA + np.absolute(Tf.data)**2 /self.ST)
        else:
            hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /self.SA)
        SNR = 4.0*Xs.df* hh
        return np.sqrt(SNR)

    def intrinsic_SNR_T(self, pGBs):
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
        Tf = (Zs_total + Ys_total + Xs_total)/np.sqrt(3.0)
        hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /np.absolute(Tf.data)**2)
        SNR = 4.0*Xs.df* hh
        return np.sqrt(SNR)

    def intrinsic_SNR_old(self, pGBs):
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
        hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /self.SA)
        SNR = 4.0*Xs.df* hh
        return np.sqrt(SNR)

    def SNR_with_rolling_mean(self, pGBs):
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
        residual_A = self.DAf - Af.data
        residual_E = self.DEf - Ef.data
        average_length = 7
        noise_rolling_mean_A = moving_average(np.abs(residual_A.values)**2, n=average_length)

        hh = np.sum((np.absolute(Af.data[int((average_length-1)/2):len(Af.data)-int((average_length-1)/2)])**2 + np.absolute(Ef.data[int((average_length-1)/2):len(Af.data)-int((average_length-1)/2)])**2) /noise_rolling_mean_A)
        SNR = 4.0*Xs.df* hh
        return np.sqrt(SNR)

    def differential_evolution_search(self, frequency_boundaries, initial_guess = None):
        bounds = []
        for signal in range(number_of_signals):
            for i in range(7):
                bounds.append((0,1))

        maxpGB = []
        self.boundaries_reduced = deepcopy(self.boundaries)
        self.boundaries_reduced['Frequency'] = frequency_boundaries
        if initial_guess != None:
            initial_guess01 = np.zeros((len(parameters)-1)*number_of_signals)
            for signal in range(number_of_signals):
                pGBstart01 = scaleto01(initial_guess[signal], self.boundaries_reduced)

                for count, parameter in enumerate(parameters_no_amplitude):
                    if pGBstart01[parameter] < 0:
                        pGBstart01[parameter] = 0
                    if pGBstart01[parameter] > 1:
                        pGBstart01[parameter] = 1
                    initial_guess01[count+(len(parameters_no_amplitude))*signal] = pGBstart01[parameter]
            start = time.time()
            res = differential_evolution(self.function_evolution, bounds=bounds, disp=False, strategy='best1exp', popsize=8,tol= 1e-8 , maxiter=1000, recombination= self.recombination, mutation=(0.5,1), x0=initial_guess01)
            print('time',time.time()-start)
        else:
            start = time.time()
            res = differential_evolution(self.function_evolution, bounds=bounds, disp=False, strategy='best1exp', popsize=8, tol= 1e-8 , maxiter=1000, recombination= self.recombination, mutation=(0.5,1))
            print('time',time.time()-start)
        for signal in range(number_of_signals):
            pGB01 = [0.5] + res.x[signal*7:signal*7+7].tolist()
            maxpGB.append(scaletooriginal(pGB01,self.boundaries_reduced))
        print(res)
        print(maxpGB)
        print('log-likelihood',self.loglikelihood(maxpGB))
        # print(pGB)
        return [maxpGB], res.nfev

    def differential_evolution_search_F(self, frequency_boundaries):
        bounds = []
        for signal in range(number_of_signals):
            for i in range(4):
                bounds.append((0,1))

        maxpGB = []
        self.boundaries_reduced = deepcopy(self.boundaries)
        self.boundaries_reduced['Frequency'] = frequency_boundaries
        start = time.time()
        res, energies = differential_evolution(self.function_evolution_F, bounds=bounds, disp=True, strategy='best1exp', popsize=5,tol= 1e-6 , maxiter=300, recombination= self.recombination, mutation=(0.5,1))
        print('time',time.time()-start)
        for signal in range(number_of_signals):
            pGB01 = [0.5] + res.x[signal*4:signal*4+4].tolist() + [0.5,0.5,0.5] 
            maxpGB.append(scaletooriginal(pGB01,self.boundaries_reduced))
        print(res)
        print(maxpGB)
        print(self.loglikelihood(maxpGB))
        # print(pGB)
        return [maxpGB], energies

    def searchCD(self):
        # np.random.seed(42)

        parameters_recorded = [None] * 10

        n_trials = 50
        number_of_evaluations = 0
        for n in range(len(parameters_recorded)):
            start = time.time()
            parameters_recorded[n] = CoordinateMC(n, self.pGBs, self.boundaries, parameters_recorded, self.SNR, n_trials=n_trials)

            print('n',n+1,'time', int(time.time()-start), np.round(parameters_recorded[n][0]['Loglikelihood'][-1],2),len(parameters_recorded[n][0]['Loglikelihood']))
            number_of_evaluations += len(parameters_recorded[n][0]['Loglikelihood']) * n_trials
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
        good_runs = loglikelihoodofruns > 0
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
                for parameter in parameters_no_amplitude + ['Loglikelihood']:
                    if parameter == 'Loglikelihood':
                        pGBmodes[-1][i][parameter] = parameters_recorded[n][0][parameter][-1]
                    else:
                        pGBmodes[-1][i][parameter] = parameters_recorded[n][i][parameter][-1]
        if len(pGBmodes) > 5:
            for signal in range(number_of_signals):
                pGBmodes = pGBmodes[:5]
        pGBmodes_opt = self.optimize_without_amplitude(pGBmodes)
        return [pGBmodes_opt], number_of_evaluations

    def optimize(self, pGBmodes, boundaries = None):
        if boundaries == None:
            boundaries = self.boundaries
        bounds = ()
        number_of_signals_optimize = len(pGBmodes[0])
        for signal in range(number_of_signals_optimize):
            bounds += ((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1))
        for i in range(len(pGBmodes)):
            maxpGB = []
            boundaries_reduced = []
            pGBs01 = []

            # print('initial loglikelihood noise matrix', self.loglikelihood_noise_matrix(pGBmodes[0]),pGBmodes[0][0]['Frequency'])
            # print('initial loglikelihood', self.loglikelihood(pGBmodes[0]),pGBmodes[0][0]['Frequency'])
            for j in range(2):
                x = []
                for signal in range(number_of_signals_optimize):
                    if j == 0:
                        maxpGB.append({})
                        boundaries_reduced.append({})
                        for parameter in parameters:
                            maxpGB[signal][parameter] = pGBmodes[i][signal][parameter]
                    # print(maxpGB)
                    # boundaries_reduced[signal] = reduce_boundaries(maxpGB[signal], boundaries,ratio=0.1)
                    if j > 0:
                        boundaries_reduced[signal] = reduce_boundaries(maxpGB[signal], boundaries,ratio=0.4)
                    if j in [0]:
                        boundaries_reduced[signal] = deepcopy(boundaries)
                    pGBs01.append({})
                    for parameter in parameters:
                        if parameter in ["EclipticLatitude"]:
                            pGBs01[signal][parameter] = (np.sin(maxpGB[signal][parameter]) - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
                        elif parameter in ["Inclination"]:
                            pGBs01[signal][parameter] = (np.cos(maxpGB[signal][parameter]) - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
                        elif parameter in parameters_log_uniform:
                            pGBs01[signal][parameter] = (np.log10(maxpGB[signal][parameter]) - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
                        else:
                            pGBs01[signal][parameter] = (maxpGB[signal][parameter] - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
                    for parameter in parameters:
                        x.append(pGBs01[signal][parameter])
                # print(loglikelihood(maxpGB))
                res = scipy.optimize.minimize(self.function, x, args=boundaries_reduced, method='SLSQP', bounds=bounds, tol=1e-5)
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
        # print('final optimized loglikelihood noise matrix', self.loglikelihood_noise_matrix(maxpGB),maxpGB[0]['Frequency'])
        # print('final optimized loglikelihood', self.loglikelihood(maxpGB),maxpGB[0]['Frequency'])
        return maxpGB

    def optimize_without_amplitude(self, pGBmodes, boundaries = None):
        if boundaries == None:
            boundaries = self.boundaries
        bounds = ()
        number_of_signals_optimize = len(pGBmodes[0])
        for signal in range(number_of_signals_optimize):
            bounds += ((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1))
        for i in range(len(pGBmodes)):
            maxpGB = []
            pGBs01 = []

            for j in range(2):
                x = []
                for signal in range(number_of_signals_optimize):
                    if j == 0:
                        maxpGB.append({})
                        self.boundaries_reduced = {}
                        for parameter in parameters_no_amplitude:
                            maxpGB[signal][parameter] = pGBmodes[i][signal][parameter]
                    # print(maxpGB)
                    self.boundaries_reduced = reduce_boundaries(maxpGB[signal], boundaries,ratio=0.1)
                    if j == 2:
                        self.boundaries_reduced = reduce_boundaries(maxpGB[signal], boundaries,ratio=0.1)
                    if j in [0,1]:
                        self.boundaries_reduced = deepcopy(boundaries)
                    pGBs01.append({})
                    for parameter in parameters_no_amplitude:
                        if parameter in ["EclipticLatitude"]:
                            pGBs01[signal][parameter] = (np.sin(maxpGB[signal][parameter]) - self.boundaries_reduced[parameter][0]) / (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])
                        elif parameter in ["Inclination"]:
                            pGBs01[signal][parameter] = (np.cos(maxpGB[signal][parameter]) - self.boundaries_reduced[parameter][0]) / (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])
                        elif parameter in parameters_log_uniform:
                            pGBs01[signal][parameter] = (np.log10(maxpGB[signal][parameter]) - self.boundaries_reduced[parameter][0]) / (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])
                        else:
                            pGBs01[signal][parameter] = (maxpGB[signal][parameter] - self.boundaries_reduced[parameter][0]) / (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])
                    for parameter in parameters_no_amplitude:
                        x.append(pGBs01[signal][parameter])
                # print(loglikelihood(maxpGB))
                res = scipy.optimize.minimize(self.function_evolution, x, method='SLSQP', bounds=bounds, tol=1e-10)
                # res = scipy.optimize.minimize(self.function, x, args=self.boundaries_reduced, method='Nelder-Mead', tol=1e-10)
                # res = scipy.optimize.least_squares(self.function, x, args=self.boundaries_reduced, bounds=bounds)
                for signal in range(number_of_signals_optimize):
                    pGB01 = [0.5] + res.x[signal*7:signal*7+7].tolist()
                    maxpGB[signal] = scaletooriginal(pGB01,self.boundaries_reduced)
                    # maxpGB[signal] = scaletooriginal(res.x[signal*7:signal*7+7],self.boundaries_reduced)
                # print('optimized loglikelihood', loglikelihood(maxpGB),maxpGB)
                # print('boundaries reduced', self.boundaries_reduced)

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

    def calculate_Amplitude(self, pGBs):
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
        logliks = 4.0*Xs.df*( SNR2 - 0.5 * hh )
        scalar_product_hh = 4.0*Xs.df* hh
        scalar_product_dh = 4.0*Xs.df* SNR2
        A = scalar_product_dh / scalar_product_hh
        return A


    def optimizeA(self, pGBmodes, boundaries = None):
        if boundaries == None:
            boundaries = self.boundaries
        bounds = ()
        number_of_signals_optimize = len(pGBmodes[0])
        for signal in range(number_of_signals_optimize):
            bounds += ((0,1))
        for i in range(len(pGBmodes)):
            maxpGB = []
            boundaries_reduced = []
            pGBs01 = []

            for j in range(2):
                x = []
                pGBx = []
                for signal in range(number_of_signals_optimize):
                    if j == 0:
                        maxpGB.append({})
                        boundaries_reduced.append({})
                        for parameter in parameters:
                            maxpGB[signal][parameter] = pGBmodes[i][signal][parameter]
                    # print(maxpGB)
                    boundaries_reduced[signal] = reduce_boundaries(maxpGB[signal], boundaries,ratio=0.1)
                    if j == 2:
                        boundaries_reduced[signal] = reduce_boundaries(maxpGB[signal], boundaries,ratio=0.1)
                    if j in [0,1]:
                        boundaries_reduced[signal] = deepcopy(boundaries)
                    pGBs01.append({})
                    for parameter in parameters:
                        if parameter in ["EclipticLatitude"]:
                            pGBs01[signal][parameter] = (np.sin(maxpGB[signal][parameter]) - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
                        elif parameter in ["Inclination"]:
                            pGBs01[signal][parameter] = (np.cos(maxpGB[signal][parameter]) - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
                        elif parameter in parameters_log_uniform:
                            pGBs01[signal][parameter] = (np.log10(maxpGB[signal][parameter]) - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
                        else:
                            pGBs01[signal][parameter] = (maxpGB[signal][parameter] - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
                    for parameter in ['Amplitude']:
                        x.append(pGBs01[signal][parameter])
                    for parameter in parameters:
                        pGBx.append(pGBs01[signal][parameter])
                self.pGBx = pGBx
                res = scipy.optimize.minimize(self.functiona, x, args=(pGBx, boundaries_reduced), method='trust-constr', bounds=[bounds], tol=1e-1)
                # res = scipy.optimize.minimize(self.function, x, args=boundaries_reduced, method='Nelder-Mead', tol=1e-10)
                # res = scipy.optimize.least_squares(self.function, x, args=boundaries_reduced, bounds=bounds)
                for signal in range(number_of_signals_optimize):
                    pGB01 = deepcopy(pGBx)
                    pGB01[signal*8] = res.x[signal]
                    maxpGB[signal] = scaletooriginal(pGB01,boundaries_reduced[signal])
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

    def functionA(self, a):
        pGBs01 = self.pGBx
        boundaries_reduced = deepcopy(self.boundaries)
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
                elif parameter in ["FrequencyDerivative"]:
                    pGBs[signal][parameter] = 10**((pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
                elif parameter in ['Amplitude']:
                    pGBs[signal][parameter] = 10**((a[0] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
                else:
                    pGBs[signal][parameter] = (pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0]
                i += 1
        
        print(pGBs)
        p = -self.loglikelihood(pGBs)
        return p#/10**4

    def functiona(self,a, pGBs01, boundaries_reduced):
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
                elif parameter in ["FrequencyDerivative"]:
                    pGBs[signal][parameter] = 10**((pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
                elif parameter in ['Amplitude']:
                    pGBs[signal][parameter] = 10**((a[0] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
                else:
                    pGBs[signal][parameter] = (pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0]
                i += 1
        p = -self.loglikelihood(pGBs)
        return p#/10**4

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
                elif parameter in parameters_log_uniform:
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
                    pGBs[signal][parameter] = np.arcsin((pGBs01[signal*7:signal*7+7][i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
                elif parameter in ["Inclination"]:
                    pGBs[signal][parameter] = np.arccos((pGBs01[signal*7:signal*7+7][i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
                elif parameter in ['Amplitude']:
                    i -= 1
                    pGBs[signal][parameter] = 10**((0.1 * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
                # elif parameter in ["FrequencyDerivative"]:
                #     pGBs[signal][parameter] = 10**((pGBs01[signal*7:signal*7+7][i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
                else:
                    pGBs[signal][parameter] = (pGBs01[signal*7:signal*7+7][i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0]
                i += 1
        p = -self.SNR(pGBs)
        return p

    def function_evolution_F(self, pGBs01):
        pGBs =  {}
        i = 0
        for parameter in intrinsic_parameters:
            if parameter in ["EclipticLatitude"]:
                pGBs[parameter] = np.arcsin((pGBs01[i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
            elif parameter in ["Inclination"]:
                pGBs[parameter] = np.arccos((pGBs01[i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
            elif parameter in ['Amplitude']:
                i -= 1
                pGBs[parameter] = 10**((0.1 * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
            elif parameter in ["FrequencyDerivative"]:
                pGBs[parameter] = 10**((pGBs01[i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
            else:
                pGBs[parameter] = (pGBs01[i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0]
            i += 1
        intrinsic_parameter_values = {}
        for parameter in intrinsic_parameters:
            intrinsic_parameter_values[parameter] = pGBs[parameter]
        p = -self.F(intrinsic_parameter_values)
        return p

    def function_evolution_8(self, pGBs01):
        pGBs = []
        for signal in range(number_of_signals):
            pGBs.append({})
            i = 0
            for parameter in parameters:
                if parameter in ["EclipticLatitude"]:
                    pGBs[signal][parameter] = np.arcsin((pGBs01[signal*8:signal*8+8][i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
                elif parameter in ["Inclination"]:
                    pGBs[signal][parameter] = np.arccos((pGBs01[signal*8:signal*8+8][i] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
                elif parameter in parameters_log_uniform:
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
        for i in range(1):
            for parameter in parameters:
                if i == 0:
                    step_size[parameter] = 1e-9
                    # if parameter == 'Frequency':
                    #     step_size[parameter] = 0.00001
                else:
                    step_size[parameter] = 0.001/np.sqrt(inner_product[parameter][parameter])
                # if step_size[parameter] > 1e-9:
                #     step_size[parameter] = 1e-9
                pGB_low = maxpGB01[parameter] #- step_size[parameter]/2
                pGB_high = maxpGB01[parameter] + step_size[parameter]
                # print(parameter, step_size[parameter],i)
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
            print(step_size['Amplitude'],inner_product['Amplitude']['Amplitude'],step_size['Frequency'],inner_product['Frequency']['Frequency'])
        return inner_product

def objective(n,tdi_fs,Tobs):
    print('index',n)
    search = Search(n,tdi_fs,Tobs)
    # search.plot(search.pGBs)
    pGBmodes =  search.search()
    maxpGB, pGB =  search.optimize(pGBmodes)
    return maxpGB, pGB

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

Radler = True
version = '1'
reduction = 2

if Radler:
    DATAPATH = grandparent+"/LDC/Radler/data/"
    SAVEPATH = grandparent+"/LDC/pictures/LDC1-4/"
    SAVEPATH = grandparent+"/LDC/Radler/LDC1-4_evaluation/"
else:
    DATAPATH = grandparent+"/LDC/Sangria/data/"
    SAVEPATH = grandparent+"/LDC/pictures/Sangria/"

if Radler:
    # sangria_fn = DATAPATH + "/LDC1-4_GB_v2.hdf5"
    sangria_fn = DATAPATH + "/LDC1-4_GB_v"+version+".hdf5"
    # sangria_fn = DATAPATH + "/LDC1-3_VGB_v2_FD_noiseless.hdf5"
    # sangria_fn = DATAPATH + "/LDC1-3_VGB_v2.hdf5"
else:
    sangria_fn = DATAPATH + "/LDC2_sangria_training_v2.h5"
fid = h5py.File(sangria_fn)


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
tdi_ts = dict([(k, TimeSeries(td[k][:int(len(td[k][:])/reduction)], dt=dt)) for k in ["X", "Y", "Z"]])
tdi_fs = xr.Dataset(dict([(k, tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
GB = fastGB.FastGB(delta_t=dt, T=Tobs)  # in seconds

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


def tdi_subtraction(tdi_fs,found_sources_mp_subtract, frequencies_search):
    #subtract the found sources from original
    tdi_fs_subtracted2 = deepcopy(tdi_fs)
    for i in range(len(found_sources_mp_subtract)):
        # for j in range(len(found_sources_to_subtract[i])):
            Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_mp_subtract[i], oversample=4, simulator="synthlisa")
            source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
            index_low = np.searchsorted(tdi_fs_subtracted2["X"].f, Xs_subtracted.f[0])
            index_high = index_low+len(Xs_subtracted)
            for k in ["X", "Y", "Z"]:
                tdi_fs_subtracted2[k].data[index_low:index_high] -= source_subtracted[k].data
    return tdi_fs_subtracted2

try:
    # cat = np.load(DATAPATH+'cat_sorted_all.npy', allow_pickle = True)
    cat = np.load(SAVEPATH+'/cat_sorted_v'+version+'.npy', allow_pickle = True)
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

    parameter_list = ['Name'] + parameters
    cat_array = []
    for parameter in parameter_list:
        cat_array.append(np.append(np.append(params_vgb[names_vgb.index(parameter)],params_igb[names_igb.index(parameter)]),params_dgb[names_dgb.index(parameter)]))
    cat = np.rec.fromarrays(cat_array, names=parameter_list)

    indexes = np.argsort(cat['Frequency'])
    cat = cat[indexes]
    np.save(SAVEPATH+'cat_sorted_all.npy',cat)

cat_sorted = cat

cat_window = cat_sorted[cat_sorted['Frequency']>0.005208333333333333]
cat_window = cat_window[cat_window['Frequency']<0.005208333333333333+512/62914560]


DAt = (tdi_ts['Z'] - tdi_ts['X'])/np.sqrt(2.0)
DEt = (tdi_ts['Z'] - 2.0*tdi_ts['Y'] + tdi_ts['X'])/np.sqrt(6.0)
DTt = (tdi_ts['Z'] + tdi_ts['Y'] + tdi_ts['X'])/np.sqrt(3.0)

DAf = (tdi_fs.Z - tdi_fs.X)/np.sqrt(2.0)
DEf = (tdi_fs.Z - 2.0*tdi_fs.Y + tdi_fs.X)/np.sqrt(6.0)
DTf = (tdi_fs.Z + tdi_fs.Y + tdi_fs.X)/np.sqrt(3.0)


# fig = plt.figure(figsize=fig_size)
# plt.semilogy(tdi_fs.f[tdi_fs.f<0.03][1000:]*10**3,np.abs(DAf[tdi_fs.f<0.03][1000:]), 'k')
# plt.xlabel('Frequency (mHz)')
# plt.ylabel('|A|')
# plt.show()

# fig = plt.figure(figsize=fig_size)
# plt.plot(tdi_ts['X'].t,DAt, 'k')
# plt.xlabel('Time (s)')
# plt.ylabel('TDI A')
# plt.show()


# LDC1-4 #####################################
frequencies = []
frequencies_even = []
frequencies_odd = []
# search_range = [0.00398, 0.0041]
# search_range = [0.0039885, 0.0040205]
# search_range = [0.0039935, 0.0039965]
f_Nyquist = 1/dt/2
search_range = [0.0003, f_Nyquist]
# search_range = [0.0002, 0.0003]
# search_range = [0.0019935, 0.0020135]
# search_range = [0.0029935, 0.0030135]
# window_length = 1*10**-7 # Hz
number_of_windows = 0
current_frequency = search_range[0]
while current_frequency < search_range[1]:
    f_smear = current_frequency *3* 10**-4
    # if current_frequency < 0.004:
    #     f_smear = current_frequency *3* 10**-4 *(1+ 3/0.004*(0.004 -current_frequency))
        # f_smear = current_frequency *3* 10**-4 *(1+ 4*np.log10(0.004 -current_frequency))
    f_deviation = frequency_derivative(current_frequency,2)*Tobs
    # window_length = np.max([f_smear, f_deviation])
    window_length = f_smear + f_deviation
    window_length += 4*32*10**-9*2
    upper_limit = current_frequency+window_length
    frequencies.append([current_frequency, upper_limit])
    current_frequency = deepcopy(upper_limit)
    number_of_windows += 1

# frequencies = frequencies[:32]
frequencies_even = frequencies[::2]
frequencies_odd = frequencies[1::2]

start_index = np.searchsorted(np.asarray(frequencies_odd)[:,0], 0.0040489)-1

# save_name = 'Sangria_1_full_cut'
# save_name = 'Radler_1_full'
save_name = 'LDC1-4_2_optimized_second' ### ETH submission
save_name = 'Montana'
# save_name = 'APC'
# save_name = 'LDC1-4_half_year'

# duration = '3932160'
# duration = '7864320'
# duration = '15728640'
duration = '31457280'
save_name = 'Montana2022_'+duration
frequencies_search = frequencies
# batch_index = int(sys.argv[1])
batch_index = 28
start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.003977)-1
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.019)
start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.004654)-1
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.007977)
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], cat_sorted[-2]['Frequency'])-1
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.0004)-1
batch_size = 2
# start_index = batch_size*batch_index
print('batch',batch_index, start_index)
# frequencies_search = frequencies_search[start_index:start_index+batch_size]
# print(i, frequencies_search[0])
### highest + padding has to be less than f Nyqist
while frequencies_search[-1][1] + (frequencies_search[-1][1] - frequencies_search[-1][0])/2 > f_Nyquist:
    frequencies_search = frequencies_search[:-1]

frequencies_min_montana = 0.00398148268199693
frequencies_max_montana = 0.004003419101309928
frequencies_min_APC = 0.003999281
frequencies_max_APC = 0.00410031292
frequencies_min = 0.00031
frequencies_max = 0.0005
start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], frequencies_min)-1
end_index = np.searchsorted(np.asarray(frequencies_search)[:,0], frequencies_max)
frequencies_search = frequencies_search[start_index:end_index]

search_range = [frequencies_search[0][0],frequencies_search[-1][1]]
# search_range = [1619472*10**-8,2689639*10**-8]
print('search range '+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))))

# i = 44
# search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
# search1.plot()

def l2_norm_match(pGB_injected, pGB_found):
    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB_found, oversample=4, simulator="synthlisa")
    Xs_injected, Ys_injected, Zs_injected = GB.get_fd_tdixyz(template=pGB_injected, oversample=4, simulator="synthlisa")
    Xs_aligned = xr.align(Xs_injected, Xs, join='left',fill_value=0)[1]
    Ys_aligned = xr.align(Ys_injected, Ys, join='left',fill_value=0)[1]
    Zs_aligned = xr.align(Zs_injected, Zs, join='left',fill_value=0)[1]
        
    fmin, fmax = float(Xs_injected.f[0]), float(Xs_injected.f[-1] + Xs_injected.attrs["df"])
    freq = np.array(Xs_injected.sel(f=slice(fmin, fmax)).f)
    Nmodel = get_noise_model(noise_model, freq)

    Af = (Zs_aligned - Xs_aligned)/np.sqrt(2.0)
    Ef = (Zs_aligned - 2.0*Ys_aligned + Xs_aligned)/np.sqrt(6.0)
    Af_injected = (Zs_injected - Xs_injected)/np.sqrt(2.0)
    Ef_injected = (Zs_injected - 2.0*Ys_injected + Xs_injected)/np.sqrt(6.0)
    frequency_T_threshold = 19.1*10**-3/2
    if Af.f[-1] > frequency_T_threshold:
        Tf = (Zs_aligned + Ys_aligned + Xs_aligned)/np.sqrt(3.0)
        Tf_injected = (Zs_injected + Ys_injected + Xs_injected)/np.sqrt(3.0)
        norm = np.sqrt(np.linalg.norm(Af_injected - Af.data)**2 + np.linalg.norm(Ef_injected - Ef.data)**2 + np.linalg.norm(Tf_injected - Tf.data)**2)
        norm_s = np.sqrt(np.linalg.norm(Af.data)**2 + np.linalg.norm(Ef.data)**2 + np.linalg.norm(Tf.data)**2)
        norm_relative = norm/norm_s
    else:
        norm = np.sqrt(np.linalg.norm(Af_injected - Af.data)**2 + np.linalg.norm(Ef_injected - Ef.data)**2)
        norm_s = np.sqrt(np.linalg.norm(Af.data)**2 + np.linalg.norm(Ef.data)**2)
        norm_relative = norm/norm_s
    return norm_relative

def SNR_match(pGB_injected, pGB_found):
    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB_found, oversample=4, simulator="synthlisa")
    Xs_injected, Ys_injected, Zs_injected = GB.get_fd_tdixyz(template=pGB_injected, oversample=4, simulator="synthlisa")
    Xs_aligned = xr.align(Xs_injected, Xs, join='left',fill_value=0)[1]
    Ys_aligned = xr.align(Ys_injected, Ys, join='left',fill_value=0)[1]
    Zs_aligned = xr.align(Zs_injected, Zs, join='left',fill_value=0)[1]
        
    fmin, fmax = float(Xs_injected.f[0]), float(Xs_injected.f[-1] + Xs_injected.attrs["df"])
    freq = np.array(Xs_injected.sel(f=slice(fmin, fmax)).f)
    Nmodel = get_noise_model(noise_model, freq)
    SA = Nmodel.psd(freq=freq, option="A")
    ST = Nmodel.psd(freq=freq, option="T")

    Af = (Zs_aligned - Xs_aligned)/np.sqrt(2.0)
    Ef = (Zs_aligned - 2.0*Ys_aligned + Xs_aligned)/np.sqrt(6.0)
    Af_injected = (Zs_injected - Xs_injected)/np.sqrt(2.0)
    Ef_injected = (Zs_injected - 2.0*Ys_injected + Xs_injected)/np.sqrt(6.0)
    frequency_T_threshold = 19.1*10**-3/2
    if Af.f[-1] > frequency_T_threshold:
        Tf = (Zs_aligned + Ys_aligned + Xs_aligned)/np.sqrt(3.0)
        Tf_injected = (Zs_injected + Ys_injected + Xs_injected)/np.sqrt(3.0)
        SNR2 = np.sum( np.real(Af_injected * np.conjugate(Af.data) + Ef_injected * np.conjugate(Ef.data))/SA + np.real(Tf_injected * np.conjugate(Tf.data))/ST)
        hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2)/SA + np.absolute(Tf.data)**2 /ST)
        ss = np.sum((np.absolute(Af_injected.data)**2 + np.absolute(Ef_injected.data)**2)/SA + np.absolute(Tf_injected.data)**2 /ST)
    else:
        SNR2 = np.sum( np.real(Af_injected * np.conjugate(Af.data) + Ef_injected * np.conjugate(Ef.data))/SA)
        hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /SA)
        ss = np.sum((np.absolute(Af_injected.data)**2 + np.absolute(Ef_injected.data)**2) /SA)
    SNR = 4.0*Xs.df* hh
    SNR2 = 4.0*Xs.df* SNR2
    SNR3 = SNR2 / (np.sqrt(SNR)*np.sqrt(4.0*Xs.df* ss))
    return SNR3.values

def SNR_match_amplitude_condsiered(pGB_injected, pGB_found):
    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB_found, oversample=4, simulator="synthlisa")
    Xs_injected, Ys_injected, Zs_injected = GB.get_fd_tdixyz(template=pGB_injected, oversample=4, simulator="synthlisa")
    Xs_aligned = xr.align(Xs_injected, Xs, join='left',fill_value=0)[1]
    Ys_aligned = xr.align(Ys_injected, Ys, join='left',fill_value=0)[1]
    Zs_aligned = xr.align(Zs_injected, Zs, join='left',fill_value=0)[1]
        
    fmin, fmax = float(Xs_injected.f[0]), float(Xs_injected.f[-1] + Xs_injected.attrs["df"])
    freq = np.array(Xs_injected.sel(f=slice(fmin, fmax)).f)
    Nmodel = get_noise_model(noise_model, freq)
    SA = Nmodel.psd(freq=freq, option="A")
    ST = Nmodel.psd(freq=freq, option="T")

    Af = (Zs_aligned - Xs_aligned)/np.sqrt(2.0)
    Ef = (Zs_aligned - 2.0*Ys_aligned + Xs_aligned)/np.sqrt(6.0)
    Af_injected = (Zs_injected - Xs_injected)/np.sqrt(2.0)
    Ef_injected = (Zs_injected - 2.0*Ys_injected + Xs_injected)/np.sqrt(6.0)
    frequency_T_threshold = 19.1*10**-3/2
    if Af.f[-1] > frequency_T_threshold:
        Tf = (Zs_aligned + Ys_aligned + Xs_aligned)/np.sqrt(3.0)
        Tf_injected = (Zs_injected + Ys_injected + Xs_injected)/np.sqrt(3.0)
        SNR2 = np.sum( np.real(Af_injected * np.conjugate(Af.data) + Ef_injected * np.conjugate(Ef.data))/SA + np.real(Tf_injected * np.conjugate(Tf.data))/ST)
        hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2)/SA + np.absolute(Tf.data)**2 /ST)
        ss = np.sum((np.absolute(Af_injected.data)**2 + np.absolute(Ef_injected.data)**2)/SA + np.absolute(Tf_injected.data)**2 /ST)
    else:
        SNR2 = np.sum( np.real(Af_injected * np.conjugate(Af.data) + Ef_injected * np.conjugate(Ef.data))/SA)
        hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /SA)
        ss = np.sum((np.absolute(Af_injected.data)**2 + np.absolute(Ef_injected.data)**2) /SA)
    SNR = 4.0*Xs.df* hh
    scalar_ss = 4.0*Xs.df* ss
    SNR2 = 4.0*Xs.df* SNR2
    amplitude_factor = np.abs(np.log10(scalar_ss/SNR))/2+1
    SNR3 = SNR2 / (np.sqrt(SNR)*np.sqrt(scalar_ss))/amplitude_factor
    return SNR3.values, amplitude_factor, SNR2 / (np.sqrt(SNR)*np.sqrt(scalar_ss))

try:
    found_sources_in_flat = np.load(SAVEPATH+'found_sources' +save_name+'_flat.npy', allow_pickle = True)
except:
    found_sources_in_flat = np.load(SAVEPATH+'found_sources' +save_name+'.npy', allow_pickle = True)
# found_sources_mp = np.load(SAVEPATH+'found_sources' +save_name+'.npy', allow_pickle = True)

# found_sources_mp_best = []
# found_sources_mp_all = []
# for i in range(len(found_sources_mp)):
#     found_sources_mp_best.append(found_sources_mp[i][0])
#     found_sources_in_window = []
#     for j in range(len(found_sources_mp[i][1])):
#         found_sources_in_window.append(found_sources_mp[i][1][j][0][0])
#     found_sources_mp_all.append(found_sources_in_window)

# found_sources_in_flat = []
# number_of_found_flat = 0
# for i in range(len(found_sources_mp)):
#     for j in range(len(found_sources_mp[i][3])):
#         found_sources_in_flat.append(found_sources_mp[i][3][j])
#         number_of_found_flat += 1
# found_sources_in_flat = np.asarray(found_sources_in_flat)

# np.save(SAVEPATH+'found_sources' +save_name+'_flat.npy', np.asarray(found_sources_in_flat))

found_sources_in_flat_frequency = []
for i in range(len(found_sources_in_flat)):
    found_sources_in_flat_frequency.append(found_sources_in_flat[i]['Frequency'])
found_sources_in_flat_frequency = np.asarray(found_sources_in_flat_frequency)
found_sources_in_flat = np.asarray(found_sources_in_flat)
indexes_in = np.argsort(found_sources_in_flat_frequency)
found_sources_in_flat_frequency = found_sources_in_flat_frequency[indexes_in]
found_sources_in_flat = found_sources_in_flat[indexes_in]

# np.save(SAVEPATH+'/found_sources_flat.npy', np.asarray(found_sources_in_flat))

found_sources_in = []
for i in range(len(frequencies_search)):
    lower_index = np.searchsorted(found_sources_in_flat_frequency,frequencies_search[i][0])
    higher_index = np.searchsorted(found_sources_in_flat_frequency,frequencies_search[i][1])
    found_sources_in.append(found_sources_in_flat[lower_index:higher_index])

# frequencies_search_with_found_only = []
# found_sources_in_with_found_only = []
# for i in range(len(frequencies_search)):
#     if len(found_sources_in[i]) > 0:
#         frequencies_search_with_found_only.append(frequencies_search[i])
#         found_sources_in_with_found_only.append(found_sources_in[i])
# frequencies_search = frequencies_search_with_found_only

new_dt = np.dtype(cat_sorted.dtype.descr + [('IntrinsicSNR','<f8')])
cat_sorted_SNR = np.zeros(cat_sorted.shape, dtype=new_dt)
for parameter in parameters:
    cat_sorted_SNR[parameter] = cat_sorted[parameter]
cat_sorted = cat_sorted_SNR

# pGB_injected_SNR = np.load(SAVEPATH+'/found_sources_pGB_injected_in_out_intrinsic_SNR_sorted'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', allow_pickle = True)
# pGB_injected_SNR = np.load(SAVEPATH+'/found_sources_pGB_injected_in_out_intrinsic_SNR_sortedsecond'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', allow_pickle = True)

# pGB_injected_no_overlap = []
# pGB_injected_flat = []
# for i in range(len(pGB_injected_SNR)):
#     pGB_injected_no_overlap.append([])
#     list_count = 0
#     for j in range(len(pGB_injected_SNR[i])):
#         if pGB_injected_SNR[i][j]['Frequency'] > frequencies_search[i][0] and pGB_injected_SNR[i][j]['Frequency'] < frequencies_search[i][1]:
#             pGB_injected_flat.append(pGB_injected_SNR[i][j])
#             pGB_injected_no_overlap[-1].append(pGB_injected_SNR[i][j])
#             list_count += 1
#         if list_count > 35:
#             break
# pGB_injected_flat_array = np.asarray(pGB_injected_flat)
# pGB_injected_SNR = np.load(SAVEPATH+'/found_sources_pGB_injected_in_out_intrinsic_SNR_sorted30000to3316929LDC1-4_half_year.npy', allow_pickle = True)

# for i in range(len(pGB_injected)):
#     pGB_injected[i] = np.asarray(pGB_injected[i])
# pGB_injected_flat = np.concatenate(np.asarray(pGB_injected))

def get_SNR(pGB_injected, lower_frequency, upper_frequency):
    search1 = Search(tdi_fs,Tobs, lower_frequency, upper_frequency)
    intrinsic_SNR_injected = []
    print(lower_frequency)
    for i in range(len(pGB_injected)):
        pGB_injected_dict = {}
        for parameter in parameters:
            try:
                pGB_injected_dict[parameter] = pGB_injected[i][parameter]
            except:
                pGB_injected_dict[parameter] = pGB_injected.iloc[i][parameter]
        intrinsic_SNR_injected.append(search1.intrinsic_SNR([pGB_injected_dict]))
        if i > 30:
            break
    return intrinsic_SNR_injected

try:
    found_sources_in = np.load(SAVEPATH+'/found_sources_in_SNR'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', allow_pickle= True)

# delete 7400 element from montana solution
# found_sources_in[5218] = np.delete(found_sources_in[5218],2)
# np.save(SAVEPATH+'/found_sources_in_SNR'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', np.asarray(found_sources_in))

except:
    #### parallel found sources
    input = []
    index = []
    for i in range(len(found_sources_in)):
    # for i in range(16):
        if len(found_sources_in[i]) > 0:
            index.append(i)
            input.append((found_sources_in[i],frequencies_search[i][0], frequencies_search[i][1]))
    start = time.time()
    pool = mp.Pool(16)
    SNR_intrinsic = pool.starmap(get_SNR, input)
    pool.close()
    pool.join()
    print('time to calculate SNR for', len(frequencies_search), 'windows: ', time.time()-start)
    for i in range(len(SNR_intrinsic)):
        for j in range(len(SNR_intrinsic[i])):
            if len(found_sources_in[index[i]]) > 0:
                found_sources_in[index[i]][j]['IntrinsicSNR'] = SNR_intrinsic[i][j]
    np.save(SAVEPATH+'/found_sources_in_SNR'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', np.asarray(found_sources_in))
# found_sources_in = np.load(SAVEPATH+'/found_sources_in_SNR30000to3312283Montana2022_31457280.npy', allow_pickle= True)

# try:
#     pGB_injected = np.load(SAVEPATH+'/found_sources_pGB_injected_no_SNR'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', allow_pickle = True)
# except:
#     get_injected_sources = pd.DataFrame(cat_sorted)
#     pGB_injected = []
#     for j in range(len(frequencies_search)):
#         if j%10 == 0:
#             print(j)
#         # padding = (frequencies_search[j][1] - frequencies_search[j][0])/2*0
#         # index_low = np.searchsorted(get_injected_sources['Frequency'], frequencies_search[j][0])
#         # index_high = np.searchsorted(get_injected_sources['Frequency'], frequencies_search[j][1])
#         # try:
#         #     if get_injected_sources['Frequency'][index_high] < frequencies_search[j][1]:
#         #         index_high -= 1
#         # except:
#         #     pass
#         # pGB_injected.append(get_injected_sources[index_low:index_high].sort_values(by='Amplitude', ascending=False))

#         pGB_injected.append(get_injected_sources[(get_injected_sources['Frequency']> frequencies_search[j][0]) & (get_injected_sources['Frequency']< frequencies_search[j][1])].sort_values(by='Amplitude', ascending=False))

#     np.save(SAVEPATH+'/found_sources_pGB_injected_no_SNR'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', np.asarray(pGB_injected))

# number_of_injected_signals= 0
# for i in range(len(pGB_injected)):
#     number_of_injected_signals += len(pGB_injected[i])

# pGB_injected_flat = np.concatenate(pGB_injected)



#### parallel
# input = []
# index = []
# for i in range(len(pGB_injected)):
# # for i in range(16):
#     if len(pGB_injected[i]) > 0:
#         index.append(i)
#         input.append((pGB_injected[i],frequencies_search[i][0], frequencies_search[i][1]))
# start = time.time()
# pool = mp.Pool(16)
# SNR_intrinsic = pool.starmap(get_SNR, input)
# pool.close()
# pool.join()
# print('time to calculate SNR for', len(frequencies_search), 'windows: ', time.time()-start)
# np.save(SAVEPATH+'/SNR_intrinsic' +save_name+'.npy', np.asarray(SNR_intrinsic))
# np.save(SAVEPATH+'/index_intrinsic' +save_name+'.npy', np.asarray(index))
# for i in range(len(SNR_intrinsic)):
#     for j in range(len(SNR_intrinsic[i])):
#         if len(pGB_injected[index[i]]) > 0:
#             pGB_injected[index[i]]['IntrinsicSNR'].iloc[j] = SNR_intrinsic[i][j]
# ## outside of normal range
# for i in range(len(pGB_injected)):
#     pGB_injected[i] = pGB_injected[i][:50]

# for i in range(len(pGB_injected)):
#     pGB_injected[i] = pGB_injected[i].to_records()



#### sequential
# for i in range(len(pGB_injected)):
#     print(i)
#     # if i != 5:
#     #     continue
#     if len(pGB_injected[i]) > 0:
#         search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
#     for j in range(len(pGB_injected[i])):
#         pGB_injected_dict = {}
#         for parameter in parameters:
#             pGB_injected_dict[parameter] = pGB_injected[i].iloc[j][parameter]
#         pGB_injected[i].iloc[j]['IntrinsicSNR'] = search1.intrinsic_SNR([pGB_injected_dict])
#         if j > 300:
#             break
        # print('SNR for noise model', noise_model, intrinsic_SNR_injected[-1], 'loglikelihood ratio',search1.loglikelihood([pGB_injected[i][j]]), 'SNR data',search1.loglikelihood_SNR([pGB_injected[i][j]]))

# pGB_injected_SNR_sorted = []
# for i in range(len(pGB_injected)):
#     indexesSNR = np.argsort(-pGB_injected[i]['IntrinsicSNR'])
#     pGB_injected_SNR_sorted.append(pGB_injected[i].iloc[indexesSNR])
# pGB_injected = pGB_injected_SNR_sorted




# np.save(SAVEPATH+'/found_sources_pGB_injectedSNR'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', np.asarray(pGB_injected))
# np.save(SAVEPATH+'/found_sources_pGB_injected_in_out_intrinsic_SNR_sortedthird'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', np.asarray(pGB_injected_SNR_sorted))
# pGB_injected = np.load(SAVEPATH+'/found_sources_pGB_injected'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', allow_pickle = True)
# pGB_injected_even = np.load(SAVEPATH+'/found_sources_pGB_injected_array2_half_year30000to3305084LDC1-4_half_even_T'+'.npy', allow_pickle = True)
# pGB_injected_odd = np.load(SAVEPATH+'/found_sources_pGB_injected_in_out_intrinsic_SNR_sorted30035to3316929LDC1-4_half_odd_T'+'.npy', allow_pickle = True)

if Radler:
    pGB_injected = np.load(DATAPATH+'/pGB_injectedthird30000to3293090LDC1-4_2_years_full.npy', allow_pickle = True)
else:
    pGB_injected = np.load(DATAPATH+'found_sources_pGB_injectedSNR30000to9261573Sangria_1_full_opt2.npy', allow_pickle = True)
    for i in range(len(pGB_injected)):
        pGB_injected[i] = pGB_injected[i].to_records()



### get injected on desired windows
pGB_injected_flat = np.concatenate(pGB_injected)
indexes_in = np.argsort(pGB_injected_flat['Frequency'])
pGB_injected_flat = pGB_injected_flat[indexes_in]

pGB_injected = []
for j in range(len(frequencies_search)):
    # padding = (frequencies_search[j][1] - frequencies_search[j][0])/2*0
    index_low = np.searchsorted(pGB_injected_flat['Frequency'], frequencies_search[j][0])
    index_high = np.searchsorted(pGB_injected_flat['Frequency'], frequencies_search[j][1])
    try:
        if pGB_injected_flat['Frequency'][index_high] < frequencies_search[j][1]:
            index_high -= 1
    except:
        pass
    pGB_injected.append(pGB_injected_flat[index_low:index_high])
pGB_injected = np.array(pGB_injected)
for i in range(len(pGB_injected)):
    # pGB_injected_amplitude = pGB_injected[i]['IntrinsicSNR']
    pGB_injected[i] = pGB_injected[i][np.argsort(-pGB_injected[i]['IntrinsicSNR'])]

pGB_injected_no_overlap = []
for i in range(len(frequencies_search)):
    pGB_injected_no_overlap.append(pGB_injected[i][(pGB_injected[i]['Frequency']> frequencies_search[i][0]) & (pGB_injected[i]['Frequency']< frequencies_search[i][1])])

### get injected windows with overlaps
# pGB_injected = np.concatenate((pGB_injected_even, pGB_injected_odd))

pGB_injected_flat = np.concatenate(pGB_injected)
indexesfrequency = np.argsort(pGB_injected_flat['Frequency'])
pGB_injected_flat = pGB_injected_flat[indexesfrequency]
get_injected_sources = pGB_injected_flat
get_pGB_injected = True
if get_pGB_injected:
    pGB_injected_overlap = []
    for j in range(len(frequencies_search)):
        padding = (frequencies_search[j][1] - frequencies_search[j][0])/2
        index_low = np.searchsorted(get_injected_sources['Frequency'], frequencies_search[j][0]-padding)
        index_high = np.searchsorted(get_injected_sources['Frequency'], frequencies_search[j][1]+padding)
        try:
            if get_injected_sources['Frequency'][index_high] < frequencies_search[j][1]:
                index_high -= 1
        except:
            pass
        # indexesA = np.argsort(-get_injected_sources[index_low:index_high]['Amplitude'])
        pGB_injected_overlap.append(get_injected_sources[index_low:index_high])

pGB_injected_flat_overlap = np.concatenate(pGB_injected_overlap)
pGB_injected_SNR_sorted_overlap = []
for i in range(len(pGB_injected_overlap)):
    # indexesSNR = np.argsort(-pGB_injected_overlap[i]['IntrinsicSNR'])
    indexesSNR = np.argsort(-pGB_injected_overlap[i]['Amplitude'])
    pGB_injected_SNR_sorted_overlap.append(pGB_injected_overlap[i][indexesSNR])

### Reduce the number of signals per window
for i in range(len(pGB_injected_overlap)):
    pGB_injected_SNR_sorted_overlap[i] = pGB_injected_SNR_sorted_overlap[i][:50]
# pGB_injected_overlap = pGB_injected_SNR_sorted_overlap
# pGB_injected_overlap_flat = np.concatenate(pGB_injected_overlap)

try:
    found_sources_in = np.load(SAVEPATH+'/found_sources_not_anticorrelated'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', allow_pickle=True)
except:
    found_sources_anitcorrelate = []
    found_sources_not_anitcorrelated = deepcopy(found_sources_in)
    number_of_matched_signals = 0
    correlation_list = []
    start = time.time()
    percentage = 0
    for i in range(len(found_sources_in)):
        found_sources_anitcorrelate.append([])
        # if i > 100:
        #     continue
        try:
            if i%int(len(found_sources_in)/100) == 0:
                percentage += 1
                print(percentage,'%')
        except:
            pass
        for j in range(len(found_sources_in[i])):
            found_match = False
            correlation_list_of_one_signal = []
            for k in range(len(found_sources_in[i])):
                pGB_injected_dict = {}
                found_dict = {}
                for parameter in parameters:
                    pGB_injected_dict[parameter] = found_sources_in[i][k][parameter]
                    found_dict[parameter] = found_sources_in[i][j][parameter]
                correlation = SNR_match(pGB_injected_dict,found_dict)
                correlation_list_of_one_signal.append(correlation)
                if k > 19:
                    print('k')
                    break
            if 0 == len(correlation_list_of_one_signal):
                break
            max_index = np.argmin(correlation_list_of_one_signal)
            if correlation_list_of_one_signal[max_index] < -0.90:
                found_match = True
            if found_match:
                found_sources_anitcorrelate[-1].append(found_sources_in[i][j])
                correlation_list.append(correlation_list_of_one_signal[max_index])
                found_sources_not_anitcorrelated[i][j] = None
                found_sources_not_anitcorrelated[i][max_index] = None
                number_of_matched_signals += 1
    print('time to match', time.time()-start)
    for i in range(len(found_sources_not_anitcorrelated)):
        found_sources_not_anitcorrelated[i] = list(filter(None, found_sources_not_anitcorrelated[i]))
    found_sources_not_anitcorrelated_flat = np.concatenate(found_sources_not_anitcorrelated)
    found_sources_anitcorrelate_flat = np.concatenate(found_sources_anitcorrelate)

    found_sources_not_anitcorrelated2 = deepcopy(found_sources_in)
    correlation_list2 = []
    for i in range(int(len(found_sources_anitcorrelate_flat)/2)):
        index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0], found_sources_anitcorrelate_flat[i*2]['Frequency'])-1
        for j in range(len(found_sources_in[index_of_interest_to_plot])):
            found_match = False
            correlation_list_of_one_signal = []
            for k in range(len(found_sources_in[index_of_interest_to_plot])):
                pGB_injected_dict = {}
                found_dict = {}
                for parameter in parameters:
                    pGB_injected_dict[parameter] = found_sources_in[index_of_interest_to_plot][k][parameter]
                    found_dict[parameter] = found_sources_in[index_of_interest_to_plot][j][parameter]
                correlation = SNR_match(pGB_injected_dict,found_dict)
                correlation_list_of_one_signal.append(correlation)
                if k > 19:
                    print('k')
                    break
            if 0 == len(correlation_list_of_one_signal):
                break
            max_index = np.argmin(correlation_list_of_one_signal)
            if correlation_list_of_one_signal[max_index] < -0.97:
                found_match = True
            if found_match:
                print('found anti',index_of_interest_to_plot,j,max_index)
                correlation_list2.append(correlation_list_of_one_signal[max_index])
                found_sources_not_anitcorrelated2[index_of_interest_to_plot][j] = None
                found_sources_not_anitcorrelated2[index_of_interest_to_plot][max_index] = None
    for i in range(len(found_sources_not_anitcorrelated2)):
        found_sources_not_anitcorrelated2[i] = list(filter(None, found_sources_not_anitcorrelated2[i]))

    ## determine index of a specific frequency
    for j in range(int(len(found_sources_anitcorrelate_flat)/2)):
        # if j != 0:
        #     continue
        index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0], found_sources_anitcorrelate_flat[j*2]['Frequency'])-1
        #### plot strains
        # for i in range(len(frequencies_search)):
        #     if i != index_of_interest_to_plot:
        #         continue
        #     search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
        #     found_extended = found_sources_not_anitcorrelated2[i]
        #     if len(pGB_injected_SNR_sorted_overlap[i]) > 20:
        #         search1.plot(found_sources_in=found_extended, pGB_injected= pGB_injected_SNR_sorted_overlap[i][:20], saving_label =SAVEPATH+'/strain added'+ str(int(np.round(frequencies_search[i][0]*10**8))) +save_name+'anitcorrelation.png') 
        #     else:
        #         search1.plot(found_sources_in=found_extended, pGB_injected= pGB_injected_SNR_sorted_overlap[i], saving_label =SAVEPATH+'/strain added'+ str(int(np.round(frequencies_search[i][0]*10**8))) +save_name+'anitcorrelation.png') 
        #         # search1.plot(found_sources_in=found_sources_mp_best[i], pGB_injected=pGB_injected[i][:10], pGB_injected_matched= matched_extended, saving_label =SAVEPATH+'/strain added'+ str(int(np.round(frequencies_search[i][0]*10**8))) +save_name+'in.png') 

    found_sources_not_anitcorrelated2_array = []
    for i in range(len(found_sources_not_anitcorrelated2)):
        found_sources_not_anitcorrelated2_array.append(np.asarray(found_sources_not_anitcorrelated2[i]))
    np.save(SAVEPATH+'/found_sources_not_anticorrelated'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', found_sources_not_anitcorrelated2)

found_sources_in_flat = np.concatenate(found_sources_in)

found_sources_in_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_in_flat]) for attribute in found_sources_in_flat[0].keys()}
found_sources_in_flat_df = pd.DataFrame(found_sources_in_flat_array)

##### save yaml file
# import yaml
# member = deepcopy(found_sources_in_flat)
# for i in range(len(member)):
#     for parameter in parameters:
#         member[i][parameter] = str(member[i][parameter])
#     member[i]['IntrinsicSNR'] = str(member[i]['IntrinsicSNR'])
# data = {}
# data["author"] = np.asarray(["Stefan Strub"])
# data["e-mail"] = np.asarray(["stefan.strub@erdw.ethz.ch"])
# data["date"] = np.asarray(["20.07.2022"])
# data["challenge"] = np.asarray(["LDC1-4"])
# data["dataset"] = np.asarray(["LDC1-4_GB_v2"])
# data["esimates"] = member
# with open(SAVEPATH+'/ETH_LDC1-4_4mHz.yaml', 'w') as file:
#     documents = yaml.dump(member, file)

# #### parallel
# input = []
# index = []
# for i in range(len(found_sources_in)):
# # for i in range(16):
#     if len(found_sources_in[i]) > 0:
#         index.append(i)
#         input.append((found_sources_in[i],frequencies_search[i][0], frequencies_search[i][1]))
# start = time.time()
# pool = mp.Pool(16)
# SNR_intrinsic = pool.starmap(get_SNR, input)
# pool.close()
# pool.join()
# print('time to calculate SNR for', len(frequencies_search), 'windows: ', time.time()-start)
# for i in range(len(SNR_intrinsic)):
#     for j in range(len(SNR_intrinsic[i])):
#         if len(SNR_intrinsic[i]) > 0:
#             found_sources_in[index[i]][j]['IntrinsicSNR'] = SNR_intrinsic[i][j]
# np.save(SAVEPATH+'/found_sources_not_anitcorrelatedsecond' +save_name+'.npy', found_sources_in)

if version == '1':
    for i in range(len(found_sources_in)):
        for j in range(len(found_sources_in[i])):
            found_sources_in[i][j]['InitialPhase'] = -(found_sources_in[i][j]['InitialPhase'])
    for i in range(len(pGB_injected_SNR_sorted_overlap)):
        for j in range(len(pGB_injected_SNR_sorted_overlap[i])):
            pGB_injected_SNR_sorted_overlap[i][j]['InitialPhase'] = -(pGB_injected_SNR_sorted_overlap[i][j]['InitialPhase'])

############## match
def match_function(found_sources_in, pGB_injected_not_matched, found_sources_not_matched):
    pGB_injected_matched = []
    found_sources_matched = []
    correlation_list = []
    for j in range(len(found_sources_in)):
        found_match = False
        correlation_list_of_one_signal = []
        found_dict = {}
        for parameter in parameters:
            found_dict[parameter] = found_sources_in[j][parameter]
        for k in range(len(pGB_injected_not_matched)):
            pGB_injected_dict = {}
            for parameter in parameters:
                pGB_injected_dict[parameter] = pGB_injected_not_matched[k][parameter]
            # correlation = l2_norm_match(pGB_injected_dict,found_dict)
            # correlation = SNR_match(pGB_injected_dict,found_dict)
            correlation, amplitude_factor, cross_correlation = SNR_match_amplitude_condsiered(pGB_injected_dict,found_dict)
            correlation_list_of_one_signal.append(correlation)
            if k > 39:
                break
        if 0 == len(correlation_list_of_one_signal):
            break
        max_index = np.nanargmax(correlation_list_of_one_signal)
        if correlation_list_of_one_signal[max_index] > 0.5:
            found_match = True
        if found_match:
            pGB_injected_matched.append(pGB_injected_not_matched[max_index])
            found_sources_matched.append(found_sources_in[j])
            # pGB_injected_not_matched = np.delete(pGB_injected_not_matched, max_index)
            found_sources_not_matched[j] = None
        correlation_list.append(correlation_list_of_one_signal[max_index])
    return found_sources_in, pGB_injected_not_matched, correlation_list, found_sources_not_matched, pGB_injected_matched, found_sources_matched
## 0.73, 0.39
# 0.93,  0.81
do_match_parallelized = True
if do_match_parallelized:
    pGB_injected_matched = []
    found_sources_matched = []
    correlation_list = []
    amplitude_factor = []
    pGB_injected_not_matched = deepcopy(pGB_injected_SNR_sorted_overlap)
    found_sources_not_matched = deepcopy(found_sources_in)
    number_of_matched_signals = 0
    input = []
    for i in range(len(found_sources_in)):
        input.append((found_sources_in[i], pGB_injected_not_matched[i], found_sources_not_matched[i]))
    start = time.time()
    pool = mp.Pool(16)
    matches = pool.starmap(match_function, input)
    pool.close()
    pool.join()
    print('time to match', time.time()-start)

    for i in range(len(matches)):
        pGB_injected_matched.append([])
        found_sources_matched.append([])
        correlation_list.append([])
        amplitude_factor.append([])
        found_sources_in[i], pGB_injected_not_matched[i], correlation_list[i], found_sources_not_matched[i], pGB_injected_matched[i], found_sources_matched[i] = matches[i]

    found_sources_matched_flat = np.concatenate(found_sources_matched)
    for i in range(len(found_sources_not_matched)):
        try:
            found_sources_not_matched[i] = list(filter(None, found_sources_not_matched[i]))
        except:
            pass

do_match_sequential = False
if do_match_sequential:
    pGB_injected_matched = []
    found_sources_matched = []
    pGB_injected_not_matched = deepcopy(pGB_injected_SNR_sorted_overlap)
    found_sources_not_matched = deepcopy(found_sources_in)
    number_of_matched_signals = 0
    correlation_list = []
    start = time.time()
    percentage = 0
    for i in range(len(found_sources_in)):
    # for i in range(3300, 3400):
        # if i != index_of_interest_to_plot:
        #     continue
        pGB_injected_matched.append([])
        found_sources_matched.append([])
        # if i < 3310 or i > 3350:
        #     continue
        try:
            if i%int(len(found_sources_in)/100) == 0:
                percentage += 1
                print(percentage,'%')
        except:
            pass
        for j in range(len(found_sources_in[i])):
            found_match = False
            correlation_list_of_one_signal = []
            correlation_list_of_one_signal_original = []
            correlation_list_of_one_signal_amplitude = []
            found_dict = {}
            for parameter in parameters:
                found_dict[parameter] = found_sources_in[i][j][parameter]
            for k in range(len(pGB_injected_not_matched[i])):
                pGB_injected_dict = {}
                for parameter in parameters:
                    pGB_injected_dict[parameter] = pGB_injected_not_matched[i][k][parameter]
                l2_norm = l2_norm_match(pGB_injected_dict,found_dict)
                correlation_apmlitude, amplitude_factor, cross_correlation = SNR_match_amplitude_condsiered(pGB_injected_dict,found_dict)
                correlation_list_of_one_signal.append(l2_norm)
                correlation_list_of_one_signal_original.append(np.array(cross_correlation))
                correlation_list_of_one_signal_amplitude.append(np.array(correlation_apmlitude))
                if k > 39:
                    break
            if 0 == len(correlation_list_of_one_signal):
                break
            max_index = np.nanargmax(correlation_list_of_one_signal)
            if correlation_list_of_one_signal[max_index] > 0.5:
                found_match = True
            if found_match:
                pGB_injected_matched[-1].append(pGB_injected_not_matched[i][max_index])
                found_sources_matched[-1].append(found_sources_in[i][j])
                pGB_injected_not_matched[i] = np.delete(pGB_injected_not_matched[i], max_index)
                correlation_list.append(correlation_list_of_one_signal[max_index])
                found_sources_not_matched[i][j] = None
                number_of_matched_signals += 1
    print('time to match', time.time()-start)
    
    for i in range(len(found_sources_not_matched)):
        found_sources_not_matched[i] = list(filter(None, found_sources_not_matched[i]))

# pGB_injected_not_matched_with_overlap_of_windows = deepcopy(pGB_injected_not_matched)
# # pGB_injected_not_matched_with_overlap_of_windows2 = deepcopy(pGB_injected_not_matched_with_overlap_of_windows)
# pGB_injected_not_matched = []
# for i in range(len(pGB_injected_not_matched_with_overlap_of_windows)):
#     pGB_injected_not_matched.append([])
#     list_count = 0
#     for j in range(len(pGB_injected_not_matched_with_overlap_of_windows[i])):
#         if pGB_injected_not_matched_with_overlap_of_windows[i][j]['Frequency'] > frequencies_search[i][0] and pGB_injected_not_matched_with_overlap_of_windows[i][j]['Frequency'] < frequencies_search[i][1]:
#             pGB_injected_not_matched[i].append(pGB_injected_not_matched_with_overlap_of_windows[i][j])
#             list_count += 1
#         if list_count > 19:
#             break

# pGB_injected_flat = []
# for i in range(len(pGB_injected)):
#     for j in range(len(pGB_injected[i])):
#         if pGB_injected[i][j]['Frequency'] > frequencies_search[i][0] and pGB_injected[i][j]['Frequency'] < frequencies_search[i][1]:
#             pGB_injected_flat.append(pGB_injected[i][j])
# number_of_injected_signals = len(pGB_injected_flat)

# for i in range(len(pGB_injected_matched)):
#     # if i != 5:
#     #     continue
#     if len(pGB_injected_matched[i]) > 0:
#         search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
#     for j in range(len(pGB_injected_matched[i])):
#         pGB_dict = {}
#         for parameter in parameters:
#             pGB_dict[parameter] = pGB_injected_matched[i][j][parameter]
#         pGB_injected_matched[i][j]['IntrinsicSNR'] = search1.intrinsic_SNR([pGB_dict])

# for i in range(len(pGB_injected_not_matched)):
#     # if i != 5:
#     #     continue
#     if len(pGB_injected_not_matched[i]) > 0:
#         search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
#     for j in range(len(pGB_injected_not_matched[i])):
#         pGB_dict = {}
#         for parameter in parameters:
#             pGB_dict[parameter] = pGB_injected_not_matched[i][j][parameter]
#         pGB_injected_not_matched[i][j]['IntrinsicSNR'] = search1.intrinsic_SNR([pGB_dict])

# count = 0
# for i in range(len(found_sources_in)):
#     count += 1
#     print(count)
#     if len(found_sources_in[i]) > 0:
#         search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
#     for j in range(len(found_sources_in[i])):
#         pGB_dict = {}
#         for parameter in parameters:
#             pGB_dict[parameter] = found_sources_in[i][j][parameter]
#         found_sources_in[i][j]['IntrinsicSNR'] = search1.intrinsic_SNR([pGB_dict])

# for i in range(len(found_sources_matched)):
#     if len(found_sources_matched[i]) > 0:
#         search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
#     for j in range(len(found_sources_matched[i])):
#         pGB_dict = {}
#         for parameter in parameters:
#             pGB_dict[parameter] = found_sources_matched[i][j][parameter]
#         found_sources_matched[i][j]['IntrinsicSNR'] = search1.intrinsic_SNR([pGB_dict])

# for i in range(len(found_sources_not_matched)):
#     if len(found_sources_not_matched[i]) > 0:
#         search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
#     for j in range(len(found_sources_not_matched[i])):
#         pGB_dict = {}
#         for parameter in parameters:
#             pGB_dict[parameter] = found_sources_not_matched[i][j][parameter]
#         found_sources_not_matched[i][j]['IntrinsicSNR'] = search1.intrinsic_SNR([pGB_dict])     


end_string = 'amplitude_considered_sqrt_03mHz'
# end_string = ''

found_sources_matched_array = []
for i in range(len(found_sources_matched)):
    found_sources_matched_array.append(np.asarray(found_sources_matched[i]))
found_sources_matched_array = np.asarray(found_sources_matched_array)
found_sources_not_matched_array = []
for i in range(len(found_sources_not_matched)):
    found_sources_not_matched_array.append(np.asarray(found_sources_not_matched[i]))
found_sources_not_matched_array = np.asarray(found_sources_not_matched_array)
np.save(SAVEPATH+'/found_sources_matched' +save_name+end_string+'.npy', found_sources_matched_array)
np.save(SAVEPATH+'/found_sources_not_matched' +save_name+end_string+'.npy', found_sources_not_matched_array)
np.save(SAVEPATH+'/injected_not_matched_windows' +save_name+end_string+'.npy', pGB_injected_not_matched)
np.save(SAVEPATH+'/injected_matched_windows' +save_name+end_string+'.npy', pGB_injected_matched)
np.save(SAVEPATH+'/correlation_list' +save_name+end_string+'.npy', correlation_list)

# found_sources_matched = np.load(SAVEPATH+'/found_sources_matched' +save_name+end_string+'.npy', allow_pickle=True)
# found_sources_not_matched = np.load(SAVEPATH+'/found_sources_not_matched' +save_name+end_string+'.npy', allow_pickle=True)
# pGB_injected_not_matched = np.load(SAVEPATH+'/injected_not_matched_windows' +save_name+end_string+'.npy', allow_pickle=True)
# pGB_injected_matched = np.load(SAVEPATH+'/injected_matched_windows' +save_name+end_string+'.npy', allow_pickle=True)

# found_sources_not_matched_high_f_old = np.concatenate(found_sources_not_matched_old)
# found_sources_not_matched_high_f = np.concatenate(found_sources_not_matched)
# found_sources_matched_high_f_old = np.concatenate(found_sources_matched_old)
# found_sources_matched_high_f = np.concatenate(found_sources_matched)
# pGB_injected_high_f = np.concatenate(pGB_injected)

# found_sources_matched_high_f_old_array = {attribute: np.asarray([x[attribute] for x in found_sources_matched_high_f_old]) for attribute in found_sources_matched_high_f_old[0].keys()}
# found_sources_matched_high_f_old_df = pd.DataFrame(found_sources_matched_high_f_old_array)

# found_sources_matched_high_f_array = {attribute: np.asarray([x[attribute] for x in found_sources_matched_high_f]) for attribute in found_sources_matched_high_f[0].keys()}
# found_sources_matched_high_f_df = pd.DataFrame(found_sources_matched_high_f_array)

# fig = plt.figure()
# plt.semilogy(found_sources_matched_high_f_old_df['Frequency'],found_sources_matched_high_f_old_df['IntrinsicSNR'], 'r.')
# plt.semilogy(found_sources_matched_high_f_df['Frequency'],found_sources_matched_high_f_df['IntrinsicSNR'], 'g.')
# plt.semilogy(pGB_injected_matched['Frequency'],pGB_injected_matched['IntrinsicSNR'], 'o.')
# plt.show()

number_of_injected_signals = 0
for i in range(len(pGB_injected)):
    number_of_injected_signals += len(pGB_injected[i])
number_of_injected_signals_SNR_high = 0
for i in range(len(pGB_injected)):
    number_of_injected_signals_SNR_high += len(pGB_injected[i][pGB_injected[i]['IntrinsicSNR']>10])
number_of_found_signals = 0
for i in range(len(found_sources_in)):
    number_of_found_signals += len(found_sources_in[i])

### 32, 128, 187, 222, 235, 257, 268, 274
### 128, 222
number_of_found_signals_not_matched = 0
for i in range(len(found_sources_not_matched)):
    number_of_found_signals_not_matched += len(found_sources_not_matched[i])
number_of_matched_signals = len(np.concatenate(found_sources_matched))
print(number_of_matched_signals ,'matched signals out of', number_of_injected_signals , 'injected signals and',number_of_found_signals, 'found signals')
print('sensitivity = matched signals/injected signals:', np.round(number_of_matched_signals/number_of_injected_signals,2))
print('number_of_injected_signals_SNR_high:', np.round(number_of_injected_signals_SNR_high,2))
print('matched signals/found signals:', np.round(number_of_matched_signals/number_of_found_signals,2))

### match frequency histogram
# found_sources_matched2 = np.load(SAVEPATH+'/found_sources_matched' +save_name+'.npy', allow_pickle=True)

# number_of_found_signals = []
# for i in range(len(frequencies_search)):
#     number_of_found_signals.append(len(found_sources_matched[i]))
# number_of_matched_signals_array2 = []
# for i in range(len(frequencies_search)):
#     number_of_matched_signals_array2.append(len(found_sources_matched2[i]))
# number_of_matched_signals_2 = 0
# for i in range(len(found_sources_in)):
#     number_of_matched_signals_2 += len(found_sources_matched2[i])
# figure = plt.figure()
# plt.semilogx(frequencies_search,number_of_matched_signals_array2, 'r.')
# plt.semilogx(frequencies_search,number_of_found_signals, '.', alpha=0.5)
# plt.show()


found_sources_matched_flat = np.concatenate(found_sources_matched)
found_sources_matched_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_matched_flat]) for attribute in found_sources_matched_flat[0].keys()}
found_sources_matched_flat_df = pd.DataFrame(found_sources_matched_flat_array)
found_sources_not_matched_flat = np.concatenate(found_sources_not_matched)
found_sources_not_matched_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_not_matched_flat]) for attribute in found_sources_not_matched_flat[0].keys()}
found_sources_not_matched_flat_df = pd.DataFrame(found_sources_not_matched_flat_array)

pGB_injected_reduced = []
for i in range(len(pGB_injected)):
    try:
        pGB_injected_reduced.append(pGB_injected[i].iloc[:40])
    except:
        pGB_injected_reduced.append(pGB_injected[i])
pGB_injected_flat_reduced = np.concatenate(np.asarray(pGB_injected_reduced))   
pGB_injected_flat_reduced_df = pd.DataFrame(np.asarray(pGB_injected_flat_reduced))
pGB_injected_flat_reduced_df = pGB_injected_flat_reduced_df.sort_values(by='Frequency')
pGB_injected_flat_highSNR_df = pGB_injected_flat_reduced_df[pGB_injected_flat_reduced_df['IntrinsicSNR']>10]

pGB_injected_matched_flat = []
for i in range(len(pGB_injected_matched)):
    for j in range(len(pGB_injected_matched[i])):
        pGB_injected_matched_flat.append(pGB_injected_matched[i][j])
# pGB_injected_matched_flat = np.asarray(pGB_injected_matched_flat)
pGB_injected_matched_flat_df = pd.DataFrame(np.asarray(pGB_injected_matched_flat))
pGB_injected_matched_flat_df = pGB_injected_matched_flat_df.sort_values(by= 'Frequency')

pGB_injected_not_matched_flat_df = pGB_injected_flat_reduced_df.append(pGB_injected_matched_flat_df)
pGB_injected_not_matched_flat_df = pGB_injected_not_matched_flat_df.drop_duplicates(keep=False)
# amplitude_considered_
found_sources_matched_flat_df.to_pickle(SAVEPATH+'/found_sources_matched' +save_name+end_string+'_df')
found_sources_not_matched_flat_df.to_pickle(SAVEPATH+'/found_sources_not_matched' +save_name+end_string+'_df')
pGB_injected_matched_flat_df.to_pickle(SAVEPATH+'/injected_matched_windows' +save_name+end_string+'_df')
pGB_injected_not_matched_flat_df.to_pickle(SAVEPATH+'/injected_not_matched_windows' +save_name+end_string+'_df')

# found_sources_matched_flat_df = pd.read_pickle(SAVEPATH+'/found_sources_matched' +save_name+'df')
# found_sources_not_matched_flat_df = pd.read_pickle(SAVEPATH+'/found_sources_not_matched' +save_name+'df')
# pGB_injected_matched_flat_df = pd.read_pickle(SAVEPATH+'/injected_matched_windows' +save_name+'df')
# pGB_injected_not_matched_flat_df = pd.read_pickle(SAVEPATH+'/injected_not_matched_windows' +save_name+'df')


found_sources_matched2_flat = np.concatenate(found_sources_matched2)
found_sources_matched2_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_matched2_flat]) for attribute in found_sources_matched2_flat[0].keys()}
found_sources_matched2_flat_df = pd.DataFrame(found_sources_matched2_flat_array)

figure = plt.figure()
# histogram on linear scale
plt.subplot(211)

hist, bins, _ = plt.hist(found_sources_matched_flat_df['IntrinsicSNR'], bins=16)
hist, bins2, _ = plt.hist(found_sources_matched2_flat_df['IntrinsicSNR'], bins=16)

# histogram on log scale. 
# Use non-equal bin sizes, such that they look equal on log scale.
logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
logbins2 = np.logspace(np.log10(bins2[0]),np.log10(bins2[-1]),len(bins2))
plt.subplot(212)
plt.hist(found_sources_matched_flat_df['IntrinsicSNR'], bins=logbins, alpha=0.5, color='blue')
plt.hist(found_sources_matched2_flat_df['IntrinsicSNR'], bins=logbins2, alpha=0.5, color='red')
plt.xscale('log')
# plt.yscale('log')
plt.show()

plt.hist(found_sources_matched_flat_df['IntrinsicSNR'])
plt.show()

# for i in range(len(pGB_injected_not_matched)):
#     # if i != 6:
#     #     continue
#     # for j in range(len(found_sources_matched[i])):
#     #     #subtract the found sources from original
#     #     tdi_fs_subtracted = deepcopy(tdi_fs)
#     #     for n in range(len( found_sources_matched[i][:j])):
#     #         Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_matched[i][n], oversample=4, simulator="synthlisa")
#     #         source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
#     #         index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
#     #         index_high = index_low+len(Xs_subtracted)
#     #         for k in ["X", "Y", "Z"]:
#     #             tdi_fs_subtracted[k].data[index_low:index_high] = tdi_fs_subtracted[k].data[index_low:index_high] - source_subtracted[k].data
#     if len(found_sources_not_matched[i]) > 0:
#         search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
#     for j in range(len(pGB_injected_not_matched[i])):
#         pGB_injected_not_matched[i][j]['IntrinsicSNR'] = search1.intrinsic_SNR([pGB_injected_not_matched[i][j]])
#         print('SNR for noise model', noise_model, search1.intrinsic_SNR([pGB_injected_not_matched[i][j]]), 'loglikelihood ratio',search1.loglikelihood([pGB_injected_not_matched[i][j]]), 'SNR data',search1.loglikelihood_SNR([pGB_injected_not_matched[i][j]]), 'frequency', pGB_injected_not_matched[i][j]['Frequency'])
#     for j in range(len(found_sources_not_matched[i])):
#         found_sources_not_matched[i][j]['IntrinsicSNR'] = search1.intrinsic_SNR([found_sources_not_matched[i][j]])
#         # print('SNR for noise model found', noise_model, search1.intrinsic_SNR([found_sources_not_matched[i][j]]), 'loglikelihood ratio',search1.loglikelihood([found_sources_not_matched[i][j]]), 'SNR data',search1.loglikelihood_SNR([found_sources_not_matched[i][j]]), 'frequency', found_sources_not_matched[i][j]['Frequency'])


found_sources_in_flat = np.concatenate(found_sources_in)
found_sources_in_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_in_flat]) for attribute in found_sources_in_flat[0].keys()}
found_sources_in_flat_df = pd.DataFrame(found_sources_in_flat_array)


found_sources_in = []
for i in range(len(frequencies_search)):
    lower_index = np.searchsorted(found_sources_in_flat_frequency,frequencies_search[i][0])
    higher_index = np.searchsorted(found_sources_in_flat_frequency,frequencies_search[i][1])
    if i == 102:
        print(lower_index, higher_index)
    found_sources_in.append(found_sources_in_flat[lower_index:higher_index])


found_sources_in_SNR_sorted = []
for i in range(len(found_sources_in)):
    listSNR = []
    for j in range(len(found_sources_in[i])):
        listSNR.append(-found_sources_in[i][j]['IntrinsicSNR'])
    indexesSNR = np.argsort(listSNR)
    found_sources_in_SNR_sorted.append(np.array(found_sources_in[i])[indexesSNR])
found_sources_in = deepcopy(found_sources_in_SNR_sorted)




### determine index of a specific frequency
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.00264612)
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.002084612)
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.003988)+1
index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.005208333333333333)
index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.004)
index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.0037305)
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.0018776)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.001851)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.002825)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.009297)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.001197)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.003988)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.008605)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.004654)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.00966)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.01022)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.01137)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.012144005913068922)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.012594655254592595)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.012911925072033257)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.012962249851778103)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.01306191115728545)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.014131028697661994)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.014880504053657114)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.015381574175293712)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.016056465377498327)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.01803821797758112)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.02364600963097699)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.04296)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.05994535926230327)-1

# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  pGB_injected_not_matched_flat_df.iloc[-1]['Frequency'])-1
# index_of_interest_to_plot = 4710
# + shift 4671
index_of_interest_to_plot = 128
# index_of_interest_to_plot = 273
# index_of_interest_to_plot = 448
# index_of_interest_to_plot = 458
# index_of_interest_to_plot = 4
# index_of_interest_to_plot = 392
#plot strains
number_of_windows = 1
for i in range(len(frequencies_search)):
    if i != index_of_interest_to_plot:
        continue
    search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i+number_of_windows-1][1])
    vertical_lines = []
    for j in range(number_of_windows):
        vertical_lines.append(frequencies_search[i+j][0])
    found_extended = np.concatenate(found_sources_matched[i:i+number_of_windows])
    found_not_matched_extended = np.concatenate(found_sources_not_matched[i:i+number_of_windows])
    # matched_extended = np.concatenate(pGB_injected_matched[i:i+number_of_windows])
    # if len(pGB_injected_SNR[i]) > 0:
    pGB_injected_dict_list = []
    matched_extended = []
    not_matched_extended = []
    for j in range(number_of_windows):
        for k in range(len(pGB_injected_SNR_sorted_overlap[i+j])):
            pGB_injected_dict_list.append({})
            for parameter in parameters:
                pGB_injected_dict_list[-1][parameter] = pGB_injected_SNR_sorted_overlap[i+j][k][parameter]
            if k > 20:
                break
    for j in range(number_of_windows+2):
        for k in range(len(pGB_injected_matched[i+j-1])):
            matched_extended.append({})
            for parameter in parameters:
                matched_extended[-1][parameter] = pGB_injected_matched[i+j-1][k][parameter]
    
    # padding = (frequencies_search[i+number_of_windows-1][1] - frequencies_search[i][0])/2
    # pGB_injected_not_matched_flat_df_red = pGB_injected_not_matched_flat_df[pGB_injected_not_matched_flat_df['Frequency'] > frequencies_search[i][0]-padding]
    # pGB_injected_not_matched_flat_df_red = pGB_injected_not_matched_flat_df_red[pGB_injected_not_matched_flat_df_red['Frequency'] < frequencies_search[i+number_of_windows-1][1]+padding]
    # for k in range(len(pGB_injected_not_matched_flat_df_red)):
    #     not_matched_extended.append({})
    #     for parameter in parameters:
    #         not_matched_extended[-1][parameter] = pGB_injected_not_matched_flat_df_red.iloc[k][parameter]

    # found_not_matched_extended = []
    # found_not_matched_extended = found_not_matched_extended[:1]
    # not_matched_extended = []
    # found_extended = deepcopy(matched_extended)

    # found_extended = []
    # matched_extended = [matched_extended[3],matched_extended[4]]
    # matched_extended = []
    save_name_path = SAVEPATH+'/strain added Amplitude'+ str(int(np.round(frequencies_search[i][0]*10**8))) +save_name+str(int(len(matched_extended)))+'kristen.png'
    if len(pGB_injected_SNR_sorted_overlap[i]) > 20:
        search1.plot(found_sources_in=found_extended, found_sources_not_matched = found_not_matched_extended, pGB_injected= not_matched_extended[:20],  pGB_injected_matched= matched_extended, vertical_lines= vertical_lines, saving_label =save_name_path) 
    else:
        search1.plot(found_sources_in=found_extended, found_sources_not_matched = found_not_matched_extended, pGB_injected= not_matched_extended,  pGB_injected_matched= matched_extended, vertical_lines= vertical_lines, saving_label =save_name_path) 
        # search1.plot(found_sources_in=found_sources_mp_best[i], pGB_injected=pGB_injected[i][:10], pGB_injected_matched= matched_extended, saving_label =SAVEPATH+'/strain added'+ str(int(np.round(frequencies_search[i][0]*10**8))) +save_name+'in.png') 

print(SNR_match(pGB_injected_dict_list[0],found_not_matched_extended[0]))
print(search1.SNR([pGB_injected_dict_list[0]]))
print(search1.SNR(pGB_injected_dict_list[:5]))
print(search1.SNR(found_extended[:5]))
print(search1.SNR(matched_extended[:5]))
print(search1.SNR([found_extended[2]]))
print(search1.SNR([found_not_matched_extended[0]]))
print(search1.SNR([found_not_matched_extended[1]]))
print(search1.loglikelihood([found_not_matched_extended[0]]))
print(search1.loglikelihood([found_not_matched_extended[1]]))
print(search1.SNR([found_not_matched_extended[1],found_not_matched_extended[0]]))
print(search1.SNR([pGB_injected_dict_list[1],pGB_injected_dict_list[2]]))
print(search1.SNR([found_not_matched_extended[0],found_not_matched_extended[1],found_extended[0],found_extended[1],found_extended[2],found_extended[3]]))
for j in range(len(found_sources_mp_all[index_of_interest_to_plot])):
    print(search1.SNR([found_sources_mp_all[index_of_interest_to_plot][j]]))
for j in range(len(found_sources_in[index_of_interest_to_plot])):
    print(search1.SNR([found_sources_in[index_of_interest_to_plot][j]]))

found_sources_mp_all_df = pd.DataFrame(found_sources_mp_all[index_of_interest_to_plot])
found_sources_best_df = pd.DataFrame(found_sources_mp_best[index_of_interest_to_plot])
found_sources_mp_all_df2 = found_sources_mp_all_df.append(found_sources_best_df)
found_sources_mp_all_df2 = found_sources_mp_all_df2.drop_duplicates(subset=['Frequency'], keep=False)
found_sources_out = found_sources_mp_best[index_of_interest_to_plot][-2:]


######### check pipeline subtraction
start_index_shift = 0
found_sources_in_optimize = []
found_sources_out_optimize = []
found_sources = []
### subtract neighbors
tdi_fs_neighbors_subtracted = deepcopy(tdi_fs)
found_sources_to_subtract = np.concatenate([found_sources_mp[start_index_shift+index_of_interest_to_plot-1][3],found_sources_mp[start_index_shift+index_of_interest_to_plot+1][3]])
for j in range(len(found_sources_to_subtract)):
    Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_to_subtract[j], oversample=4, simulator="synthlisa")
    source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
    index_low = np.searchsorted(tdi_fs_neighbors_subtracted["X"].f, Xs_subtracted.f[0])
    index_high = index_low+len(Xs_subtracted)
    for k in ["X", "Y", "Z"]:
        tdi_fs_neighbors_subtracted[k].data[index_low:index_high] = tdi_fs_neighbors_subtracted[k].data[index_low:index_high] - source_subtracted[k].data
tdi_fs_subtracted = deepcopy(tdi_fs_neighbors_subtracted)
plot_subraction = True
if plot_subraction:
    search1 = Search(tdi_fs,Tobs, frequencies_search[index_of_interest_to_plot][0], frequencies_search[index_of_interest_to_plot][1])
    search1.plot(second_data= tdi_fs_subtracted, found_sources_not_matched= found_sources_mp[start_index_shift+index_of_interest_to_plot][3])
# for i in range(len(found_sources_in[index_of_interest_to_plot])):
for i in range(5):
    search_out_subtracted = Search(tdi_fs_subtracted, Tobs, frequencies_search[index_of_interest_to_plot][0], frequencies_search[index_of_interest_to_plot][1])

    lower_frequency = frequencies_search[index_of_interest_to_plot][0]
    upper_frequency = frequencies_search[index_of_interest_to_plot][1]
    # create two sets of found sources. found_sources_in with signals inside the boundary and founce_sources_out with outside sources
    found_sources_in_optimize = []
    found_sources_out_optimize = []
    for j in range(len(found_sources)):
        if found_sources[j]['Frequency'] > lower_frequency and found_sources[j]['Frequency'] < upper_frequency:
            found_sources_in_optimize.append(found_sources[j])
        else:
            found_sources_out_optimize.append(found_sources[j])

    found_sources_out_optimize = [] # ignore signals outside
    try:
        print('injected',search_out_subtracted.SNR([pGB_injected_dict_list[i]]))
        print('injected original',search1.SNR([pGB_injected_dict_list[i]]))
    except:
        pass
    SNR_list = []
    for j in range(3):
        if found_sources_mp_all[start_index_shift+index_of_interest_to_plot][i*3+j]['Amplitude'] < 0:
            print('neg amplitude', found_sources_mp_all[start_index_shift+index_of_interest_to_plot][i*3+j]['Amplitude'] )
            found_sources_mp_all[start_index_shift+index_of_interest_to_plot][i*3+j]['Amplitude'] *= -1
            print('neg amplitude', found_sources_mp_all[start_index_shift+index_of_interest_to_plot][i*3+j]['Amplitude'] )
        SNR_list.append(float(search_out_subtracted.SNR([found_sources_mp_all[start_index_shift+index_of_interest_to_plot][i*3+j]])))
    max_index = np.nanargmax(SNR_list)
    found_sources_in_optimize.append(found_sources_mp_all[start_index_shift+index_of_interest_to_plot][i*3+max_index]) ## should exclude outsiders

    A_optimized = search1.calculate_Amplitude([found_sources_in_optimize[-1]])
    found_sources_in_optimize[-1]['Amplitude'] *= A_optimized.values
    print('found',search_out_subtracted.SNR([found_sources_in_optimize[-1]]))
    print('found log l',search_out_subtracted.loglikelihood([found_sources_in_optimize[-1]]))
    print('found original',search1.SNR([found_sources_in_optimize[-1]]))

    tdi_fs_subtracted = deepcopy(tdi_fs_neighbors_subtracted)
    for j in range(len(found_sources_out_optimize)):
        Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_out_optimize[j], oversample=4, simulator="synthlisa")
        source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
        index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
        index_high = index_low+len(Xs_subtracted)
        for k in ["X", "Y", "Z"]:
            tdi_fs_subtracted[k].data[index_low:index_high] = tdi_fs_subtracted[k].data[index_low:index_high] - source_subtracted[k].data

    tdi_fs_subtracted = deepcopy(tdi_fs_neighbors_subtracted)
    search_out_subtracted = Search(tdi_fs_subtracted, Tobs, frequencies_search[index_of_interest_to_plot][0], frequencies_search[index_of_interest_to_plot][1])

    print('injected all',search_out_subtracted.SNR(pGB_injected_dict_list[:i+1]))
    print('injected all original',search1.SNR(pGB_injected_dict_list[:i+1]))
    print('found all',search_out_subtracted.SNR(found_sources_in_optimize))
    print('found all log l',search_out_subtracted.loglikelihood(found_sources_in_optimize))
    print('found all original',search1.SNR(found_sources_in_optimize))
    total_boundaries = deepcopy(search1.boundaries)
    amplitudes = []
    for k in range(len(found_sources_in_optimize)):
        amplitudes.append(found_sources_in_optimize[k]['Amplitude'])
    total_boundaries['Amplitude'] = [np.min(amplitudes),np.max(amplitudes)]
    amplitudes_length = search1.boundaries['Amplitude'][1] - search1.boundaries['Amplitude'][0]
    total_boundaries['Amplitude'] = [np.log10(total_boundaries['Amplitude'][0]), np.log10(total_boundaries['Amplitude'][1])]
    total_boundaries['Amplitude'] = [total_boundaries['Amplitude'][0] - amplitudes_length/10,total_boundaries['Amplitude'][1] + amplitudes_length/10,]
                
    start = time.time()
    found_sources_in_optimize = search_out_subtracted.optimize([found_sources_in_optimize], boundaries= total_boundaries)
    print('found optimized',search_out_subtracted.SNR(found_sources_in_optimize))
    print('global optimization time', time.time()-start)

    found_sources = found_sources_in_optimize + found_sources_out_optimize
    # found_sources = [pGB_injected_dict_list[0]]
    #subtract the found sources from original
    found_sources = pGB_injected_dict_list
    tdi_fs_subtracted = deepcopy(tdi_fs_neighbors_subtracted)
    for j in range(len(found_sources)):
        Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources[j], oversample=4, simulator="synthlisa")
        source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
        index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
        index_high = index_low+len(Xs_subtracted)
        for k in ["X", "Y", "Z"]:
            tdi_fs_subtracted[k].data[index_low:index_high] = tdi_fs_subtracted[k].data[index_low:index_high] - source_subtracted[k].data
    search_out_subtracted = Search(tdi_fs_subtracted, Tobs, frequencies_search[index_of_interest_to_plot][0], frequencies_search[index_of_interest_to_plot][1])
    plot_subraction = True
    if plot_subraction:
        search1 = Search(tdi_fs,Tobs, frequencies_search[index_of_interest_to_plot][0], frequencies_search[index_of_interest_to_plot][1])
        search1.plot(second_data= tdi_fs_subtracted, found_sources_not_matched = found_sources)

##### final global optimization
def tdi_subtraction(tdi_fs,found_sources_to_subtract):
    #subtract the found sources from original
    tdi_fs_subtracted2 = deepcopy(tdi_fs)
    for i in range(len(found_sources_to_subtract)):
        Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_to_subtract[i], oversample=4, simulator="synthlisa")
        source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
        index_low = np.searchsorted(tdi_fs_subtracted2["X"].f, Xs_subtracted.f[0])
        index_high = index_low+len(Xs_subtracted)
        for k in ["X", "Y", "Z"]:
            tdi_fs_subtracted2[k].data[index_low:index_high] -= source_subtracted[k].data
    return tdi_fs_subtracted2

i = deepcopy(index_of_interest_to_plot)
found_sources_to_subtract = np.concatenate([found_sources_in[i-1],found_sources_in[i+1]])
tdi_fs_subtracted = tdi_subtraction(tdi_fs,found_sources_to_subtract)
plot_subraction = True
if plot_subraction:
    search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
    search1.plot(second_data= tdi_fs_subtracted)

correlation_list2 = []
for i in range(len(frequencies_search)):
    if i != index_of_interest_to_plot:
        continue        
    for j in range(len(found_sources_in[i])):
            found_match = False
            # if j != 1:
            #     continue
            # print('i', i, 'j',j)
            correlation_list_of_one_signal = []
            for k in range(len(pGB_injected_SNR_sorted_overlap[i])):
                pGB_injected_dict = {}
                found_dict = {}
                for parameter in parameters:
                    pGB_injected_dict[parameter] = pGB_injected_SNR_sorted_overlap[i].iloc[k][parameter]
                    found_dict[parameter] = found_sources_in[i][j][parameter]
                # print('SNR', SNR_match(pGB_injected_not_matched[i][k],found_sources_in[i][j]),'parameter comparison:',pGB_injected_not_matched[i][k]['EclipticLatitude'],found_sources_in[i][j]['EclipticLatitude'],eclipticlongitude, found_sources_in[i][j]['EclipticLongitude'])
                correlation = SNR_match(pGB_injected_dict_list[0],found_dict)
                correlation_list_of_one_signal.append(correlation)
                if k > 19:
                    break
            correlation_list2.append(correlation_list_of_one_signal)

plt.imshow(correlation_list2)

#### plot SNR - frequency
markersize = 3
alpha = 0.5
parameter_to_plot = 'IntrinsicSNR'
fig = plt.figure()
# plt.plot(pGB_injected_flat_df['Frequency']*10**3,pGB_injected_flat_df['IntrinsicSNR'], '.', color= colors[0], label = 'Injected', markersize= markersize, alpha = alpha)
# plt.plot(pGB_injected_matched_flat_df['Frequency']*10**3,pGB_injected_matched_flat_df['IntrinsicSNR'], '.', color= colors[1], label = 'Injected matched', markersize= markersize, alpha = alpha)
# plt.plot(pGB_injected_flat_df_high_SNR['Frequency']*10**3,pGB_injected_flat_df_high_SNR['IntrinsicSNR'],'.', color= colors[1], markersize= markersize, label = 'Injected SNR > 10', alpha = alpha)
plt.plot(found_sources_matched_flat_df['Frequency']*10**3,found_sources_matched_flat_df['IntrinsicSNR'],'g.', label = 'Found matched', markersize= markersize, alpha = alpha)
# plt.plot(pGB_injected_not_matched_flat_df['Frequency']*10**3,pGB_injected_not_matched_flat_df['IntrinsicSNR'], '+', color = 'r', label = 'Injected not matched', markersize= markersize, alpha = alpha)
plt.plot(found_sources_not_matched_flat_df['Frequency']*10**3,found_sources_not_matched_flat_df['IntrinsicSNR'],'+',color= 'blue',  markerfacecolor='None', markersize= markersize, label = 'Found not matched', alpha = alpha)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('f (mHz)')
plt.xlim(0.3,100)
plt.ylim(0.1,2000)
if parameter_to_plot == 'IntrinsicSNR':
    plt.ylabel('Intrinsic SNR')
else:
    plt.ylabel(parameter_to_plot)    
plt.legend(markerscale=4, loc = 'lower right')
plt.savefig(SAVEPATH+'/Evaluation/'+parameter_to_plot+save_name+'injected_not_matched_found_matched_found_not_matched',dpi=300,bbox_inches='tight')
plt.show()

### plot match histo gramm - frequency
pGB_injected_matched_frequencies = []
for i in range(len(pGB_injected_matched)):
    for j in range(len(pGB_injected_matched[i])):
        pGB_injected_matched_frequencies.append(pGB_injected_matched[i][j]['Frequency']*10**3)
pGB_injected_not_matched_frequencies = []
for i in range(len(pGB_injected_not_matched)):
    for j in range(len(pGB_injected_not_matched[i])):
        pGB_injected_not_matched_frequencies.append(pGB_injected_not_matched[i][j]['Frequency']*10**3)
fig = plt.figure()
plt.hist(pGB_injected_matched_frequencies, 20, color = 'green')
plt.hist(pGB_injected_not_matched_frequencies, 20, color = 'red')
plt.xlabel('f (mHz)')
plt.legend()
plt.savefig(SAVEPATH+'/Evaluation/_SNR_histo'+save_name,dpi=300,bbox_inches='tight')
plt.show()


# ######### get errors
# error = []
# for i in range(len(pGB_injected_matched)):
#     error.append([])
#     for j in range(len(pGB_injected_matched[i])):
#         error[i].append({})
#         for parameter in parameters:
#             if parameter == 'EclipticLongitude':
#                 if pGB_injected_matched[i][j][parameter] > np.pi:
#                     pGB_injected_matched[i][j][parameter] -= 2*np.pi
#             error[i][j][parameter] = pGB_injected_matched[i][j][parameter] - found_sources_matched[i][j][parameter]
#             if parameter in ['EclipticLongitude', 'EclipticLatitude', 'Inclination', 'InitialPhase', 'Polarization']:
#                 error[i][j][parameter] = np.abs(np.arcsin(np.sin(pGB_injected_matched[i][j][parameter] - found_sources_matched[i][j][parameter])))
#             if parameter == 'Amplitude':
#                 # error[i][j][parameter] = np.log10(pGB_injected_matched[i][j][parameter]) - np.log10(found_sources_matched[i][j][parameter])
#                 error[i][j][parameter] = np.abs(pGB_injected_matched[i][j][parameter] - found_sources_matched[i][j][parameter])/pGB_injected_matched[i][j][parameter]
#             found_sources_matched[i][j][parameter+'Error'] = error[i][j][parameter]
# found_sources_matched_flat = np.concatenate(found_sources_matched)
# found_sources_matched_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_matched_flat]) for attribute in found_sources_matched_flat[0].keys()}
# found_sources_matched_flat_df = pd.DataFrame(found_sources_matched_flat_array)


######### get errors
def get_errors(pGB_injected_matched_flat_df, found_sources_matched_flat_df):
    error = {}
    for parameter in parameters:
        error[parameter] = []
        # found_sources_matched_flat_df[parameter+'Error'] = []
        for i in range(len(pGB_injected_matched_flat_df)):
            if parameter == 'EclipticLongitude':
                if pGB_injected_matched_flat_df[parameter][i] > np.pi:
                    pGB_injected_matched_flat_df[parameter][i] -= 2*np.pi
            if parameter in ['EclipticLongitude', 'EclipticLatitude', 'Inclination', 'InitialPhase', 'Polarization']:
                error[parameter].append(np.abs(np.arcsin(np.sin(pGB_injected_matched_flat_df[parameter][i] - found_sources_matched_flat_df[parameter][i]))))
            elif parameter == 'Amplitude':
                # error[parameter].append(np.log10(pGB_injected_matched_flat_df[parameter][i]) - np.log10(found_sources_matched_flat_df[parameter][i]))
                error[parameter].append(np.abs(pGB_injected_matched_flat_df[parameter][i] - found_sources_matched_flat_df[parameter][i])/pGB_injected_matched_flat_df[parameter][i])
            # found_sources_matched_flat_df[parameter+'Error'].append(error[parameter][i])
            else:
                error[parameter].append(pGB_injected_matched_flat_df[parameter][i] - found_sources_matched_flat_df[parameter][i])
        found_sources_matched_flat_df[parameter+'Error'] = error[parameter]
    return found_sources_matched_flat_df

found_sources_matched_flat_df = get_errors(pGB_injected_matched_flat_df, found_sources_matched_flat_df)

#### get galactic coordinates
coordinates_found = coord.SkyCoord(found_sources_matched_flat_df['EclipticLongitude'], found_sources_matched_flat_df['EclipticLatitude'], unit='rad', frame='barycentricmeanecliptic')
found_sources_matched_flat_df['GalacticLongitude'] = coordinates_found.galactic.l.value
found_sources_matched_flat_df['GalacticLatitude'] = coordinates_found.galactic.b.value

coordinates_injected = coord.SkyCoord(pGB_injected_matched_flat_df['EclipticLongitude'], pGB_injected_matched_flat_df['EclipticLatitude'], unit='rad', frame='barycentricmeanecliptic')
pGB_injected_matched_flat_df['GalacticLongitude'] = coordinates_injected.galactic.l.value
pGB_injected_matched_flat_df['GalacticLatitude'] = coordinates_injected.galactic.b.value

found_sources_matched_flat_df['SkylocationError'] = coordinates_injected.separation(coordinates_found)

error_flat = {}
for parameter in parameters:
    error_flat[parameter] = []
    for i in range(len(error)):
        for j in range(len(error[i])):
            error_flat[parameter].append(error[i][j][parameter] )
error_flat['Skylocation'] = found_sources_matched_flat_df['SkylocationError']

mean_error = {}
std_table = {}
for parameter in parameters:
    mean_error[parameter] = []
    std_table[parameter] = []
    for i in range(len(error)):
        mean = 0
        std = 0
        for j in range(len(error[i])):
            mean += error[i][j][parameter]
        for j in range(len(error[i])):
            std += (mean - error[i][j][parameter])**2
        std = np.sqrt(std)
        mean_error[parameter].append(mean)
        std_table[parameter].append(std)
mean_error = pd.DataFrame.from_dict(mean_error)
std_table = pd.DataFrame.from_dict(std_table)

###### plot errors
frequencies_search = np.asarray(frequencies_search)
for parameter in parameters:
    fig = plt.figure()
    plt.errorbar(frequencies_search[:,0],mean_error[parameter],yerr=std_table[parameter])
    plt.xlabel(parameter)
    plt.ylabel('Error')
    plt.savefig(SAVEPATH+'/Evaluation/'+parameter+'_error_bar'+save_name,dpi=300,bbox_inches='tight')
    plt.show()

##### plot SNR to error
# for parameter in  parameters +['Skylocation']:
# for parameter in  ['SkylocationError']:
x_parameter = 'Frequency'
x_parameter = 'IntrinsicSNR'
for parameter in  ['EclipticLatitude']:
    # x_parameter = parameter
    markersize = 4
    alpha = 0.5
    fig = plt.figure()
    # plt.plot(found_sources_matched_flat_df[parameter+'Error']/found_sources_matched_flat_df[parameter],found_sources_matched_flat_df[y_parameter],'.', label = 'found', markersize=markersize)
    if parameter in ['Frequency', 'FrequencyDerivative']:
        plt.scatter(found_sources_matched_flat_df[x_parameter],found_sources_matched_flat_df[parameter+'Error']/found_sources_matched_flat_df['Frequency'], c=found_sources_matched_flat_df['Frequency'])
    else:
        plt.plot(found_sources_matched_flat_df[x_parameter],found_sources_matched_flat_df[parameter+'Error'],'.', alpha =0.4, markersize=6)
        # plt.scatter(found_sources_matched_flat_df[x_parameter],found_sources_matched_flat_df[parameter], c=found_sources_matched_flat_df['Frequency'])
    plt.ylabel(parameter+' Error')
    if parameter in ['Frequency', 'FrequencyDerivative']:
        plt.ylabel('Relative '+parameter+' Error')
    plt.xlabel(x_parameter)
    plt.xscale('log')
    plt.yscale('log')
    # plt.legend(markerscale=4, loc = 'upper right')
    plt.savefig(SAVEPATH+'/Evaluation/SNRto'+parameter+'Error'+save_name)
    plt.show()

##### plot SNR to frequency, color error
parameter = 'SkylocationError'
fig = plt.figure()
fig.set_size_inches(8,10)
plt.scatter(found_sources_matched_flat_df['Frequency'],found_sources_matched_flat_df['IntrinsicSNR'], alpha =0.2)
# fig.colorbar(im, ax=ax0)
# ax0.set_title(parameter)
plt.xscale('log')
fig.tight_layout()
plt.savefig(SAVEPATH+'/Evaluation/SNRtoFrequencycolor'+parameter+save_name)
plt.show() 

### histogram
# for parameter in parameters + ['Skylocation']:
for parameter in  ['Frequency']:
    fig = plt.figure()
    if parameter == 'Skylocation':
        plt.hist(np.abs(error_flat[parameter]), bins= np.logspace(-2,2, 100))
    elif parameter == 'Frequency':
        # plt.hist(np.abs(error_flat[parameter]), bins= np.logspace(-14,-7, 100))
        
        plt.hist(np.abs(found_sources_matched_flat_df[parameter+'Error']/pGB_injected_matched_flat_df[parameter]), bins= np.logspace(-9,-4, 100), log=True)
        # plt.hist(error_flat[parameter], bins= np.linspace(-10**-7,10**-7, 100), log=True, density=True)
    else:
        plt.hist(found_sources_matched_flat_df[parameter+'Error'], bins=100, density=True)
    plt.xlabel('Error of '+parameter)
    plt.ylabel('count')
    if parameter == 'Skylocation':
        plt.xlabel('Error of '+parameter+' (deg)')
        plt.xscale('log')
    if parameter == 'Frequency':
        plt.xscale('log')
        plt.ylim(0,10**3)
        plt.xlabel('Relative error of '+parameter)
    # plt.yscale('log')
    plt.savefig(SAVEPATH+'/Evaluation/'+parameter+'_error_histogram'+save_name,dpi=300,bbox_inches='tight')
    plt.show()


#### get distance
def get_distance2(amplitude, frequency, frequency_derivative):
    M_chirp = (96/5*np.pi**(8/3)*frequency**(11/3)/frequency_derivative)**(-3/5)
    G = 6.674*10**(-11)
    c = 3*10**8
    M_c= M_chirp / (2*10**30)
    Mc_s = M_c/G*c**3
    print(M_chirp)
    distance = 2*M_chirp**(5/3)*np.pi*(2/3)*frequency**(2/3)/amplitude
    print('Mc',Mc_s)
    return distance
def get_distance(amplitude, frequency, frequency_derivative):
    c = 3*10**8
    distance = 2/(96/5*np.pi**(8/3)*frequency**(11/3)/frequency_derivative)*np.pi*(2/3)*frequency**(2/3)/amplitude*c /3.086e+19 #to kpc
    return distance

def get_distance_for_dataframe(dataframe):
    dataframe['Distance'] = np.empty(len(dataframe['Amplitude']))
    postitive_fd_mask = dataframe['FrequencyDerivative'] >= 0
    distance = get_distance(dataframe['Amplitude'][postitive_fd_mask], dataframe['Frequency'][postitive_fd_mask],  dataframe['FrequencyDerivative'][postitive_fd_mask])
    dataframe['Distance'][postitive_fd_mask] = distance
    return dataframe

# found_sources_matched_flat_df['Amplitude'] = pGB_injected_matched_flat_df['Amplitude']
# found_sources_matched_flat_df['FrequencyDerivative'] = pGB_injected_matched_flat_df['FrequencyDerivative']
found_sources_matched_flat_df = get_distance_for_dataframe(found_sources_matched_flat_df)
pGB_injected_matched_flat_df = get_distance_for_dataframe(pGB_injected_matched_flat_df)
pGB_injected_not_matched_flat_df = get_distance_for_dataframe(pGB_injected_not_matched_flat_df)
pGB_injected_flat_reduced_df = get_distance_for_dataframe(pGB_injected_flat_reduced_df) ### chaned to reduced
# pGB_injected_flat_highSNR_df = get_distance_for_dataframe(pGB_injected_flat_highSNR_df)

found_sources_matched_flat_positive_fd_df = found_sources_matched_flat_df[found_sources_matched_flat_df['FrequencyDerivative']>0]
pGB_injected_not_matched_flat_positive_fd_df = pGB_injected_not_matched_flat_df[pGB_injected_not_matched_flat_df['FrequencyDerivative']>0]
pGB_injected_flat_positive_fd_df = pGB_injected_flat_reduced_df[pGB_injected_flat_reduced_df['FrequencyDerivative']>0]  ### chaned to reduced
# pGB_injected_flat_highSNR_positive_fd_df = pGB_injected_flat_highSNR_df[pGB_injected_flat_highSNR_df['FrequencyDerivative']>0]

i = 1001
search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
boundaries = deepcopy(search1.boundaries)
boundaries['GalacticLongitude'] = [0,360]
boundaries['GalacticLatitude'] = [-90,90]

#### get galactic coordinates
coordinates_found = coord.SkyCoord(found_sources_matched_flat_df['EclipticLongitude'], found_sources_matched_flat_df['EclipticLatitude'], unit='rad', frame='barycentricmeanecliptic')
found_sources_matched_flat_df['GalacticLongitude'] = coordinates_found.galactic.l.value
found_sources_matched_flat_df['GalacticLatitude'] = coordinates_found.galactic.b.value

coordinates_injected = coord.SkyCoord(pGB_injected_matched_flat_df['EclipticLongitude'], pGB_injected_matched_flat_df['EclipticLatitude'], unit='rad', frame='barycentricmeanecliptic')
pGB_injected_matched_flat_df['GalacticLongitude'] = coordinates_injected.galactic.l.value
pGB_injected_matched_flat_df['GalacticLatitude'] = coordinates_injected.galactic.b.value

found_sources_matched_flat_df['SkylocationError'] = coordinates_injected.separation(coordinates_found)

# i = 1001
# N = 1
# search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
# start = time.time()
# for j in range(N):
#     loglikelihood = search1.loglikelihood([found_sources_mp[i][3][0]])
# print('time loglikelihood', time.time() - start)
# start = time.time()
# for j in range(N):
#     Xs, Ys, Zs = search1.GB.get_fd_tdixyz(template=found_sources_mp[i][3][0], oversample=4, simulator="synthlisa")
# print('time tdi', time.time() - start)

#### get galactic coordinates with distance
coordinates_found = coord.SkyCoord(found_sources_matched_flat_df['EclipticLongitude'], found_sources_matched_flat_df['EclipticLatitude'], found_sources_matched_flat_df['Distance'], unit='rad', frame='barycentricmeanecliptic')
found_sources_matched_flat_df['GalacticLongitude'] = coordinates_found.galactic.l.value
found_sources_matched_flat_df['GalacticLatitude'] = coordinates_found.galactic.b.value
found_sources_matched_flat_df['GalacticDistance'] = coordinates_found.galactic.distance.value

coordinates_injected = coord.SkyCoord(pGB_injected_matched_flat_df['EclipticLongitude'], pGB_injected_matched_flat_df['EclipticLatitude'], pGB_injected_matched_flat_df['Distance'], unit='rad', frame='barycentricmeanecliptic')
pGB_injected_matched_flat_df['GalacticLongitude'] = coordinates_injected.galactic.l.value
pGB_injected_matched_flat_df['GalacticLatitude'] = coordinates_injected.galactic.b.value
pGB_injected_matched_flat_df['GalacticDistance'] = coordinates_injected.galactic.distance.value

coordinates = coord.SkyCoord(pGB_injected_not_matched_flat_df['EclipticLongitude'], pGB_injected_not_matched_flat_df['EclipticLatitude'], pGB_injected_not_matched_flat_df['Distance'], unit='rad', frame='barycentricmeanecliptic')
pGB_injected_not_matched_flat_df['GalacticLongitude'] = coordinates.galactic.l.value
pGB_injected_not_matched_flat_df['GalacticLatitude'] = coordinates.galactic.b.value
pGB_injected_not_matched_flat_df['GalacticDistance'] = coordinates.galactic.distance.value

# coordinates = coord.SkyCoord(pGB_injected_flat_df['EclipticLongitude'], pGB_injected_flat_df['EclipticLatitude'], pGB_injected_flat_df['Distance'], unit='rad', frame='barycentricmeanecliptic')
# pGB_injected_flat_df['GalacticLongitude'] = coordinates.galactic.l.value
# pGB_injected_flat_df['GalacticLatitude'] = coordinates.galactic.b.value
# pGB_injected_flat_df['GalacticDistance'] = coordinates.galactic.distance.value

coordinates = coord.SkyCoord(pGB_injected_flat_highSNR_df['EclipticLongitude'], pGB_injected_flat_highSNR_df['EclipticLatitude'], pGB_injected_flat_highSNR_df['Distance'], unit='rad', frame='barycentricmeanecliptic')
pGB_injected_flat_highSNR_df['GalacticLongitude'] = coordinates.galactic.l.value
pGB_injected_flat_highSNR_df['GalacticLatitude'] = coordinates.galactic.b.value
pGB_injected_flat_highSNR_df['GalacticDistance'] = coordinates.galactic.distance.value

boundaries = deepcopy(search1.boundaries)
boundaries['GalacticLongitude'] = [0,360]
boundaries['GalacticLatitude'] = [-90,90]
boundaries['GalacticDistance'] = [np.nanmin(found_sources_matched_flat_df['GalacticDistance']),np.nanmax(found_sources_matched_flat_df['GalacticDistance'])]


#### get Galactocentric coordinates with distance
coordinates = coord.SkyCoord(found_sources_matched_flat_positive_fd_df['EclipticLongitude']* u.rad, found_sources_matched_flat_positive_fd_df['EclipticLatitude']* u.rad, found_sources_matched_flat_positive_fd_df['Distance'] * u.kpc, frame='barycentricmeanecliptic')
coordinates = coordinates.transform_to(coord.Galactocentric) 
found_sources_matched_flat_positive_fd_df['GalactocentricX'] = coordinates.galactocentric.x.value
found_sources_matched_flat_positive_fd_df['GalactocentricY'] = coordinates.galactocentric.y.value
found_sources_matched_flat_positive_fd_df['GalactocentricZ'] = coordinates.galactocentric.z.value

coordinates = coord.SkyCoord(pGB_injected_not_matched_flat_positive_fd_df['EclipticLongitude']* u.rad, pGB_injected_not_matched_flat_positive_fd_df['EclipticLatitude']* u.rad, pGB_injected_not_matched_flat_positive_fd_df['Distance'] * u.kpc, frame='barycentricmeanecliptic')
coordinates = coordinates.transform_to(coord.Galactocentric) 
pGB_injected_not_matched_flat_positive_fd_df['GalactocentricX'] = coordinates.galactocentric.x.value
pGB_injected_not_matched_flat_positive_fd_df['GalactocentricY'] = coordinates.galactocentric.y.value
pGB_injected_not_matched_flat_positive_fd_df['GalactocentricZ'] = coordinates.galactocentric.z.value

# coordinates = coord.SkyCoord(np.asarray(pGB_injected_flat_positive_fd_df['EclipticLongitude'])* u.rad, np.asarray(pGB_injected_flat_positive_fd_df['EclipticLatitude'])* u.rad, np.asarray(pGB_injected_flat_positive_fd_df['Distance']) * u.kpc, frame='barycentricmeanecliptic')
# coordinates = coordinates.transform_to(coord.Galactocentric) 
# pGB_injected_flat_positive_fd_df['GalactocentricX'] = coordinates.galactocentric.x.value
# pGB_injected_flat_positive_fd_df['GalactocentricY'] = coordinates.galactocentric.y.value
# pGB_injected_flat_positive_fd_df['GalactocentricZ'] = coordinates.galactocentric.z.value

coordinates = coord.SkyCoord(np.asarray(pGB_injected_flat_highSNR_positive_fd_df['EclipticLongitude'])* u.rad, np.asarray(pGB_injected_flat_highSNR_positive_fd_df['EclipticLatitude'])* u.rad, np.asarray(pGB_injected_flat_highSNR_positive_fd_df['Distance']) * u.kpc, frame='barycentricmeanecliptic')
coordinates = coordinates.transform_to(coord.Galactocentric) 
pGB_injected_flat_highSNR_positive_fd_df['GalactocentricX'] = coordinates.galactocentric.x.value
pGB_injected_flat_highSNR_positive_fd_df['GalactocentricY'] = coordinates.galactocentric.y.value
pGB_injected_flat_highSNR_positive_fd_df['GalactocentricZ'] = coordinates.galactocentric.z.value

boundaries = deepcopy(search1.boundaries)
boundaries['GalactocentricLongitude'] = [0,360]
boundaries['GalactocentricLatitude'] = [-90,90]

##### plot distance
markersize = 5
alpha = 0.5
fig = plt.figure()
postitive_fd_mask = found_sources_matched_flat_df['FrequencyDerivative'] >= 0
plt.plot(found_sources_matched_flat_df['GalacticLatitude'][postitive_fd_mask],found_sources_matched_flat_df['Distance'][postitive_fd_mask],'.', label = 'found', markersize=1)
# postitive_fd_mask = pGB_injected_flat_df['FrequencyDerivative'] >= 0
# plt.plot(pGB_injected_flat_df['GalacticLatitude'][postitive_fd_mask],pGB_injected_flat_df['Distance'][postitive_fd_mask],'g.', label = 'Injected', markersize= 1, zorder=1)
# postitive_fd_mask = pGB_injected_not_matched_flat_df['FrequencyDerivative'] >= 0
# plt.plot(pGB_injected_not_matched_flat_df['GalacticLatitude'][postitive_fd_mask],pGB_injected_not_matched_flat_df['Distance'][postitive_fd_mask],'+', label = 'not matched', color = 'r', markersize=2, zorder= 1)
postitive_fd_mask = pGB_injected_flat_highSNR_df['FrequencyDerivative'] >= 0
plt.plot(pGB_injected_flat_highSNR_df['GalacticLatitude'][postitive_fd_mask],pGB_injected_flat_highSNR_df['Distance'][postitive_fd_mask],'+', label = 'injected SNR>10', color = 'r', markersize=2, zorder= 4)
plt.xlabel('GalacticLatitude')
plt.ylabel('Distance [kpc]')
plt.legend(markerscale=4, loc = 'upper right')
plt.savefig(SAVEPATH+'/galacticLatitudeDistancehighSNR'+save_name)
plt.show()

fig = plt.figure()
postitive_fd_mask = found_sources_matched_flat_df['FrequencyDerivative'] >= 0
plt.plot(found_sources_matched_flat_df['GalacticLongitude'][postitive_fd_mask],found_sources_matched_flat_df['Distance'][postitive_fd_mask],'.', label = 'found', markersize=1)
# postitive_fd_mask = pGB_injected_flat_df['FrequencyDerivative'] >= 0
# plt.plot(pGB_injected_flat_df['GalacticLongitude'][postitive_fd_mask],pGB_injected_flat_df['Distance'][postitive_fd_mask],'g.', label = 'Injected', markersize= 1, zorder= 1)
postitive_fd_mask = pGB_injected_not_matched_flat_df['FrequencyDerivative'] >= 0
plt.plot(pGB_injected_not_matched_flat_df['GalacticLongitude'][postitive_fd_mask],pGB_injected_not_matched_flat_df['Distance'][postitive_fd_mask],'+', color = 'r', label = 'not matched', markersize=2, zorder= 1)
postitive_fd_mask = pGB_injected_flat_highSNR_df['FrequencyDerivative'] >= 0
plt.plot(pGB_injected_flat_highSNR_df['GalacticLongitude'][postitive_fd_mask],pGB_injected_flat_highSNR_df['Distance'][postitive_fd_mask],'+', label = 'injected SNR>10', color = 'r', markersize=2, zorder= 4)
plt.xlabel('GalacticLongitude')
plt.ylabel('Distance [kpc]')
plt.legend(markerscale=4, loc = 'upper right')
plt.savefig(SAVEPATH+'/galacticLongitudeDistancehighSNR'+save_name)
plt.show()

markersize = 5
alpha = 0.5
parameter_x = 'GalactocentricY'
parameter_y = 'GalactocentricZ'
fig = plt.figure()
plt.plot(found_sources_matched_flat_positive_fd_df[parameter_x],found_sources_matched_flat_positive_fd_df[parameter_y],'.', label = 'found', markersize=1, zorder=5)
# plt.plot(pGB_injected_flat_positive_fd_df[parameter_x],pGB_injected_flat_positive_fd_df[parameter_y],'g.', label = 'Injected', markersize= 1, zorder=1)
# plt.plot(pGB_injected_not_matched_flat_df[parameter_x],pGB_injected_not_matched_flat_df[parameter_y],'+', label = 'not matched', color = 'r', markersize=2, zorder= 1)
plt.plot(pGB_injected_flat_highSNR_positive_fd_df[parameter_x],pGB_injected_flat_highSNR_positive_fd_df[parameter_y],'+', label = 'injected SNR>10', color = 'r', markersize=2, zorder= 6)
plt.xlabel(parameter_x)
plt.ylabel(parameter_y)
plt.legend(markerscale=4, loc = 'upper right')
plt.savefig(SAVEPATH+'/'+parameter_x+parameter_y+'highSNR'+save_name)
plt.show()

##### error histogram
search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
boundaries = deepcopy(search1.boundaries)
boundaries['GalactocentricLongitude'] = [0,360]
boundaries['GalactocentricLatitude'] = [-90,90]
boundaries['Frequency'] = [frequencies_search[0][0],frequencies_search[-1][1]]
boundaries['IntrinsicSNR'] = [np.min(found_sources_matched_flat_df['IntrinsicSNR']),np.max(found_sources_matched_flat_df['IntrinsicSNR'])]
# parameter_x = 'GalacticLongitude'
# parameter_y = 'GalacticLatitude'
# parameter_x = 'EclipticLongitude'
# parameter_y = 'EclipticLatitude'
x_scale = 'log'
parameter_x = 'Frequency'
parameter_y = 'IntrinsicSNR'
n_bins = 50
x_coordinates = []
y_coordinates = []
for i in range(n_bins+1):
    length = (boundaries[parameter_x][1] - boundaries[parameter_x][0])/n_bins
    x_coordinates.append(boundaries[parameter_x][0]+length*i)
    length = (boundaries[parameter_y][1] - boundaries[parameter_y][0])/n_bins
    y_coordinates.append(boundaries[parameter_y][0]+length*i)

error = []
std = []
count = []
parameter_to_plot = 'SkylocationError'
parameter_to_plot = 'AmplitudeError'
found_sources_matched_flat_df_parameter_x_sorted = found_sources_matched_flat_df.sort_values(by=parameter_x)
# bin_boundaries_x = np.linspace(boundaries[parameter_x][0], boundaries[parameter_x][1],n_bins+1)
# bin_boundaries_x = np.log10(bin_boundaries_x)
bin_boundaries_y = np.linspace(boundaries[parameter_y][0], boundaries[parameter_y][1],n_bins+1)
# bin_boundaries_y = np.log10(bin_boundaries_y)
bin_boundaries_x = np.logspace(np.log10(boundaries[parameter_x][0]), np.log10(boundaries[parameter_x][1]),n_bins+1)
bin_boundaries_y = np.logspace(np.log10(boundaries[parameter_y][0]), np.log10(boundaries[parameter_y][1]),n_bins+1)
for i in range(n_bins):
    error.append([])
    std.append([])
    count.append([])
    # length = (boundaries[parameter_x][1] - boundaries[parameter_x][0])/n_bins
    start_index = np.searchsorted(found_sources_matched_flat_df_parameter_x_sorted[parameter_x],bin_boundaries_x[i], side='left')
    end_index = np.searchsorted(found_sources_matched_flat_df_parameter_x_sorted[parameter_x],bin_boundaries_x[i+1], side='left')
    section = found_sources_matched_flat_df_parameter_x_sorted[start_index:end_index]
    for j in range(n_bins):
        section_parameter_y_sorted = section.sort_values(by=parameter_y)
        # length = (boundaries[parameter_y][1] - boundaries[parameter_y][0])/n_bins
        start_index = np.searchsorted(section_parameter_y_sorted[parameter_y],bin_boundaries_y[j], side='left')
        end_index = np.searchsorted(section_parameter_y_sorted[parameter_y],bin_boundaries_y[j+1], side='left')
        field = section[start_index:end_index]
        change_to_deg = False
        if parameter_to_plot == 'SkylocationError':
            change_to_deg = False
        if change_to_deg:
            field[parameter_to_plot] *= 180/np.pi
        error[-1].append(np.mean(np.abs(field[parameter_to_plot])))
        std[-1].append(np.std(field[parameter_to_plot]))
        count[-1].append(len(field[parameter_to_plot]))

for i in range(len(count)):
    for j in range(len(count[i])):
        if count[i][j] == 0:
            count[i][j] = np.nan

fig, (ax0, ax1, ax2) = plt.subplots(nrows=3)
fig.set_size_inches(8,10)
im = ax0.pcolormesh(bin_boundaries_x,bin_boundaries_y, np.array(error).T)
im.set_clim(0,5)
fig.colorbar(im, ax=ax0)
ax0.set_title('mean '+parameter_to_plot)
im1 = ax1.pcolormesh(bin_boundaries_x,bin_boundaries_y, np.array(std).T)
im1.set_clim(0,5)
fig.colorbar(im1, ax=ax1)
ax1.set_title('standard deviation')
im2 = ax2.pcolormesh(bin_boundaries_x,bin_boundaries_y, np.array(count).T)
fig.colorbar(im2, ax=ax2)
ax2.set_title('number of signals')
ax0.set_yscale('log')
ax1.set_yscale('log')
ax2.set_yscale('log')
ax0.set_xscale('log')
ax1.set_xscale('log')
ax2.set_xscale('log')
ax2.set_xlabel(parameter_x)
ax0.set_ylabel(parameter_y)
ax1.set_ylabel(parameter_y)
ax2.set_ylabel(parameter_y)
fig.tight_layout()
plt.savefig(SAVEPATH+'/Evaluation/'+parameter_x+parameter_y+parameter_to_plot+save_name)
plt.show()


fig = plt.figure()
plt.title('Error')
plt.pcolormesh(x_coordinates,y_coordinates, error)
plt.xlabel(parameter_x)
plt.ylabel(parameter_y)
plt.show()


fig = plt.figure()
plt.pcolormesh(std)
plt.show()

fig = plt.figure()
plt.pcolormesh(count)
plt.show()


from fast_histogram import histogram2d
error_flat = np.concatenate(error)
error_flat_array = {attribute: np.asarray([x[attribute] for x in error_flat]) for attribute in error_flat[0].keys()}
error_flat_df = pd.DataFrame(error_flat_array)

###### plot errors histogramm
parameter_x = 'Frequency'
parameter_y = 'Amplitude'
bounds = [[found_sources_matched_flat_df[parameter_x].min(), found_sources_matched_flat_df[parameter_x].max()], [error_flat_df[parameter_y].min(),error_flat_df[parameter_y].max()]]
h = histogram2d(found_sources_matched_flat_df[parameter_x], error_flat_df[parameter_y], range=bounds, bins=100)
fig = plt.figure()
plt.imshow(h)
plt.axis('off')
plt.show()

###### plot frequency - SNR histogramm
import matplotlib.colors as colors_o
pGB_injected_flat = np.concatenate(pGB_injected)
parameter_x = 'Frequency'
parameter_y = 'Amplitude'
bounds = [[pGB_injected_flat[parameter_x].min(), pGB_injected_flat[parameter_x].max()], [pGB_injected_flat[parameter_y].min(),pGB_injected_flat[parameter_y].max()]]
h = histogram2d(pGB_injected_flat[parameter_x], pGB_injected_flat[parameter_y], range=bounds, bins=100)
fig = plt.figure()
plt.imshow(h, norm=colors_o.LogNorm(vmin=1, vmax=h.max()))
plt.axis('off')
plt.show()

bounds = [[pGB_injected_flat[parameter_x].min(), pGB_injected_flat[parameter_x].max()], [pGB_injected_flat[parameter_y].min(),pGB_injected_flat[parameter_y].max()]]
h = plt.hist2d(pGB_injected_flat[parameter_x], pGB_injected_flat[parameter_y], bins=100)
fig = plt.figure()
plt.imshow(h)
plt.axis('off')
plt.show()

###### plot errors
for parameter in parameters:
# if parameter in ['EclipticLongitude', 'EclipticLongitude', 'Inclination', 'InitialPhase', 'Polarization']:
    fig = plt.figure()
    plt.plot(found_sources_matched_flat_df['Frequency']*10**3,np.abs(found_sources_matched_flat_df[parameter+'Error']),'.', alpha= 0.1, color = 'green')
    plt.xlabel('f (mHz)')
    plt.ylabel(parameter+' error')   
    if parameter == 'Amplitude':
        plt.ylabel('Log '+parameter+' Error')    
    if parameter in ['EclipticLongitude','EclipticLatitude','Inclination','Polarization','InitialPhase']:
        plt.ylabel(parameter+' Error (radians)')    
    if parameter == 'Frequency':
        plt.ylabel(parameter+' Error (Hz)')    
    if parameter == 'FrequencyDerivative':
        plt.ylabel(parameter+' Error (Hz/s)')    
    plt.xscale('log')
    # plt.yscale('log')
    plt.savefig(SAVEPATH+'/Evaluation/'+parameter+'_error_loglog_shaded_closest_angle'+save_name,dpi=300,bbox_inches='tight')
    plt.show()

##### plot correlation histogramm
correlation_list= np.asarray(correlation_list)
correlation_list_flat = np.concatenate(correlation_list)
fig = plt.figure()
plt.hist(correlation_list_flat,50)
plt.xlabel('Correlation')
plt.ylabel('Count')
# plt.yscale('log')
plt.ylim(0,2000)
plt.xlim(0,1)
plt.savefig(SAVEPATH+'/Evaluation/correlation'+save_name+end_string,dpi=300,bbox_inches='tight')
plt.show()
