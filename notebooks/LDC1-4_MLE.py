#%%
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
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

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

    def plot(self, maxpGBs=None, pGBadded=None, second_data = None , found_sources_in= [], pGB_injected = [], pGB_injected_matched = [], added_label='Injection2', saving_label =None):
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
            ax1.plot(Af.f*10**3,Af,'k--',zorder= 1, linewidth = 2, label = 'Data subtracted')
            ax2.plot(Ef.f*10**3,np.abs(Ef),'k--',zorder= 1, linewidth = 2, label = 'Data subtracted')

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
            Zs = xr.align(self.dataZ, Zs, join='left',fill_value=0)[1]
            Af = (Zs - Xs)/np.sqrt(2.0)
            Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)
            ax1.plot(Af.f* 1000, Af.data,'--', color= colors[j%10], linewidth = 1.6)
            ax2.plot(Ef.f* 1000, np.abs(Ef.data),'--', color= colors[j%10], linewidth = 1.6)

        # ax1.plot(Xs_added2.f * 1000, Xs_added2.values.real, label="VGB2", marker=".", zorder=5)
        ax1.axvline(self.lower_frequency* 1000, color= 'red', label='Boundaries')
        ax1.axvline(self.upper_frequency* 1000, color= 'red')
        ax2.axvline(self.lower_frequency* 1000, color= 'red')
        ax2.axvline(self.upper_frequency* 1000, color= 'red')
        # ax2.axvline(self.lower_frequency* 1000- 4*32*10**-6, color= 'green')
        # ax2.axvline(self.upper_frequency* 1000+ 4*32*10**-6, color= 'green')
        # if self.reduced_frequency_boundaries != None:
        #     ax1.axvline(self.reduced_frequency_boundaries[0]* 1000, color= 'green', label='Reduced Boundaries')
        #     ax1.axvline(self.reduced_frequency_boundaries[1]* 1000, color= 'green')

        # ax1.plot(Xs.f * 1000, dataX.values.real - Xs.values.real, label="residual", alpha=0.8, color="red", marker=".")
        plt.xlabel('f (mHz)')
        ax1.set_ylabel('A')    
        ax2.set_ylabel('|E|') 
        # ax1.set_yscale('log')  
        ax2.set_yscale('log')   
        ax1.set_xlim((self.lower_frequency-self.padding)*10**3, (self.upper_frequency+self.padding)*10**3)
        ax2.set_xlim((self.lower_frequency-self.padding)*10**3, (self.upper_frequency+self.padding)*10**3)
        ax1.xaxis.set_major_locator(plt.MaxNLocator(4))
        ax2.xaxis.set_major_locator(plt.MaxNLocator(4))
        # plt.legend()
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
        if self.use_T_component:
            Tf = (Zs_total + Ys_total + Xs_total)/np.sqrt(3.0)
            hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2)/self.SA + np.absolute(Tf.data)**2 /self.ST)
            SNR2 = np.sum( np.real(self.DAf * np.conjugate(Af.data) + self.DEf * np.conjugate(Ef.data))/self.SA + np.real(self.DTf * np.conjugate(Tf.data))/self.ST )
        else:
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
        logliks = 4.0*Xs.df*( SNR2 - 0.5 * hh )
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
                res = scipy.optimize.minimize(self.function, x, args=boundaries_reduced, method='SLSQP', bounds=bounds, tol=1e-5, options= {'maxiter':100})
                # res = scipy.optimize.minimize(self.function, x, args=boundaries_reduced, method='Nelder-Mead', tol=1e-6)
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
pGBadded7['Amplitude'] = 1.45e-22*0.1
pGBadded7['EclipticLatitude'] = -0.2
pGBadded7['EclipticLongitude'] = 1.4
pGBadded7['Frequency'] = 0.00210457
pGBadded7['FrequencyDerivative'] = 1e-17
pGBadded7['Inclination'] = 0.5
pGBadded7['InitialPhase'] = 3
pGBadded7['Polarization'] = 2
pGBadded8 = {}
pGBadded8['Amplitude'] = 1.45e-22*0.1
pGBadded8['EclipticLatitude'] = -0.2
pGBadded8['EclipticLongitude'] = 1.4
pGBadded8['Frequency'] = 0.00180457
pGBadded8['FrequencyDerivative'] = 1e-17
pGBadded8['Inclination'] = 0.5
pGBadded8['InitialPhase'] = 3
pGBadded8['Polarization'] = 2
pGBadded9 = {} # low SNR
pGBadded9['Amplitude'] = 1.45e-22*0.1
pGBadded9['EclipticLatitude'] = 0.2
pGBadded9['EclipticLongitude'] = 1.5
pGBadded9['Frequency'] = 0.00190457
pGBadded9['FrequencyDerivative'] = 1e-17
pGBadded9['Inclination'] = 0.5
pGBadded9['InitialPhase'] = 3
pGBadded9['Polarization'] = 2
pGBadded10 = {}
pGBadded10['Amplitude'] = 1.36368e-22*0.1
pGBadded10['EclipticLatitude'] = -0.4
pGBadded10['EclipticLongitude'] = 1.4
pGBadded10['Frequency'] = 0.00240457
pGBadded10['FrequencyDerivative'] = 1e-17
pGBadded10['Inclination'] = 0.5
pGBadded10['InitialPhase'] = 3
pGBadded10['Polarization'] = 2
pGBadded11 = {}
pGBadded11['Amplitude'] = 1.36368e-22
pGBadded11['EclipticLatitude'] = -0.2
pGBadded11['EclipticLongitude'] = 1.4
pGBadded11['Frequency'] = 0.00201457
pGBadded11['FrequencyDerivative'] = 1e-17
pGBadded11['Inclination'] = 0.8
pGBadded11['InitialPhase'] = 2
pGBadded11['Polarization'] = 1
pGBadded12 = {}
pGBadded12['Amplitude'] = 1.36368e-22
pGBadded12['EclipticLatitude'] = 0.4
pGBadded12['EclipticLongitude'] = -1
pGBadded12['Frequency'] = 0.00201457
pGBadded12['FrequencyDerivative'] = 1e-17
pGBadded12['Inclination'] = 0.8
pGBadded12['InitialPhase'] = 2
pGBadded12['Polarization'] = 1
pGBadded13 = deepcopy(pGBadded11)
pGBadded13['Frequency'] = 0.00191457
pGBadded14 = deepcopy(pGBadded13)
pGBadded14['EclipticLatitude'] = -0.3
pGBadded14['EclipticLongitude'] = 1.5
pGBadded15 = deepcopy(pGBadded11)
pGBadded15['Frequency'] = 0.00181457
pGBadded16 = deepcopy(pGBadded15)
pGBadded16['EclipticLatitude'] = -0.4
pGBadded16['EclipticLongitude'] = 1.6
pGBadded17 = {}# low SNR
pGBadded17['Amplitude'] = 4.15e-23
pGBadded17['EclipticLatitude'] = -0.2
pGBadded17['EclipticLongitude'] = 1.4
pGBadded17['Frequency'] = 0.001104517
pGBadded17['FrequencyDerivative'] = 8*1e-18
pGBadded17['Inclination'] = 0.5
pGBadded17['InitialPhase'] = 3
pGBadded17['Polarization'] = 2
pGBadded18 = {}
pGBadded18['Amplitude'] = 4.15e-23
pGBadded18['EclipticLatitude'] = -0.2
pGBadded18['EclipticLongitude'] = 1.4
pGBadded18['Frequency'] = 0.00120457
pGBadded18['FrequencyDerivative'] = 5*1e-18
pGBadded18['Inclination'] = 0.5
pGBadded18['InitialPhase'] = 3
pGBadded18['Polarization'] = 2
pGBadded19 = {} 
pGBadded19['Amplitude'] = 4.15e-23
pGBadded19['EclipticLatitude'] = 0.2
pGBadded19['EclipticLongitude'] = 1.5
pGBadded19['Frequency'] = 0.00130457
pGBadded19['FrequencyDerivative'] = 5*1e-18
pGBadded19['Inclination'] = 1.2
pGBadded19['InitialPhase'] = 3
pGBadded19['Polarization'] = 2
pGBadded20 = {}
pGBadded20['Amplitude'] = 4.55e-23
pGBadded20['EclipticLatitude'] = 0.2
pGBadded20['EclipticLongitude'] = 1.5
pGBadded20['Frequency'] = 0.00140457
pGBadded20['FrequencyDerivative'] = 5*1e-18
pGBadded20['Inclination'] = 1.2
pGBadded20['InitialPhase'] = 3
pGBadded20['Polarization'] = 2
pGBadded21 = {}
pGBadded21['Amplitude'] = 1.45e-22*0.2
pGBadded21['EclipticLatitude'] = -0.2
pGBadded21['EclipticLongitude'] = 1.5
pGBadded21['Frequency'] = 0.00150457
pGBadded21['FrequencyDerivative'] = 5*1e-18
pGBadded21['Inclination'] = 0.9
pGBadded21['InitialPhase'] = 3
pGBadded21['Polarization'] = 2
pGBadded22 = {}
pGBadded22['Amplitude'] = 4.15e-23
pGBadded22['EclipticLatitude'] = -0.2
pGBadded22['EclipticLongitude'] = 1.5
pGBadded22['Frequency'] = 0.00160457
pGBadded22['FrequencyDerivative'] = 5*1e-18
pGBadded22['Inclination'] = 1.2
pGBadded22['InitialPhase'] = 3
pGBadded22['Polarization'] = 2
pGBadded23 = {}
pGBadded23['Amplitude'] = 4.15e-23
pGBadded23['EclipticLatitude'] = -0.2
pGBadded23['EclipticLongitude'] = 1.5
pGBadded23['Frequency'] = 0.00170457
pGBadded23['FrequencyDerivative'] = 1e-17
pGBadded23['Inclination'] = 1.2
pGBadded23['InitialPhase'] = 3
pGBadded23['Polarization'] = 2

pGBadded24 = deepcopy(pGBadded11)
pGBadded24['Amplitude'] *= 2
pGBadded24['EclipticLatitude'] = -0.3
pGBadded24['EclipticLongitude'] = 1.5
# pGBadded13 = {'Amplitude': 1.36368e-22, 'EclipticLatitude': -0.529009, 'EclipticLongitude': -2.51031, 'Frequency': 0.00125313, 'FrequencyDerivative': 9.159587298288947e-19, 'Inclination': 0.244346, 'InitialPhase': 2.64414439, 'Polarization': 2.229426357}
# pGBadded13['Frequency'] -= 2*10**-6
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

reduction = 2

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

print(pGBadded7["FrequencyDerivative"] * Tobs)
print('smear f', pGBadded7["Frequency"] *2* 10**-4)

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


class MLP_search():
    def __init__(self,tdi_fs, Tobs, signals_per_window, found_sources_previous = None, strategy = 'DE'):
        self.tdi_fs = tdi_fs
        self.Tobs = Tobs
        self.signals_per_window = signals_per_window
        self.strategy = strategy
        self.found_sources_previous = found_sources_previous

    def search(self, lower_frequency, upper_frequency):
        found_sources = []
        tdi_fs_search = deepcopy(self.tdi_fs)

        initial_guess = []
        if len(self.found_sources_previous) > 0:
            if do_subtract:
                padding_of_initial_guess_range = 0
            else:
                padding_of_initial_guess_range = (upper_frequency - lower_frequency)/2
            found_sources_previous_in_range = self.found_sources_previous[self.found_sources_previous['Frequency'] > lower_frequency-padding_of_initial_guess_range]
            found_sources_previous_in_range = found_sources_previous_in_range[found_sources_previous_in_range['Frequency'] < upper_frequency+padding_of_initial_guess_range]
            indexesA = np.argsort(-found_sources_previous_in_range['Amplitude'])
            pGB_stacked = {}
            for parameter in parameters:
                pGB_stacked[parameter] = found_sources_previous_in_range[parameter][indexesA]
            for i in range(len(found_sources_previous_in_range['Amplitude'])):
                pGBs = {}
                for parameter in parameters:
                    pGBs[parameter] = pGB_stacked[parameter][i]
                initial_guess.append(pGBs)
            
            search1 = Search(tdi_fs_search,self.Tobs, lower_frequency, upper_frequency)
            ### sort the initial guesses such that the highest loglikelihhod guess comes first
            SNR_guesses = []
            for i in range(len(initial_guess)):
                SNR_guesses.append(search1.SNR([initial_guess[i]]))
            indexes = np.argsort(SNR_guesses)[::-1]
            initial_guess = [initial_guess[i] for i in indexes]
        # indexes = np.argsort(p.get('Frequency'))
        # index_low = np.searchsorted(p.get('Frequency')[indexes], lower_frequency)
        # index_high = np.searchsorted(p.get('Frequency')[indexes], upper_frequency)
        # pGB_injected = []
        # for i in range(len(p.get('Amplitude')[indexes][index_low:index_high])):
        #     pGBs = {}
        #     for parameter in parameters:
        #         pGBs[parameter] = p.get(parameter)[indexes][index_low:index_high][i]
        #     pGB_injected.append(pGBs)
        # previous_found_sources = [{'Amplitude': 4.084935966774485e-22, 'EclipticLatitude': 0.8719934546490874, 'EclipticLongitude': 0.48611009683797857, 'Frequency': 0.003995221087430858, 'FrequencyDerivative': 1.0704703957490903e-16, 'Inclination': 1.0245091695238984, 'InitialPhase': 2.320136113624083, 'Polarization': 2.65883774239409}, {'Amplitude': 1.170377953453263e-22, 'EclipticLatitude': -1.1827019140449202, 'EclipticLongitude': -2.6708716710257203, 'Frequency': 0.003994619937260686, 'FrequencyDerivative': 9.604827167870394e-17, 'Inclination': 1.9399867466326164, 'InitialPhase': 2.468693959968005, 'Polarization': 2.5128702009090644}]
        found_sources_all = []
        number_of_evaluations_all = []
        found_sources_in = []
        current_SNR = 100
        SNR_threshold = 10
        loglikelihood_ratio_threshold = 50
        f_transfer = 19.1*10**-3
        if lower_frequency > f_transfer:
            loglikelihood_ratio_threshold = 200
        current_loglikelihood_ratio = 1000
        ind = 0
        # while current_SNR > SNR_threshold and ind < self.signals_per_window:
        while current_loglikelihood_ratio > loglikelihood_ratio_threshold and ind < self.signals_per_window:
            ind += 1
            
            search1 = Search(tdi_fs_search,self.Tobs, lower_frequency, upper_frequency)
            # N_frequency = 5
            # F_stat, frequencies_F_stat, eclipticlatitude_F_stat, eclipticlongitude_F_stat =  search1.f_statistic(N_frequency,5)
            # ind = np.unravel_index(np.argmax(F_stat, axis=None), F_stat.shape)
            # if ind[0]>1:
            #     lower_index = ind[0]-2
            # else:
            #     lower_index = ind[0]
            # if ind[0]<N_frequency-2:
            #     upper_index = ind[0]+2
            # else:
            #     upper_index = ind[0]
            # print(frequencies_F_stat[ind[0]])
            # search1.reduced_frequency_boundaries = [frequencies_F_stat[lower_index],frequencies_F_stat[upper_index]]
            start = time.time()

            # print('SNR ',np.round(search1.SNR([search1.pGB])))
            # print('SNR2', np.round(search1.loglikelihood([search1.pGB])))
            # print('SNR2', np.round(search1.loglikelihoodsdf([search1.pGB])))
            # print('SNR', np.round(search1.SNR([search1.pGB]),3))
            # print('SNRflat', np.round(search1.loglikelihoodflat([search1.pGB])))
            # search1.plot()#pGBadded=pGBadded5)
            # print(pGBadded7["FrequencyDerivative"] * self.Tobs)
            # print('smear f', 300*pGBadded7["Frequency"] * 10**3 / 10**9)
            # print(search1.reduced_frequency_boundaries)
            if ind <= len(initial_guess):
                search_repetitions = 3
            else:
                search_repetitions = 3
            for i in range(search_repetitions):
                # if i > 0:
                #     search1.recombination = 0.75
                #     maxpGBsearch_new, energies =  search1.differential_evolution_search(search1.boundaries['Frequency'], initial_guess=maxpGBsearch[0])
                # else:
                if self.strategy == 'DE':
                    if ind <= len(initial_guess) and i == 0:
                        maxpGBsearch_new, number_of_evaluations =  search1.differential_evolution_search(search1.boundaries['Frequency'], initial_guess = [initial_guess[ind-1]])
                    else:
                        maxpGBsearch_new, number_of_evaluations =  search1.differential_evolution_search(search1.boundaries['Frequency'])
                if self.strategy == 'CD':
                    maxpGBsearch_new, number_of_evaluations =  search1.searchCD()

                found_sources_all.append(maxpGBsearch_new)
                number_of_evaluations_all.append(number_of_evaluations)
                new_SNR = search1.SNR(maxpGBsearch_new[0])
                print('SNR of found signal', np.round(new_SNR,3))
                print('which signal per window', ind,'and repetition:', i)
                if i == 0:
                    current_SNR = deepcopy(new_SNR)
                    maxpGBsearch = deepcopy(maxpGBsearch_new)
                if new_SNR >= current_SNR:
                    current_SNR = deepcopy(new_SNR)
                    maxpGBsearch = deepcopy(maxpGBsearch_new)
                print('current SNR', current_SNR)
                found_sources_all[-1] = maxpGBsearch_new

            print('to optimize Amplitude', maxpGBsearch[0])
            for j in range(len(maxpGBsearch[0])):
                A_optimized = search1.calculate_Amplitude([maxpGBsearch[0][j]])
                # A_optimized.values = -1
                if A_optimized.values > 0:
                    maxpGBsearch[0][j]['Amplitude'] *= A_optimized.values
                else:
                    for i in range(30):
                        print('pGB error!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                    print('Amplitude optimization failed with parameters:', maxpGBsearch[0][j])
                    print('switch to optimize with scipy minimize trust-constr')
                    maxpGBsearch[0][j] = search1.optimizeA([[maxpGBsearch[0][j]]])[0]
                    print('Optimized with parameters:', maxpGBsearch[0][j])
                print('loglikelihood optimized amplitude',search1.loglikelihood([maxpGBsearch[0][j]]))
            print('in range', maxpGBsearch[0][0]['Frequency'] > lower_frequency and maxpGBsearch[0][0]['Frequency'] < upper_frequency)
            # new_SNR = search1.SNR(maxpGBsearch[0])

            current_loglikelihood_ratio = search1.loglikelihood(maxpGBsearch[0])
            print('current loglikelihood ratio', current_loglikelihood_ratio)


            if current_loglikelihood_ratio < loglikelihood_ratio_threshold:
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
            if len(found_sources_in) > 0:
                tdi_fs_subtracted = deepcopy(self.tdi_fs)
                for i in range(len(found_sources_out)):
                    Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_out[i], oversample=4, simulator="synthlisa")
                    source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
                    index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
                    index_high = index_low+len(Xs_subtracted)
                    for k in ["X", "Y", "Z"]:
                        tdi_fs_subtracted[k].data[index_low:index_high] = tdi_fs_subtracted[k].data[index_low:index_high] - source_subtracted[k].data

                search_out_subtracted = Search(tdi_fs_subtracted,self.Tobs, lower_frequency, upper_frequency)

                total_boundaries = deepcopy(search_out_subtracted.boundaries)
                # amplitudes = []
                # for i in range(len(found_sources_in)):
                #     amplitudes.append(found_sources_in[i]['Amplitude'])
                # total_boundaries['Amplitude'] = [np.min(amplitudes),np.max(amplitudes)]
                # amplitudes_length = search1.boundaries['Amplitude'][1] - search1.boundaries['Amplitude'][0]
                # total_boundaries['Amplitude'] = [np.log10(total_boundaries['Amplitude'][0]), np.log10(total_boundaries['Amplitude'][1])]
                # total_boundaries['Amplitude'] = [total_boundaries['Amplitude'][0] - amplitudes_length/10,total_boundaries['Amplitude'][1] + amplitudes_length/10,]
                
                start = time.time()
                found_sources_in_opt = search_out_subtracted.optimize([found_sources_in], boundaries= total_boundaries)
                print('global optimization time', time.time()-start)

                #### after optimization a signal inside window could lay outside. Therefore new selection is required
                if search_out_subtracted.loglikelihood(found_sources_in_opt) > search_out_subtracted.loglikelihood(found_sources_in):
                    found_sources_in = []
                    for i in range(len(found_sources_in_opt)):
                        if found_sources_in_opt[i]['Frequency'] > lower_frequency and found_sources_in_opt[i]['Frequency'] < upper_frequency:
                            found_sources_in.append(found_sources_in_opt[i])
                        else:
                            found_sources_out.append(found_sources_in_opt[i])
                else:
                    for i in range(30):
                        print('optimization failed: ', 'new loglikelihood', search_out_subtracted.loglikelihood(found_sources_in_opt), 'old loglikelihood', search_out_subtracted.loglikelihood(found_sources_in))


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
        return found_sources, found_sources_all, number_of_evaluations_all, found_sources_in, [lower_frequency, upper_frequency]

class Global_optimizer():
    def __init__(self,tdi_fs, Tobs):
        self.tdi_fs = tdi_fs
        self.Tobs = Tobs

    def optimize(self, lower_frequency, upper_frequency, found_sources):
        # create two sets of found sources. found_sources_in with signals inside the boundary and founce_sources_out with outside sources
        found_sources_in = []
        found_sources_out = []
        for i in range(len(found_sources)):
            if found_sources[i]['Frequency'] > lower_frequency and found_sources[i]['Frequency'] < upper_frequency:
                found_sources_in.append(found_sources[i])
            else:
                found_sources_out.append(found_sources[i])

        #global optimization
        if len(found_sources_in) > 0:
            tdi_fs_subtracted = deepcopy(self.tdi_fs)
            search_out_subtracted = Search(tdi_fs_subtracted,self.Tobs, lower_frequency, upper_frequency)

            total_boundaries = deepcopy(search_out_subtracted.boundaries)
            start = time.time()
            start_loglikelihood = search_out_subtracted.loglikelihood(found_sources_in)
            found_sources_in_new = search_out_subtracted.optimize([found_sources_in], boundaries= total_boundaries)
            optimized_loglikelihood = search_out_subtracted.loglikelihood(found_sources_in_new)
            if optimized_loglikelihood > start_loglikelihood:
                found_sources_in = found_sources_in_new
            print('global optimization time', np.round(time.time()-start), 'initial loglikelihood', np.round(start_loglikelihood,5), 'optimized_loglikelihood', np.round(optimized_loglikelihood,5), 'difference loglikelihood', np.round(optimized_loglikelihood-start_loglikelihood,5), 'frequency', lower_frequency )

            found_sources = found_sources_in + found_sources_out

        return found_sources,[],[], found_sources_in, [lower_frequency, upper_frequency]

def tdi_subtraction(tdi_fs,found_sources_mp_subtract, frequencies_search):

    # found_sources_mp_best = []
    # for i in range(len(found_sources_mp_subtract)):
    #     found_sources_mp_best.append(found_sources_mp_subtract[i][0])

    # frequencies_search = np.asarray(frequencies_search)
    # found_sources_to_subtract = []
    # for i in range(len(found_sources_mp_best)):
    #     found_sources_to_subtract.append([])
    #     for j in range(len(found_sources_mp_best[i])):        
    #         # find closest frequency window
    #         frequency_window_index = np.searchsorted(frequencies_search[:,0], found_sources_mp_best[i][j]['Frequency'])-1
    #         if frequency_window_index < 0:
    #             found_sources_to_subtract[i].append(found_sources_mp_best[i][j])
    #         elif found_sources_mp_best[i][j]['Frequency'] > frequencies_search[frequency_window_index][1]:
    #             found_sources_to_subtract[i].append(found_sources_mp_best[i][j])

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

# try:
#     cat = np.load(SAVEPATH+'cat_sorted.npy', allow_pickle = True)
#     print('cat sorted loaded')
# except:
#     # get the source parameters
#     # Radler
#     if Radler:
#         names = np.array(fid['H5LISA/GWSources/GalBinaries']) # radler
#         params = [fid['H5LISA/GWSources/GalBinaries'][k] for k in names]
#         reduced_names = []
#         i = 0
#         for p in params:
#             i += 1
#             if p.shape:
#                 reduced_names.append(names[i-1])
#         params = [np.array(p) for p in params if p.shape]
#         names = reduced_names
#     # Sangria
#     else:
#         names_dgb = fid["sky/dgb/cat"].dtype.names # Sangria
#         params_dgb = [np.array(fid["sky/dgb/cat"][k]).squeeze() for k in names_dgb]
#         names_igb = fid["sky/igb/cat"].dtype.names # Sangria
#         params_igb = [np.array(fid["sky/igb/cat"][k]).squeeze() for k in names_igb]
#         names_vgb = fid["sky/vgb/cat"].dtype.names # Sangria
#         params_vgb = [np.array(fid["sky/vgb/cat"][k]).squeeze() for k in names_vgb]

#     cat = np.rec.fromarrays(params_dgb, names=list(names_dgb))
#     indexes = np.argsort(cat['Frequency'])
#     cat = cat[indexes]
#     np.save(SAVEPATH+'cat_sorted.npy',cat)

# LDC1-4 #####################################
frequencies = []
frequencies_even = []
frequencies_odd = []
# search_range = [0.00398, 0.0041]
# search_range = [0.0039885, 0.0040205]
# search_range = [0.0039935, 0.0039965]
f_Nyquist = 1/dt/2
search_range = [0.0003, f_Nyquist]
search_range = [0.0001, 0.11]
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

frequencies_half_shifted = []
for i in range(len(frequencies)-1):
    frequencies_half_shifted.append([(frequencies[i][1]-frequencies[i][0])/2 +frequencies[i][0],(frequencies[i+1][1]-frequencies[i+1][0])/2 +frequencies[i+1][0]])
# frequencies = frequencies_half_shifted #### if shifted
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

frequencies_search = np.asarray(frequencies)
# figure = plt.figure()
# plt.loglog(frequencies_search[:,1],counts, '.')
# plt.xlabel('Frequency [Hz]')
# plt.ylabel('Number of signals')
# plt.show()

figure = plt.figure()
plt.loglog(frequencies_search[:,0],frequencies_search[:,1]-frequencies_search[:,0],  linewidth= 4, label= '$B$')
plt.loglog(frequencies_search[:,0],frequency_derivative(frequencies_search[:,0],2)*Tobs, label= '$B_{F}$')
plt.loglog(frequencies_search[:,0],frequencies_search[:,0]*3* 10**-4, label= '$3 \cdot B_{O}$')
plt.loglog(frequencies_search[:,0],np.ones(len(frequencies_search[:,0]))*4*32*10**-9*2, label= '$2 \cdot B_{C}$')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Frequency window witdh [Hz]')
plt.xlim(search_range[0],0.1)
plt.ylim(bottom=(frequencies_search[0,1]-frequencies_search[0,0])/10**1)
plt.legend()
plt.show()
plt.savefig(SAVEPATH+'bandwidth.png')


save_name = 'Radler_1_even10'
# for i in range(65):
frequencies_search = frequencies_even
frequencies_search_full = deepcopy(frequencies_search)
# batch_index = int(sys.argv[1])
batch_index = 0
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
batch_size = 64
start_index = batch_size*batch_index
print('batch',batch_index, start_index)
frequencies_search = frequencies_search[start_index:start_index+batch_size]
# frequencies_search = frequencies_search[int(6757/2):]

# print(i, frequencies_search[0])
### highest + padding has to be less than f Nyqist
while frequencies_search[-1][1] + (frequencies_search[-1][1] - frequencies_search[-1][0])/2 > f_Nyquist:
    frequencies_search = frequencies_search[:-1]
# frequencies_search = frequencies_search[70:80]
# frequencies_search = frequencies_search[25:]

change_to_optimized_values = False
if change_to_optimized_values:
    ####### change to optimized signals with neighboors removed and inside singals only
    # found_sources_mp_o = np.load(SAVEPATH+'found_sourcesLDC1-4_2_years_full.npy', allow_pickle = True)
    # found_sources_mp_o = np.load(SAVEPATH+'found_sourcesLDC1-4_2_even_optimized.npy', allow_pickle = True)
    # found_sources_mp_o = np.load(SAVEPATH+'found_sourcesLDC1-4_2_optimized.npy', allow_pickle = True)
    # found_sources_mp_o = np.load(SAVEPATH+'found_sourcesSangria_1_full.npy', allow_pickle = True)
    found_sources_mp_o = np.load(SAVEPATH+'found_sourcesSangria_1_full_opt2_even.npy', allow_pickle = True)
    # found_sources_mp_o = np.load(SAVEPATH+'found_sourcesSangria_1_even.npy', allow_pickle = True)
    found_sources_mp_o = np.load(SAVEPATH+'found_sources' +save_name+'.npy', allow_pickle = True)

    frequencies_mp = []
    for i in range(len(found_sources_mp_o)):
        frequencies_mp.append(found_sources_mp_o[i][4])
    # found_sources_mp_o = np.load(SAVEPATH+'/found_sources387654to408404LDC1-4_2year_even10noise_matrix.npy', allow_pickle = True)
    # found_sources_mp = np.load(SAVEPATH+'found_sources387812to408573LDC1-4_v1final_optimized.npy', allow_pickle = True) #odd
    # found_sources_mp = np.load(SAVEPATH+'found_sources387654to408404LDC1-4_v1final_optimized.npy', allow_pickle = True) #even
    # found_sources_mp = np.load(SAVEPATH+'found_sources387812to408573LDC1-4_2_oddfinal_optimized.npy', allow_pickle = True) #odd
    # found_sources_mp = np.load(SAVEPATH+'found_sources387654to408404LDC1-4_2_evenfinal_optimized.npy', allow_pickle = True) #even
    # found_sources_mp = np.load(SAVEPATH+'found_sourcesLDC1-4_2_even_optimized.npy', allow_pickle = True) #even
    # found_sources_mp = np.load(SAVEPATH+'found_sourcesLDC1-4_2_optimized_odd_second.npy', allow_pickle = True) #odd
    # found_sources_mp = np.load(SAVEPATH+'found_sourcesSangria_1_odd_opt2.npy', allow_pickle = True) #odd
    found_sources_mp = np.load(SAVEPATH+'found_sourcesSangria_1_full_opt2high_frequency_opt2_odd.npy', allow_pickle = True) #odd
    
    found_sources_mp_2 = deepcopy(found_sources_mp_o)
    # extra_out_signals = []
    # for i in range(len(found_sources_mp)):
    #     start_index = np.searchsorted(frequencies_mp, found_sources_mp[i][4][0])
    #     found_in = []
    #     for j in range(len(found_sources_mp_o[start_index][0])):
    #         if found_sources_mp_o[start_index][0][j]['Frequency'] > found_sources_mp_o[start_index][4][0] and found_sources_mp_o[start_index][0][j]['Frequency'] < found_sources_mp_o[start_index][4][1]:
    #             found_in.append(found_sources_mp_o[start_index][0][j])
    #     if len(found_in) != len(found_sources_mp[i][3]):
    #         print(found_sources_mp_o[start_index][4][0])
    #         print(len(found_sources_mp_o[start_index][3]), len(found_sources_mp[i][3]))
    #         # if found_sources_mp_o[start_index][4][0] > 0.01630 and found_sources_mp_o[start_index][4][0] < 0.01631:
    #         print(start_index)
            # search1 = Search(tdi_fs,Tobs, found_sources_mp_o[start_index][4][0], found_sources_mp_o[start_index][4][1])
            # if len(found_sources_mp_o[start_index][3]) > 0:
            #     print(search1.SNR(found_sources_mp_o[start_index][3]))
            # print(search1.SNR(found_sources_mp_o[start_index][0]))
            # # create two sets of found sources. found_sources_in with signals inside the boundary and founce_sources_out with outside sources
            # found_sources_in = []
            # found_sources_out = []
            # for k in range(len(found_sources_mp_o[start_index][0])):
            #     if found_sources_mp_o[start_index][0][k]['Frequency'] > found_sources_mp_o[start_index][4][0] and found_sources_mp_o[start_index][0][k]['Frequency'] < found_sources_mp_o[start_index][4][1]:
            #         found_sources_in.append(found_sources_mp_o[start_index][0][k])
            #     else:
            #         found_sources_out.append(found_sources_mp_o[start_index][0][k])
            # for k in range(len(found_sources_out)):
            #     print(search1.SNR([found_sources_out[k]]),search1.loglikelihood([found_sources_out[k]]))
            #     if search1.loglikelihood([found_sources_out[k]]) > 50:
            #         extra_out_signals.append([found_sources_out[k], search1.loglikelihood([found_sources_out[k]])])
            #         if found_sources_out[k]['Frequency'] < found_sources_mp_o[start_index][4][0]:
            #             found_sources_mp_2[start_index-1][0].append(found_sources_out[k])
            #             found_sources_mp_2[start_index-1][3].append(found_sources_out[k])
            #         if found_sources_out[k]['Frequency'] > found_sources_mp_o[start_index][4][1]:
            #             found_sources_mp_2[start_index+1][3].append(found_sources_out[k])
            #             found_sources_mp_2[start_index+1][3].append(found_sources_out[k])

        # if len(found_sources_mp_2[start_index][3]) != len(found_sources_mp[i][3]):
        #     print(start_index, frequencies_mp[start_index],found_sources_mp[i][4][0], len(found_sources_mp_2[start_index][3]), len(found_sources_mp[i][3]), len(found_sources_mp_2[start_index][0]), len(found_sources_mp[i][0]))
        #     print(found_sources_mp_2[start_index][3], found_sources_mp[i][3])
        #     print(found_sources_mp_2[start_index][3], found_sources_mp[i][3])
            # for j in range(len(found_sources_mp[i][0])):
            #     if found_sources_mp[i][0][j] < found_sources_mp_2[start_index][4][1]:


    # found_sources_in_flat = []
    # found_sources_in_flat_frequency = []
    # for i in range(len(found_sources_mp)):
    #     for j in range(len(found_sources_mp[i][3])):
    #         found_sources_in_flat.append(found_sources_mp[i][3][j])
    #         found_sources_in_flat_frequency.append(found_sources_in_flat[-1]['Frequency'])
    # found_sources_in_flat_frequency = np.asarray(found_sources_in_flat_frequency)
    # found_sources_in_flat = np.asarray(found_sources_in_flat)
    # indexes_in = np.argsort(found_sources_in_flat_frequency)
    # found_sources_in_flat_frequency = found_sources_in_flat_frequency[indexes_in]
    # found_sources_in_flat = found_sources_in_flat[indexes_in]

    found_sources_new_flat = []
    for i in range(len(found_sources_mp)):
        for j in range(len(found_sources_mp[i][3])):
            found_sources_new_flat.append(found_sources_mp[i][3][j])
    found_sources_new_flat = np.asarray(found_sources_new_flat)
    found_sources_new_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_new_flat]) for attribute in found_sources_new_flat[0].keys()}
    found_sources_new_flat_df = pd.DataFrame(found_sources_new_flat_array)
    found_sources_new_flat_df = found_sources_new_flat_df.sort_values('Frequency')

    found_sources_flat = []
    for i in range(len(found_sources_mp_o)):
        for j in range(len(found_sources_mp_o[i][3])):
            found_sources_flat.append(found_sources_mp_o[i][3][j])
    found_sources_flat = np.asarray(found_sources_flat)
    found_sources_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_flat]) for attribute in found_sources_flat[0].keys()}
    found_sources_flat_df = pd.DataFrame(found_sources_flat_array)
    found_sources_flat_df = found_sources_flat_df.sort_values('Frequency')
    for i in range(len(found_sources_mp)):
        found_sources_flat_df = found_sources_flat_df[(found_sources_flat_df['Frequency']< found_sources_mp[i][4][0]) | (found_sources_flat_df['Frequency']> found_sources_mp[i][4][1])]
    found_sources_combined_flat_df = found_sources_flat_df.append(found_sources_new_flat_df, ignore_index=True)
    found_sources_combined_flat_df = found_sources_combined_flat_df.sort_values('Frequency')

    for i in range(len(found_sources_mp_2)):
        found_sources_mp_2[i][3] = found_sources_combined_flat_df[(found_sources_combined_flat_df['Frequency']> frequencies_mp[i][0]) & (found_sources_combined_flat_df['Frequency']< frequencies_mp[i][1])].to_dict(orient='records')

    np.save(SAVEPATH+'found_sourcesSangria_1_full_opt2high_frequency_opt2.npy', found_sources_mp_2)

# target_frequencies = []
# index_low = np.searchsorted(cat_sorted['Frequency'], search_range[0])
# for i in range(10):
#     target_frequencies.append(cat_sorted[-10+i-1]['Frequency'])
#     # target_frequencies.append(cat_sorted[index_low+i*100]['Frequency'])
# # target_frequencies = cat_sorted[-17:-1]['Frequency']
# frequencies_search = []
# for i in range(len(target_frequencies)):
#     current_frequency = target_frequencies[i]
#     f_smear = current_frequency *3* 10**-4
#     f_deviation = frequency_derivative(current_frequency,M_chirp_upper_boundary)*Tobs
#     print(current_frequency,frequency_derivative(current_frequency,M_chirp_upper_boundary))
#     window_length = f_smear + f_deviation
#     window_length += 4*32*10**-9*2
#     window_shift = ((np.random.random(1)-0.5)*window_length*0.5)[0]
#     frequencies_search.append([target_frequencies[i]-window_length/2+window_shift,target_frequencies[i]+window_length/2+window_shift])

search_range = [frequencies_search[0][0],frequencies_search[-1][1]]
# search_range = [1619472*10**-8,2689639*10**-8]
print('search range '+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))))

# i = 0
# lower_frequency = frequencies_search[i][0]
# upper_frequency = frequencies_search[i][1]
# search1 = Search(tdi_fs,Tobs, lower_frequency, upper_frequency)
# search1.plot(second_data= tdi_fs_dgb)


do_subtract = True
if do_subtract:
    start = time.time()
    # save_name_previous = 'found_sources397769to400619LDC1-4_4mHz_half_year_even3'
    # save_name_previous = 'found_sources397919to400770LDC1-4_4mHz_half_year_odd'
    # save_name_previous = 'found_sources397769to400619LDC1-4_4mHz_2_year_initial_half_even3'
    # save_name_previous = 'found_sources397956to401074LDC1-4_4mHz_2_year_initial_half_even_SNR10'
    # save_name_previous = 'found_sources397793to400909LDC1-4_4mHz_2_year_initial_half_odd_SNR10'
    # save_name_previous = 'found_sources397956to401074LDC1-4_4mHz_loglikelihood_ratio_threshold_even3'
    # save_name_previous = 'found_sources397793to400909LDC1-4_4mHz_loglikelihood_ratio_threshold_odd'
    # save_name_previous = 'found_sources40709to41430LDC1-4_04mHz_loglikelihood_ratio_threshold_even3'
    # save_name_previous = 'found_sources40747to41468LDC1-4_04mHz_loglikelihood_ratio_threshold_odd'
    # save_name_previous = 'LDC1-4 odd'
    # save_name_previous = 'found_sourcesLDC1-4_half_even'
    # save_name_previous = 'found_sourcesLDC1-4_half_even_T'
    # save_name_previous = 'found_sourcesLDC1-4_half_odd'
    # save_name_previous = 'found_sourcesLDC1-4_2_even10'
    # save_name_previous = 'found_sourcesLDC1-4_2_odd'
    # save_name_previous = 'found_sourcesLDC1-4_2_years_full'
    # save_name_previous = 'found_sourcesLDC1-4_2_odd_optimized'
    # save_name_previous = 'found_sourcesLDC1-4_2_even_optimized'
    # save_name_previous = 'found_sourcesLDC1-4_2_optimized'
    # Sangria
    # save_name_previous = 'found_sourcesSangria_half_even3'
    # save_name_previous = 'found_sourcesSangria_half_odd'
    # save_name_previous = 'found_sourcesSangria_1_full'
    # save_name_previous = 'found_signals_1_even/found_sources1630803to2001429Sangria_1_even'
    # save_name_previous = 'found_sourcesSangria_odd'
    # save_name_previous = 'found_sourcesRadler_1_even3'
    save_name_previous = 'found_sourcesRadler_1_odd'
    # save_name_previous = 'found_sourcesSangria_1_even_opt1'
    # save_name_previous = 'found_sources'+save_name
    found_sources_mp_subtract = np.load(SAVEPATH+save_name_previous+'.npy', allow_pickle = True)

    found_sources_flat = []
    for i in range(len(found_sources_mp_subtract)):
        for j in range(len(found_sources_mp_subtract[i][3])):
            found_sources_flat.append(found_sources_mp_subtract[i][3][j])
    found_sources_flat = np.asarray(found_sources_flat)
    found_sources_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_flat]) for attribute in found_sources_flat[0].keys()}
    found_sources_flat_df = pd.DataFrame(found_sources_flat_array)
    found_sources_flat_df = found_sources_flat_df.sort_values('Frequency')
    for i in range(len(frequencies_search_full)):
        found_sources_flat_df = found_sources_flat_df[(found_sources_flat_df['Frequency']< frequencies_search_full[i][0]) | (found_sources_flat_df['Frequency']> frequencies_search_full[i][1])]
    found_sources_flat_df = found_sources_flat_df.sort_values('Frequency')
    found_sources_out_flat = found_sources_flat_df.to_dict(orient='records')
    tdi_fs_subtracted = tdi_subtraction(tdi_fs,found_sources_out_flat, frequencies_search_full)

    print('subtraction time', time.time()-start)
    plot_subtraction = True
    if plot_subtraction:
        i = 1
        lower_frequency = frequencies_search[i][0]
        upper_frequency = frequencies_search[i][1]
        search1 = Search(tdi_fs,Tobs, lower_frequency, upper_frequency)
        source = [{'Amplitude': 4.500916389929765e-20, 'EclipticLatitude': 0.8528320149942861, 'EclipticLongitude': -0.9418744765040503, 'Frequency': frequencies_search[i][0]+(frequencies_search[i][1]-frequencies_search[i][0])/2, 'FrequencyDerivative': 2.1688352300259018e-22, 'Inclination': 1.343872907043714, 'InitialPhase': 3.583816929574315, 'Polarization': 2.69557290704741}]
        # search1.plot(second_data= tdi_fs_subtracted, found_sources_in=found_sources_out_flat)
        search1.plot(found_sources_in=source)
        # search1.plot(second_data= tdi_fs_subtracted, found_sources_in=found_sources_mp_o[start_index][0])
        
    tdi_fs = deepcopy(tdi_fs_subtracted)

# i = 1718
# search1 = Search(tdi_fs,Tobs, found_sources_mp_subtract[i][4][0], found_sources_mp_subtract[i][4][1])
# search1.plot(second_data= tdi_fs_subtracted, found_sources_in=found_sources_mp_subtract[i][0])
# search1.plot(second_data= tdi_fs_subtracted, found_sources_in=found_sources_mp_subtract[i][3])

# save_name_previous = 'found_sourcesSangria_1_even'
# found_sources_mp_subtract = np.load(SAVEPATH+save_name_previous+'.npy', allow_pickle = True)
# print(search1.SNR(found_sources_mp_subtract[17*64+8][3]))
# print(search1.SNR(found_sources_mp[17*64+8][3]))
# search1.plot(second_data= tdi_fs_subtracted, found_sources_in=found_sources_in_flat)

# do_subtract = True
# if do_subtract:
#     start = time.time()
#     # save_name_previous = 'found_sources397769to400619LDC1-4_4mHz_half_year_even3'
#     # save_name_previous = 'found_sources397919to400770LDC1-4_4mHz_half_year_odd'
#     # save_name_previous = 'found_sources397956to401074LDC1-4_4mHz_2_year_initial_half_even3'
#     # save_name_previous = 'found_sources397769to400619LDC1-4_4mHz_2_year_initial_half_even3'
#     # save_name_previous = 'found_sources397956to401074LDC1-4_4mHz_2_year_initial_half_even_SNR10'
#     save_name_previous = 'found_sources397793to400909LDC1-4_4mHz_2_year_initial_half_odd_SNR10'
#     # save_name_previous = 'found_sourcesLDC1-4_half_odd'
#     found_sources_mp_subtract = np.load(SAVEPATH+'/'+save_name_previous+'.npy', allow_pickle = True)

#     found_sources_in_flat = []
#     found_sources_in_flat_frequency = []
#     for i in range(len(found_sources_mp_subtract)):
#         for j in range(len(found_sources_mp_subtract[i][3])):
#             found_sources_in_flat.append(found_sources_mp_subtract[i][3][j])
#             found_sources_in_flat_frequency.append(found_sources_in_flat[-1]['Frequency'])
#     found_sources_in_flat_frequency = np.asarray(found_sources_in_flat_frequency)
#     found_sources_in_flat = np.asarray(found_sources_in_flat)
#     indexes_in = np.argsort(found_sources_in_flat_frequency)
#     found_sources_in_flat_frequency = found_sources_in_flat_frequency[indexes_in]
#     found_sources_in_flat = found_sources_in_flat[indexes_in]

#     found_sources_to_subtract = []
#     for i in range(len(frequencies_search)):
#         lower_index = np.searchsorted(found_sources_in_flat_frequency,frequencies_search[i][0])
#         higher_index = np.searchsorted(found_sources_in_flat_frequency,frequencies_search[i][1])
#         found_sources_to_subtract.append(found_sources_in_flat[lower_index:higher_index])

#     #subtract the found sources from original
#     tdi_fs_subtracted = deepcopy(tdi_fs)
#     for i in range(len(found_sources_to_subtract)):
#         for j in range(len(found_sources_to_subtract[i])):
#             Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_to_subtract[i][j], oversample=4, simulator="synthlisa")
#             source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
#             index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
#             index_high = index_low+len(Xs_subtracted)
#             for k in ["X", "Y", "Z"]:
#                 tdi_fs_subtracted[k].data[index_low:index_high] -= source_subtracted[k].data
#     print('subtraction time', time.time()-start)
#     i = 5
#     lower_frequency = frequencies_search[i][0]
#     upper_frequency = frequencies_search[i][1]
#     search1 = Search(tdi_fs,Tobs, lower_frequency, upper_frequency)
#     search1.plot()
#     search1 = Search(tdi_fs_subtracted,Tobs, lower_frequency, upper_frequency)
#     search1.plot()
#     tdi_fs = deepcopy(tdi_fs_subtracted)


# length = 100
# for i in range(length):
#     print(i,found_sources_mp[-i][0])
#     print(i,found_sources_mp_even[-i][0])
# found_sources_mp = found_sources_mp_even[-100:]
# frequencies_search = frequencies_even[-100:]

found_sources_sorted = []
use_initial_guess = True
if use_initial_guess:
    # save_name_found_sources_previous = 'found_sources397769to400619LDC1-4_4mHz_half_year_even10'
    # save_name_found_sources_previous = 'found_sources397919to400770LDC1-4_4mHz_half_year_odd'
    # save_name_found_sources_previous = 'found_sources2537595to3305084LDC1-4_4mHz_half_year_even'
    save_name_found_sources_previous = 'found_sourcesLDC1-4_half_even10'
    # save_name_found_sources_previous = 'found_sourcesLDC1-4_half_odd'
    found_sources_mp_subtract = np.load(SAVEPATH+save_name_found_sources_previous+'.npy', allow_pickle = True)

    found_sources_flat = []
    for i in range(len(found_sources_mp_subtract)):
        for j in range(len(found_sources_mp_subtract[i][0])):
            found_sources_flat.append(found_sources_mp_subtract[i][0][j])
    found_sources_flat = np.asarray(found_sources_flat)
    found_sources_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_flat]) for attribute in found_sources_flat[0].keys()}
    found_sources_flat_df = pd.DataFrame(found_sources_flat_array)
    found_sources_flat_df = found_sources_flat_df.sort_values('Frequency')
    found_sources_sorted = found_sources_flat_df.to_dict(orient='records')

    found_sources_loaded = []
    found_sources_loaded.append(np.load(SAVEPATH+save_name_found_sources_previous+'.npy', allow_pickle = True))
    found_sources_previous = []
    for i in range(len(found_sources_loaded)):
        for j in range(len(found_sources_loaded[i])):
            for k in range(len(found_sources_loaded[i][j][0])):
                found_sources_previous.append(found_sources_loaded[i][j][0][k])

    found_sources_array = np.zeros((len(found_sources_previous),len(parameters)))
    for i in range(len(found_sources_previous)):
        for j, parameter in enumerate(parameters):
            found_sources_array[i,j] = found_sources_previous[i][parameter]
    found_sources_tuples = []
    for i in range(len(found_sources_array)):
        found_sources_tuples.append(tuple(found_sources_array[i]))
    found_sources_sorted = np.array(found_sources_tuples, dtype=[('Amplitude', '<f8'), ('EclipticLatitude', '<f8'), ('EclipticLongitude', '<f8'), ('Frequency', '<f8'), ('FrequencyDerivative', '<f8'), ('Inclination', '<f8'), ('InitialPhase', '<f8'), ('Polarization', '<f8')])
    indexes = np.argsort(found_sources_sorted['Frequency'])
    found_sources_sorted = found_sources_sorted[indexes]

# i = 4
# search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
# search1.plot(pGB_injected=pGB_injected[i][:10])
# pGBadded['Frequency'] = frequencies_search[0][0]
# start = time.time()
# for i in range(10**3):
#     search1.SNR([pGBadded])
# print('full time', time.time()-start)
# start = time.time()
# for i in range(10**3):
#     Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBadded, oversample=4, simulator="synthlisa")
# print('full time', time.time()-start)

# cat_negative_frequency_derivative = cat_sorted['Frequency'][cat_sorted['FrequencyDerivative']<0]
# np.max(cat_negative_frequency_derivative)


# from sources import *
# MLP = MLP_search(tdi_fs, Tobs, signals_per_window = 1, found_sources_previous = found_sources_sorted,  strategy = 'DE')
# found_sources_mp = MLP.search(frequencies_search[1][0], frequencies_search[1][1])
# found_sources_mp = [found_sources_mp]
# frequencies_search = [frequencies_search[7]]

# search1 = Search(tdi_fs,Tobs, lower_frequency, upper_frequency)
# search1.plot(found_sources_in= np.concatenate(found_sources_mp[index:index+2])[0], pGB_injected=np.concatenate(pGB_injected[index:index+2]))

# search1.SNR(found_sources_mp[index][0])
# search1.SNR(pGB_injected[index])

do_search = True
if do_search:
    MLP = MLP_search(tdi_fs, Tobs, signals_per_window = 10, found_sources_previous = found_sources_sorted, strategy = 'DE')
    start = time.time()
    cpu_cores = 16
    pool = mp.Pool(cpu_cores)
    found_sources_mp = pool.starmap(MLP.search, frequencies_search)
    pool.close()
    pool.join()
    print('time to search ', len(frequencies_search), 'windows: ', time.time()-start)
    np.save(SAVEPATH+'found_signals/found_sources'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', found_sources_mp)

final_optimization = False
if final_optimization:
    # found_sources_mp = np.load(SAVEPATH+'found_sources'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'noise_matrix.npy', allow_pickle = True)
    # found_sources_mp = np.load(SAVEPATH+'found_sourcesLDC1-4_2_years_full.npy', allow_pickle = True)
    # found_sources_mp = np.load(SAVEPATH+'found_sources387812to408573LDC1-4_2year_oddnoise_matrix.npy', allow_pickle = True)
    # found_sources_mp = np.load(SAVEPATH+'found_sourcesLDC1-4_2_odd_optimized.npy', allow_pickle = True)
    # found_sources_mp = np.load(SAVEPATH+'found_sourcesLDC1-4_2_even_optimized.npy', allow_pickle = True)
    # found_sources_mp = np.load(SAVEPATH+'found_sourcesLDC1-4_2_optimized.npy', allow_pickle = True)
    # found_sources_mp = np.load(SAVEPATH+'found_signals/found_sources'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', allow_pickle = True)
    found_sources_mp = np.load(SAVEPATH+'found_sources'+save_name+'.npy', allow_pickle = True)
    optimizer = Global_optimizer(tdi_fs, Tobs)
    input = []

    # found_sources_in_flat = []
    # found_sources_in_flat_frequency = []
    # for i in range(len(found_sources_mp)):
    #     for j in range(len(found_sources_mp[i][3])):
    #         found_sources_in_flat.append(found_sources_mp[i][3][j])
    #         found_sources_in_flat_frequency.append(found_sources_in_flat[-1]['Frequency'])
    # found_sources_in_flat_frequency = np.asarray(found_sources_in_flat_frequency)
    # found_sources_in_flat = np.asarray(found_sources_in_flat)
    # indexes_in = np.argsort(found_sources_in_flat_frequency)
    # found_sources_in_flat_frequency = found_sources_in_flat_frequency[indexes_in]
    # found_sources_in_flat = found_sources_in_flat[indexes_in]
    # found_sources_in = [] 
    # ##### error
    # for i in range(len(frequencies_search)):
    #     lower_index = np.searchsorted(found_sources_in_flat_frequency,frequencies_search[i][0])
    #     higher_index = np.searchsorted(found_sources_in_flat_frequency,frequencies_search[i][1])
    #     found_sources_in.append(found_sources_in_flat[lower_index:higher_index])
    #     input.append([frequencies_search[i][0],frequencies_search[i][1],found_sources_in_flat[lower_index:higher_index]])

    found_sources_flat = []
    for i in range(len(found_sources_mp)):
        for j in range(len(found_sources_mp[i][3])):
            found_sources_flat.append(found_sources_mp[i][3][j])
    found_sources_flat = np.asarray(found_sources_flat)
    found_sources_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_flat]) for attribute in found_sources_flat[0].keys()}
    found_sources_flat_df = pd.DataFrame(found_sources_flat_array)
    found_sources_flat_df = found_sources_flat_df.sort_values('Frequency')
    for i in range(len(frequencies_search)):
        input_dict = found_sources_flat_df[(found_sources_flat_df['Frequency'] > frequencies_search[i][0]) & (found_sources_flat_df['Frequency'] < frequencies_search[i][1])].to_dict(orient='records')
        input.append([frequencies_search[i][0],frequencies_search[i][1],input_dict])

    start = time.time()
    cpu_cores = 16
    pool = mp.Pool(cpu_cores)
    found_sources_mp = pool.starmap(optimizer.optimize, input)
    pool.close()
    pool.join()
    print('time to optimize', len(frequencies_search), 'windows: ', time.time()-start)
    np.save(SAVEPATH+'optimized/found_sources'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'_opt1_even.npy', found_sources_mp)

# better_pGB = {'Amplitude': 3.7080776510756e-23, 'EclipticLatitude': -0.0864329194471405, 'EclipticLongitude': -1.5608489415225566, 'Frequency': 0.011097063538503463, 'FrequencyDerivative': 3.795997584356877e-15, 'Inclination': 1.3544536642993756, 'InitialPhase': 3.802341846303522, 'Polarization': 3.0450807858161113}
# search_out_subtracted = Search(tdi_fs,Tobs, input[0][0], input[0][1])
# print(search_out_subtracted.loglikelihood(found_sources_mp[0][0]))
# print(search_out_subtracted.loglikelihood([better_pGB]))

# i = 53
# search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
# search1.plot(found_sources_in=found_sources_mp_o[i][0][:-1])
# search1.plot(found_sources_in=found_sources_mp[i][0][:-1])
# print(search1.SNR(found_sources_mp_o[i][0][:-1]))
# print(search1.SNR(found_sources_mp[i][0][:-1]))
# print(search1.SNR(found_sources_mp_o[i][0]))
# print(search1.SNR(found_sources_mp[i][0]))

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

    Af = (Zs_aligned - Xs_aligned)/np.sqrt(2.0)
    Ef = (Zs_aligned - 2.0*Ys_aligned + Xs_aligned)/np.sqrt(6.0)
    Af_injected = (Zs_injected - Xs_injected)/np.sqrt(2.0)
    Ef_injected = (Zs_injected - 2.0*Ys_injected + Xs_injected)/np.sqrt(6.0)
    SNR2 = np.sum( np.real(Af_injected * np.conjugate(Af.data) + Ef_injected * np.conjugate(Ef.data))/SA )
    hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /SA)
    ss = np.sum((np.absolute(Af_injected.data)**2 + np.absolute(Ef_injected.data)**2) /SA)
    do_plot = False
    if do_plot:
        fig, ax = plt.subplots(nrows=2, sharex=True) 
        ax[0].semilogy(Af.f, np.abs(Af_injected))
        ax[0].semilogy(Af.f, np.abs(Af.data))
        
        ax[1].semilogy(Af.f, np.abs(Ef_injected))
        ax[1].semilogy(Af.f, np.abs(Ef.data))
        plt.show()
        
    SNR = 4.0*Xs.df* hh
    SNR2 = 4.0*Xs.df* SNR2
    SNR3 = SNR2 / (np.sqrt(SNR)*np.sqrt(4.0*Xs.df* ss))
    return SNR3.values

def SNR_match_XYZ(pGB_injected, pGB_found):
    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB_found, oversample=4, simulator="synthlisa")
    Xs_injected, Ys_injected, Zs_injected = GB.get_fd_tdixyz(template=pGB_injected, oversample=4, simulator="synthlisa")
    Xs_aligned = xr.align(Xs_injected, Xs, join='left',fill_value=0)[1]
    Ys_aligned = xr.align(Ys_injected, Ys, join='left',fill_value=0)[1]
    Zs_aligned = xr.align(Zs_injected, Zs, join='left',fill_value=0)[1]
        
    fmin, fmax = float(Xs_injected.f[0]), float(Xs_injected.f[-1] + Xs_injected.attrs["df"])
    freq = np.array(Xs_injected.sel(f=slice(fmin, fmax)).f)
    Nmodel = get_noise_model(noise_model, freq)
    SA = Nmodel.psd(freq=freq, option="A")

    Af = (Zs_aligned - Xs_aligned)/np.sqrt(2.0)
    Ef = (Zs_aligned - 2.0*Ys_aligned + Xs_aligned)/np.sqrt(6.0)
    Af_injected = (Zs_injected - Xs_injected)/np.sqrt(2.0)
    Ef_injected = (Zs_injected - 2.0*Ys_injected + Xs_injected)/np.sqrt(6.0)
    SNR2 = np.sum( np.real(Af_injected * np.conjugate(Af.data) + Ef_injected * np.conjugate(Ef.data))/SA )
    hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /SA)
    ss = np.sum((np.absolute(Af_injected.data)**2 + np.absolute(Ef_injected.data)**2) /SA)
    do_plot = False
    if do_plot:
        fig, ax = plt.subplots(nrows=2, sharex=True) 
        ax[0].semilogy(Af.f, np.abs(Af_injected))
        ax[0].semilogy(Af.f, np.abs(Af.data))
        
        ax[1].semilogy(Af.f, np.abs(Ef_injected))
        ax[1].semilogy(Af.f, np.abs(Ef.data))
        plt.show()
        
    SNR = 4.0*Xs.df* hh
    SNR2 = 4.0*Xs.df* SNR2
    SNR3 = SNR2 / (np.sqrt(SNR)*np.sqrt(4.0*Xs.df* ss))
    return SNR3.values

do_print = False
if do_print:
    # found_sources_mp = np.load(SAVEPATH+'/LDC1-4/Found_signals_half_year_odd/found_sources'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', allow_pickle = True)
    # found_sources_mp = np.load(SAVEPATH+'/found_sources' +save_name+'.npy', allow_pickle = True)
    found_sources_mp = np.load(SAVEPATH+'/found_signals/found_sources'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', allow_pickle = True)
    
    # found_sources_mp = np.load(SAVEPATH+'/found_sources397956to401074LDC1-4_4mHz_loglikelihood_ratio_threshold_even10.npy', allow_pickle = True)
    

    found_sources_mp_best = []
    found_sources_mp_all = []
    for i in range(len(found_sources_mp)):
        found_sources_mp_best.append(found_sources_mp[i][0])
        found_sources_in_window = []
        for j in range(len(found_sources_mp[i][1])):
            found_sources_in_window.append(found_sources_mp[i][1][j][0][0])
        found_sources_mp_all.append(found_sources_in_window)

    found_sources_in_flat = []
    found_sources_in_flat_frequency = []
    for i in range(len(found_sources_mp)):
        for j in range(len(found_sources_mp[i][3])):
            found_sources_in_flat.append(found_sources_mp[i][3][j])
            found_sources_in_flat_frequency.append(found_sources_in_flat[-1]['Frequency'])
    found_sources_in_flat_frequency = np.asarray(found_sources_in_flat_frequency)
    found_sources_in_flat = np.asarray(found_sources_in_flat)
    indexes_in = np.argsort(found_sources_in_flat_frequency)
    found_sources_in_flat_frequency = found_sources_in_flat_frequency[indexes_in]
    found_sources_in_flat = found_sources_in_flat[indexes_in]

    found_sources_in = []
    for i in range(len(frequencies_search)):
        lower_index = np.searchsorted(found_sources_in_flat_frequency,frequencies_search[i][0])
        higher_index = np.searchsorted(found_sources_in_flat_frequency,frequencies_search[i][1])
        found_sources_in.append(found_sources_in_flat[lower_index:higher_index])

    for index in range(10):
        j = 0
        best_first_signal = []
        optimized_first_signal = []
        search1 = Search(tdi_fs,Tobs, frequencies_search[index][0], frequencies_search[index][1])
        for source in [found_sources_mp]:
            j += 1
            SNR_list = []
            for i in range(3):
                print(j)
                print('intrinsic SNR A E T',search1.intrinsic_SNR(source[index][1][i][0]))
                print('intrinsic SNR A E',search1.intrinsic_SNR_old(source[index][1][i][0]))
                print('intrinsic SNR T',search1.intrinsic_SNR_T(source[index][1][i][0]))
                print('SNR A E T',search1.SNR(source[index][1][i][0]))
                print('SNR A E',search1.SNR_AE(source[index][1][i][0]))
                print('SNR XYZ',search1.SNR_XYZ(source[index][1][i][0]))
                SNR_list.append(search1.SNR_noise_matrix(source[index][1][i][0]))
                print('SNR noise matrix',search1.SNR_noise_matrix(source[index][1][i][0]))
                print(search1.SNR(source[index][1][i][0]))
                print(source[index][1][i][0][0])
            print('best index',np.argmax(SNR_list))
            best_first_signal.append(source[index][1][np.argmax(SNR_list)][0])

    found_sources_out = []
    index = 0
    for found_sources in found_sources_mp_best:
        # create two sets of found sources. found_sources_in with signals inside the boundary and founce_sources_out with outside sources
        found_sources_out.append([])
        for i in range(len(found_sources)):
            if found_sources[i]['Frequency'] > search1.lower_frequency and found_sources[i]['Frequency'] < search1.upper_frequency:
                pass
            else:
                found_sources_out[index].append(found_sources[i])
        index += 1

    #### global optimization
    optimized_signals = []
    for i in range(len(found_sources_mp_best)):
        if i != 0:
            continue
        search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
        optimized_signals.append([])
        #global optimization
        tdi_fs_subtracted = deepcopy(tdi_fs)
        for j in range(len(found_sources_out[i])):
            Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_out[i][j], oversample=4, simulator="synthlisa")
            source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
            index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
            index_high = index_low+len(Xs_subtracted)
            for k in ["X", "Y", "Z"]:
                tdi_fs_subtracted[k].data[index_low:index_high] = tdi_fs_subtracted[k].data[index_low:index_high] - source_subtracted[k].data
        search_out_subtracted = Search(tdi_fs_subtracted,Tobs, search1.lower_frequency, search1.upper_frequency)
        lower_frequency = frequencies_search[i][0]
        upper_frequency = frequencies_search[i][1]
        search2 = Search(tdi_fs,Tobs, lower_frequency, upper_frequency)
        search2.plot(found_sources_in=found_sources_out[i])
        search2 = Search(tdi_fs_subtracted,Tobs, lower_frequency, upper_frequency)
        search2.plot(found_sources_in=found_sources_out[i])

        total_boundaries = deepcopy(search1.boundaries)
        start = time.time()
        found_sources_in[i] = search_out_subtracted.optimize([found_sources_in[i]], boundaries= total_boundaries)
        optimized_signals[i].append(found_sources_in[i])
        print('global optimization time', time.time()-start)

    for parameter in parameters:
        print(parameter, (optimized_signals[index][0][parameter]-optimized_signals[index][0][parameter])/optimized_signals[index][0][parameter])
        # print(parameter, (best_first_signal[0][0][parameter]-best_first_signal[1][0][parameter])/best_first_signal[0][0][parameter])

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

    i = 3
    search2 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
    search2.plot(found_sources_in=found_sources_in[i], pGB_injected=pGB_injected[i])

    fig = plt.figure()
    parameter1 = 'EclipticLongitude'
    parameter2 = 'EclipticLatitude'
    plt.scatter(best_first_signal[0][0][parameter1],best_first_signal[0][0][parameter2], color='blue')
    plt.scatter(best_first_signal[1][0][parameter1],best_first_signal[1][0][parameter2], color='red')
    plt.scatter(optimized_first_signal[0][0][parameter1],optimized_first_signal[0][0][parameter2], color='orange')
    plt.scatter(optimized_first_signal[1][0][parameter1],optimized_first_signal[1][0][parameter2], color='green')
    # plt.scatter(pGB_injected[1][0][parameter1]-np.pi*2,pGB_injected[1][0][parameter2], color='black')
    plt.show()

    # list_best_SNR_index = []
    # list_best_unoptimized_parameters = []
    # for i in range(len(found_sources_mp_all)-9):
    #     list_best_SNR_index.append([])
    #     list_best_unoptimized_parameters.append([])
    #     search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
    #     for j in range(int(len(found_sources_mp_all[i])/3)):
    #         list_SNR = []
    #         for k in range(3):
    #             list_SNR.append(search1.SNR([found_sources_mp_all[i][j+k]]))
    #         list_best_SNR_index[-1].append(np.argmax(list_SNR))
    #         list_best_unoptimized_parameters[-1].append(found_sources_mp_all[i][j+list_best_SNR_index[-1][0]])

    # unoptimized_in = []
    # unoptimized_out = []
    # for i in range(len(list_best_unoptimized_parameters)):
    #     if list_best_unoptimized_parameters[i]['Frequency'] < frequencies_search[i][1] and list_best_unoptimized_parameters[i]['Frequency'] > frequencies_search[i][0]:
    #         unoptimized_in.append(list_best_unoptimized_parameters[i])
    #     else:
    #         unoptimized_out.append(list_best_unoptimized_parameters[i])


    #extend found sources
    # found_sources_mp2 = []
    # for i in range(len(found_sources_mp)):
    #     found_sources_mp2.append([])
    #     list_found_sources = list(found_sources_mp[i])
    #     list_found_sources.append(found_sources_in[i])
    #     found_sources_mp2[i] = list_found_sources
    # found_sources_mp = found_sources_mp2


    # index_low = 100000
    # index_high = 102000
    # fig = plt.figure()
    # plt.plot(tdi_fs_subtracted['X'].f[index_low:index_high],tdi_fs_subtracted['X'][index_low:index_high].values)
    # plt.plot(tdi_fs['X'].f[index_low:index_high],tdi_fs['X'][index_low:index_high].values)
    # plt.savefig(SAVEPATH+'/subtracted.png')
    # frequencies_search = frequencies_odd[-100:]

    index = 0
    search1 = Search(tdi_fs,Tobs, frequencies_search[index][0], frequencies_search[index][1])
    print('intrinsic SNR A E T',search1.intrinsic_SNR([pGB_injected[index][0]]))
    print('intrinsic SNR A E',search1.intrinsic_SNR_old([pGB_injected[index][0]]))
    print('intrinsic SNR T',search1.intrinsic_SNR_T([pGB_injected[index][0]]))
    print('SNR A E T',search1.SNR([pGB_injected[index][0]]))
    print('SNR A E T compute',search1.SNR_AET_compute([pGB_injected[index][0]]))
    print('SNR A E',search1.SNR_AE([pGB_injected[index][0]]))
    print('SNR XYZ',search1.SNR_XYZ([pGB_injected[index][0]]))
    print('SNR XYZ Sa',search1.SNR_XYZ_Sa([pGB_injected[index][0]]))
    print('SNR noise matrix',search1.SNR_noise_matrix([pGB_injected[index][0]]))
    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB_injected[index][0], oversample=4)
    tdi = dict({"X":Xs, "Y":Ys, "Z":Zs})
    print(np.sqrt(compute_tdi_snr(tdi, Nmodel)["tot2"]))

    print(search1.SNR_with_rolling_mean([pGB_injected[index][0]]))

    start = time.time()
    for i in range(100):
        search1.SNR([pGB_injected[index][0]])
    print('time SNR',time.time()-start)
    start = time.time()
    for i in range(100):
        search1.SNR_AE([pGB_injected[index][0]])
    print('time AE',time.time()-start)
    start = time.time()
    for i in range(100):
        search1.SNR_T([pGB_injected[index][0]])
    print('time T',time.time()-start)
    start = time.time()
    for i in range(100):
        search1.SNR_AET_compute([pGB_injected[index][0]])
    print('time AET compute',time.time()-start)
    start = time.time()
    for i in range(100):
        search1.SNR_XYZ_Sa([pGB_injected[index][0]])
    print('time XYZ',time.time()-start)
    start = time.time()
    for i in range(100):
        search1.SNR_noise_matrix([pGB_injected[index][0]])
    print('time noise matrix XYZ',time.time()-start)

    start = time.time()
    for i in range(100):
        Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB_injected[index][0], oversample=4)
        tdi = dict({"X":Xs, "Y":Ys, "Z":Zs})
        np.sqrt(compute_tdi_snr(tdi, Nmodel)["tot2"])
    print('time',time.time()-start)

    do_match = True
    pGB_injected_matched = []
    pGB_injected_not_matched = deepcopy(pGB_injected)
    number_of_matched_signals = 0
    if do_match:
        start = time.time()
        for i in range(len(found_sources_in)):
            pGB_injected_matched.append([])
            # if i != 3:
            #     continue
            for j in range(len(found_sources_in[i])):
                found_match = False
                # if j != 1:
                #     continue
                # print('i', i, 'j',j)
                for k in range(len(pGB_injected_not_matched[i])):
                    eclipticlongitude = pGB_injected_not_matched[i][k]['EclipticLongitude']
                    if pGB_injected_not_matched[i][k]['EclipticLongitude'] > np.pi:
                        eclipticlongitude -= np.pi*2
                    # print('SNR', SNR_match(pGB_injected_not_matched[i][k],found_sources_in[i][j]),'parameter comparison:',pGB_injected_not_matched[i][k]['EclipticLatitude'],found_sources_in[i][j]['EclipticLatitude'],eclipticlongitude, found_sources_in[i][j]['EclipticLongitude'])
                    if SNR_match(pGB_injected_not_matched[i][k],found_sources_in[i][j]) > 0.5:
                        found_match = True
                    if found_match:
                        pGB_injected_matched[-1].append(pGB_injected_not_matched[i][k])
                        pGB_injected_not_matched[i].pop(k)
                        number_of_matched_signals += 1
                        break
        print('time to match', time.time()-start)
                    # if k == len(pGB_injected_not_matched[i])-1:
                    #     if not(found_match):
                    #         pGB_injected_not_matched[-1].append(found_sources_in[i][j])
                            
    intrinsic_SNR = []
    for i in range(len(pGB_injected)):
        if i != 5:
            continue
        search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
        for j in range(len(pGB_injected[i])):
            intrinsic_SNR.append(search1.intrinsic_SNR([pGB_injected[i][j]]))
            print('SNR for noise model', noise_model, intrinsic_SNR[-1],'rolling mean SNR', search1.SNR_with_rolling_mean([pGB_injected[i][j]]), 'loglikelihood ratio',search1.loglikelihood([pGB_injected[i][j]]), 'SNR data',search1.SNR([pGB_injected[i][j]]))
                            
    intrinsic_SNR = []
    for i in range(len(pGB_injected_not_matched)):
        # if i != 5:
        #     continue
        search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
        for j in range(len(pGB_injected_not_matched[i])):
            intrinsic_SNR.append(search1.intrinsic_SNR([pGB_injected_not_matched[i][j]]))
            print('SNR for noise model', noise_model, intrinsic_SNR[-1], 'loglikelihood ratio',search1.loglikelihood([pGB_injected_not_matched[i][j]]), 'SNR data',search1.SNR([pGB_injected_not_matched[i][j]]), 'frequency', pGB_injected_not_matched[i][j]['Frequency'])

    number_of_injected_signals = 0
    for i in range(len(pGB_injected)):
        for j in range(len(pGB_injected[i])):
            number_of_injected_signals += 1
    number_of_found_signals = 0
    for i in range(len(found_sources_in)):
        for j in range(len(found_sources_in[i])):
            number_of_found_signals += 1
    print(number_of_matched_signals ,'matched signals out of', number_of_injected_signals , 'injected signals and',number_of_found_signals, 'found signals')
    print('sensitivity = matched signals/injected signals:', number_of_matched_signals/number_of_injected_signals)
    # pGB_injected = pGB_injected_matched

    #plot strains
    for i in range(len(frequencies_search)):
        if i != 6:
            continue
        lower_frequency = frequencies_search[i][0]
        upper_frequency = frequencies_search[i][1]
        search1 = Search(tdi_fs,Tobs, lower_frequency, upper_frequency)
        found_extended = found_sources_in[i]#+found_sources_in[i+1]
        injected_extended = pGB_injected[i]#+pGB_injected[i+1]
        matched_extended = pGB_injected_matched[i]#+pGB_injected_matched[i+1]
        if len(pGB_injected[i]) > 0:
            search1.plot(found_sources_in=found_extended, pGB_injected= injected_extended, pGB_injected_matched= matched_extended, saving_label =SAVEPATH+'/strain added'+ str(int(np.round(lower_frequency*10**8))) +save_name+'.png') 
            # search1.plot(found_sources_in=found_sources_mp_best[i], pGB_injected=pGB_injected[i][:10], pGB_injected_matched= matched_extended, saving_label =SAVEPATH+'/strain added'+ str(int(np.round(lower_frequency*10**8))) +save_name+'in.png') 



    SNR_threshold = 10
    number_of_found_signals = 0
    for i in range(int(len(found_sources_in))):
        for j in range(len(found_sources_in[i])):
            number_of_found_signals += 1
    number_of_injected_signals = 0
    number_of_injected_in_window = {}
    for i in range(6):
        number_of_injected_in_window[str(i+1)] = 0
    for i in range(int(len(pGB_injected))):
        search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
        if 1 == len(pGB_injected[i]):
            number_of_injected_in_window['1'] += 1
        if 2 == len(pGB_injected[i]):
            number_of_injected_in_window['2'] += 1
        if 3 == len(pGB_injected[i]):
            number_of_injected_in_window['3'] += 1
        if 4 == len(pGB_injected[i]):
            number_of_injected_in_window['4'] += 1
        if 5 ==  len(pGB_injected[i]):
            number_of_injected_in_window['5'] += 1
        if 5 < len(pGB_injected[i]):
            number_of_injected_in_window['6'] += 1
        for j in range(len(pGB_injected[i])):
            number_of_injected_signals += 1
            pGB_injected[i][j]['SNR'] = float(search1.SNR([pGB_injected[i][j]])[2].values)
    number_of_injected_signals_high_SNR = 0
    for i in range(int(len(pGB_injected))):
        for j in range(len(pGB_injected[i])):
            if pGB_injected[i][j]['SNR'] > SNR_threshold:
                number_of_injected_signals_high_SNR += 1  
   

    
    found_sources_in_all = []
    for i in range(len(found_sources_mp_all)):
        found_sources_in_all.append([])
        for j in range(len(found_sources_mp_all[i])):
            if found_sources_mp_all[i][j]['Frequency'] > frequencies_search[i][0] and found_sources_mp_all[i][j]['Frequency'] < frequencies_search[i][1]:
                found_sources_in_all[i].append(found_sources_mp_all[i][j])

    found_sources_in_all = []
    number_of_evaluations = []
    for i in range(len(found_sources_mp)):
        number_of_evaluations.append([])
        for j in range(len(found_sources_mp[i][1])):
            number_of_evaluations[i].append(found_sources_mp[i][2][j])


    #check loglikelihood
    higherSNR = 0
    search_results = {}
    search_results['Frequency'] = []
    search_results['success rate'] = []
    search_results['nfe'] = []
    for i in range(len(found_sources_mp_all)):
        higher_loglikelihood = 0
        search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
        # for j in range(len( pGB_injected[i])):
            # print(frequencies_search[i], search1.loglikelihood([pGB_injected[i][j]]))
        for j in range(len(found_sources_mp_all[i])):
            # print('found', search1.loglikelihood([found_sources_mp_all[i][j]]))
            if search1.loglikelihood([pGB_injected[i][0]]) < search1.loglikelihood([found_sources_mp_all[i][j]]):
                higherSNR += 1
                higher_loglikelihood += 1    
        search_results['success rate'].append(higher_loglikelihood/len(found_sources_mp_all[i]))
        search_results['nfe'].append(int(np.mean(number_of_evaluations[i])))
        search_results['Frequency'].append(pGB_injected[i][0]['Frequency']*1000)
        print('higherloglikelihood ',higher_loglikelihood, 'number of evaluations', np.mean(number_of_evaluations[i]))
    df_search_results = pd.DataFrame(data=search_results)
    print(df_search_results.to_latex(index=False))
    # padding_of_initial_guess_range = 0
    # found_sources_in = []
    # for i in range(len(frequencies_search)):
    #     found_sources_previous_in_range = found_sources_sorted[found_sources_sorted['Frequency'] > frequencies_search[i][0]-padding_of_initial_guess_range]
    #     found_sources_in.append(found_sources_previous_in_range[found_sources_previous_in_range['Frequency'] < frequencies_search[i][1]+padding_of_initial_guess_range])

    #check SNR
    for i in range(len(found_sources_in)):
        # if i != 3:
        #     continue
        print('frequency range', frequencies_search[i][0],frequencies_search[i][1])
        search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
        for j in range(len( pGB_injected[i][:10])):
            #subtract the found sources from original
            tdi_fs_subtracted = deepcopy(tdi_fs)
            for n in range(len( pGB_injected[i][:j])):
                Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=pGB_injected[i][n], oversample=4, simulator="synthlisa")
                source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
                index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
                index_high = index_low+len(Xs_subtracted)
                for k in ["X", "Y", "Z"]:
                    tdi_fs_subtracted[k].data[index_low:index_high] = tdi_fs_subtracted[k].data[index_low:index_high] - source_subtracted[k].data
            search_subtracted = Search(tdi_fs_subtracted,Tobs, frequencies_search[i][0], frequencies_search[i][1])
            print('true subtracted',np.round(search_subtracted.SNR([pGB_injected[i][j]]).values,2), 'original data', np.round(search1.SNR([pGB_injected[i][j]]).values,2))
            print('true subtracted ratio',np.round(search_subtracted.loglikelihood([pGB_injected[i][j]]),2), 'original data ratio', np.round(search1.loglikelihood([pGB_injected[i][j]]),2))
        for j in range(len(found_sources_in[i])):
            #subtract the found sources from original
            tdi_fs_subtracted = deepcopy(tdi_fs)
            for n in range(len( found_sources_in[i][:j])):
                Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_in[i][n], oversample=4, simulator="synthlisa")
                source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
                index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
                index_high = index_low+len(Xs_subtracted)
                for k in ["X", "Y", "Z"]:
                    tdi_fs_subtracted[k].data[index_low:index_high] = tdi_fs_subtracted[k].data[index_low:index_high] - source_subtracted[k].data
            search_subtracted = Search(tdi_fs_subtracted,Tobs, frequencies_search[i][0], frequencies_search[i][1])
            print('found subtracted',np.round(search_subtracted.SNR([found_sources_in[i][j]]).values,2), 'original data', np.round(search1.SNR([found_sources_in[i][j]]).values,2))
            print('found subtracted ratio',np.round(search_subtracted.loglikelihood([found_sources_in[i][j]]),2), 'original data ratio', np.round(search1.loglikelihood([found_sources_in[i][j]]),2))
            # print('found', search1.SNR([found_sources_mp_even_all[i][j]]))
    # for j in range(len(found_sources_mp_all[i])):
    #     print(np.round(search1.loglikelihood([found_sources_mp_all[i][j]]),2))

    # #check loglikelihood all
    # higherSNR = 0
    # total_searches = 0
    # for i in range(len(found_sources_mp_even_all)):
    #     if len(found_sources_mp_even_all[i]) == 0:
    #         pass
    #     else:
    #         search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
    #         for j in range(len( pGB_injected[i][:10])):
    #             print('true',search1.loglikelihood([pGB_injected[i][j]]))
    #         for j in range(len(found_sources_mp_even_all[i])):
    #             print('found', search1.loglikelihood([found_sources_mp_even_all[i][j]]))
    #             if search1.loglikelihood([pGB_injected[i][0]]) < search1.loglikelihood([found_sources_mp_even_all[i][j]]):
    #                 higherSNR += 1
    #             total_searches += 1
    # print('Number of higher SNR signals',higherSNR, 'out of', total_searches, 'searches')

    #plot strains
    for i in range(len(frequencies_search)):
        if i != 0:
            continue
        lower_frequency = frequencies_search[i][0]
        upper_frequency = frequencies_search[i][1]
        search1 = Search(tdi_fs,Tobs, lower_frequency, upper_frequency)

        A_optimized = search1.calculate_Amplitude([found_sources_mp_all[i][6]])
        found_sources_mp_all[i][6]['Amplitude'] *= A_optimized.values
        if len(pGB_injected[i]) > 0:
            search1.plot(found_sources_in=found_sources_mp_best[i], pGB_injected=pGB_injected[i], saving_label =SAVEPATH+'/strain added'+ str(int(np.round(lower_frequency*10**8))) +save_name+'.png') 
            # search1.plot(pGB_injected=pGB_injected[i], saving_label =SAVEPATH+'/strain added'+ str(int(np.round(lower_frequency*10**8))) +save_name+'.png') 
            # search1.plot(found_sources_in=found_sources_in[i], pGB_injected=pGB_injected[i][:10], saving_label =SAVEPATH+'/strain added'+ str(int(np.round(lower_frequency*10**8))) +save_name+'in.png') 
        correlation = SNR_match(found_sources_mp_best[i][1],found_sources_mp_best[i][3])
#     #subtract the found sources from original
#     tdi_fs_subtracted = deepcopy(tdi_fs)
#     for i in range(len(found_sources_in)):
#         for j in range(len(found_sources_in[i])):
#             Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_in[i][j], oversample=4, simulator="synthlisa")
#             source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
#             index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
#             index_high = index_low+len(Xs_subtracted)
#             for k in ["X", "Y", "Z"]:
#                 tdi_fs_subtracted[k].data[index_low:index_high] = tdi_fs_subtracted[k].data[index_low:index_high] - source_subtracted[k].data

#     MLP = MLP_search(tdi_fs_subtracted, Tobs, signals_per_window = 10)
#     start = time.time()
#     pool = mp.Pool(mp.cpu_count())
#     found_sources_mp_odd = pool.starmap(MLP.search, frequencies_odd)
#     pool.close()
#     pool.join()
#     print('time to search ', number_of_windows, 'windows: ', time.time()-start)

#     found_sources_in = []
#     found_sources_out = []
#     for i in range(len(found_sources_mp_odd)):
#         found_sources_in.append([])
#         found_sources_out.append([])
#         for j in range(len(found_sources_mp_odd[i])):
#             if found_sources_mp_odd[i][j]['Frequency'] > frequencies_odd[i][0] and found_sources_mp_odd[i][j]['Frequency'] < frequencies_odd[i][1]:
#                 found_sources_in[i].append(found_sources_mp_odd[i][j])
#             else:
#                 found_sources_out[i].append(found_sources_mp_odd[i][j])

#     #subtract the found sources from original
#     tdi_fs_subtracted = deepcopy(tdi_fs)
#     for i in range(len(found_sources_in)):
#         for j in range(len(found_sources_in[i])):
#             Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_in[i][j], oversample=4, simulator="synthlisa")
#             source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
#             index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
#             index_high = index_low+len(Xs_subtracted)
#             for k in ["X", "Y", "Z"]:
#                 tdi_fs_subtracted[k].data[index_low:index_high] = tdi_fs_subtracted[k].data[index_low:index_high] - source_subtracted[k].data

#     MLP = MLP_search(tdi_fs_subtracted, Tobs, signals_per_window = 10)
#     start = time.time()
#     pool = mp.Pool(mp.cpu_count())
#     found_sources_mp_even = pool.starmap(MLP.search, frequencies_even)
#     pool.close()
#     pool.join()
#     print('time to search ', number_of_windows, 'windows: ', time.time()-start)

#     found_sources_mp =[]
#     for i in range(number_of_windows):
#         ind = int(i/2)
#         if i % 2 == 0:
#             found_sources_mp.append(found_sources_mp_even[ind])
#         else:
#             found_sources_mp.append(found_sources_mp_odd[ind])

#     np.save(SAVEPATH+'/found_sources_'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', found_sources_mp)
# else:
#     found_sources_mp = np.load(SAVEPATH+'/found_sources_'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', allow_pickle= True)

# nperseg = 5 * 1.0/ dt / 1e-6
# nperseg = len(tdi_fs["X"])
# noise_model = "SciRDv1"
# f, psd_x_noisy = scipy.signal.welch(tdi_ts["X"], fs=1.0/dt, window='hanning', nperseg=nperseg)
# fmin, fmax = 0.00001, 0.1
# freq = np.array(tdi_fs['X'].sel(f=slice(fmin, fmax)).f)
# freq = f[f>0]
# Nmodel = get_noise_model(noise_model, freq)
# npsd = Nmodel.psd(option='X') # could be A, E, XY
# noise_model = "MRDv1"
# Nmodel = get_noise_model(noise_model, freq)
# npsd_mrd = Nmodel.psd(option='X') # could be A, E, XY

# plt.figure(figsize=(15,6))
# plt.loglog(f, np.sqrt(psd_x_noisy), label='Signal + noise', color='orange')
# plt.loglog(freq, np.sqrt(npsd), label='Noise PSD', color='green')
# plt.loglog(freq, np.sqrt(npsd_mrd), label='Noise PSD mrd')
# # plt.axis([1e-5, 1/dt/2, 1e-24, 1e-18])
# plt.ylabel("TDI X")
# plt.xlabel("Freq [Hz]")
# plt.legend(loc='upper left')
# plt.show()

if False:
    f_line = np.logspace(-4,-1, num=20)
    # f_line = f_line[10:12]
    # f_line = [0.0003,0.0005]
    above_indexes = np.searchsorted(cat['Frequency'],0.003977)
    cat_reduced = cat[above_indexes:]
    indexes_below = np.searchsorted(cat_reduced['Frequency'],0.00401)
    cat_reduced = cat_reduced[:indexes_below]
    negative_fd_indexes = cat_reduced['FrequencyDerivative'] < 0
    print(cat_reduced[negative_fd_indexes])
    print(len(cat_reduced[negative_fd_indexes])/len(cat_reduced))


    pGB_injected = []
    frequencies_plot = []
    for i in range(len(f_line)-1):
        frequencies_plot.append([f_line[i],f_line[i+1]])
    # frequencies_plot.append([0.00950403-0.0001,0.00950403+0.0001])
    # frequencies_plot = [frequencies_plot[5]]
    cat_plot = []
    for j in range(len(frequencies_plot)):
        index_low = np.searchsorted(cat_reduced['Frequency'], frequencies_plot[j][0])
        index_high = np.searchsorted(cat_reduced['Frequency'], frequencies_plot[j][1])
        try:
            cat_plot.append(cat_reduced[index_low:index_low+50])
        except:
            cat_plot.append(cat_reduced[index_low:])
    for j in range(len(frequencies_plot)):
        index_low = np.searchsorted(cat_reduced[negative_fd_indexes]['Frequency'], frequencies_plot[j][0])
        index_high = np.searchsorted(cat_reduced[negative_fd_indexes]['Frequency'], frequencies_plot[j][1])
        try:
            cat_plot.append(cat_reduced[negative_fd_indexes][index_low:index_low+50])
        except:
            cat_plot.append(cat_reduced[negative_fd_indexes][index_low:])

    fig = plt.figure()
    parameter_x = 'Frequency'
    parameter_y = 'FrequencyDerivative'
    for i in range(len(cat_plot)):
        for j in range(len(cat_plot[i])):
            plt.scatter(cat_plot[i][parameter_x][j],cat_plot[i][parameter_y][j])
    plt.plot(f_line, frequency_derivative(f_line,0.1))
    plt.plot(f_line, frequency_derivative(f_line,0.5))
    plt.plot(f_line, frequency_derivative(f_line,M_chirp_upper_boundary))
    # plt.plot(f_line, frequency_derivative(f_line,100))
    plt.plot(f_line, frequency_derivative_tyson(f_line))
    plt.plot(f_line, frequency_derivative_tyson_lower(f_line))
    plt.hlines(0.01/Tobs**2, xmin=f_line[0], xmax=f_line[-1], linestyles='--')
    # plt.hlines(0.01/(Tobs/2)**2, xmin=f_line[0], xmax=f_line[-1], linestyles='--')
    plt.hlines(-0.01/(Tobs)**2, xmin=f_line[0], xmax=f_line[-1], linestyles='--')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.xlabel('$\log $f $[$Hz$]$')
    plt.ylabel('$\log  \dot{f} [s^{-2}]$')
    # plt.savefig(SAVEPATH+'/found_sources_'+save_name+'f-fd.png')