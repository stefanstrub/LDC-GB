#%%
from re import A
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
import yaml

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
        # elif parameter == "FrequencyDerivative":
        #     boundaries_reduced[parameter] = [boundaries[parameter][0],-16.5]
        #     boundaries_reduced[parameter] = [
        #         np.log10(maxpGB[parameter]) - length * ratio / 2*2,
        #         np.log10(maxpGB[parameter]) + length * ratio / 2*2,
        #     ]
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

from search_class import Search

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

DATAPATH = "/home/stefan/LDC/Radler/data"
DATAPATH = grandparent+"/LDC/Radler/data"
SAVEPATH = grandparent+"/LDC/Radler/LDC1-4_evaluation"

# sangria_fn = DATAPATH + "/dgb-tdi.h5"
# sangria_fn = DATAPATH + "/LDC1-3_VGB_v2.hdf5"
sangria_fn = DATAPATH + "/LDC1-4_GB_v2.hdf5"
# sangria_fn = DATAPATH + "/LDC1-3_VGB_v2_FD_noiseless.hdf5"
fid = h5py.File(sangria_fn)
# get the source parameters
names = np.array(fid['H5LISA/GWSources/GalBinaries'])
params = [fid['H5LISA/GWSources/GalBinaries'][k] for k in names]
reduced_names = []
i = 0
for p in params:
    i += 1
    if p.shape:
        reduced_names.append(names[i-1])
params = [np.array(p) for p in params if p.shape]
cat = np.rec.fromarrays(params, names=list(reduced_names))

# get TDI 
td = np.array(fid["H5LISA/PreProcess/TDIdata"])
td = np.rec.fromarrays(list(td.T), names=["t", "X", "Y", "Z"])
del_t = float(np.array(fid['H5LISA/GWSources/GalBinaries']['Cadence']))
reduction = 1
Tobs = float(int(np.array(fid['H5LISA/GWSources/GalBinaries']['ObservationDuration']))/reduction)

dt = del_t
# Build timeseries and frequencyseries object for X,Y,Z
tdi_ts = dict([(k, TimeSeries(td[k][:int(len(td[k][:])/reduction)], dt=dt)) for k in ["X", "Y", "Z"]])
# tdi_ts = xr.Dataset(dict([(k, TimeSeries(td[k][:int(len(td[k][:])/reduction)], dt=dt)) for k in ["X", "Y", "Z"]]))
# tdi_ts = xr.Dataset(dict([(k,TimeSeries(tdi_ts[k][:,1], dt=dt)) for k in ["X", "Y", "Z"]]))
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

# data_file_name = 'Montana'
data_file_name = 'APC'
data_file_name = 'ETH'
save_name = 'LDC1-4_4mHz'+data_file_name
try:
    cat = np.load(SAVEPATH+'/cat_sorted.npy', allow_pickle = True)
    print('cat sorted loaded')
except:
    indexes = np.argsort(cat['Frequency'])
    cat = cat[indexes]
    np.save(SAVEPATH+'/cat_sorted.npy',cat)
cat_sorted = cat
# LDC1-4 #####################################
frequencies = []
frequencies_even = []
frequencies_odd = []
# search_range = [0.00398, 0.0041]
# search_range = [0.0039885, 0.0040205]
# search_range = [0.0039935, 0.0039965]
f_Nyquist = 1/dt/2
search_range = [0.0003, f_Nyquist]
# search_range = [0.0003, 0.03]
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

# found_sources = np.load(SAVEPATH+'/ETH.npy', allow_pickle = True)
# found_sources = np.load(SAVEPATH+'/Montana.npy', allow_pickle = True)
found_sources = np.load(SAVEPATH+'/'+data_file_name+'.npy', allow_pickle = True)
# found_sources = np.load(SAVEPATH+'/found_sources_not_matched.npy', allow_pickle = True)
# found_sources = np.concatenate(found_sources)
# found_sources = np.load(SAVEPATH+'/APC.npy', allow_pickle = True)

found_sources_in_flat_frequency = []
for i in range(len(found_sources)):
    found_sources_in_flat_frequency.append(found_sources[i]['Frequency'])
found_sources_in_flat_frequency = np.asarray(found_sources_in_flat_frequency)
frequency_sort_index = np.argsort(found_sources_in_flat_frequency)
found_sources_in_flat_frequency = found_sources_in_flat_frequency[frequency_sort_index]
found_sources = found_sources[frequency_sort_index]


frequencies_search = frequencies
# batch_index = int(sys.argv[1])
# batch_index = 65
start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], found_sources_in_flat_frequency.min())-1
end_index = np.searchsorted(np.asarray(frequencies_search)[:,0], found_sources_in_flat_frequency.max())
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.00040707)
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.007977)
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], cat_sorted[-2]['Frequency'])-1
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.0004)-1
# batch_size = 20
# start_index = batch_size*batch_index
# print('batch',batch_index, start_index)
frequencies_search = frequencies_search[start_index:end_index]
# print(i, frequencies_search[0])
### highest + padding has to be less than f Nyqist
while frequencies_search[-1][1] + (frequencies_search[-1][1] - frequencies_search[-1][0])/2 > f_Nyquist:
    frequencies_search = frequencies_search[:-1]
# frequencies_search = frequencies_search[70:80]

search_range = [frequencies_search[0][0],frequencies_search[-1][1]]
# search_range = [1619472*10**-8,2689639*10**-8]
print('search range '+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))))

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
        # SNR2 = np.sum( np.real(Af_injected * np.conjugate(Af.data) + Ef_injected * np.conjugate(Ef.data))/SA)
        SNR2 = np.sum( np.real(np.abs(Af_injected) * np.abs(Af.data) + np.abs(Ef_injected) * np.abs(Ef.data))/SA)
        hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /SA)
        ss = np.sum((np.absolute(Af_injected.data)**2 + np.absolute(Ef_injected.data)**2) /SA)
    SNR = 4.0*Xs.df* hh
    SNR2 = 4.0*Xs.df* SNR2
    SNR3 = SNR2 / (np.sqrt(SNR)*np.sqrt(4.0*Xs.df* ss))
    return SNR3.values


found_sources_in = []
for i in range(len(frequencies_search)):
    lower_index = np.searchsorted(found_sources_in_flat_frequency,frequencies_search[i][0])
    higher_index = np.searchsorted(found_sources_in_flat_frequency,frequencies_search[i][1])
    found_sources_in.append(found_sources[lower_index:higher_index])



# for i in range(len(found_sources_in)):
#     if len(found_sources_in[i]) > 0:
#         search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
#         for j in range(len(found_sources_in[i])):
#             pGB_dict = {}
#             for parameter in parameters:
#                 pGB_dict[parameter] = found_sources_in[i][j][parameter]
#             found_sources_in[i][j]['IntrinsicSNR'] = search1.intrinsic_SNR([pGB_dict])

# found_sources_in_flat = np.concatenate(found_sources_in)

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
# with open(SAVEPATH+'ETH_LDC1-4_4mHz_full.yaml', 'w') as file:
#     documents = yaml.dump(member, file)



frequencies_search = np.asarray(frequencies_search)
xx = np.arange(stop=len(frequencies_search))
fig = plt.figure()
for i in range(len(frequencies_search)):
    plt.plot(frequencies_search[i,0],0,'r.')
    plt.plot(frequencies_search[i,1],0,'g.')
for i in range(2604, 2604+len(found_sources_in_flat_frequency[2604:2634])):
    plt.plot(found_sources_in_flat_frequency[i],0,'b.')
plt.show()

new_dt = np.dtype(cat_sorted.dtype.descr + [('IntrinsicSNR','<f8')])
cat_sorted_SNR = np.zeros(cat_sorted.shape, dtype=new_dt)
for parameter in parameters:
    cat_sorted_SNR[parameter] = cat_sorted[parameter]
cat_sorted = cat_sorted_SNR

get_pGB_injected = True
if get_pGB_injected:
    pGB_injected = []
    for j in range(len(frequencies_search)):
        padding = (frequencies_search[j][1] - frequencies_search[j][0])/2 *0
        index_low = np.searchsorted(cat_sorted['Frequency'], frequencies_search[j][0]-padding)
        index_high = np.searchsorted(cat_sorted['Frequency'], frequencies_search[j][1]+padding)
        try:
            if cat_sorted['Frequency'][index_high] < frequencies_search[j][1]:
                index_high -= 1
        except:
            pass
        indexesA = np.argsort(-cat_sorted[index_low:index_high]['Amplitude'])
        pGB_injected.append(cat_sorted[index_low:index_high][indexesA])

def get_SNR(pGB_injected, lower_frequency, upper_frequency):
    search1 = Search(tdi_fs,Tobs, lower_frequency, upper_frequency)
    intrinsic_SNR_injected = []
    # print(lower_frequency)
    for i in range(len(pGB_injected)):
        pGB_injected_dict = {}
        for parameter in parameters:
            pGB_injected_dict[parameter] = pGB_injected[i][parameter]
        intrinsic_SNR_injected.append(search1.intrinsic_SNR([pGB_injected_dict]))
        if i > 99:
            break
    return intrinsic_SNR_injected

##### parallel
intrinsic_SNR_injected = []
input = []
index = []
for i in range(len(pGB_injected)):
# for i in range(16):
    intrinsic_SNR_injected.append([])
    if len(pGB_injected[i]) > 0:
        index.append(i)
        input.append((pGB_injected[i],frequencies_search[i][0], frequencies_search[i][1]))
start = time.time()
pool = mp.Pool(16)
SNR_intrinsic = pool.starmap(get_SNR, input)
pool.close()
pool.join()
print('time to calculate SNR for', len(frequencies_search), 'windows: ', time.time()-start)
for i in range(len(SNR_intrinsic)):
    for j in range(len(SNR_intrinsic[i])):
        if len(pGB_injected[i]) > 0:
            pGB_injected[index[i]][j]['IntrinsicSNR'] = SNR_intrinsic[i][j]

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
#             pGB_injected_dict[parameter] = pGB_injected[i][j][parameter]
#         pGB_injected[i][j]['IntrinsicSNR'] = search1.intrinsic_SNR([pGB_injected_dict])
#         if j > 300:
#             break
        # print('SNR for noise model', noise_model, intrinsic_SNR_injected[-1], 'loglikelihood ratio',search1.loglikelihood([pGB_injected[i][j]]), 'SNR data',search1.loglikelihood_SNR([pGB_injected[i][j]]))

pGB_injected_SNR_sorted = []
for i in range(len(pGB_injected)):
    indexesSNR = np.argsort(-pGB_injected[i]['IntrinsicSNR'])
    pGB_injected_SNR_sorted.append(pGB_injected[i][indexesSNR])
pGB_injected = pGB_injected_SNR_sorted

# np.save(SAVEPATH+'/found_sources_pGB_injected'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', np.asarray(pGB_injected))
# pGB_injected = np.load(SAVEPATH+'/found_sources_pGB_injected'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', allow_pickle = True)

pGB_injected_SNR_sorted = pGB_injected


do_match = True
pGB_injected_matched = []
found_sources_matched = []
pGB_injected_not_matched = deepcopy(pGB_injected_SNR_sorted)
found_sources_not_matched = deepcopy(found_sources_in)
number_of_matched_signals = 0
correlation_list = []
if do_match:
    start = time.time()
    percentage = 0
    for i in range(len(found_sources_in)):
    # for i in range(2000,2200):
        # found_sources_not_matched[i] = found_sources_not_matched[i].tolist()
        pGB_injected_matched.append([])
        found_sources_matched.append([])
        try:
            if i%int(len(found_sources_in)/100) == 0:
                percentage += 1
                print(percentage,'%')
        except:
            pass
        # if i != 3:
        #     continue
        for j in range(len(found_sources_in[i])):
            found_match = False
            # if j != 1:
            #     continue
            # print('i', i, 'j',j)
            correlation_list_of_one_signal = []
            for k in range(len(pGB_injected_not_matched[i])):
                pGB_injected_dict = {}
                found_dict = {}
                for parameter in parameters:
                    pGB_injected_dict[parameter] = pGB_injected_not_matched[i][k][parameter]
                    found_dict[parameter] = found_sources_in[i][j][parameter]
                # print('SNR', SNR_match(pGB_injected_not_matched[i][k],found_sources_in[i][j]),'parameter comparison:',pGB_injected_not_matched[i][k]['EclipticLatitude'],found_sources_in[i][j]['EclipticLatitude'],eclipticlongitude, found_sources_in[i][j]['EclipticLongitude'])
                correlation = SNR_match(pGB_injected_dict,found_dict)
                correlation_list_of_one_signal.append(correlation)
                if k > 19:
                    break
                # correlation = SNR_match(pGB_injected_dict,found_dict)
            if 0 == len(correlation_list_of_one_signal):
                break
            max_index = np.argmax(correlation_list_of_one_signal)
            if correlation_list_of_one_signal[max_index] > 0.5:
                found_match = True
            if found_match:
                pGB_injected_matched[-1].append(pGB_injected_not_matched[i][max_index])
                found_sources_matched[-1].append(found_sources_in[i][j])
                # pGB_injected_not_matched[i] = np.delete(pGB_injected_not_matched[i], max_index)
                correlation_list.append(correlation_list_of_one_signal[max_index])
                found_sources_not_matched[i][j] = None
                number_of_matched_signals += 1
    print('time to match', time.time()-start)
    
    for i in range(len(found_sources_not_matched)):
        found_sources_not_matched[i] = list(filter(None, found_sources_not_matched[i]))

pGB_injected_not_matched_with_overlap_of_windows = deepcopy(pGB_injected_not_matched)
pGB_injected_not_matched = []
for i in range(len(pGB_injected_not_matched_with_overlap_of_windows)):
    for j in range(len(pGB_injected_not_matched_with_overlap_of_windows[i])):
        if pGB_injected_not_matched_with_overlap_of_windows[i][j]['Frequency'] > frequencies_search[i][0] and pGB_injected_not_matched_with_overlap_of_windows[i][j]['Frequency'] < frequencies_search[i][1]:
            pGB_injected_not_matched.append(pGB_injected_not_matched_with_overlap_of_windows[i][j])
pGB_injected_flat = []
for i in range(len(pGB_injected)):
    for j in range(len(pGB_injected[i])):
        if pGB_injected[i][j]['Frequency'] > frequencies_search[i][0] and pGB_injected[i][j]['Frequency'] < frequencies_search[i][1]:
            pGB_injected_flat.append(pGB_injected[i][j])


# number_of_injected_signals = 0
# for i in range(len(pGB_injected)):
#     for j in range(len(pGB_injected[i])):
#         number_of_injected_signals += 1
number_of_injected_signals = len(pGB_injected_flat)
number_of_found_signals = 0
for i in range(len(found_sources_in)):
    for j in range(len(found_sources_in[i])):
        number_of_found_signals += 1
print(number_of_matched_signals ,'matched signals out of', number_of_injected_signals , 'injected signals and',number_of_found_signals, 'found signals')
print('sensitivity = matched signals/injected signals:', np.round(number_of_matched_signals/number_of_injected_signals,2))
print('matched signals/found signals:', np.round(number_of_matched_signals/number_of_found_signals,2))
# pGB_injected = pGB_injected_matched


pGB_injected_matched_flat = []
for i in range(len(pGB_injected_matched)):
    for j in range(len(pGB_injected_matched[i])):
        pGB_injected_matched_flat.append(pGB_injected_matched[i][j])

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

for i in range(len(found_sources_matched)):
    if len(found_sources_matched[i]) > 0:
        search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
    for j in range(len(found_sources_matched[i])):
        pGB_dict = {}
        for parameter in parameters:
            pGB_dict[parameter] = found_sources_matched[i][j][parameter]
        found_sources_matched[i][j]['IntrinsicSNR'] = search1.intrinsic_SNR([pGB_dict])

for i in range(len(found_sources_not_matched)):
    if len(found_sources_not_matched[i]) > 0:
        search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
    for j in range(len(found_sources_not_matched[i])):
        pGB_dict = {}
        for parameter in parameters:
            pGB_dict[parameter] = found_sources_not_matched[i][j][parameter]
        found_sources_not_matched[i][j]['IntrinsicSNR'] = search1.intrinsic_SNR([pGB_dict])
        
# np.save(SAVEPATH+'/found_sources_not_matched.npy', found_sources_not_matched)

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


### determine index of a specific frequency
index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.003988)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.0040225)-1
#plot strains
for i in range(len(frequencies_search)):
    if i != index_of_interest_to_plot:
        continue
    search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
    found_extended = found_sources_matched[i]#+found_sources_not_matched[i]
    matched_extended = pGB_injected_matched[i]#+pGB_injected_matched[i+1]
    if len(pGB_injected_SNR_sorted[i]) > 0:
        if len(pGB_injected_SNR_sorted[i]) > 20:
            search1.plot(found_sources_in=found_extended, pGB_injected= pGB_injected_SNR_sorted[i][:20], found_sources_not_matched=found_sources_not_matched[i], pGB_injected_matched= matched_extended, saving_label =SAVEPATH+'/strain added'+ str(int(np.round(frequencies_search[i][0]*10**8))) +save_name+'.png') 
        else:
            search1.plot(found_sources_in=found_extended, pGB_injected= pGB_injected_SNR_sorted[i], found_sources_not_matched=found_sources_not_matched[i], pGB_injected_matched= matched_extended, saving_label =SAVEPATH+'/strain added'+ str(int(np.round(frequencies_search[i][0]*10**8))) +save_name+'.png') 
        # search1.plot(found_sources_in=found_sources_mp_best[i], pGB_injected=pGB_injected[i][:10], pGB_injected_matched= matched_extended, saving_label =SAVEPATH+'/strain added'+ str(int(np.round(frequencies_search[i][0]*10**8))) +save_name+'in.png') 

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
            for k in range(len(pGB_injected_SNR_sorted[i])):
                pGB_injected_dict = {}
                found_dict = {}
                for parameter in parameters:
                    pGB_injected_dict[parameter] = pGB_injected_SNR_sorted[i][k][parameter]
                    found_dict[parameter] = found_sources_in[i][j][parameter]
                    if parameter in ['InitialPhase','Polarization']:
                        found_dict[parameter] = pGB_injected_SNR_sorted[i][k][parameter]
                # print('SNR', SNR_match(pGB_injected_not_matched[i][k],found_sources_in[i][j]),'parameter comparison:',pGB_injected_not_matched[i][k]['EclipticLatitude'],found_sources_in[i][j]['EclipticLatitude'],eclipticlongitude, found_sources_in[i][j]['EclipticLongitude'])
                correlation = SNR_match(pGB_injected_dict,found_dict)
                print(search1.SNR([pGB_injected_dict]))
                print('found',search1.SNR([found_dict]))
                print(pGB_injected_dict,found_dict)
                print(correlation)
                correlation_list_of_one_signal.append(correlation)
                if k > 19:
                    break
            correlation_list2.append(correlation_list_of_one_signal)

#### plot SNR - frequency
parameter_to_plot = 'IntrinsicSNR'
fig = plt.figure()
# is_labeled = False
# for i in range(len(pGB_injected_matched)):
#     for j in range(len(pGB_injected_matched[i])):
#         if not(is_labeled):
#             plt.scatter(pGB_injected_matched[i][j]['Frequency']*10**3,pGB_injected_matched[i][j]['IntrinsicSNR'], color = 'green', label = 'Match')
#             is_labeled = True
#         else:
#             plt.scatter(pGB_injected_matched[i][j]['Frequency']*10**3,pGB_injected_matched[i][j]['IntrinsicSNR'], color = 'green')
is_labeled = False
for i in range(len(found_sources_matched)):
    for j in range(len(found_sources_matched[i])):
        if not(is_labeled):
            plt.scatter(found_sources_matched[i][j]['Frequency']*10**3,found_sources_matched[i][j][parameter_to_plot], color = 'orange', alpha= 1, label = 'Match')
            is_labeled = True
        else:
            plt.scatter(found_sources_matched[i][j]['Frequency']*10**3,found_sources_matched[i][j][parameter_to_plot], color = 'orange', alpha= 1)
is_labeled = False
for i in range(len(found_sources_not_matched)):
    for j in range(len(found_sources_not_matched[i])):
        if not(is_labeled):
            plt.scatter(found_sources_not_matched[i][j]['Frequency']*10**3,found_sources_not_matched[i][j][parameter_to_plot], s=80, facecolors = 'none', edgecolors = 'blue', alpha= 1, label = 'found not matched')
            is_labeled = True
        else:
            plt.scatter(found_sources_not_matched[i][j]['Frequency']*10**3,found_sources_not_matched[i][j][parameter_to_plot], s=80, facecolors = 'none', edgecolors = 'blue', alpha= 1)
# is_labeled = False
# for i in range(3500,len(pGB_injected)):
#     for j in range(len(pGB_injected[i])):
#         if not(is_labeled):
#             plt.scatter(pGB_injected[i][j]['Frequency']*10**3,pGB_injected[i][j]['IntrinsicSNR'], color = 'purple', label = 'injected')
#             is_labeled = True
#         else:
#             plt.scatter(pGB_injected[i][j]['Frequency']*10**3,pGB_injected[i][j]['IntrinsicSNR'], color = 'purple')
# is_labeled = False
# for i in range(len(pGB_injected_not_matched)):
#     for j in range(len(pGB_injected_not_matched[i])):
#         if not(is_labeled):
#             plt.scatter(pGB_injected_not_matched[i][j]['Frequency']*10**3,pGB_injected_not_matched[i][j]['IntrinsicSNR'],marker='x', color = 'red', alpha= 1, label = 'injected not matched')
#             is_labeled = True
#         else:
#             plt.scatter(pGB_injected_not_matched[i][j]['Frequency']*10**3,pGB_injected_not_matched[i][j]['IntrinsicSNR'],marker='x', color = 'red', alpha= 1)
#         if j > 3:
#             break
is_labeled = False
for i in range(len(pGB_injected_not_matched)):
    if not(is_labeled):
        plt.scatter(pGB_injected_not_matched[i]['Frequency']*10**3,pGB_injected_not_matched[i]['IntrinsicSNR'],marker='x', color = 'red', alpha= 1, label = 'injected not matched')
        is_labeled = True
    else:
        plt.scatter(pGB_injected_not_matched[i]['Frequency']*10**3,pGB_injected_not_matched[i]['IntrinsicSNR'],marker='x', color = 'red', alpha= 1)
plt.yscale('log')
# plt.xscale('log')
plt.xlabel('f (mHz)')
plt.ylabel(parameter_to_plot)    
plt.legend()
plt.savefig(SAVEPATH+'/Evaluation/'+parameter_to_plot+save_name,dpi=300,bbox_inches='tight')
plt.show()

found_sources_in_flat_frequency

### plot match histo gramm - frequency
pGB_injected_matched_frequencies = []
for i in range(len(pGB_injected_matched)):
    for j in range(len(pGB_injected_matched[i])):
        pGB_injected_matched_frequencies.append(pGB_injected_matched[i][j]['Frequency']*10**3)
pGB_injected_not_matched_frequencies = []
for i in range(len(pGB_injected_not_matched)):
    pGB_injected_not_matched_frequencies.append(pGB_injected_not_matched[i]['Frequency']*10**3)
fig = plt.figure()
plt.hist(pGB_injected_matched_frequencies, 20, color = 'green')
plt.hist(pGB_injected_not_matched_frequencies, 20, color = 'red')
plt.xlabel('f (mHz)')
plt.legend()
plt.savefig(SAVEPATH+'/Evaluation/_SNR_histo'+save_name,dpi=300,bbox_inches='tight')
plt.show()



######### get errors
error = []
for i in range(len(pGB_injected_matched)):
    error.append([])
    for j in range(len(pGB_injected_matched[i])):
        error[i].append({})
        for parameter in parameters:
            if parameter == 'EclipticLongitude':
                if pGB_injected_matched[i][j][parameter] > np.pi:
                    pGB_injected_matched[i][j][parameter] -= 2*np.pi
            error[i][j][parameter] = pGB_injected_matched[i][j][parameter] - found_sources_matched[i][j][parameter]
            if parameter in ['EclipticLongitude', 'EclipticLongitude', 'Inclination', 'InitialPhase', 'Polarization']:
                error[i][j][parameter] = np.abs(np.arcsin(np.sin(pGB_injected_matched[i][j][parameter] - found_sources_matched[i][j][parameter])))
            if parameter == 'Amplitude':
                error[i][j][parameter] = np.log10(pGB_injected_matched[i][j][parameter]) - np.log10(found_sources_matched[i][j][parameter])


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
parameter = 'Amplitude'
for parameter in parameters:
    fig = plt.figure()
    plt.errorbar(frequencies_search[:,0],mean_error[parameter],yerr=std_table[parameter])
    plt.xlabel(parameter)
    plt.ylabel('Error')
    plt.savefig(SAVEPATH+'/Evaluation/'+parameter+'_error_histo'+save_name,dpi=300,bbox_inches='tight')
    plt.show()

###### plot errors scatter
for parameter in parameters:
    fig = plt.figure()
    for i in range(len(pGB_injected_matched)):
        for j in range(len(pGB_injected_matched[i])):
            plt.scatter(pGB_injected_matched[i][j]['Frequency']*10**3,error[i][j][parameter], color = 'green')
    # plt.yscale('log')
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
    plt.savefig(SAVEPATH+'/Evaluation/'+parameter+'_error'+save_name,dpi=300,bbox_inches='tight')
    plt.show()

##### plot correlation histogramm
correlation_list= np.asarray(correlation_list)
fig = plt.figure()
plt.hist(correlation_list,50)
plt.xlabel('Correlation')
plt.ylabel('Count')
plt.yscale('log')
plt.savefig(SAVEPATH+'/Evaluation/correlation'+save_name,dpi=300,bbox_inches='tight')
plt.show()
