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
import pickle
import sys
sys.path.append('/cluster/home/sstrub/Repositories/LDC/lib/lib64/python3.8/site-packages/ldc-0.1-py3.8-linux-x86_64.egg')

from astropy import units as u
import astropy.coordinates as coord

from ldc.lisa.noise import get_noise_model
from ldc.common.series import TimeSeries
import ldc.waveform.fastGB as fastGB
from ldc.common.tools import compute_tdi_snr, window

from sources2 import *

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
version = '2'
reduction = 1

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
tdi_ts = dict([(k, TimeSeries(td[k][:int(len(td[k][:])/reduction)], dt=dt, t0=0)) for k in ["X", "Y", "Z"]])
tdi_fs = xr.Dataset(dict([(k, tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
GB = fastGB.FastGB(delta_t=dt, T=Tobs)  # in seconds

noise_model = "SciRDv1"
# noise_model = "sangria"
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
            Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_mp_subtract[i], oversample=4)
            source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
            index_low = np.searchsorted(tdi_fs_subtracted2["X"].f, Xs_subtracted.f[0])
            index_high = index_low+len(Xs_subtracted)
            for k in ["X", "Y", "Z"]:
                tdi_fs_subtracted2[k].data[index_low:index_high] -= source_subtracted[k].data
    return tdi_fs_subtracted2

try:
    cat = np.load(SAVEPATH+'cat_sorted.npy', allow_pickle = True)
    # cat = np.load(SAVEPATH+'/cat_sorted_v'+version+'.npy', allow_pickle = True)
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

frequencies = create_frequency_windows(search_range, Tobs)


# frequencies = frequencies[:32]
frequencies_even = frequencies[::2]
frequencies_odd = frequencies[1::2]

start_index = np.searchsorted(np.asarray(frequencies_odd)[:,0], 0.0040489)-1

save_name = 'Sangria_1_full_cut'
save_name = 'Sangria_12m'
# save_name = 'Radler_24m'
# save_name = 'Radler_half_even_dynamic_noise'
# save_name = 'LDC1-4_2_optimized_second' ### ETH submission
# save_name = 'Montana'
# save_name = 'APC'
# save_name = 'LDC1-4_half_year'

# duration = '3932160'
# duration = '7864320'
# duration = '15728640'
duration = '31457280'
# save_name = 'Montana2022_'+duration
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
frequencies_min = 0.0095
frequencies_max = 0.00955
start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], frequencies_min_APC)-1
end_index = np.searchsorted(np.asarray(frequencies_search)[:,0], frequencies_max_APC)
# frequencies_search = frequencies_search[start_index:end_index]

search_range = [frequencies_search[0][0],frequencies_search[-1][1]]
# search_range = [1619472*10**-8,2689639*10**-8]
print('search range '+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))))

# i = 44
# search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
# search1.plot()

def l2_norm_match(pGB_injected, pGB_found):
    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB_found, oversample=4)
    Xs_injected, Ys_injected, Zs_injected = GB.get_fd_tdixyz(template=pGB_injected, oversample=4)
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

def correlation_match(pGB_injected, pGB_found):
    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB_found, oversample=4)
    Xs_injected, Ys_injected, Zs_injected = GB.get_fd_tdixyz(template=pGB_injected, oversample=4)
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

def SNR_match_scaled(pGB_injected, pGB_found):
    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB_found, oversample=4)
    Xs_injected, Ys_injected, Zs_injected = GB.get_fd_tdixyz(template=pGB_injected, oversample=4)
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
        SNR2 = np.sum((np.absolute(Af_injected-Af.data)**2 + np.absolute(Ef_injected-Ef.data)**2)/SA + np.absolute(Tf_injected-Tf.data)**2 /ST)
        # SNR2_complex_mult = np.sum( (np.real((Af_injected-Af.data) * np.conjugate((Af_injected-Af.data))) + np.real((Ef_injected-Ef.data) * np.conjugate((Ef_injected-Ef.data))))/SA + np.real((Tf_injected-Tf.data) * np.conjugate(Tf_injected-Tf.data))/ST)
        # hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2)/SA + np.absolute(Tf.data)**2 /ST)
        ss = np.sum((np.absolute(Af_injected.data)**2 + np.absolute(Ef_injected.data)**2)/SA + np.absolute(Tf_injected.data)**2 /ST)
    else:
        # SNR2 = np.sum( np.real((Af_injected-Af.data) * np.conjugate((Af_injected-Af.data)) + (Ef_injected-Ef.data) * np.conjugate((Ef_injected-Ef.data)))/SA)
        SNR2 = np.sum((np.absolute(Af_injected-Af.data)**2 + np.absolute(Ef_injected-Ef.data)**2)/SA)
        # hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /SA)
        ss = np.sum((np.absolute(Af_injected.data)**2 + np.absolute(Ef_injected.data)**2) /SA)
    # SNR = 4.0*Xs.df* hh
    SNR3 = SNR2/ss
    return SNR3.values

def SNR_match_amplitude_condsiered(pGB_injected, pGB_found):
    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB_found, oversample=4)
    Xs_injected, Ys_injected, Zs_injected = GB.get_fd_tdixyz(template=pGB_injected, oversample=4)
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
    found_sources_in_flat = np.load(SAVEPATH+'found_sources_' +save_name+'_flat.npy', allow_pickle = True)
except:
    # found_sources_in_flat = np.load(SAVEPATH+'found_sources' +save_name+'.npy', allow_pickle = True)
    found_sources_mp = np.load(SAVEPATH+'found_sources_' +save_name+'.npy', allow_pickle = True)

    # found_sources_mp_best = []
    # found_sources_mp_all = []
    # for i in range(len(found_sources_mp)):
    #     found_sources_mp_best.append(found_sources_mp[i][0])
    #     found_sources_in_window = []
    #     for j in range(len(found_sources_mp[i][1])):
    #         found_sources_in_window.append(found_sources_mp[i][1][j][0][0])
    #     found_sources_mp_all.append(found_sources_in_window)

    found_sources_in_flat = []
    for i in range(len(found_sources_mp)):
        for j in range(len(found_sources_mp[i][3])):
            found_sources_in_flat.append(found_sources_mp[i][3][j])
    found_sources_in_flat = np.asarray(found_sources_in_flat)

    np.save(SAVEPATH+'found_sources_' +save_name+'_flat.npy', np.asarray(found_sources_in_flat))

found_sources_in_flat_frequency = []
for i in range(len(found_sources_in_flat)):
    found_sources_in_flat_frequency.append(found_sources_in_flat[i]['Frequency'])
found_sources_in_flat_frequency = np.asarray(found_sources_in_flat_frequency)
found_sources_in_flat = np.asarray(found_sources_in_flat)
indexes_in = np.argsort(found_sources_in_flat_frequency)
found_sources_in_flat_frequency = found_sources_in_flat_frequency[indexes_in]
found_sources_in_flat = found_sources_in_flat[indexes_in]

# np.save(SAVEPATH+'/found_sources_flat.npy', np.asarray(found_sources_in_flat))

# found_sources_in = []
# for i in range(len(frequencies_search)):
#     lower_index = np.searchsorted(found_sources_in_flat_frequency,frequencies_search[i][0])
#     higher_index = np.searchsorted(found_sources_in_flat_frequency,frequencies_search[i][1])
#     found_sources_in.append(found_sources_in_flat[lower_index:higher_index])

found_sources_in_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_in_flat]) for attribute in found_sources_in_flat[0].keys()}
found_sources_in_flat_df = pd.DataFrame(found_sources_in_flat_array)
found_sources_in_flat_df = found_sources_in_flat_df.sort_values('Frequency')
found_sources_in = []
for i in range(len(frequencies_search)):
    found_sources_in.append(found_sources_in_flat_df[(found_sources_in_flat_df['Frequency'] > frequencies_search[i][0]) & (found_sources_in_flat_df['Frequency'] < frequencies_search[i][1])].to_dict(orient='records'))

# frequencies_search_with_found_only = []
# found_sources_in_with_found_only = []
# for i in range(len(frequencies_search)):
#     if len(found_sources_in[i]) > 0:
#         frequencies_search_with_found_only.append(frequencies_search[i])
#         found_sources_in_with_found_only.append(found_sources_in[i])
# frequencies_search = frequencies_search_with_found_only
# found_sources_in = []
# for i in range(len(frequencies_search)):
#     lower_index = np.searchsorted(found_sources_in_flat_frequency,frequencies_search[i][0])
#     higher_index = np.searchsorted(found_sources_in_flat_frequency,frequencies_search[i][1])
#     found_sources_in.append(found_sources_in_flat[lower_index:higher_index])

new_dt = np.dtype(cat_sorted.dtype.descr + [('IntrinsicSNR','<f8')])
cat_sorted_SNR = np.zeros(cat_sorted.shape, dtype=cat_sorted.dtype.descr + [('IntrinsicSNR','<f8')])
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
    search1 = Search(tdi_fs,Tobs, lower_frequency, upper_frequency, update_noise=False)
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
    # found_sources_in = np.load(SAVEPATH+'/found_sources_in_SNR_'+save_name+'.npy', allow_pickle= True)
    
    fn = SAVEPATH+'/found_sources_in_SNR_'+save_name+'.pkl'
    found_sources_in = pickle.load(open(fn, 'rb'))

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

    fn = SAVEPATH+'/found_sources_in_SNR_'+save_name+'.pkl'
    pickle.dump(found_sources_in, open(fn, "wb"))
    # np.save(SAVEPATH+'/found_sources_in_SNR'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy',found_sources_in)

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
# # outside of normal range
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


# np.save(SAVEPATH+'/pGB_injected_intrinsic_SNR'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', np.asarray(pGB_injected))
# np.save(SAVEPATH+'/pGB_injected_intrinsic_SNR_sorted'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', np.asarray(pGB_injected_SNR_sorted))
# pGB_injected = np.load(SAVEPATH+'/found_sources_pGB_injected'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', allow_pickle = True)
# pGB_injected_even = np.load(SAVEPATH+'/found_sources_pGB_injected_array2_half_year30000to3305084LDC1-4_half_even_T'+'.npy', allow_pickle = True)
# pGB_injected_odd = np.load(SAVEPATH+'/found_sources_pGB_injected_in_out_intrinsic_SNR_sorted30035to3316929LDC1-4_half_odd_T'+'.npy', allow_pickle = True)

if Radler:
    if reduction == 1:
        # pGB_injected = np.load(DATAPATH+'/pGB_injectedthird30000to3293090LDC1-4_2_years_full.npy', allow_pickle = True)
        pGB_injected = np.load(SAVEPATH+'/pGB_injected_intrinsic_SNR_sorted30000to3293090LDC1-4_2_optimized_second.npy', allow_pickle = True)
    elif reduction == 2:
        pGB_injected = np.load(SAVEPATH+'/pGB_injected_intrinsic_SNR_sorted30000to3312283Montana2022_31457280.npy', allow_pickle = True)
    elif reduction == 4:
        pGB_injected = np.load(SAVEPATH+'/pGB_injected_intrinsic_SNR_sorted30000to3316929LDC1-4_half_year.npy', allow_pickle = True)
    for i in range(len(pGB_injected)):
        pGB_injected[i] = pGB_injected[i].to_records()
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
# pGB_injected = np.array(pGB_injected)
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
    indexesSNR = np.argsort(-pGB_injected_overlap[i]['IntrinsicSNR'])
    # indexesSNR = np.argsort(-pGB_injected_overlap[i]['Amplitude'])
    pGB_injected_SNR_sorted_overlap.append(pGB_injected_overlap[i][indexesSNR])

### Reduce the number of signals per window
for i in range(len(pGB_injected_overlap)):
    pGB_injected_SNR_sorted_overlap[i] = pGB_injected_SNR_sorted_overlap[i][:50]

### just take the SNR > threshold signals
for i in range(len(pGB_injected_overlap)):
    pGB_injected_SNR_sorted_overlap[i] = pGB_injected_SNR_sorted_overlap[i][pGB_injected_SNR_sorted_overlap[i]['IntrinsicSNR'] > 5]

# i = 6837
# pGB_injected_overlap = pGB_injected_SNR_sorted_overlap
# pGB_injected_overlap_flat = np.concatenate(pGB_injected_overlap)

try:
    fn = SAVEPATH+'/found_sources_not_anticorrelated_'+save_name+'.pkl'
    found_sources_in = pickle.load(open(fn, 'rb'))
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
                correlation = correlation_match(pGB_injected_dict,found_dict)
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
                correlation = correlation_match(pGB_injected_dict,found_dict)
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

    fn = SAVEPATH+'/found_sources_not_anticorrelated_'+save_name+'.pkl'
    pickle.dump(found_sources_in, open(fn, "wb"))
    # np.save(SAVEPATH+'/found_sources_not_anticorrelated_'+save_name+'.npy', found_sources_not_anitcorrelated2)
    found_sources_in = found_sources_not_anitcorrelated2

found_sources_in_flat = np.concatenate(found_sources_in)

found_sources_in_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_in_flat]) for attribute in found_sources_in_flat[0].keys()}
found_sources_in_flat_df = pd.DataFrame(found_sources_in_flat_array)

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
    match_list = []
    pGB_best_list = []
    match_best_list = []
    for j in range(len(found_sources_in)):
        found_match = False
        match_list_one_found_signal = []
        found_dict = {}
        for parameter in parameters:
            found_dict[parameter] = found_sources_in[j][parameter]
        for k in range(len(pGB_injected_not_matched)):
            pGB_injected_dict = {}
            for parameter in parameters:
                pGB_injected_dict[parameter] = pGB_injected_not_matched[k][parameter]
            # correlation = l2_norm_match(pGB_injected_dict,found_dict)
            # correlation = correlation_match(pGB_injected_dict,found_dict)
            SNR_scaled = SNR_match_scaled(pGB_injected_dict,found_dict)
            # correlation, amplitude_factor, cross_correlation = SNR_match_amplitude_condsiered(pGB_injected_dict,found_dict)
            match_list_one_found_signal.append(SNR_scaled)
            if k > 39:
                break
        if 0 == len(match_list_one_found_signal):
            break
        try:
            best_index = np.nanargmin(match_list_one_found_signal)
        except:
            print('all NAN:', match_list_one_found_signal, found_sources_in[0], pGB_injected_not_matched, found_sources_in)
            break
        if match_list_one_found_signal[best_index] < 0.3:
            found_match = True
        if found_match:
            pGB_injected_matched.append(pGB_injected_not_matched[best_index])
            found_sources_matched.append(found_sources_in[j])
            # pGB_injected_not_matched = np.delete(pGB_injected_not_matched, best_index)
            found_sources_not_matched[j] = None
        match_list.append(match_list_one_found_signal[best_index])
        pGB_best_list.append(pGB_injected_not_matched[best_index])
        match_best_list.append(found_sources_in[j])
    return found_sources_in, pGB_injected_not_matched, match_list, pGB_best_list, match_best_list, found_sources_not_matched, pGB_injected_matched, found_sources_matched

do_match_parallelized = True
if do_match_parallelized:
    pGB_injected_matched = []
    found_sources_matched = []
    match_list = []
    pGB_best_list = []
    match_best_list = []
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
        match_list.append([])
        pGB_best_list.append([])
        match_best_list.append([])
        amplitude_factor.append([])
        found_sources_in[i], pGB_injected_not_matched[i], match_list[i],  pGB_best_list[i], match_best_list[i], found_sources_not_matched[i], pGB_injected_matched[i], found_sources_matched[i] = matches[i]

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
    match_list = []
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
            match_list_one_found_signal = []
            match_list_one_found_signal_original = []
            match_list_one_found_signal_amplitude = []
            found_dict = {}
            for parameter in parameters:
                found_dict[parameter] = found_sources_in[i][j][parameter]
            for k in range(len(pGB_injected_not_matched[i])):
                pGB_injected_dict = {}
                for parameter in parameters:
                    pGB_injected_dict[parameter] = pGB_injected_not_matched[i][k][parameter]
                l2_norm = l2_norm_match(pGB_injected_dict,found_dict)
                SNR_scaled = SNR_match_scaled(pGB_injected_dict,found_dict)
                correlation_apmlitude, amplitude_factor, cross_correlation = SNR_match_amplitude_condsiered(pGB_injected_dict,found_dict)
                match_list_one_found_signal.append(SNR_scaled)
                match_list_one_found_signal_original.append(np.array(cross_correlation))
                match_list_one_found_signal_amplitude.append(np.array(correlation_apmlitude))
                if k > 39:
                    break
            if 0 == len(match_list_one_found_signal):
                break
            best_index = np.nanargmin(match_list_one_found_signal)
            if match_list_one_found_signal[best_index] < 1:
                found_match = True
            if found_match:
                pGB_injected_matched[-1].append(pGB_injected_not_matched[i][best_index])
                found_sources_matched[-1].append(found_sources_in[i][j])
                pGB_injected_not_matched[i] = np.delete(pGB_injected_not_matched[i], best_index)
                match_list.append(match_list_one_found_signal[best_index])
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


end_string = '_SNR_scaled_03_injected_snr5'
# end_string = '_correlation_injected_09_snr7'
# end_string = ''

try:
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
    np.save(SAVEPATH+'/match_list' +save_name+end_string+'.npy', match_list)
    np.save(SAVEPATH+'/pGB_best_list' +save_name+end_string+'.npy', pGB_best_list)
    np.save(SAVEPATH+'/match_best_list' +save_name+end_string+'.npy', match_best_list)
except:
    found_sources_matched = np.load(SAVEPATH+'/found_sources_matched' +save_name+end_string+'.npy', allow_pickle=True)
    found_sources_not_matched = np.load(SAVEPATH+'/found_sources_not_matched' +save_name+end_string+'.npy', allow_pickle=True)
    pGB_injected_not_matched = np.load(SAVEPATH+'/injected_not_matched_windows' +save_name+end_string+'.npy', allow_pickle=True)
    pGB_injected_matched = np.load(SAVEPATH+'/injected_matched_windows' +save_name+end_string+'.npy', allow_pickle=True)
    match_list = np.load(SAVEPATH+'/match_list' +save_name+end_string+'.npy', allow_pickle=True)
    pGB_best_list = np.load(SAVEPATH+'/pGB_best_list' +save_name+end_string+'.npy', allow_pickle=True)
    match_best_list = np.load(SAVEPATH+'/match_best_list' +save_name+end_string+'.npy', allow_pickle=True)

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


# metric_threshold = 0.5
# end_string = '_SNR_scaled_'+str(int(metric_threshold))+'_injected_snr5'

# match_list= np.asarray(match_list)
# for i in range(len(match_list)):
#     match_list[i] = np.asarray(match_list[i])

# match_list_flat = np.concatenate(match_list)
# print(len(match_list_flat))

# for i in range(len(match_list)):
#     if len(match_list[i]) > 0:
#         # indexes = np.asarray(match_list[i])<metric_threshold
#         found_sources_matched[i] = found_sources_matched[i][np.asarray(match_list[i])<metric_threshold]
#         # match_list[i] = match_list[i][np.asarray(match_list[i])<metric_threshold]

# match_list_flat = np.concatenate(match_list)
# print(len(match_list_flat))

number_of_injected_signals = 0
for i in range(len(pGB_injected)):
    number_of_injected_signals += len(pGB_injected[i])
number_of_injected_signals_SNR_high = 0
for i in range(len(pGB_injected)):
    number_of_injected_signals_SNR_high += len(pGB_injected[i][pGB_injected[i]['IntrinsicSNR']>10])
number_of_injected_signals_SNR_high2 = len(pGB_injected_flat[pGB_injected_flat['IntrinsicSNR']>10])
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

if version == '1':
    for i in range(len(found_sources_matched)):
        for j in range(len(found_sources_matched[i])):
            a = deepcopy(found_sources_matched[i][j])
            a['InitialPhase'] = -(a['InitialPhase'])
            found_sources_matched[i][j] = deepcopy(a)
    for i in range(len(found_sources_not_matched)):
        for j in range(len(found_sources_not_matched[i])):
            a = deepcopy(found_sources_not_matched[i][j])
            a['InitialPhase'] = -(a['InitialPhase'])
            found_sources_not_matched[i][j] = deepcopy(a)
    for i in range(len(pGB_injected_matched)):
        for j in range(len(pGB_injected_matched[i])):
            a = deepcopy(pGB_injected_matched[i][j])
            a['InitialPhase'] = -(a['InitialPhase'])
            pGB_injected_matched[i][j] = deepcopy(a)


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
# pGB_injected_flat_highSNR_df = pGB_injected_flat_reduced_df[pGB_injected_flat_reduced_df['IntrinsicSNR']>10]

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


index =433
search1 = Search(tdi_fs,Tobs, frequencies_search[index][0], frequencies_search[index][1], update_noise=False)
print(search1.SNR(found_sources_not_matched[index]))
search1.plot(found_sources_in= found_sources_not_matched[index], pGB_injected= pGB_injected_not_matched[index])

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
#     #         Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_matched[i][n], oversample=4)
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
index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.004)+15
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.01)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.0037305)
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
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.0268)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.04296)-1
# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.05994535926230327)-1

# index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  pGB_injected_not_matched_flat_df.iloc[-1]['Frequency'])-1
# index_of_interest_to_plot = 4710
# + shift 4671
# index_of_interest_to_plot = 3
# index_of_interest_to_plot = 128
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
    search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i+number_of_windows-1][1], update_noise=False)
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
    for j in range(number_of_windows):
        for k in range(len(pGB_injected_matched[i+j])):
            matched_extended.append({})
            for parameter in parameters:
                matched_extended[-1][parameter] = pGB_injected_matched[i+j][k][parameter]
    
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
        search1.plot(found_sources_in=found_extended, found_sources_not_matched = found_not_matched_extended, pGB_injected= pGB_injected_dict_list[:20],  pGB_injected_matched= matched_extended, vertical_lines= vertical_lines, saving_label =save_name_path) 
    else:
        search1.plot(found_sources_in=found_extended, found_sources_not_matched = found_not_matched_extended, pGB_injected= pGB_injected_dict_list,  pGB_injected_matched= matched_extended, vertical_lines= vertical_lines, saving_label =save_name_path) 
        # search1.plot(found_sources_in=found_sources_mp_best[i], pGB_injected=pGB_injected[i][:10], pGB_injected_matched= matched_extended, saving_label =SAVEPATH+'/strain added'+ str(int(np.round(frequencies_search[i][0]*10**8))) +save_name+'in.png') 

print(correlation_match(pGB_injected_dict_list[0],found_not_matched_extended[0]))
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
    Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_to_subtract[j], oversample=4)
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
        Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_out_optimize[j], oversample=4)
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
        Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources[j], oversample=4)
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
        Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_to_subtract[i], oversample=4)
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
                # print('SNR', correlation_match(pGB_injected_not_matched[i][k],found_sources_in[i][j]),'parameter comparison:',pGB_injected_not_matched[i][k]['EclipticLatitude'],found_sources_in[i][j]['EclipticLatitude'],eclipticlongitude, found_sources_in[i][j]['EclipticLongitude'])
                correlation = correlation_match(pGB_injected_dict_list[0],found_dict)
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
#     Xs, Ys, Zs = search1.GB.get_fd_tdixyz(template=found_sources_mp[i][3][0], oversample=4)
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
postitive_fd_mask = pGB_injected_matched_flat_df['FrequencyDerivative'] >= 0
plt.plot(pGB_injected_matched_flat_df['GalacticLatitude'][postitive_fd_mask],pGB_injected_matched_flat_df['Distance'][postitive_fd_mask],'g.', label = 'Injected', markersize= 1, zorder=1)
postitive_fd_mask = pGB_injected_not_matched_flat_df['FrequencyDerivative'] >= 0
plt.plot(pGB_injected_not_matched_flat_df['GalacticLatitude'][postitive_fd_mask],pGB_injected_not_matched_flat_df['Distance'][postitive_fd_mask],'g.', label = 'Injected', markersize= 1, zorder=1)
postitive_fd_mask = pGB_injected_not_matched_flat_df['FrequencyDerivative'] >= 0
plt.plot(pGB_injected_not_matched_flat_df['GalacticLatitude'][postitive_fd_mask],pGB_injected_not_matched_flat_df['Distance'][postitive_fd_mask],'+', label = 'not matched', color = 'r', markersize=2, zorder= 1)
# postitive_fd_mask = pGB_injected_flat_highSNR_df['FrequencyDerivative'] >= 0
# plt.plot(pGB_injected_flat_highSNR_df['GalacticLatitude'][postitive_fd_mask],pGB_injected_flat_highSNR_df['Distance'][postitive_fd_mask],'+', label = 'injected SNR>10', color = 'r', markersize=2, zorder= 4)
plt.xlabel('GalacticLatitude')
plt.ylabel('Distance [kpc]')
plt.legend(markerscale=4, loc = 'upper right')
plt.savefig(SAVEPATH+'/galacticLatitudeDistance'+save_name)
plt.show()

fig = plt.figure()
postitive_fd_mask = found_sources_matched_flat_df['FrequencyDerivative'] >= 0
plt.plot(found_sources_matched_flat_df['GalacticLongitude'][postitive_fd_mask],found_sources_matched_flat_df['Distance'][postitive_fd_mask],'.', label = 'found', markersize=1)
postitive_fd_mask = pGB_injected_flat_df['FrequencyDerivative'] >= 0
plt.plot(pGB_injected_flat_df['GalacticLongitude'][postitive_fd_mask],pGB_injected_flat_df['Distance'][postitive_fd_mask],'g.', label = 'Injected', markersize= 1, zorder= 1)
postitive_fd_mask = pGB_injected_not_matched_flat_df['FrequencyDerivative'] >= 0
plt.plot(pGB_injected_not_matched_flat_df['GalacticLongitude'][postitive_fd_mask],pGB_injected_not_matched_flat_df['Distance'][postitive_fd_mask],'+', color = 'r', label = 'not matched', markersize=2, zorder= 1)
# postitive_fd_mask = pGB_injected_flat_highSNR_df['FrequencyDerivative'] >= 0
# plt.plot(pGB_injected_flat_highSNR_df['GalacticLongitude'][postitive_fd_mask],pGB_injected_flat_highSNR_df['Distance'][postitive_fd_mask],'+', label = 'injected SNR>10', color = 'r', markersize=2, zorder= 4)
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
n_bins= 50
match_list= np.asarray(match_list)
match_list_flat = np.concatenate(match_list)
fig = plt.figure()
plt.hist(match_list_flat, bins= np.linspace(0,0.3,n_bins), histtype='step',)
plt.xlabel('match metric')
plt.ylabel('Count')
# plt.yscale('log')
plt.ylim(0,2000)
# plt.xlim(0,1)
plt.savefig(SAVEPATH+'/Evaluation/match'+save_name+end_string,dpi=300,bbox_inches='tight')
plt.show()
