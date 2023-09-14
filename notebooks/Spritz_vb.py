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
import pickle
sys.path.append('/cluster/home/sstrub/Repositories/LDC/lib/lib64/python3.8/site-packages/ldc-0.1-py3.8-linux-x86_64.egg')


from ldc.lisa.noise import get_noise_model
from ldc.common.series import TimeSeries
import ldc.waveform.fastGB as fastGB
from ldc.common.tools import compute_tdi_snr

from sources import *

# customized settings
plot_parameter = {  # 'backend': 'ps',
    "font.family": "DeJavu Serif",
    # "font.serif": "Times",
    "font.serif" : ["Computer Modern Serif"],
    "font.size": 16,
    "mathtext.fontset": "cm",
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

# dataset = 'Radler'
# dataset = 'Sangria'
dataset = 'Spritz'
if dataset == 'Radler':
    DATAPATH = grandparent+"/LDC/Radler/data"
    SAVEPATH = grandparent+"/LDC/pictures/LDC1-4/"
elif dataset == 'Sangria':
    DATAPATH = grandparent+"/LDC/Sangria/data"
    SAVEPATH = grandparent+"/LDC/pictures/Sangria/"
elif dataset == 'Spritz':
    DATAPATH = grandparent+"/LDC/Spritz/data"
    SAVEPATH = grandparent+"/LDC/Spritz/evaluation"

if dataset == 'Radler':
    data_fn = DATAPATH + "/LDC1-4_GB_v2.hdf5"
    # data_fn = DATAPATH + "/LDC1-3_VGB_v2_FD_noiseless.hdf5"
    # data_fn = DATAPATH + "/LDC1-3_VGB_v2.hdf5"
elif dataset == 'Sangria':
    data_fn = DATAPATH + "/LDC2_sangria_training_v2.h5"
elif dataset == 'Spritz':
    data_fn = DATAPATH + "/LDC2_spritz_vgb_training_v2.h5"
fid = h5py.File(data_fn)

def print_attrs(name, obj):
    shift = name.count('/') * '    '
    print(shift + name)
    for key, val in obj.attrs.items():
        print(shift + '    ' + f"{key}: {val}")

# fid.visititems(print_attrs)

reduction = 1

# get TDI 
if dataset == 'Radler':
    td = np.array(fid["H5LISA/PreProcess/TDIdata"])
    td = np.rec.fromarrays(list(td.T), names=["t", "X", "Y", "Z"])
    dt = float(np.array(fid['H5LISA/GWSources/GalBinaries']['Cadence']))
    Tobs = float(int(np.array(fid['H5LISA/GWSources/GalBinaries']['ObservationDuration']))/reduction)
elif dataset == 'Sangria':
    td = fid["obs/tdi"][()]
    td = np.rec.fromarrays(list(td.T), names=["t", "X", "Y", "Z"])
    td = td['t']
    dt = td["t"][1]-td["t"][0]
    
    td_mbhb = fid["sky/mbhb/tdi"][()]
    # cat_mbhb = fid["sky/mbhb/cat"]
    td_mbhb  = np.rec.fromarrays(list(td_mbhb .T), names=["t", "X", "Y", "Z"])
    td_mbhb  = td_mbhb ['t']

    Tobs = float(int(np.array(fid['obs/config/t_max']))/reduction)
    Tobs = float(td['t'][-1]/reduction) + dt
    for k in ["X", "Y", "Z"]:
        td[k] = td[k] - td_mbhb[k]

elif dataset == 'Spritz':
    names = fid["sky/cat"].dtype.names
    cat_vgb = dict(zip(names, [fid["sky/cat"][name] for name in names]))
    cat = []
    for i in range(len(cat_vgb['Frequency'])):
        cat.append({})
        for parameter in parameters:
            cat[i][parameter] = cat_vgb[parameter][i][0]

    # for i in range(len(cat_vgb['Frequency'])):
    #     cat[i]['InitialPhase'] *= -1

    # print(cat_vgb)
    td = fid["obs/tdi"][()]
    td = np.rec.fromarrays(list(td.T), names=["t", "X", "Y", "Z"])
    td = td['t']
    td_clean = fid["clean/tdi"][()]
    td_clean = np.rec.fromarrays(list(td_clean.T), names=["t", "X", "Y", "Z"])
    td_clean = td_clean['t']
    td_sky = fid["sky/tdi"][()]
    td_sky = np.rec.fromarrays(list(td_sky.T), names=["t", "X", "Y", "Z"])
    td_sky = td_sky['t']
    td_noisefree = fid["noisefree/tdi"][()]
    td_noisefree = np.rec.fromarrays(list(td_noisefree.T), names=["t", "X", "Y", "Z"])
    td_noisefree = td_noisefree['t']
    dt = td["t"][1]-td["t"][0]
    Tobs = float(td['t'][-1]/reduction) + dt
    # Tobs = float(int(np.array(fid['obs/config/t_max']))/reduction) * 31536000.0


# Build timeseries and frequencyseries object for X,Y,Z
tdi_ts = dict([(k, TimeSeries(td[k][:int(len(td[k][:])/reduction)], dt=dt, t0=td.t[0])) for k in ["X", "Y", "Z"]])
gaps = np.isnan(tdi_ts['X'])
tdi_ts['X'][gaps] = 0
tdi_fs = xr.Dataset(dict([(k, tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))

tdi_ts_clean = dict([(k, TimeSeries(td_clean[k][:int(len(td_clean[k][:])/reduction)], dt=dt, t0=td_clean.t[0])) for k in ["X", "Y", "Z"]])
tdi_fs_clean = xr.Dataset(dict([(k, tdi_ts_clean[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
tdi_ts_noisefree = dict([(k, TimeSeries(td_noisefree[k][:int(len(td_noisefree[k][:])/reduction)], dt=dt, t0=td_noisefree.t[0])) for k in ["X", "Y", "Z"]])
tdi_fs_noisefree = xr.Dataset(dict([(k, tdi_ts_noisefree[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
tdi_ts_sky = dict([(k, TimeSeries(td_sky[k][:int(len(td_sky[k][:])/reduction)], dt=dt, t0=td_sky.t[0])) for k in ["X", "Y", "Z"]])
tdi_fs_sky = xr.Dataset(dict([(k, tdi_ts_sky[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))

GB = fastGB.FastGB(delta_t=dt, T=Tobs)  # in seconds

start = 200
gaps = {}
tdi_ts_clean_gaps = deepcopy(tdi_ts_clean)
if dataset == 'Spritz':
    for k in ["X", "Y", "Z"]:
        gap = np.isnan(tdi_ts[k])
        tdi_ts[k][gap] = 0
        gaps[k] = tdi_ts[k] == 0
        tdi_ts_clean_gaps[k][gaps[k]] = 0
        # gaps = np.isnan(tdi_ts_noisefree[k])
        tdi_ts_noisefree[k][gaps[k]] = 0

        mad = scipy.stats.median_abs_deviation(tdi_ts[k])
        peaks, properties = scipy.signal.find_peaks(np.abs(tdi_ts[k]), height=10*mad, threshold=None, distance=1)
        # Turning glitches into gaps
        for pk in peaks:
            tdi_ts[k][pk-10:pk+10] = 0.0

        mad = scipy.stats.median_abs_deviation(tdi_ts_noisefree[k])
        peaks, properties = scipy.signal.find_peaks(np.abs(tdi_ts_noisefree[k]), height=10*mad, threshold=None, distance=1)
        # Turning glitches into gaps
        for pk in peaks:
            tdi_ts_noisefree[k][pk-10:pk+10] = 0.0

    tdi_ts_glitches = deepcopy(tdi_ts)
    tdi_ts_signals = deepcopy(tdi_ts)
    for k in ["X", "Y", "Z"]:
        tdi_ts_glitches[k].values = tdi_ts_sky[k].values - tdi_ts_noisefree[k].values
        # tdi_ts_signals[k].values = tdi_ts_noisefree[k].values - tdi_ts_glitches[k].values

    
    plt.figure(figsize=(16, 6))
    plt.subplot(211)
    # plt.plot(tdi_ts['X'].t[start:], tdi_ts['X'][start:])
    # plt.plot(tdi_ts_clean['X'].t[start:], tdi_ts_clean["X"][start:])
    # plt.plot(tdi_ts_clean_gaps['X'].t[start:], tdi_ts_clean_gaps["X"][start:])
    # plt.plot(tdi_ts_noisefree['X'].t[start:], tdi_ts_noisefree["X"][start:])
    plt.plot(tdi_ts_sky['X'].t[start:], tdi_ts_sky["X"][start:])
    plt.plot(tdi_ts_glitches['X'].t[start:], tdi_ts_glitches["X"][start:])
    # plt.plot(tdi_ts_signals['X'].t[start:], tdi_ts_signals["X"][start:])
    plt.ylabel("TDI X")
    plt.xlabel("Time [s]")
    plt.xlim([1.8e7, 2.2e7])
    # plt.ylim([-1e-18, 1e-18])
    plt.show()
    
    
    tdi_fs = xr.Dataset(dict([(k, tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
    tdi_fs_noisefree = xr.Dataset(dict([(k, tdi_ts_noisefree[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
    tdi_fs_sky = xr.Dataset(dict([(k, tdi_ts_sky[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
    tdi_fs_glitches = xr.Dataset(dict([(k, tdi_ts_glitches[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
    # tdi_fs_signals = xr.Dataset(dict([(k, tdi_ts_signals[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))

# cat = [{'Amplitude': 6.058847755160548e-23, 'EclipticLatitude': 1.0725776721181566, 'EclipticLongitude': 5.204804002612463, 'Frequency': 0.0018421295017039697, 'FrequencyDerivative': 3.856131474449337e-18, 'Inclination': 0.2617993877991494, 'InitialPhase': 3.1429922509263313, 'Polarization': 3.741917727527864}]

i = 0
lower_frequency = cat[i]['Frequency']-0.0000005
upper_frequency = cat[i]['Frequency']+0.0000005
search = Search(tdi_fs_sky,Tobs, lower_frequency, upper_frequency)

search.plot(found_sources_in= [cat[i]])
print(search.SNR([cat[i]]))
print(search.SNR_scaled([cat[i]]))


Xs, Ys, Zs = GB.get_fd_tdixyz(template=cat[i], oversample=4)
Af = (Zs - Xs)/np.sqrt(2.0)
Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)   
plt.figure()
plt.semilogx(Xs.f*1000,np.real(Xs.values)/2.5)
plt.semilogx(tdi_fs['X'].f.values*1000,np.real(tdi_fs['X'].values))
plt.semilogx(tdi_fs_clean['X'].f.values*1000,np.real(tdi_fs_clean['X'].values))
plt.semilogx(tdi_fs_noisefree['X'].f.values*1000,np.real(tdi_fs_noisefree['X'].values))
plt.semilogx(tdi_fs_sky['X'].f.values*1000,np.real(tdi_fs_sky['X'].values))
plt.xlim(lower_frequency*1000,upper_frequency*1000)
# plt.semilogx(Ef.f*1000,Ef.values)
plt.show()

start = 200
tdix_wo_nan = tdi_ts["X"].copy()[start:]
gaps = np.isnan(tdi_ts['X'])[start:]
tdix_wo_nan[gaps] = 0 # set to 0 to compute PSD


plt.figure(figsize=(16, 6))
plt.subplot(211)
plt.plot(tdi_ts_clean['X'].t[start:], tdi_ts_clean["X"][start:])
plt.plot(tdi_ts['X'].t[start:], tdi_ts["X"][start:])
plt.ylabel("TDI X")
plt.xlabel("Time [s]")
plt.xlim([1.8e7, 2.2e7])
plt.ylim([-1e-18, 1e-18])
plt.show()

plt.figure()
# plt.subplot(212)
# plt.figure(figsize=(16, 6))
f, psdX =  scipy.signal.welch(tdix_wo_nan, fs=1.0/dt, window='hann', nperseg=int(1e6 / dt))
plt.loglog(f, np.sqrt(psdX))
plt.loglog(tdi_fs['X'].f, np.abs(tdi_fs['X'].values))
plt.ylabel("$\sqrt{PSD X}$")
plt.xlabel("Freq [Hz]")
plt.show()

# Define a low-pass filter function
def low_pass_filter(y, n=5, fc=0.05, fs=10.0):

    b, a = scipy.signal.butter(n, fc, 'low', analog=False, fs=fs)
    y = scipy.signal.filtfilt(b, a, y)

    return y


 # Low-pass filter the data to better see the glitches
# tdi_x_filtered = low_pass_filter(tdix_wo_nan, n=5, fc=0.05, fs=1/dt)
tdi_x_filtered = deepcopy(tdix_wo_nan)
# Re-add the gaps
tdi_x_filtered[gaps] = 0
tdi_t = tdi_ts['X'].t[start:]

plt.figure(figsize=(16, 6))
i1 = 1000
i2 = 1000000
plt.plot(tdi_t[i1:i2], tdix_wo_nan[i1:i2], label='Unfiltered data')
plt.plot(tdi_t[i1:i2], tdi_x_filtered[i1:i2], label='Filtered data')
plt.ylabel("TDI X")
plt.xlabel("Time [s]")
# plt.xlim([1.8e7, 2.2e7])
# plt.ylim([-1e-18, 1e-18])
# plt.legend()
plt.show()

plt.figure(figsize=(16, 6))
i1 = 1000
i2 = 10000000
f, psdX =  scipy.signal.welch(tdi_x_filtered, fs=1.0/dt, window='hann', nperseg=int(1e6 / dt))
plt.loglog(f, np.sqrt(psdX), label='Filtered data')
f, psdX =  scipy.signal.welch(tdix_wo_nan, fs=1.0/dt, window='hann', nperseg=int(1e6 / dt))
plt.loglog(f, np.sqrt(psdX), label='Unfiltered data')
# plt.plot(tdi_t[i1:i2], tdi_x_filtered[i1:i2], label='Filtered data')
# plt.plot(tdi_t[i1:i2], tdix_wo_nan[i1:i2], label='Unfiltered data')
plt.ylabel("TDI X")
plt.xlabel("Time [s]")
# plt.xlim([1.8e7, 2.2e7])
# plt.ylim([-1e-18, 1e-18])
# plt.legend()
plt.show()

 # Detecting outliers by sigma-clipping
mad = scipy.stats.median_abs_deviation(tdi_x_filtered[i1:6000])
median = np.median(tdi_x_filtered[i1:6000])
maxval = np.max(np.abs(tdi_x_filtered[i1:i2]))
peaks, properties = scipy.signal.find_peaks(np.abs(tdi_x_filtered), height=10*mad, threshold=None, distance=1)

 
plt.figure(figsize=(8, 6))
plt.plot(tdi_t[i1:i2], np.abs(tdi_x_filtered[i1:i2]), label='Filtered data')
plt.plot(tdi_ts_clean['X'].t[start:], tdi_ts_clean["X"][start:]-tdi_ts["X"][start:])
plt.vlines(tdi_t[peaks], ymin=-1.1*maxval, ymax=0, color='red', linestyle='dashed', label='Detected outliers')
plt.ylabel("TDI X")
plt.xlabel("Time [s]")
plt.xlim([tdi_t[i1], tdi_t[i2]])
# plt.ylim([-1e-18, 1e-18])
plt.legend(loc='upper right')
plt.show()

# Turning glitches into gaps
tdix_wo_glitches = np.copy(tdix_wo_nan)
for pk in peaks:
    tdix_wo_glitches[pk-10:pk+10] = 0.0
f, psdX_wo_glitches =  scipy.signal.welch(tdix_wo_glitches, fs=1.0/dt, window='hann', nperseg=int(1e6 / dt))
tdix_wo_glitches = TimeSeries(tdix_wo_glitches[:int(len(tdix_wo_glitches[:])/reduction)], dt=dt, t0=td.t[0])


plt.figure(figsize=(8, 6))
plt.subplot(211)
plt.plot(tdi_t, tdix_wo_glitches, label='Observed data, glitches set to zero', rasterized=True)
# plt.plot(tdi_t, tdix, label='VGB signal', rasterized=True)
plt.ylabel("TDI X")
plt.xlabel("Time [s]")
# plt.xlim([1.8e7, 2.2e7])
# plt.ylim([-1e-18, 1e-18])
plt.legend(loc='upper right')




tdi_noiseless = fid["sky/tdi"][()].squeeze()
tdi_x_noiseless = tdi_noiseless["X"].copy()[start:]

plt.figure(figsize=(16, 6))
plt.subplot(211)
plt.plot(tdi_noiseless['t'][start:], tdi_x_noiseless, color='orange')
plt.ylabel("TDI X")
plt.xlabel("Time [s]")
plt.xlim([1.8e7, 2.2e7])
# plt.ylim([-1e-18, 1e-18])
plt.title("Noiseless time series")

plt.subplot(212)
f, psdX_noiseless =  scipy.signal.welch(tdi_x_noiseless, fs=1.0/dt, window='hann', nperseg=int(1e6 / dt))
plt.loglog(f, np.sqrt(psdX_noiseless), color='orange')
plt.ylabel("$\sqrt{PSD X}$")
plt.xlabel("Freq [Hz]")
plt.show()

tdi_x_noiseless = TimeSeries(tdi_x_noiseless[:int(len(tdi_x_noiseless[:])/reduction)], dt=dt, t0=td.t[0])
tdi_fs_wo_glitches = tdix_wo_glitches.ts.fft(win=window)
tdi_fs_noiseless = tdi_x_noiseless.ts.fft(win=window)

plt.figure()

# plt.subplot(212)
plt.loglog(f, np.sqrt(psdX_wo_glitches), label='Observed data, glitches set to zero', rasterized=True)
plt.loglog(f, np.sqrt(psdX))
plt.loglog(tdi_fs_wo_glitches.f.values, np.abs(tdi_fs_wo_glitches.values))
plt.loglog(tdi_fs_noiseless.f.values, np.abs(tdi_fs_noiseless.values))
f, psdX_noiseless =  scipy.signal.welch(tdi_x_noiseless, fs=1.0/dt, window='hann', nperseg=int(1e6 / dt))
plt.loglog(f, np.sqrt(psdX_noiseless), color='orange')
# plt.plot(tdi_noiseless['t'][start:], np.abs(tdi_x_noiseless), color='orange')
# plt.loglog(f, np.sqrt(psdX_noiseless), label='VGB signal', rasterized=True)
# plt.loglog(fxx[fxx>0], np.sqrt(npsd), label='Theoretical noise PSD (without artefacts)', rasterized=True)
plt.ylabel("$\sqrt{PSD X}$")
plt.xlabel("Freq [Hz]")
plt.xlim([1e-4, 1e-1])
# plt.legend(loc='lower left')
plt.show()


tdi_X_diffg = tdi["X"]-tdi_gf["X"] # with - without galaxy -> galaxy only
tdi_X_diffn = tdi["X"]-tdi_nf["X"] # with - without noise -> noise only
tdi_X_g = tdi_gf_sky["X"]-tdi_gf_nf["X"] # sky - noisefree -> glitch only

plt.figure(figsize=(16,8))
plt.subplot(311)
plt.plot(tdi['t'][start:], tdi["X"][start:], label='total')
plt.axis([None, None, -0.5e-18, 0.5e-18])
plt.ylabel("TDI X")
plt.xlabel("Time [s]")
plt.xlim([1.8e7, 2.2e7])
plt.ylim([-1e-19, 1e-19])
plt.legend()
plt.subplot(312)
plt.plot(tdi['t'][start:], tdi_X_g[start:], label='glitch only')
plt.axis([None, None, -0.5e-18, 0.5e-18])
plt.ylabel("TDI X")
plt.xlabel("Time [s]")
plt.xlim([1.8e7, 2.2e7])
plt.ylim([-2e-31, 2e-31])
plt.legend()
plt.subplot(313)
plt.plot(tdi['t'][start:], tdi_X_diffg[start:], label='galaxy only')
plt.axis([None, None, -0.5e-18, 0.5e-18])
plt.ylabel("TDI X")
plt.xlabel("Time [s]")
plt.xlim([1.8e7, 2.2e7])
plt.ylim([-2e-22, 2e-22])
plt.legend()

print('end')