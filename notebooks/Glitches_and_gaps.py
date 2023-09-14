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
from ldc.common.tools import compute_tdi_snr, window

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

dataset = 'Radler'
# dataset = 'Sangria'
# dataset = 'Spritz'
if dataset == 'Radler':
    DATAPATH = grandparent+"/LDC/Radler/data"
    SAVEPATH = grandparent+"/LDC/gaps/"
elif dataset == 'Sangria':
    DATAPATH = grandparent+"/LDC/Sangria/data"
    SAVEPATH = grandparent+"/LDC/pictures/Sangria/"
elif dataset == 'Spritz':
    DATAPATH = grandparent+"/LDC/Spritz/data"
    SAVEPATH = grandparent+"/LDC/Spritz/evaluation"

if dataset == 'Radler':
    # data_fn = DATAPATH + "/LDC1-4_GB_v2.hdf5"
    # data_fn = DATAPATH + "/LDC1-3_VGB_v2_FD_noiseless.hdf5"
    data_fn = DATAPATH + "/LDC1-3_VGB_v2.hdf5"
elif dataset == 'Sangria':
    data_fn = DATAPATH + "/LDC2_sangria_training_v2.h5"
elif dataset == 'Spritz':
    data_fn = DATAPATH + "/LDC2_spritz_vgb_training_v2.h5"
fid = h5py.File(data_fn)

DATAPATH_spritz = grandparent+"/LDC/Spritz/data"
data_fn_spritz = DATAPATH_spritz + "/LDC2_spritz_vgb_training_v2.h5"
fid_spritz = h5py.File(data_fn_spritz)

def print_attrs(name, obj):
    shift = name.count('/') * '    '
    print(shift + name)
    for key, val in obj.attrs.items():
        print(shift + '    ' + f"{key}: {val}")

# fid_spritz.visititems(print_attrs)

reduction = 1

# get TDI 
if dataset == 'Radler':
    td = np.array(fid["H5LISA/PreProcess/TDIdata"])
    td = np.rec.fromarrays(list(td.T), names=["t", "X", "Y", "Z"])
    dt = float(np.array(fid['H5LISA/GWSources/GalBinaries']['Cadence']))
    Tobs = float(int(np.array(fid['H5LISA/GWSources/GalBinaries']['ObservationDuration']))/reduction)
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
    dt = td["t"][1]-td["t"][0]
    Tobs = float(td['t'][-1]/reduction) + dt

td_obs = fid_spritz["obs/tdi"][()]
td_obs = np.rec.fromarrays(list(td_obs.T), names=["t", "X", "Y", "Z"])
td_obs = td_obs['t']
td_clean = fid_spritz["clean/tdi"][()]
td_clean = np.rec.fromarrays(list(td_clean.T), names=["t", "X", "Y", "Z"])
td_clean = td_clean['t']
td_sky = fid_spritz["sky/tdi"][()]
td_sky = np.rec.fromarrays(list(td_sky.T), names=["t", "X", "Y", "Z"])
td_sky = td_sky['t']
td_noisefree = fid_spritz["noisefree/tdi"][()]
td_noisefree = np.rec.fromarrays(list(td_noisefree.T), names=["t", "X", "Y", "Z"])
td_noisefree = td_noisefree['t']
td_galaxy = fid_spritz["gal/tdi"][()]
td_galaxy = np.rec.fromarrays(list(td_galaxy.T), names=["t", "X", "Y", "Z"])
td_galaxy = td_galaxy['t']

tdi = fid_spritz["obs/tdi"][()].squeeze()
tdi_nf = fid_spritz["noisefree/tdi"][()].squeeze()
dt_spritz = tdi['t'][1]-tdi['t'][0]

# Build timeseries and frequencyseries object for X,Y,Z
tdi_ts = dict([(k, TimeSeries(td[k][:int(len(td[k][:])/reduction)], dt=dt, t0=td.t[0])) for k in ["X", "Y", "Z"]])
tdi_fs = xr.Dataset(dict([(k, tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))

tdi_ts_obs = dict([(k, TimeSeries(td_obs[k][:int(len(td_obs[k][:])/reduction)], dt=dt_spritz, t0=td_obs.t[0])) for k in ["X", "Y", "Z"]])
tdi_fs_obs = xr.Dataset(dict([(k, tdi_ts_obs[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
tdi_ts_clean = dict([(k, TimeSeries(td_clean[k][:int(len(td_clean[k][:])/reduction)], dt=dt_spritz, t0=td_clean.t[0])) for k in ["X", "Y", "Z"]])
tdi_fs_clean = xr.Dataset(dict([(k, tdi_ts_clean[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
# tdi_ts_noisefree = dict([(k, TimeSeries(td_noisefree[k][:int(len(td_noisefree[k][:])/reduction)], dt=dt_spritz, t0=td_noisefree.t[0])) for k in ["X", "Y", "Z"]])
# tdi_fs_noisefree = xr.Dataset(dict([(k, tdi_ts_noisefree[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
tdi_ts_sky = dict([(k, TimeSeries(td_sky[k][:int(len(td_sky[k][:])/reduction)], dt=dt_spritz, t0=td_sky.t[0])) for k in ["X", "Y", "Z"]])
tdi_fs_sky = xr.Dataset(dict([(k, tdi_ts_sky[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
tdi_ts_galaxy = dict([(k, TimeSeries(td_galaxy[k][:int(len(td_galaxy[k][:])/reduction)], dt=dt_spritz, t0=td_galaxy.t[0])) for k in ["X", "Y", "Z"]])
tdi_fs_galaxy = xr.Dataset(dict([(k, tdi_ts_galaxy[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))

tdi_ts_glitches = deepcopy(tdi_ts_obs)
for k in ["X", "Y", "Z"]:
    # tdi_ts_glitches[k].values = tdi_ts_sky[k].values - tdi_ts_noisefree[k].values - tdi_ts_galaxy[k].values
    tdi_ts_glitches[k].values = tdi_ts_obs[k].values - tdi_ts_clean[k].values - tdi_ts_galaxy[k].values


## add glitches to tdi
tdi_ts_with_glitches = deepcopy(tdi_ts)
if dataset != 'Spritz':
    for k in ["X", "Y", "Z"]:
        tdi_ts_with_glitches[k].values = tdi_ts[k].values + tdi_ts_glitches[k].values[:len(tdi_ts[k].values)]

plt.figure()
plt.plot(tdi_ts['X'].t, tdi_ts['X'])
# plt.plot(tdi_ts_obs['X'].t, tdi_ts_obs['X'])
# plt.plot(tdi_ts_clean['X'].t, tdi_ts_clean["X"])
# plt.plot(tdi_ts_noisefree['X'].t, tdi_ts_noisefree["X"])
# plt.plot(tdi_ts_sky['X'].t, tdi_ts_sky["X"])
# plt.plot(tdi_ts_galaxy['X'].t, tdi_ts_galaxy["X"])
plt.plot(tdi_ts_glitches['X'].t, tdi_ts_glitches["X"])
plt.plot(tdi_ts_with_glitches['X'].t, tdi_ts_with_glitches["X"])
plt.show()

GB = fastGB.FastGB(delta_t=dt, T=Tobs)  # in seconds

plt.figure()
plt.plot(tdi_ts_with_glitches['X'].t, tdi_ts_with_glitches['X'])
plt.plot(tdi_ts_with_glitches['Y'].t, tdi_ts_with_glitches['Y'])
plt.plot(tdi_ts_with_glitches['Z'].t, tdi_ts_with_glitches['Z'])
# plt.plot(tdi_ts['X'].t, tdi_ts['X'])
plt.show()

start = 200

gaps = {}
for k in ["X", "Y", "Z"]:
    gap = np.isnan(tdi_ts_with_glitches[k])
    tdi_ts_with_glitches[k][gap] = 0
    gaps[k] = tdi_ts_with_glitches[k] == 0
    # gaps = np.isnan(tdi_ts_with_glitches[k])
    tdi_ts_with_glitches[k][gaps[k]] = 0

    mad = scipy.stats.median_abs_deviation(tdi_ts_with_glitches[k])
    peaks, properties = scipy.signal.find_peaks(np.abs(tdi_ts_with_glitches[k]), height=10*mad, threshold=None, distance=1)
    # Turning glitches into gaps
    for pk in peaks:
        tdi_ts_with_glitches[k][pk-10:pk+10] = 0.0

    mad = scipy.stats.median_abs_deviation(tdi_ts_with_glitches[k])
    peaks, properties = scipy.signal.find_peaks(np.abs(tdi_ts_with_glitches[k]), height=10*mad, threshold=None, distance=1)
    # Turning glitches into gaps
    for pk in peaks:
        tdi_ts_with_glitches[k][pk-10:pk+10] = 0.0

    mad = scipy.stats.median_abs_deviation(tdi_ts_with_glitches[k])
    peaks, properties = scipy.signal.find_peaks(np.abs(tdi_ts_with_glitches[k]), height=10*mad, threshold=None, distance=1)
    # Turning glitches into gaps
    for pk in peaks:
        tdi_ts_with_glitches[k][pk-10:pk+10] = 0.0


for k in ["X", "Y", "Z"]:
    tdi_ts_with_glitches[k][:300] = 0

plt.figure()
plt.plot(tdi_ts_with_glitches['X'].t, tdi_ts_with_glitches['X'])
plt.plot(tdi_ts_with_glitches['Y'].t, tdi_ts_with_glitches['Y'])
plt.plot(tdi_ts_with_glitches['Z'].t, tdi_ts_with_glitches['Z'])
# plt.plot(tdi_ts['X'].t, tdi_ts['X'])
plt.show()

### fill gaps
for k in ["X", "Y", "Z"]:
    groups = []
    gaps = tdi_ts_with_glitches[k] == 0
    gaps = tdi_ts_with_glitches[k][gaps]
    start_points = []
    if gaps.t[0].values == tdi_ts_with_glitches[k].t[0].values:
        start_points = [0]
    differences = gaps.t[1:].values - gaps.t[:-1].values
    jumps = differences > 15
    end_points = list(gaps[:-1][jumps].t.values) + list([gaps.t[-1].values])
    start_points = list(start_points) + list(gaps[1:][jumps].t.values)

    for i in range(len(start_points)):
        index_start = int((start_points[i]-tdi_ts_with_glitches[k].t[0].values)/dt)
        index_end = int((end_points[i]-tdi_ts_with_glitches[k].t[0])/dt)
        length_gap = len(tdi_ts_with_glitches[k][index_start:index_end])+1
        tdi_ts_with_glitches[k][index_start:index_start+length_gap] = tdi_ts_with_glitches[k][index_end+1:index_end+1+length_gap].values
        # else:
        #     tdi_ts_with_glitches[k][index_start:index_start+int(length_gap/2)] = tdi_ts_with_glitches[k][index_start-int(length_gap/2):index_start].values
        #     tdi_ts_with_glitches[k][index_end-int(length_gap/2):index_end] = tdi_ts_with_glitches[k][index_end:index_end+int(length_gap/2)].values
    
plt.figure()
plt.plot(tdi_ts_with_glitches['X'].t, tdi_ts_with_glitches['X'])
plt.plot(tdi_ts_with_glitches['Y'].t, tdi_ts_with_glitches['Y'])
plt.plot(tdi_ts_with_glitches['Z'].t, tdi_ts_with_glitches['Z'])
# plt.plot(tdi_ts['X'].t, tdi_ts['X'])
plt.show()

tdi_ts_signals = deepcopy(tdi_ts)
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
# tdi_fs_noisefree = xr.Dataset(dict([(k, tdi_ts_noisefree[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
tdi_fs_sky = xr.Dataset(dict([(k, tdi_ts_sky[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
tdi_fs_glitches = xr.Dataset(dict([(k, tdi_ts_glitches[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
tdi_fs_with_glitches = xr.Dataset(dict([(k, tdi_ts_with_glitches[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
# tdi_fs_signals = xr.Dataset(dict([(k, tdi_ts_signals[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))


i = 0
lower_frequency = cat[i]['Frequency']-0.0000005
upper_frequency = cat[i]['Frequency']+0.0000005
search = Search(tdi_fs,Tobs, lower_frequency, upper_frequency)

search.plot(found_sources_in= [cat[i]])
print(search.SNR([cat[i]]))
print(search.SNR_scaled([cat[i]]))

Xs, Ys, Zs = GB.get_fd_tdixyz(template=cat[i], oversample=4)
Af = (Zs - Xs)/np.sqrt(2.0)
Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)   
plt.figure()
plt.semilogx(Xs.f*1000,np.real(Xs.values))
plt.semilogx(tdi_fs_with_glitches['X'].f.values*1000,np.real(tdi_fs_with_glitches['X'].values))
plt.semilogx(tdi_fs['X'].f.values*1000,np.real(tdi_fs['X'].values))
# plt.semilogx(tdi_fs_clean['X'].f.values*1000,np.real(tdi_fs_clean['X'].values))
# plt.semilogx(tdi_fs_noisefree['X'].f.values*1000,np.real(tdi_fs_noisefree['X'].values))
plt.semilogx(tdi_fs_sky['X'].f.values*1000,np.real(tdi_fs_sky['X'].values))
plt.xlim(lower_frequency*1000,upper_frequency*1000)
# plt.semilogx(Ef.f*1000,Ef.values)
plt.show()

plt.figure()
plt.plot(tdi_ts_with_glitches['X'].t, tdi_ts_with_glitches['X'])
plt.plot(tdi_ts['X'].t, tdi_ts['X'])
plt.show()

print('end')