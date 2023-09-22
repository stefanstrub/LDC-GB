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
from ldc.common.series import TimeSeries, FrequencySeries
import ldc.waveform.fastGB as fastGB

from ldc.common.tools import compute_tdi_snr, window
try:
    import cupy as xp

    gpu_available = True

except (ImportError, ModuleNotFoundError) as e:
    import numpy as xp

    gpu_available = False

from gbgpu.gbgpu import GBGPU

from gbgpu.utils.constants import *

from sources import *

# gpu_available = False
gb_gpu = GBGPU(use_gpu=gpu_available)

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

fill_gaps = False
dataset = 'Radler'
# dataset = 'Sangria'
dataset = 'Spritz'
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
# plt.plot(tdi_ts_with_glitches['X'].t, tdi_ts_with_glitches["X"])
plt.plot(tdi_ts['X'].t, tdi_ts['X'])
# plt.plot(tdi_ts_obs['X'].t, tdi_ts_obs['X'])
# plt.plot(tdi_ts_clean['X'].t, tdi_ts_clean["X"])
# plt.plot(tdi_ts_noisefree['X'].t, tdi_ts_noisefree["X"])
# plt.plot(tdi_ts_sky['X'].t, tdi_ts_sky["X"])
# plt.plot(tdi_ts_galaxy['X'].t, tdi_ts_galaxy["X"])
# plt.plot(tdi_ts_glitches['X'].t, tdi_ts_glitches["X"])
plt.xlabel('time')
plt.ylabel('TDI X')
plt.show()

GB = fastGB.FastGB(delta_t=dt, T=Tobs)  # in seconds

plt.figure()
plt.plot(tdi_ts_with_glitches['X'].t, tdi_ts_with_glitches['X'])
plt.plot(tdi_ts_with_glitches['Y'].t, tdi_ts_with_glitches['Y'])
plt.plot(tdi_ts_with_glitches['Z'].t, tdi_ts_with_glitches['Z'])
# plt.plot(tdi_ts['X'].t, tdi_ts['X'])
plt.xlabel('time')
plt.ylabel('TDI')
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
        tdi_ts_with_glitches['X'][pk-10:pk+10] = 0.0
        tdi_ts_with_glitches['Y'][pk-10:pk+10] = 0.0
        tdi_ts_with_glitches['Z'][pk-10:pk+10] = 0.0

    mad = scipy.stats.median_abs_deviation(tdi_ts_with_glitches[k])
    peaks, properties = scipy.signal.find_peaks(np.abs(tdi_ts_with_glitches[k]), height=10*mad, threshold=None, distance=1)
    # Turning glitches into gaps
    for pk in peaks:
        tdi_ts_with_glitches['X'][pk-10:pk+10] = 0.0
        tdi_ts_with_glitches['Y'][pk-10:pk+10] = 0.0
        tdi_ts_with_glitches['Z'][pk-10:pk+10] = 0.0

    mad = scipy.stats.median_abs_deviation(tdi_ts_with_glitches[k])
    peaks, properties = scipy.signal.find_peaks(np.abs(tdi_ts_with_glitches[k]), height=10*mad, threshold=None, distance=1)
    # Turning glitches into gaps
    for pk in peaks:
        tdi_ts_with_glitches['X'][pk-10:pk+10] = 0.0
        tdi_ts_with_glitches['Y'][pk-10:pk+10] = 0.0
        tdi_ts_with_glitches['Z'][pk-10:pk+10] = 0.0


for k in ["X", "Y", "Z"]:
    tdi_ts_with_glitches[k][:300] = 0

plt.figure()
plt.plot(tdi_ts_with_glitches['X'].t, tdi_ts_with_glitches['X'])
plt.plot(tdi_ts_with_glitches['Y'].t, tdi_ts_with_glitches['Y'])
plt.plot(tdi_ts_with_glitches['Z'].t, tdi_ts_with_glitches['Z'])
# plt.plot(tdi_ts['X'].t, tdi_ts['X'])
plt.xlabel('time')
plt.ylabel('TDI')
plt.show()

if fill_gaps:
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


data_fd= {}
gap_less_duration = {}
for k in ["X", "Y", "Z"]:
    data_fd[k] = []
    gap_less_duration[k] = []
    groups = []
    gaps = tdi_ts_with_glitches[k] == 0
    gaps = tdi_ts_with_glitches[k][gaps]
    start_points = []
    if gaps.t[0].values != tdi_ts_with_glitches[k].t[0].values:
        start_points = [0]
        end_points = gaps.t[0].values
    differences = gaps.t[1:].values - gaps.t[:-1].values
    jumps = differences > 15
    start_points = list(start_points) + list(gaps[:-1][jumps].t.values) + list([gaps.t[-1].values])
    end_points =  list(gaps[1:][jumps].t.values)
    for i in range(len(end_points)):
        index_start = int((start_points[i]-tdi_ts_with_glitches[k].t[0].values)/dt)
        index_end = int((end_points[i]-tdi_ts_with_glitches[k].t[0])/dt)
        data_fd[k].append(tdi_ts_with_glitches[k][index_start:index_end].ts.fft(win=window))
        gap_less_duration[k].append(tdi_ts_with_glitches[k][index_start:index_end].t[-1].values- tdi_ts_with_glitches[k][index_start:index_end].t[0].values)


data_fd['A'] = []
data_fd['E'] = []
for i in range(len(end_points)):
    data_fd['A'].append((data_fd['Z'][i] - data_fd['X'][i])/np.sqrt(2.0))
    data_fd['E'].append((data_fd['Z'][i] - 2.0*data_fd['Y'][i] + data_fd['X'][i])/np.sqrt(6.0))


# tdi_fs = xr.Dataset(dict([(k, tdi_ts[k][index_start:index_end].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
# # tdi_fs_noisefree = xr.Dataset(dict([(k, tdi_ts_noisefree[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
# tdi_fs_sky = xr.Dataset(dict([(k, tdi_ts_sky[k][index_start:index_end].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
# tdi_fs_glitches = xr.Dataset(dict([(k, tdi_ts_glitches[k][index_start:index_end].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
tdi_fs_with_glitches = xr.Dataset(dict([(k, tdi_ts_with_glitches[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
# tdi_fs_signals = xr.Dataset(dict([(k, tdi_ts_signals[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))

data_glitches = {}
data_glitches['A'] = (tdi_fs_with_glitches['Z'] - tdi_fs_with_glitches['X'])/np.sqrt(2.0)
data_glitches['E'] = (tdi_fs_with_glitches['Z'] - 2.0*tdi_fs_with_glitches['Y'] + tdi_fs_with_glitches['X'])/np.sqrt(6.0)

cat_index = 0
lower_frequency = cat[cat_index]['Frequency']-0.0000005
upper_frequency = cat[cat_index]['Frequency']+0.0000005
# search = Search(tdi_fs_with_glitches,Tobs, lower_frequency, upper_frequency)

# search.plot(found_sources_in= [cat[i]])
# print(search.SNR([cat[i]]))
# print(search.SNR_scaled([cat[i]]))
parameters_gpgpu = ['Amplitude','Frequency','FrequencyDerivative','Frequency2Derivative','InitialPhase','Inclination','Polarization','EclipticLongitude','EclipticLatitude']
start = time.time()
Xs, Ys, Zs = GB.get_fd_tdixyz(template=cat[cat_index], oversample=4, tdi2=True)
print('time', time.time()-start)

start = time.time()
X_td, Ys_td, Zs_td = GB.get_td_tdixyz(template=cat[cat_index], oversample=4, tdi2=True)
print('time', time.time()-start)

simulated_td = { 'X': X_td, 'Y': Ys_td, 'Z': Zs_td}
simulated_fd = {}
for k in ["X", "Y", "Z"]:
    simulated_fd[k] = []
    for i in range(len(end_points)):
        index_start = int((start_points[i]-tdi_ts_with_glitches[k].t[0].values)/dt)
        index_end = int((end_points[i]-tdi_ts_with_glitches[k].t[0])/dt)
        simulated_fd[k].append(simulated_td[k][index_start:index_end].ts.fft(win=window))

pGB = {}
for parameter in parameters:
    pGB[parameter] = cat[cat_index][parameter]
pGB['Frequency2Derivative'] = 0
pGB2 = deepcopy(pGB)
pGB2['Amplitude'] = 1e-18
num_bin = 100
params = np.array([np.full(num_bin, pGB[parameter]) for parameter in parameters_gpgpu])
# params = (np.zeros((100,len(params))) + params).T
start = time.time()
gb_gpu.run_wave(*params,
        N=None,
        T=end_points[i]-start_points[i],
        dt=dt,
        t_start=start_points[i],
        oversample=4,
        tdi2=True)
print('time', time.time()-start)

gb_gpu.run_wave(*params,
        N=None,
        T=Tobs,
        dt=dt,
        t_start=0,
        oversample=4,
        tdi2=True)

noise_model =  "SciRDv1"
data_gpu = []
PSD_GPU = []
for i in range(len(data_fd['A'])):
    data_gpu.append(xp.array([data_fd['A'][i].values[1:], data_fd['E'][i].values[1:]]))
    Nmodel = get_noise_model(noise_model, data_fd['A'][i].f[1:])
    SA_full_f = Nmodel.psd(freq=data_fd['A'][i].f[1:], option="A")
    SE_full_f = Nmodel.psd(freq=data_fd['E'][i].f[1:], option="E")
    PSD_GPU.append(xp.array([SA_full_f, SE_full_f]))
gb_gpu.d_d = 0.0
loglikelihood = []
for i in range(len(data_fd['A'])):
    loglikelihood.append(gb_gpu.get_ll(params, data_gpu[i], PSD_GPU[i], N=None, oversample=4, dt=dt, T=end_points[i]-start_points[i], start_freq_ind=0, tdi2=True, t_start=start_points[i]))


start_index = 10**5
end_index = 1
t_start = float(tdi_ts_sky[k].t[start_index].values)
t_end = float(tdi_ts_sky[k].t[-end_index].values)
tdi_fs_sky = xr.Dataset(dict([(k, tdi_ts_sky[k][start_index:-end_index].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))

data_sky = {}
data_sky['A'] = (tdi_fs_sky['Z'] - tdi_fs_sky['X'])/np.sqrt(2.0)
data_sky['E'] = (tdi_fs_sky['Z'] - 2.0*tdi_fs_sky['Y'] + tdi_fs_sky['X'])/np.sqrt(6.0)

N = 256
data_gpu_full = [xp.array(data_sky['A'].values), xp.array(data_sky['E'].values)]
Nmodel = get_noise_model(noise_model, data_sky['A'].f)
SA_full_f = Nmodel.psd(freq=data_sky['A'].f, option="A")
SE_full_f = Nmodel.psd(freq=data_sky['E'].f, option="E")
PSD_GPU = xp.array([SA_full_f, SE_full_f])
loglikelihood = gb_gpu.get_ll(params, data_gpu_full, PSD_GPU, N=N, oversample=4, dt=dt, T=t_end-t_start, start_freq_ind=0, tdi2=True, t_start=t_start)

search = Search(tdi_fs_sky,t_end-t_start, lower_frequency, upper_frequency, tdi2=True, gb_gpu=gb_gpu, t_start=t_start)
search.plot(found_sources_in= [cat[cat_index]])
loglikelihood = search.loglikelihood_gpu(params)


X_fd = X_td[start_index:-end_index].ts.fft(win=window)

gb_gpu.run_wave(*params,
        N=None,
        T=t_end-t_start,
        dt=dt,
        t_start=t_start,
        oversample=4,
        tdi2=True)

t0 = tdi_ts_with_glitches[k].t[0].values
target_length = index_end - index_start
# xr.Dataset(dict([(k, tdi_ts[k][index_start:index_end].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
# X_fd = FrequencySeries(gb_gpu.X.get()[0], fs= gb_gpu.freqs.get()[0])
# X_td = TimeSeries(X_fd.ifft() / dt, t0=t0, dt=dt)

# Y_td = TimeSeries(np.fft.irfft(gb_gpu.Y.get()[0] / dt, n=int(Tobs)), t0=t0, dt=dt)
# Z_td = TimeSeries(np.fft.irfft(gb_gpu.Z.get()[0] / dt, n=int(Tobs)), t0=t0, dt=dt)
# gb_gpu_td = xr.Dataset({'X': X_td, 'Y': Y_td, 'Z': Z_td})
A_td = TimeSeries(np.fft.irfft(gb_gpu.A.get()[0] / dt, n=int(Tobs)), t0=t0, dt=dt)
E_td = TimeSeries(np.fft.irfft(gb_gpu.E.get()[0] / dt, n=int(Tobs)), t0=t0, dt=dt)
gb_gpu_td = xr.Dataset({'A': A_td, 'E': E_td})

plt.figure()
plt.plot(X_td.t, X_td)
plt.plot(X_td[index_start:index_end].t, X_td[index_start:index_end])
plt.xlabel('time')
plt.ylabel('TDI X')
plt.show()

# X_fd = X_td[index_start:index_end].ts.fft(win=window)

tdi_fs = xr.Dataset(dict([(k, tdi_ts[k][index_start:index_end].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))

Af = (Zs - Xs)/np.sqrt(2.0)
Ef = (Zs - 2.0*Ys + Xs)/np.sqrt(6.0)   
plt.figure()
plt.semilogx(Xs.f*1000,np.real(Xs.values), label='fast GB')
plt.plot(np.array(gb_gpu.freqs.get()[0])*1000,np.real(gb_gpu.X.get()[0]), '.', label='GBGPU')
plt.plot(np.array(X_fd.f)*1000,np.real(X_fd), label='fast GB part')
plt.semilogx(tdi_fs_with_glitches['X'].f.values*1000,np.real(tdi_fs_with_glitches['X'].values),label='data')
plt.semilogx(tdi_fs['X'].f.values*1000,np.real(tdi_fs['X'].values), label='original data')
# plt.semilogx(tdi_fs_clean['X'].f.values*1000,np.real(tdi_fs_clean['X'].values))
# plt.semilogx(tdi_fs_noisefree['X'].f.values*1000,np.real(tdi_fs_noisefree['X'].values))
plt.semilogx(tdi_fs_sky['X'].f.values*1000,np.real(tdi_fs_sky['X'].values),'.',label='sky')
plt.xlim(lower_frequency*1000,upper_frequency*1000)
# plt.semilogx(Ef.f*1000,Ef.values)
plt.legend()
plt.show()

print('end')