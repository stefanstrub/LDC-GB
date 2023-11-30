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

gpu_available = False
import numpy as xp

from gbgpu.gbgpu import GBGPU

from gbgpu.utils.constants import *

from sources import *

gpu_available = False
import numpy as xp
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
taper_gaps = False
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
    SAVEPATH = grandparent+"/LDC/Spritz/"

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
    noise_model =  "spritz"
    names = fid["sky/cat"].dtype.names
    params = [np.array(fid["sky/cat"][k]).squeeze() for k in names]
    cat_vgb = np.rec.fromarrays(params, names=list(names))
    indexes = np.argsort(cat_vgb['Frequency'])
    cat_vgb = cat_vgb[indexes]
    cat = []
    for i in range(len(cat_vgb['Frequency'])):
        cat.append({})
        for parameter in parameters:
            cat[i][parameter] = cat_vgb[parameter][i]

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
# td_noisefree = fid_spritz["noisefree/tdi"][()]
# td_noisefree = np.rec.fromarrays(list(td_noisefree.T), names=["t", "X", "Y", "Z"])
# td_noisefree = td_noisefree['t']
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


# add gaps to tdi
tdi_ts_with_glitches = deepcopy(tdi_ts)
if dataset != 'Spritz':
    for k in ["X", "Y", "Z"]:
        tdi_ts_with_glitches[k].values = tdi_ts[k].values + tdi_ts_glitches[k].values- tdi_ts_glitches[k].values
if dataset == 'Spritz':
    for k in ["X", "Y", "Z"]:
        # tdi_ts_with_glitches[k].values = tdi_ts_clean[k].values + tdi_ts_glitches[k].values - tdi_ts_glitches[k].values + tdi_ts_galaxy[k].values
        tdi_ts_with_glitches[k].values = tdi_ts[k].values - tdi_ts_glitches[k].values
tdi_ts_with_glitches['A'] = TimeSeries((tdi_ts_with_glitches['Z'] - tdi_ts_with_glitches['X'])/np.sqrt(2.0), dt=dt, t0=td.t[0])
tdi_ts_with_glitches['E'] = TimeSeries((tdi_ts_with_glitches['Z']- 2.0*tdi_ts_with_glitches['Y'] + tdi_ts_with_glitches['X'])/np.sqrt(6.0), dt=dt, t0=td.t[0])

tdi_ts_clean_galaxy = deepcopy(tdi_ts)
for k in ["X", "Y", "Z"]:
    tdi_ts_clean_galaxy[k].values = tdi_ts_clean[k].values + tdi_ts_galaxy[k].values
tdi_ts_clean_galaxy['A'] = TimeSeries((tdi_ts_clean_galaxy['Z'] - tdi_ts_clean_galaxy['X'])/np.sqrt(2.0), dt=dt, t0=td.t[0])
tdi_ts_clean_galaxy['E'] = TimeSeries((tdi_ts_clean_galaxy['Z']- 2.0*tdi_ts_clean_galaxy['Y'] + tdi_ts_clean_galaxy['X'])/np.sqrt(6.0), dt=dt, t0=td.t[0])


glitchesX = np.abs(tdi_ts_glitches["X"]) > 10**-22
glitchesY = np.abs(tdi_ts_glitches["Y"]) > 10**-22
glitchesZ = np.abs(tdi_ts_glitches["Z"]) > 10**-22
plt.figure()
# plt.plot(tdi_ts_with_glitches['X'].t, tdi_ts_with_glitches["X"])
# plt.plot(tdi_ts['X'].t, tdi_ts['X'])
# plt.plot(tdi_ts_obs['X'].t, tdi_ts_obs['X'])
# plt.plot(tdi_ts_clean['X'].t, tdi_ts_clean["X"])
# plt.plot(tdi_ts_noisefree['X'].t, tdi_ts_noisefree["X"])
# plt.plot(tdi_ts_sky['X'].t, tdi_ts_sky["X"])
# plt.plot(tdi_ts_galaxy['X'].t, tdi_ts_galaxy["X"])
plt.plot(tdi_ts_glitches['X'].t, tdi_ts_glitches["X"])
plt.plot(tdi_ts_glitches['X'][glitchesX].t, tdi_ts_glitches["X"][glitchesX],'.')
plt.xlabel('time')
plt.ylabel('TDI X')
plt.show()

GB = fastGB.FastGB(delta_t=dt, T=Tobs)  # in seconds

# plt.figure()
# plt.plot(tdi_ts_with_glitches['X'].t, tdi_ts_with_glitches['X'])
# plt.plot(tdi_ts_with_glitches['Y'].t, tdi_ts_with_glitches['Y'])
# plt.plot(tdi_ts_with_glitches['Z'].t, tdi_ts_with_glitches['Z'])
# # plt.plot(tdi_ts['X'].t, tdi_ts['X'])
# plt.xlabel('time')
# plt.ylabel('TDI')
# plt.show()

gaps = {}
for k in ["X", "Y", "Z", 'A', 'E']:
    gap = np.isnan(tdi_ts_with_glitches[k])
    tdi_ts_with_glitches[k][gap] = 0
    gaps[k] = tdi_ts_with_glitches[k] == 0
    # gaps = np.isnan(tdi_ts_with_glitches[k])
    tdi_ts_with_glitches[k][gaps[k]] = 0


tdi_ts_with_glitches = deepcopy(tdi_ts_clean_galaxy)
# add gap
gap_duration = 48*3600
gap_location = 60*24*3600
gap_locations = [70*24*3600,140*24*3600,220*24*3600,300*24*3600]
for gap_location in gap_locations:
    gap_start_index = np.searchsorted(tdi_ts_with_glitches[k].t, gap_location)
    gap_end_index = np.searchsorted(tdi_ts_with_glitches[k].t, gap_location+gap_duration)
    for k in ["X", "Y", "Z", 'A', 'E']:
        tdi_ts_with_glitches[k][gap_start_index:gap_end_index] = 0


# for k in ["X", "Y", "Z"]:
#     for i in range(1):
#         mad = scipy.stats.median_abs_deviation(tdi_ts_with_glitches[k])
#         peaks, properties = scipy.signal.find_peaks(np.abs(tdi_ts_with_glitches[k]), height=10*mad, threshold=None, distance=1)
#         # Turning glitches into gaps
#         for pk in peaks:
#             tdi_ts_with_glitches['X'][pk-20:pk+20] = 0
#             tdi_ts_with_glitches['Y'][pk-20:pk+20] = 0
#             tdi_ts_with_glitches['Z'][pk-20:pk+20] = 0

plt.figure()
plt.plot(tdi_ts_with_glitches['X'].t, tdi_ts_with_glitches['X'])
# plt.plot(tdi_ts_with_glitches['Y'].t, tdi_ts_with_glitches['Y'])
# plt.plot(tdi_ts_with_glitches['Z'].t, tdi_ts_with_glitches['Z'])
# plt.plot(tdi_ts_glitches['X'].t, tdi_ts_glitches["X"])
# plt.plot(tdi_ts_glitches['X'][glitchesX].t, tdi_ts_glitches["X"][glitchesX],'.')
# plt.plot(tdi_ts_glitches['Y'][glitchesY].t, tdi_ts_glitches["Y"][glitchesY],'.')
# plt.plot(tdi_ts_glitches['Z'][glitchesZ].t, tdi_ts_glitches["Z"][glitchesZ],'.')
# plt.plot(tdi_ts['X'].t, tdi_ts['X'])
plt.xlabel('time')
plt.ylabel('TDI')
plt.show()

tdi_fs_clean_galaxy = xr.Dataset(dict([(k, tdi_ts_clean_galaxy[k].ts.fft()) for k in ["X", "Y", "Z", 'A', 'E']]))
# tdi_fs_clean_galaxy = tdi_fs_clean_galaxy.assign(A=(tdi_fs_clean_galaxy.Z - tdi_fs_clean_galaxy.X)/np.sqrt(2.0))
# tdi_fs_clean_galaxy = tdi_fs_clean_galaxy.assign(E=(tdi_fs_clean_galaxy.Z - 2.0*tdi_fs_clean_galaxy.Y + tdi_fs_clean_galaxy.X)/np.sqrt(6.0))


N = 256
data_gpu = []
PSD_GPU = []
data_gpu.append([xp.array(tdi_fs_clean_galaxy.X.values), xp.array(tdi_fs_clean_galaxy.X.values)])
Nmodel = get_noise_model(noise_model, tdi_fs_clean_galaxy.X.f)
SA_f = Nmodel.psd(freq=tdi_fs_clean_galaxy.X.f, option="A", tdi2=True)

df = tdi_fs_clean_galaxy.f[1] - tdi_fs_clean_galaxy.f[0]


f_Nyquist = 1/dt/2
search_range = [0.0003, f_Nyquist]
frequencies = create_frequency_windows(search_range, Tobs)
frequencies_even = frequencies[::2]
frequencies_odd = frequencies[1::2]

frequencies_search = []
for i in range(len(cat)):
    start_index = np.searchsorted(np.asarray(frequencies)[:,0], cat[i]['Frequency'])-1
    frequencies_search.append(frequencies[start_index])

for cat_index in range(len(cat)):
    # cat[cat_index]['FrequencyDerivative'] *= 10
    # lower_frequency = cat[cat_index]['Frequency']-0.0000005
    # upper_frequency = cat[cat_index]['Frequency']+0.0000005
    lower_frequency = frequencies_search[cat_index][0]
    upper_frequency = frequencies_search[cat_index][1]
    search = Search(tdi_fs_clean_galaxy,Tobs, lower_frequency, upper_frequency, tdi2=True)
    # search.plot(found_sources_in= [cat[cat_index]])
    print(cat[cat_index]['Frequency'],search.SNR([cat[cat_index]]))

cat_index =35
lower_frequency = frequencies_search[cat_index][0]
upper_frequency = frequencies_search[cat_index][1]
search = Search(tdi_fs_clean_galaxy,Tobs, lower_frequency, upper_frequency, tdi2=True)
lower_index = np.searchsorted(tdi_fs_clean_galaxy['A'].f, search.frequencyrange[0])
upper_index = np.searchsorted(tdi_fs_clean_galaxy['A'].f, search.frequencyrange[1])
lower_index = np.searchsorted(tdi_fs_clean_galaxy['A'].f, 0.0003)
upper_index = np.searchsorted(tdi_fs_clean_galaxy['A'].f, 0.02)

def is_odd(num):
    return num & 0x1


if fill_gaps:
    tdi_ts_filled = deepcopy(tdi_ts_with_glitches)
    for k in ["X", "Y", "Z"]:
        gaps = tdi_ts_filled[k] == 0
        gaps = tdi_ts_filled[k][gaps]
        start_points = []
        if gaps.t[0].values != tdi_ts_filled['X'].t[0].values:
            start_points = [0]
            end_points = np.array([gaps.t[0].values]-dt)
        differences = gaps.t[1:].values - gaps.t[:-1].values
        jumps = differences > 2*dt

        start_points_data = list(start_points) + list(np.array(gaps[:-1][jumps].t.values)+dt) + list([gaps.t[-1].values]+dt)
        end_points_data =  list(end_points) + list(np.array(gaps[1:][jumps].t.values)-dt) + list([float(tdi_ts_filled['X'].t[-1].values)])

        start_points_gaps = end_points_data[:-1]
        end_points_gaps = start_points_data[1:]
            
        for i in range(len(start_points)):
            index_start = int((start_points_gaps[i]-tdi_ts_with_glitches[k].t[0].values)/dt)
            index_end = int((end_points_gaps[i]-tdi_ts_with_glitches[k].t[0])/dt)
            length_gap = len(tdi_ts_with_glitches[k][index_start:index_end])+1
            tdi_ts_with_glitches[k][index_start:index_start+length_gap].values = tdi_ts_with_glitches[k][index_end+1:index_end+1+length_gap].values
            # else:
            #     tdi_ts_with_glitches[k][index_start:index_start+int(length_gap/2)] = tdi_ts_with_glitches[k][index_start-int(length_gap/2):index_start].values
            #     tdi_ts_with_glitches[k][index_end-int(length_gap/2):index_end] = tdi_ts_with_glitches[k][index_end:index_end+int(length_gap/2)].values
    tdi_fs_filled = xr.Dataset(dict([(k, tdi_ts_filled[k].ts.fft()) for k in ['A', 'E']]))

    diff_filled_clean = tdi_fs_clean_galaxy - tdi_fs_filled
    fd_error = -1/2*(4.0*df*np.sum(((diff_filled_clean.A.data*np.conjugate(diff_filled_clean.A.data) + diff_filled_clean.E.data*np.conjugate(diff_filled_clean.E.data)) /SA_f)[lower_index:upper_index]).data)


    plt.figure()
    plt.loglog(tdi_fs_filled.f,np.abs(tdi_fs_filled.A))
    plt.loglog(tdi_fs_clean_galaxy.f,np.abs(tdi_fs_clean_galaxy.A))
    plt.loglog(diff_filled_clean.A.f[lower_index:upper_index],np.abs(diff_filled_clean.A)[lower_index:upper_index])
    plt.loglog(SA_f.f,np.sqrt(SA_f))
    plt.show()

    tdi_ts_with_glitches = deepcopy(tdi_ts_filled)


# taper_lengths = [0,10,20,30,40,50,60,70,80,90,100,110,120,200,300,400,500,600]
taper_lengths = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,30,40,50,60,70,80,90,100,120,130,140,150,200]
# taper_lengths = [70]
if taper_gaps:
    fd_error = []
    plt.figure()
    for length_taper_index in taper_lengths:
        # length_taper = 90
        # length_taper_index = 1+2*int(length_taper/2/dt)
        middle_window = int((length_taper_index)/2)
        window_hann = scipy.signal.windows.barthann(length_taper_index)
        window_hann_right = window_hann[middle_window:]
        window_hann_left = window_hann[:middle_window+is_odd(length_taper_index)]

        # plt.figure()
        # plt.plot(window_hann)
        # plt.plot(window_hann_right)
        # plt.plot(window_hann_left)
        # plt.title("Hann window")
        # plt.ylabel("Amplitude")
        # plt.xlabel("Sample")
        # plt.show()

        # A = np.fft.fft(window_hann, 2048) / (len(window_hann)/2.0)
        # freq = np.linspace(-0.5, 0.5, len(A))
        # response = np.abs(np.fft.fftshift(A / np.max(np.abs(A))))
        # response = 20 * np.log10(np.maximum(response, 1e-10))
        # plt.plot(freq, response, label=length_taper)

        tdi_ts_tapered = deepcopy(tdi_ts_with_glitches)
        for k in ['A', 'E']:
        # for k in ["X", "Y", "Z"]:
            gaps = tdi_ts_tapered[k] == 0
            gaps = tdi_ts_tapered[k][gaps]
            start_points = []
            if gaps.t[0].values != tdi_ts_tapered['X'].t[0].values:
                start_points = [0]
                end_points = np.array([gaps.t[0].values]-dt)
            differences = gaps.t[1:].values - gaps.t[:-1].values
            jumps = differences > 2*dt

            start_points_data = list(start_points) + list(np.array(gaps[:-1][jumps].t.values)+dt) + list([gaps.t[-1].values]+dt)
            end_points_data =  list(end_points) + list(np.array(gaps[1:][jumps].t.values)-dt) + list([float(tdi_ts_tapered['X'].t[-1].values)])

            start_points_gaps = end_points_data[:-1]
            end_points_gaps = start_points_data[1:]

            for i in range(len(start_points_gaps)):
                index_start = int((start_points_gaps[i]-tdi_ts_tapered[k].t[0].values)/dt)+1
                index_end = int((end_points_gaps[i]-tdi_ts_tapered[k].t[0])/dt)
                tdi_ts_tapered[k][index_start-len(window_hann_right):index_start] *= window_hann_right
                tdi_ts_tapered[k][index_end:index_end+len(window_hann_left)] *= window_hann_left

        tdi_fs_tapered = xr.Dataset(dict([(k, tdi_ts_tapered[k].ts.fft()) for k in ['A', 'E']]))

        diff_tapered_clean = tdi_fs_clean_galaxy - tdi_fs_tapered
        # diff_tapered_clean = diff_tapered_clean.assign(A=diff_tapered_clean.Z - diff_tapered_clean.X)/np.sqrt(2.0)
        # diff_tapered_clean = diff_tapered_clean.assign(E=diff_tapered_clean.Z - 2.0*diff_tapered_clean.Y + diff_tapered_clean.X)/np.sqrt(6.0)
        # dh = np.sum( np.real(self.DAf * np.conjugate(tdi_fs_tapered.A.data) + self.DEf * np.conjugate(Ef.data))/self.SA )
        # hh = np.sum((np.absolute(Af.data)**2 + np.absolute(Ef.data)**2) /self.SA)
        # dh = 4.0*df*np.sum((np.real(tdi_fs_clean_galaxy.A.data * np.conjugate(tdi_fs_tapered.A.data) + tdi_fs_clean_galaxy.E.data * np.conjugate(tdi_fs_tapered.E.data)) /SA_f)[lower_index:upper_index]).data
        # hh = 4.0*df*np.sum(((np.absolute(tdi_fs_tapered.A.data)**2 + np.absolute(tdi_fs_tapered.E.data)**2) /SA_f)[lower_index:upper_index]).data
        # dd = 4.0*df*np.sum(((np.absolute(tdi_fs_clean_galaxy.A.data)**2 + np.absolute(tdi_fs_clean_galaxy.E.data)**2) /SA_f)[lower_index:upper_index]).data
        # fd_error.append(-1/2*(4.0*df*np.sum(((np.absolute(diff_tapered_clean.A.data)**2 + np.absolute(diff_tapered_clean.E.data)**2) /SA_f)[lower_index:upper_index])).data)
        fd_error.append(-1/2*(4.0*df*np.sum(((diff_tapered_clean.A.data*np.conjugate(diff_tapered_clean.A.data) + diff_tapered_clean.E.data*np.conjugate(diff_tapered_clean.E.data)) /SA_f)[lower_index:upper_index])).data)
        # fd_error.append(-1/2*(4.0*df*np.sum(((np.absolute(tdi_fs_clean_galaxy.A.data-tdi_fs_tapered.A.data)**2 + np.absolute(tdi_fs_clean_galaxy.E.data-tdi_fs_tapered.E.data)**2) /SA_f)[lower_index:upper_index])).data)
        # fd_error.append((dh/np.sqrt(hh)).data)
        # fd_error.append((dh-hh/2).data)
        # fd_error.append((dh-hh/2-dd/2).data)

    # plt.title("Frequency response of the Hann window")
    # plt.ylabel("Normalized magnitude [dB]")
    # plt.xlabel("Normalized frequency [cycles per sample]")
    # plt.legend()
    # plt.show()

    max_index = np.argmax(fd_error)
    plt.figure()
    plt.plot(taper_lengths, fd_error)
    plt.plot(taper_lengths[max_index],fd_error[max_index],'.')
    plt.xlabel('taper length (index)')
    plt.ylabel('loglikelihood of tapered signal')
    plt.show()

    plt.figure()
    plt.loglog(tdi_fs_tapered.f,np.abs(tdi_fs_tapered.A))
    plt.loglog(tdi_fs_clean_galaxy.f,np.abs(tdi_fs_clean_galaxy.A))
    plt.loglog(diff_tapered_clean.A.f[lower_index:upper_index],np.abs(diff_tapered_clean.A)[lower_index:upper_index])
    plt.loglog(SA_f.f,np.sqrt(SA_f))
    plt.show()

    tdi_ts_with_glitches = deepcopy(tdi_ts_tapered)
    # plt.figure()
    # # plt.plot(gaps.t[:-1],differences)
    # plt.plot(tdi_ts_with_glitches['X'].t, tdi_ts_with_glitches['X'])
    # # plt.plot(gaps[:-1][jumps].t.values, float(np.max(tdi_ts_with_glitches['X']))*np.ones_like(gaps[:-1][jumps].t.values), '.')
    # # plt.plot(gaps[1:][jumps].t.values, float(np.max(tdi_ts_with_glitches['X']))*np.ones_like(gaps[1:][jumps].t.values), '.')
    # plt.plot(start_points_data, float(np.max(tdi_ts_with_glitches['X']))*np.ones_like(start_points_data), '.')
    # plt.plot(end_points_data,float(np.max(tdi_ts_with_glitches['X']))*np.ones_like(end_points_data), '.')
    # plt.show()

    plt.figure()
    # plt.plot(tdi_ts_clean_galaxy['A'].t, tdi_ts_clean_galaxy['A'])
    plt.plot(tdi_ts_with_glitches['A'].t, tdi_ts_with_glitches['A'], label='Data')
    plt.plot(tdi_ts_tapered['A'].t, tdi_ts_tapered['A'], label='Tapered data')
    # plt.plot(tdi_ts_with_glitches['Y'].t, tdi_ts_with_glitches['Y'])
    # plt.plot(tdi_ts_with_glitches['Z'].t, tdi_ts_with_glitches['Z'])
    # plt.plot(tdi_ts['X'].t, tdi_ts['X'])
    plt.xlabel('time (s)')
    plt.ylabel('TDI A')
    plt.legend(loc='upper right')
    plt.show()

start = 200
tdi_ts_signals = deepcopy(tdi_ts)
plt.figure(figsize=(16, 6))
plt.subplot(211)
# plt.plot(tdi_ts['X'].t[start:], tdi_ts['X'][start:])
# plt.plot(tdi_ts_clean['X'].t[start:], tdi_ts_clean["X"][start:])
# plt.plot(tdi_ts_clean_gaps['X'].t[start:], tdi_ts_clean_gaps["X"][start:])
# plt.plot(tdi_ts_noisefree['X'].t[start:], tdi_ts_noisefree["X"][start:])
plt.plot(tdi_ts_sky['X'].t[start:], tdi_ts_sky["X"][start:])
# plt.plot(tdi_ts_glitches['X'].t[start:], tdi_ts_glitches["X"][start:])
# plt.plot(tdi_ts_signals['X'].t[start:], tdi_ts_signals["X"][start:])
plt.ylabel("TDI X")
plt.xlabel("Time [s]")
plt.xlim([1.8e7, 2.2e7])
# plt.ylim([-1e-18, 1e-18])
plt.show()


# data_fd= {}
# gap_less_duration = {}
# for k in ["X", "Y", "Z"]:
#     data_fd[k] = []
#     gap_less_duration[k] = []
#     groups = []
#     gaps = tdi_ts_with_glitches[k] == 0
#     gaps = tdi_ts_with_glitches[k][gaps]
#     start_points = []
#     if gaps.t[0].values != tdi_ts_with_glitches[k].t[0].values:
#         start_points = [0]
#         end_points = gaps.t[0].values
#     differences = gaps.t[1:].values - gaps.t[:-1].values
#     jumps = differences > 15
#     start_points = list(start_points) + list(gaps[:-1][jumps].t.values) + list([gaps.t[-1].values])
#     end_points =  list(gaps[1:][jumps].t.values)
#     for i in range(len(end_points)):
#         index_start = int((start_points[i]-tdi_ts_with_glitches[k].t[0].values)/dt)
#         index_end = int((end_points[i]-tdi_ts_with_glitches[k].t[0])/dt)
#         data_fd[k].append(tdi_ts_with_glitches[k][index_start:index_end].ts.fft(win=window))
#         gap_less_duration[k].append(tdi_ts_with_glitches[k][index_start:index_end].t[-1].values- tdi_ts_with_glitches[k][index_start:index_end].t[0].values)

data_fd = []
gap_less_duration = []
groups = []

gaps = tdi_ts_with_glitches[k] == 0
gaps = tdi_ts_with_glitches[k][gaps]
start_points = []
if gaps.t[0].values != tdi_ts_with_glitches['X'].t[0].values:
    start_points = [0]
    end_points = np.array([gaps.t[0].values]-dt)
differences = gaps.t[1:].values - gaps.t[:-1].values
jumps = differences > 2*dt

start_points = list(start_points) + list(np.array(gaps[:-1][jumps].t.values)+dt) + list([gaps.t[-1].values]+dt)
end_points =  list(end_points) + list(np.array(gaps[1:][jumps].t.values)-dt) + list([float(tdi_ts_with_glitches['X'].t[-1].values)])

start_points_gaps = end_points[:-1]
end_points_gaps = start_points[1:]

## no gaps
# tdi_ts_with_glitches = deepcopy(tdi_ts)
# tdi_ts_with_glitches = deepcopy(tdi_ts_clean_galaxy)
# tdi_ts_with_glitches = tdi_ts_clean
# start_points = np.array([0])
# end_points = np.array([tdi_ts_with_glitches['X'].t[-1].values])

Tobs_list = np.array(end_points) - np.array(start_points)
# remove short observation times
i = 0
while i < len(Tobs_list):
    if Tobs_list[i] < 300000:
        start_points.pop(i)
        end_points.pop(i)
        Tobs_list = np.array(end_points) - np.array(start_points)
    i += 1
for i in range(len(end_points)):
    index_start = int((start_points[i]-tdi_ts_with_glitches['X'].t[0].values)/dt)
    index_end = int((end_points[i]-tdi_ts_with_glitches['X'].t[0])/dt)

    tdi_ts_sky_part = {}
    for k in ["X", "Y", "Z"]:
        # tdi_ts_sky[k]['t'] = tdi_ts_sky[k].t[start_index:-end_index or None]
        tdi_ts_sky_part[k] = tdi_ts_with_glitches[k][index_start:index_end or None]
        # tdi_ts_sky_part[k] = tdi_ts_sky[k][index_start:index_end or None]
        # tdi_ts_sky_part[k]['t0'] = t_start
    data_fd.append(xr.Dataset(dict([(k, tdi_ts_sky_part[k].ts.fft()) for k in ["X", "Y", "Z"]])))
    # data_fd.append(xr.Dataset(dict([(k, np.fft.rfft(tdi_ts_sky[k][index_start:index_end])*dt) for k in ["X", "Y", "Z"]])))
    # data_fd[k].append(tdi_ts_with_glitches[k][index_start:index_end].ts.fft(win=window))
    gap_less_duration.append(tdi_ts_with_glitches['X'][index_start:index_end].t[-1].values- tdi_ts_with_glitches['X'][index_start:index_end].t[0].values)

for i in range(len(end_points)):
    data_fd[i]['A'] = (data_fd[i]['Z'] - data_fd[i]['X'])/np.sqrt(2.0)
    data_fd[i]['E'] = (data_fd[i]['Z'] - 2.0*data_fd[i]['Y'] + data_fd[i]['X'])/np.sqrt(6.0)

plt.figure()
plt.plot(tdi_ts_with_glitches['X'].t, tdi_ts_with_glitches['X'])
plt.show()

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

# print(search.SNR([cat[i]]))
# print(search.SNR_scaled([cat[i]]))
parameters_gpgpu = ['Amplitude','Frequency','FrequencyDerivative','Frequency2Derivative','InitialPhase','Inclination','Polarization','EclipticLongitude','EclipticLatitude']
start = time.time()
Xs, Ys, Zs = GB.get_fd_tdixyz(template=cat[cat_index], oversample=4, tdi2=True)
print('time', time.time()-start)

start = time.time()
X_td, Ys_td, Zs_td = GB.get_td_tdixyz(template=cat[cat_index], oversample=4, tdi2=True)
print('time', time.time()-start)

# simulated_td = { 'X': X_td, 'Y': Ys_td, 'Z': Zs_td}
# simulated_fd = {}
# for k in ["X", "Y", "Z"]:
#     simulated_fd[k] = []
#     for i in range(len(end_points)):
#         index_start = int((start_points[i]-tdi_ts_with_glitches[k].t[0].values)/dt)
#         index_end = int((end_points[i]-tdi_ts_with_glitches[k].t[0])/dt)
#         simulated_fd[k].append(simulated_td[k][index_start:index_end].ts.fft(win=window))



start_index = 1000000
# start_index = 0
end_index = 6000000
start_index = index_start
end_index = index_end
k = 'X'
t_start = float(tdi_ts_sky[k].t[start_index].values)
t_end = float(tdi_ts_sky[k].t[end_index].values)
tdi_ts_sky_part = {}
tdi_ts_data_part = {}
tdi_ts_clean_part = {}
tdi_ts_clean_galaxy_part = {}
for k in ["X", "Y", "Z"]:
    # tdi_ts_sky[k]['t'] = tdi_ts_sky[k].t[start_index:-end_index or None]
    tdi_ts_data_part[k] = tdi_ts_with_glitches[k][start_index:end_index or None]
    tdi_ts_clean_galaxy_part[k] = tdi_ts_clean_galaxy[k][start_index:end_index or None]
    tdi_ts_clean_part[k] = tdi_ts_clean[k][start_index:end_index or None]
    tdi_ts_sky_part[k] = tdi_ts_sky[k][start_index:end_index or None]
    # tdi_ts_sky_part[k]['t0'] = t_start
tdi_fs_data = xr.Dataset(dict([(k, tdi_ts_data_part[k].ts.fft()) for k in ["X", "Y", "Z"]]))
tdi_fs_clean = xr.Dataset(dict([(k, tdi_ts_clean_part[k].ts.fft()) for k in ["X", "Y", "Z"]]))
tdi_fs_clean_galaxy = xr.Dataset(dict([(k, tdi_ts_clean_galaxy_part[k].ts.fft()) for k in ["X", "Y", "Z"]]))
tdi_fs_sky = xr.Dataset(dict([(k, tdi_ts_sky_part[k].ts.fft()) for k in ["X", "Y", "Z"]]))

data_glitches = {}
data_glitches['A'] = (tdi_fs_data['Z'] - tdi_fs_data['X'])/np.sqrt(2.0)
data_glitches['E'] = (tdi_fs_data['Z'] - 2.0*tdi_fs_data['Y'] + tdi_fs_data['X'])/np.sqrt(6.0)

plt.figure()
plt.plot(tdi_ts_data_part[k].t,tdi_ts_data_part[k])
# plt.plot(tdi_ts_glitches[k].t,tdi_ts_glitches[k])
plt.plot(tdi_ts_data_part[k].t,tdi_ts_data_part[k]-tdi_ts_clean_part[k])
# plt.plot(tdi_ts_clean[k].t,tdi_ts_clean[k])
plt.show()

i =0
plt.figure()
plt.plot(data_fd[i].f.values,data_fd[i]['X'].values)
plt.plot(tdi_fs_clean.f,tdi_fs_clean['X'])
plt.plot(tdi_fs_clean_galaxy.f,tdi_fs_clean_galaxy['X'])
plt.plot(Xs.f,Xs)
plt.xlim( frequencies_search[cat_index])
plt.show()

pGB = {}
for parameter in parameters:
    pGB[parameter] = cat[cat_index][parameter]
pGB['Frequency2Derivative'] = 0
if pGB['EclipticLongitude'] > np.pi:
    pGB['EclipticLongitude'] -= 2*np.pi
# pGB['FrequencyDerivative'] = 10**-16
num_bin = 1
pGBs = deepcopy(pGB)
# pGBs['InitialPhase'] = cat[cat_index]['InitialPhase'] - pGB['Frequency']*t_start * np.pi*2 - pGB['FrequencyDerivative']*t_start**2 * np.pi - pGB['Frequency2Derivative']*t_start**3 * np.pi/3
params = np.array([np.full(num_bin, pGBs[parameter]) for parameter in parameters_gpgpu])

N = 256
noise_model =  "SciRDv1"
data_gpu = []
PSD_GPU = []
for i in range(len(data_fd)):
    data_gpu.append([xp.array(data_fd[i]['A'].values), xp.array(data_fd[i]['E'].values)])
    Nmodel = get_noise_model(noise_model, data_fd[i]['A'].f)
    SA_full_f = Nmodel.psd(freq=data_fd[i]['A'].f, option="A")
    SE_full_f = Nmodel.psd(freq=data_fd[i]['E'].f, option="E")
    PSD_GPU.append(xp.array([SA_full_f, SE_full_f]))
gb_gpu.d_d = 0.0
SNR_list = []
pGBs = deepcopy(pGB)

params = np.array([np.full(num_bin, pGBs[parameter]) for parameter in parameters_gpgpu])
for i in range(len(data_fd)):
    # pGBs['InitialPhase'] = cat[cat_index]['InitialPhase'] - pGB['Frequency']*t_start * np.pi*2 - pGB['FrequencyDerivative']*t_start**2 * np.pi - pGB['Frequency2Derivative']*t_start**3 * np.pi/3
    SNR_list.append(gb_gpu.get_ll(params, data_gpu[i], PSD_GPU[i], N=N, oversample=4, dt=dt, T=end_points[i]-start_points[i], start_freq_ind=0, tdi2=True, t_start=start_points[i], get_SNR=True))

total_gapless_time = np.sum(gap_less_duration)
SNR_average = 0
for i in range(len(SNR_list)):
    SNR_average += SNR_list[i] * gap_less_duration[i]
SNR_average /= total_gapless_time

class GBSearch_with_gaps():
    def __init__(self, data_gpu, PSD_GPU, dt, Tobs, lower_frequency, upper_frequency, start_points, end_points, gap_less_duration, parameters, recombination=0.75, get_SNR=True):
        self.data_gpu = data_gpu
        self.PSD_GPU = PSD_GPU
        self.dt = dt
        self.Tobs = Tobs
        self.lower_frequency = lower_frequency
        self.upper_frequency = upper_frequency
        self.start_points = start_points
        self.end_points = end_points
        self.gap_less_duration = gap_less_duration
        self.parameters = parameters
        self.search = Search(tdi_fs_with_glitches,Tobs, lower_frequency, upper_frequency, tdi2=True, gb_gpu=gb_gpu)
        self.search.initialize_for_gaps(end_points, start_points, gap_less_duration, data_gaps_gpu=data_gpu, PSD_gaps_gpu=PSD_GPU, get_SNR=get_SNR)
        self.parameters_log_uniform = self.search.parameters_log_uniform
        self.boundaries = self.search.boundaries
        self.recombination = recombination
        self.get_SNR = get_SNR
        self.total_gapless_time = np.sum(self.gap_less_duration)

    # def loglikelihood_gaps_gpu(self, params, get_SNR=None):
    #     if get_SNR == None:
    #         get_SNR = self.get_SNR
    #     self.loglikelihood_list = []
    #     if len(np.shape(params)) == 2:
    #         pGB_gpu = np.copy([params])
    #     else:
    #         pGB_gpu = np.copy(params)
    #     # initial_phase = np.copy(pGB_gpu[:,4])
    #     for i in range(len(self.data_gpu)):
    #         # pGB_gpu[:,4] = initial_phase - pGB_gpu[:,1]*start_points[i] * np.pi*2 - pGB_gpu[:,2]*t_start**2 * np.pi - pGB_gpu[:,3]*t_start**3 * np.pi/3
    #         # pGB_gpu = np.array([np.full(num_bin, pGB_gpu[parameter]) for parameter in pGB_gpu_gpgpu])
    #         self.loglikelihood_list.append(gb_gpu.get_ll(pGB_gpu, self.data_gpu[i], self.PSD_GPU[i], N=N, oversample=4, dt=self.dt, T=self.end_points[i]-self.start_points[i], start_freq_ind=0, tdi2=True, t_start=self.start_points[i], get_SNR=get_SNR))

    #     loglikelihood_average = 0
    #     for i in range(len(self.loglikelihood_list)):
    #         loglikelihood_average += self.loglikelihood_list[i] * self.gap_less_duration[i]
    #     loglikelihood_average /= self.total_gapless_time
    #     return loglikelihood_average

    def differential_evolution_search(self, initial_guess = None, number_of_signals = 1, get_SNR=True):
        bounds = []
        for signal in range(number_of_signals):
            for i in range(7):
                bounds.append((0,1))

        maxpGB = []
        if initial_guess != None:
            initial_guess01 = np.zeros((len(self.parameters)-1)*number_of_signals)
            for signal in range(number_of_signals):
                pGBstart01 = scaleto01(initial_guess[signal], self.boundaries, self.parameters, self.parameters_log_uniform)

                for count, parameter in enumerate(self.parameters_no_amplitude):
                    if pGBstart01[parameter] < 0:
                        pGBstart01[parameter] = 0
                    if pGBstart01[parameter] > 1:
                        pGBstart01[parameter] = 1
                    initial_guess01[count+(len(self.parameters_no_amplitude))*signal] = pGBstart01[parameter]
            start = time.time()
            res = differential_evolution(self.function_evolution_gb_gpu, bounds=bounds, disp=True, strategy='best1exp', popsize=8,tol= 1e-8 , maxiter=1000, recombination= self.recombination, mutation=(0.5,1), x0=initial_guess01)
            print('time',time.time()-start)
        else:
            start = time.time()
            res = differential_evolution(self.function_evolution_gb_gpu, bounds=bounds, disp=True, strategy='best1exp', popsize=8, tol= 1e-8 , maxiter=1000, recombination= self.recombination, mutation=(0.5,1))
            print('time differential evolution',time.time()-start)
        for signal in range(number_of_signals):
            pGB01 = [0.5] + res.x[signal*7:signal*7+7].tolist()
            maxpGB.append(scaletooriginal(pGB01,self.boundaries, self.parameters, self.parameters_log_uniform))
            maxpGB[-1]['Frequency2Derivative'] = 0
        print(res)
        print(maxpGB)
        # print('log-likelihood',self.loglikelihood(maxpGB))
        # print(pGB)
        return maxpGB, res.nfev

    def function_evolution_gb_gpu(self, pGBs01, number_of_signals = 1):
        pGBs_gpu = []
        for signal in range(number_of_signals):
            pGBs = {}
            pGBs['Frequency2Derivative'] = 0
            i = 0
            for parameter in self.parameters:
                if parameter in ["EclipticLatitude"]:
                    pGBs[parameter] = np.arcsin((pGBs01[signal*7:signal*7+7][i] * (self.boundaries[parameter][1] - self.boundaries[parameter][0])) + self.boundaries[parameter][0])
                elif parameter in ["Inclination"]:
                    pGBs[parameter] = np.arccos((pGBs01[signal*7:signal*7+7][i] * (self.boundaries[parameter][1] - self.boundaries[parameter][0])) + self.boundaries[parameter][0])
                elif parameter in ['Amplitude']:
                    i -= 1
                    pGBs[parameter] = 10**((0.1 * (self.boundaries[parameter][1] - self.boundaries[parameter][0])) + self.boundaries[parameter][0])
                # elif parameter in ["FrequencyDerivative"]:
                #     pGBs[parameter] = 10**((pGBs01[signal*7:signal*7+7][i] * (self.boundaries[parameter][1] - self.boundaries[parameter][0])) + self.boundaries[parameter][0])
                else:
                    pGBs[parameter] = (pGBs01[signal*7:signal*7+7][i] * (self.boundaries[parameter][1] - self.boundaries[parameter][0])) + self.boundaries[parameter][0]
                i += 1
            pGBs_gpu.append([np.array([pGBs[parameter]]) for parameter in parameters_gpgpu])
        p = -self.search.loglikelihood_gaps_gpu(pGBs_gpu, get_SNR=self.get_SNR)
        return p
    
params = np.array([np.full(num_bin, pGB[parameter]) for parameter in parameters_gpgpu])

save_name = '_data_gaps_ignored'
for cat_index in range(22,len(frequencies_search)):
    pGB = {}
    for parameter in parameters:
        pGB[parameter] = cat[cat_index][parameter]
    pGB['Frequency2Derivative'] = 0
    if pGB['EclipticLongitude'] > np.pi:
        pGB['EclipticLongitude'] -= 2*np.pi
    params = np.array([np.full(num_bin, pGB[parameter]) for parameter in parameters_gpgpu])
    gap_search = GBSearch_with_gaps(data_gpu, PSD_GPU, dt, Tobs, frequencies_search[cat_index][0], frequencies_search[cat_index][1], start_points=start_points, end_points=end_points, gap_less_duration=gap_less_duration, parameters=parameters, get_SNR=True)
    print(gap_search.search.loglikelihood_gaps_gpu(params, get_SNR=True))
    amplitude = gap_search.search.calculate_Amplitude_gaps(params)
    params01 = np.ones((7))/2
    loglikelihood_function = gap_search.function_evolution_gb_gpu(params01, number_of_signals=1)
    
    amplitude_list = []
    for i in range(len(data_fd)):
        amplitude_list.append(gb_gpu.get_ll(params, data_gpu[i], PSD_GPU[i], N=N, oversample=4, dt=dt, T=end_points[i]-start_points[i], start_freq_ind=0, tdi2=True, t_start=start_points[i], get_SNR=False, get_dh_hh_ratio=True))

    total_gapless_time = np.sum(gap_less_duration)
    amplitude_average = 0
    for i in range(len(amplitude_list)):
        amplitude_average += amplitude_list[i] * gap_less_duration[i]
    amplitude_average /= total_gapless_time
    amplitude = params[0]*amplitude_average

    amplitude_factor = gap_search.search.calculate_Amplitude([pGB])

    maxpGB, nfev = gap_search.differential_evolution_search(initial_guess = None, number_of_signals = 1, get_SNR=True)
    params = np.array([np.full(num_bin, maxpGB[0][parameter]) for parameter in parameters_gpgpu])
    amplitude_factor = gap_search.search.calculate_Amplitude_gaps(params)
    maxpGB[0]['Amplitude'] *= amplitude_factor[0]
    fn = SAVEPATH+'/found_signals/found_sources'+ str(int(np.round(frequencies_search[cat_index][0]*10**9)))+'nHz_to'+ str(int(np.round(frequencies_search[cat_index][1]*10**9)))+'nHz_' +save_name+'.pkl'
    pickle.dump(maxpGB, open(fn, "wb"))

found_clean = list(pickle.load(open(SAVEPATH+'/found_sources_Spritz_clean_galaxy_flat.pkl', "rb"))[21:])
found_gaps_ignored = list(pickle.load(open(SAVEPATH+'/found_sources_Spritz_4gaps_ignored_flat.pkl', "rb")))
found_tapered = list(pickle.load(open(SAVEPATH+'/found_sources_Spritz_4gaps_tapered_hann70_flat.pkl', "rb")))
pGB_injected = cat[22:]

found_clean.pop(8)
found_gaps_ignored.pop(8)
found_tapered.pop(8)
pGB_injected.pop(8)

signal_names = ['clean', 'gaps_ignored', 'tapered_hann70']
errors = {}
frequencies_found = {}
for j, signals in enumerate([found_clean,found_gaps_ignored, found_tapered]):
    errors[signal_names[j]] = {}
    frequencies_found[signal_names[j]] = {}
    for parameter in parameters:
        errors[signal_names[j]][parameter] = []
        frequencies_found[signal_names[j]] = []
        for i in range(len(signals)):
            if parameter in ['EclipticLatitude', 'Inclination', 'Polarization', 'InintialPhase', 'EclipticLongitude']:
                errors[signal_names[j]][parameter].append(np.abs(np.arcsin(np.sin(signals[i][parameter]-pGB_injected[i][parameter]))))
            else:
                errors[signal_names[j]][parameter].append(np.abs(signals[i][parameter]-pGB_injected[i][parameter]))
            frequencies_found[signal_names[j]].append(signals[i]['Frequency'])
        frequencies_found[signal_names[j]] = np.array(frequencies_found[signal_names[j]])

plt.figure()
# parameter = 'Frequency'
parameter = 'EclipticLongitude'
# parameter = 'Amplitude'
symbol_list = ['o', '.','+']
for j, signals in enumerate([found_clean,found_gaps_ignored, found_tapered]):
    plt.plot(errors[signal_names[j]][parameter], symbol_list[j], label=signal_names[j], markersize=12)
# plt.hlines([0.1], 0,35)
plt.ylabel(parameter+ ' Error')
plt.xlabel('Signal with frequency (mHz)')
plt.xticks(ticks=np.arange(len(frequencies_found[signal_names[j]])),labels=np.round(frequencies_found[signal_names[j]]*1000,2), rotation='vertical' )
plt.legend()
plt.tight_layout()
plt.savefig(SAVEPATH+parameter+'_Error_4gaps.png')
plt.show()


# found_clean[19]['Frequency2Derivative'] = 0
# params_clean = np.array([np.full(num_bin, found_clean[19][parameter]) for parameter in parameters_gpgpu])
# loglikelihood_gpu_clean = gap_search.loglikelihood_gaps_gpu(params_clean, get_SNR=True)
# [{'Amplitude': 2.9338406741898537e-22, 'EclipticLatitude': 0.08595294406518945, 'EclipticLongitude': -2.1680570099157785, 'Frequency': 0.0018137298821796902, 'FrequencyDerivative': -6.590573749734849e-18, 'Inclination': 0.9896873937025192, 'InitialPhase': 3.843957120120863, 'Polarization': 2.2832363172020513}]
# tdi_fs_sky = np.fft.rfft(tdi_ts_sky['X'][start_index:-end_index or None])*dt #* np.exp(-2j*np.pi*tdi_fs_sky_f*t_start)

# def ft(samples, Fs, t0):
#     """Approximate the Fourier Transform of a time-limited signal 
#     by means of the discrete Fourier Transform.
    
#     samples: signal values sampled at the positions t0 + n/Fs
#     Fs: Sampling frequency of the signal
#     t0: starting time of the sampling of the signal
#     """
#     f = np.linspace(-Fs/2, Fs/2, len(samples), endpoint=False)
#     return np.fft.fftshift(np.fft.fft(samples)/Fs * np.exp(-2j*np.pi*f*t0))

# tdi_fs_sky = ft(tdi_ts_sky['X'][start_index:-end_index or None], 1/dt, t_start)
# tdi_fs_sky = tdi_fs_sky[int(len(tdi_fs_sky)/2):]

data_gpu_full = [xp.array(data_glitches['A'].values), xp.array(data_glitches['E'].values)]
Nmodel = get_noise_model(noise_model, data_glitches['A'].f)
SA_full_f = Nmodel.psd(freq=data_glitches['A'].f, option="A")
SE_full_f = Nmodel.psd(freq=data_glitches['E'].f, option="E")
PSD_GPU = [xp.array(SA_full_f), xp.array(SE_full_f)]
loglikelihood = gb_gpu.get_ll(params, data_gpu_full, PSD_GPU, N=N, oversample=4, dt=dt, T=t_end-t_start, start_freq_ind=0, tdi2=True, t_start=t_start, get_SNR=True)

search = Search(tdi_fs_data,t_end-t_start, lower_frequency, upper_frequency, tdi2=True, gb_gpu=gb_gpu,use_gpu=gpu_available, t_start=t_start)
search.plot(found_sources_in= [cat[cat_index]])
start = time.time()
for i in range(100):
    loglikelihood = search.loglikelihood_gpu(params, get_SNR=True)
    # loglikelihood = gb_gpu.get_ll(params, data_gpu_full, PSD_GPU, N=N, oversample=4, dt=dt, T=t_end-t_start, start_freq_ind=0, tdi2=True, t_start=t_start, get_SNR=True)
print('time', time.time()-start)


plt.figure()
# plt.plot(np.array(tdi_ts_sky['X'].t), tdi_ts_sky['X']/dt)
plt.plot(X_td.t, X_td*dt)
# plt.plot(np.array(tdi_ts_with_glitches['X'].t), tdi_ts_with_glitches['X'])
for i in range(len(data_fd)):

    pGB['InitialPhase'] = cat[cat_index]['InitialPhase'] - pGB['Frequency']*start_points[i] * np.pi*2
    params = np.array([np.full(num_bin, pGB[parameter]) for parameter in parameters_gpgpu])

    index_start = int((start_points[i]-tdi_ts_with_glitches['X'].t[0].values)/dt)
    index_end = int((end_points[i]-tdi_ts_with_glitches['X'].t[0])/dt)
    tdi_fs_sky_f = np.fft.rfftfreq(len(tdi_ts_sky['X'][index_start:index_end or None]), d=dt)
    gb_gpu.run_wave(*params,
            N=N,
            T=end_points[i]-start_points[i],
            dt=dt,
            t_start=start_points[i],
            oversample=4,
            tdi2=True)
    
    df = tdi_fs_sky_f[1]
    target_length = int(1 / df / dt)

    gb_gpu_padded = xp.zeros(len(tdi_fs_sky), dtype=xp.complex128)
    lower_index_freq = np.searchsorted(tdi_fs_sky_f,float(gb_gpu.freqs[0][0]))
    # gb_gpu_padded[lower_index_freq:lower_index_freq+len(gb_gpu.freqs[0])] = gb_gpu_padded[lower_index_freq:lower_index_freq+len(gb_gpu.freqs[0])] + gb_gpu.X[0]

    gb_gpu_padded = xp.pad(gb_gpu.X[0], (lower_index_freq, 0))
    gb_gpu_padded_ts = xp.fft.irfft(gb_gpu_padded, n=int(target_length))

    tdi_ts_sky_return_t = xp.linspace(start_points[i], start_points[i]+target_length*dt, num=target_length, endpoint=False)
    plt.plot(tdi_ts_sky_return_t.get(), gb_gpu_padded_ts.get())
    # plt.semilogx(data_fd[i]['A'].f, xp.abs(data_fd[i]['A']), label='A')
    # plt.plot(data_fd[i]['E'].f, xp.abs(data_fd[i]['E']), label='E')
plt.show()



pGBs = deepcopy(pGB)
# pGBs['InitialPhase'] = cat[cat_index]['InitialPhase'] - pGB['Frequency']*t_start * np.pi*2 - pGB['FrequencyDerivative']*t_start**2 * np.pi - pGB['Frequency2Derivative']*t_start**3 * np.pi/3
params = np.array([np.full(num_bin, pGBs[parameter]) for parameter in parameters_gpgpu])

X_fd = X_td[start_index:end_index or None].ts.fft(win=window)

tdi_fs_sky_f = np.fft.rfftfreq(len(tdi_ts_sky['X'][start_index:end_index or None]), d=dt)
gb_gpu.run_wave(*params,
        N=N,
        T=t_end-t_start,
        dt=dt,
        t_start=t_start,
        oversample=4,
        tdi2=True)


plt.figure()
plt.plot(tdi_fs_sky.f*1000, search.data_GPU[0].get())
plt.semilogx(np.array(gb_gpu.freqs[0])*1000,np.real(gb_gpu.A[0]), '.', label='GBGPU')
plt.xlim([lower_frequency*1000, upper_frequency*1000])
plt.show()

tm = gb_gpu.xp.linspace(0, Tobs-5, num=N, endpoint=False)
Ps0 = gb_gpu._spacecraft(tm)
Ps = gb_gpu._spacecraft(gb_gpu.tm)



t0 = tdi_ts_with_glitches[k].t[0].values
target_length = index_end - index_start
# xr.Dataset(dict([(k, tdi_ts[k][index_start:index_end].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
# X_fd = FrequencySeries(gb_gpu.X.get()[0], fs= gb_gpu.freqs.get()[0])
# X_td = TimeSeries(X_fd.ifft() / dt, t0=t0, dt=dt)

# Y_td = TimeSeries(np.fft.irfft(gb_gpu.Y.get()[0] / dt, n=int(Tobs)), t0=t0, dt=dt)
# Z_td = TimeSeries(np.fft.irfft(gb_gpu.Z.get()[0] / dt, n=int(Tobs)), t0=t0, dt=dt)
# gb_gpu_td = xr.Dataset({'X': X_td, 'Y': Y_td, 'Z': Z_td})
# A_td = TimeSeries(np.fft.irfft(gb_gpu.A.get()[0] / dt, n=int(Tobs)), t0=t0, dt=dt)
# E_td = TimeSeries(np.fft.irfft(gb_gpu.E.get()[0] / dt, n=int(Tobs)), t0=t0, dt=dt)
# gb_gpu_td = xr.Dataset({'A': A_td, 'E': E_td})

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
plt.plot(Xs.f*1000,np.abs(Xs.values), label='fast GB')
plt.plot(np.array(gb_gpu.freqs[0])*1000,np.abs(gb_gpu.X[0]), '.', label='GBGPU')
plt.semilogx(np.array(X_fd.f)*1000,np.abs(X_fd), label='fast GB part')
plt.semilogx(tdi_fs_with_glitches['X'].f.values*1000,np.abs(tdi_fs_with_glitches['X'].values),label='data')
# plt.semilogx(tdi_fs['X'].f.values*1000,np.abs(tdi_fs['X'].values), label='original data')
# plt.semilogx(tdi_fs_clean['X'].f.values*1000,np.abs(tdi_fs_clean['X'].values))
# plt.semilogx(tdi_fs_noisefree['X'].f.values*1000,np.abs(tdi_fs_noisefree['X'].values))
# plt.semilogx(tdi_fs_sky['X'].f.values*1000,np.abs(tdi_fs_sky['X'].values),'.',label='sky')
plt.semilogx(tdi_fs_sky_f*1000,np.abs(tdi_fs_data),'.',label='sky')
plt.xlim(lower_frequency*1000,upper_frequency*1000)
# plt.semilogx(Ef.f*1000,Ef.values)
plt.legend()
plt.show()

plt.figure()
plt.semilogx(Xs.f*1000,np.real(Xs.values), label='fast GB')
plt.semilogx(np.array(gb_gpu.freqs[0])*1000,np.real(gb_gpu.X[0]), '.', label='GBGPU')
# plt.semilogx(np.array(X_fd.f)*1000,np.real(X_fd), label='fast GB part')
# plt.semilogx(tdi_fs_sky_f*1000,np.real(tdi_fs_sky),'.',label='sky')
plt.semilogx(tdi_fs_sky['X'].f.values*1000,np.real(tdi_fs_data['X'].values),'.',label='data')
plt.semilogx(tdi_fs_sky['X'].f.values*1000,np.real(tdi_fs_sky['X'].values),'.',label='sky')
# plt.semilogx(tdi_fs_sky['X'].f.values*1000,np.real(tdi_fs_clean['X'].values),'.',label='clean')
plt.xlim(lower_frequency*1000,upper_frequency*1000)
# plt.semilogx(Ef.f*1000,Ef.values)
plt.legend()
plt.show()

gb_gpu_padded = xp.zeros(len(tdi_fs_sky['X']), dtype=np.complex128)
lower_index_freq = np.searchsorted(tdi_fs_sky_f,gb_gpu.freqs[0][0])
gb_gpu_padded[lower_index_freq:lower_index_freq+len(gb_gpu.freqs[0])] = gb_gpu_padded[lower_index_freq:lower_index_freq+len(gb_gpu.freqs[0])] + gb_gpu.X[0]

# Xs_padded = np.zeros(len(tdi_fs_clean['X']), dtype=np.complex128)
# lower_index_freq = np.searchsorted(tdi_fs_clean['X'].f,Xs.f[0])
# Xs_padded[lower_index_freq:lower_index_freq+len(gb_gpu.freqs[0])] = Xs_padded[lower_index_freq:lower_index_freq+len(gb_gpu.freqs[0])] + Xs

# gb_gpu_padded = np.zeros(len(tdi_fs_sky), dtype=np.complex128)
# lower_index_freq = np.searchsorted(tdi_fs_sky_f,gb_gpu.freqs[0][0])
# gb_gpu_padded = np.pad(gb_gpu.X[0], (lower_index_freq, 0))
# Xs_padded = np.pad(Xs, (Xs.kmin, 0))

df = tdi_fs_sky_f[1]
target_length = int(1 / df / dt)
tdi_ts_sky_return = np.fft.irfft(tdi_fs_sky['X']/dt, n=int(target_length))
tdi_ts_sky_return_t = np.linspace(t_start, t_start+target_length*dt, num=target_length, endpoint=False)
tdi_ts_sky_return_t0 = np.linspace(0, 0+target_length*dt, num=target_length, endpoint=False)
gb_gpu_padded_ts = np.fft.irfft(gb_gpu_padded, n=int(target_length))


# target_length_full = len(X_td.t)
# Xs_padded_ts = np.fft.irfft(Xs_padded / dt, n=target_length_full)

# Xs_td = np.fft.irfft(Xs/dt, n=int(Tobs/dt))
# Xs_td_t = np.linspace(0, Tobs, num=int(Tobs/dt), endpoint=False)

plt.figure()
plt.plot(X_td.t, X_td*dt)
# plt.plot(X_td.t, Xs_padded_ts*dt)
# plt.plot(tdi_ts_sky['X'].t, tdi_ts_sky['X'])
# plt.plot(tdi_ts_sky_return_t, tdi_ts_sky_return)
plt.plot(tdi_ts_sky_return_t, gb_gpu_padded_ts)
plt.show()


df = 1/Tobs
f0 = pGB["Frequency"]
kmin = int(np.round(f0/df))
FrequencySeries(gb_gpu.X[0], df=df, kmin=kmin, t0=0, name="X")
X_td = np.fft.ifft(gb_gpu.X[0])
gp_gpu_t = np.fft.ifft(gb_gpu.X[0])


print('end')