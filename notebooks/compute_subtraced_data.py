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
import pickle
# sys.path.append('/cluster/home/sstrub/Repositories/LDC/lib/lib64/python3.8/site-packages/ldc-0.1-py3.8-linux-x86_64.egg')
sys.path.append('/cluster/home/sstrub/python/lib64/python3.8/site-packages/ldc-0.1-py3.8-linux-x86_64.egg')
sys.path.append('/cluster/home/sstrub/.local/lib/python3.9/site-packages')

from ldc.lisa.noise import get_noise_model
from ldc.common.series import TimeSeries
import ldc.waveform.fastGB as fastGB
from ldc.common.tools import compute_tdi_snr
try:
    from ldc.common.series import window ### manual install of  ldc
except:
    from ldc.common.tools import window ### pip install of ldc

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
if Radler:
    DATAPATH = grandparent+"/LDC/Radler/data"
    SAVEPATH = grandparent+"/LDC/pictures/LDC1-4/"
else:
    DATAPATH = grandparent+"/LDC/Sangria/data"
    SAVEPATH = grandparent+"/LDC/pictures/Sangria/"
    MBHBPATH = grandparent+"/LDC/MBHB/"

if Radler:
    sangria_fn = DATAPATH + "/LDC1-4_GB_v2.hdf5"
    # sangria_fn = DATAPATH + "/LDC1-3_VGB_v2_FD_noiseless.hdf5"
    # sangria_fn = DATAPATH + "/LDC1-3_VGB_v2.hdf5"
else:
    sangria_fn = DATAPATH + "/LDC2_sangria_training_v2.h5"
fid = h5py.File(sangria_fn)

seed = int(sys.argv[2])
reduction = 12
weeks = 13
Tobs = float((weeks)*7*24*3600)
# Tobs = float((weeks+1)*7*24*3600+24*3600)
# Tobs = float(2628000)
SNR_threshold = 9
HM = False
mbhbs_removed = bool(int(sys.argv[3]))

which_run = str(sys.argv[4])
# get TDI
if Radler:
    td = np.array(fid["H5LISA/PreProcess/TDIdata"])
    td = np.rec.fromarrays(list(td.T), names=["t", "X", "Y", "Z"])
    dt = float(np.array(fid['H5LISA/GWSources/GalBinaries']['Cadence']))
    # Tobs = float(int(np.array(fid['H5LISA/GWSources/GalBinaries']['ObservationDuration']))/reduction)
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
        
    # Tobs = float(int(np.array(fid['obs/config/t_max'])))
    if mbhbs_removed:
        for k in ["X", "Y", "Z"]:
            td[k] = td[k] - td_mbhb[k]
            # td_injected[k] -= td_injected[k]
    else:
        td_original = deepcopy(td)
        # if reduction == 2:
        #     # wave = pickle.load(open(MBHBPATH+dataset+"_mbhbh_found_6months.pkl", "rb"))
        #     wave = pickle.load(open(MBHBPATH+'Sangria_mbhbh_found_6months_seed'+str(seed)+'.pkl', "rb"))
        #     # wave = pickle.load(open(MBHBPATH+"Sangria_mbhbh_found_6months.pkl", "rb"))
        # else:
        #     wave = pickle.load(open(MBHBPATH+'Sangria_mbhbh_found_12months_seed'+str(seed)+'.pkl', "rb"))
        # if HM:
        #     wave = pickle.load(open(MBHBPATH+"Sangria_mbhbh_HM_found.pkl", "rb"))

        wave = pickle.load(open(MBHBPATH+'Sangria_mbhbh_found_'+str(weeks)+'w_seed'+str(seed)+'.pkl', "rb"))
        for i, k in enumerate(["X", "Y", "Z"]):
            # td[k] = td_mbhb[k]
            td[k] -= wave[k] 

# Build timeseries and frequencyseries object for X,Y,Z
t_max_index = np.searchsorted(td['t'], Tobs)
tdi_ts = dict([(k, TimeSeries(td[k][:t_max_index], dt=dt, t0=td.t[0])) for k in ["X", "Y", "Z"]])
# tdi_ts = dict([(k, TimeSeries(td[k][:int(len(td[k][:])/reduction)], dt=dt, t0=td.t[0])) for k in ["X", "Y", "Z"]])
tdi_fs = xr.Dataset(dict([(k, tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
GB = fastGB.FastGB(delta_t=dt, T=Tobs)  # in seconds


# plt.figure()
# plt.plot(td_original.t, td_original['X'])
# plt.plot(tdi_ts['X'].t, tdi_ts['X'])
# plt.show()

pGBadded20 = {}
pGBadded20['Amplitude'] = 5e-21
pGBadded20['EclipticLatitude'] = 0.2
pGBadded20['EclipticLongitude'] = 1.5
pGBadded20['Frequency'] = 0.00031
pGBadded20['FrequencyDerivative'] = 5*1e-20
pGBadded20['Inclination'] = 1.2
pGBadded20['InitialPhase'] = 3
pGBadded20['Polarization'] = 2


noise_model = "SciRDv1"
Nmodel = get_noise_model(noise_model, np.logspace(-5, -1, 100))

pGB = {}
ind = 0
found_sources = []
target_sources = []
first_start = time.time()
np.random.seed(seed) #40
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

chandrasekhar_limit = 1.4
M_chirp_upper_boundary = (chandrasekhar_limit**2)**(3/5)/(2*chandrasekhar_limit)**(1/5)

def tdi_subtraction(tdi_fs,found_sources_mp_subtract, frequencies_search):

    #subtract the found sources from original
    tdi_fs_subtracted2 = deepcopy(tdi_fs)
    for i in range(len(found_sources_mp_subtract)):
        # for j in range(len(found_sources_to_subtract[i])):
        freqs = None
        if GB.T < 365.26*24*3600/4:
            freqs = np.linspace(found_sources_mp_subtract[i]['Frequency'], found_sources_mp_subtract[i]['Frequency']+0.00001, 128)
        # freqs = np.array(search1.dataX.f)
        # freqs = None
        # if len(freqs) < 16:
        #     freqs = np.linspace(found_sources_mp_subtract[i]['Frequency'], found_sources_mp_subtract[i]['Frequency']+0.00001, 16)
        Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_mp_subtract[i], oversample=4, freqs=freqs)
        source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
        index_low = np.searchsorted(tdi_fs_subtracted2["X"].f, Xs_subtracted.f[0])
        index_high = index_low+len(Xs_subtracted)
        for k in ["X", "Y", "Z"]:
            tdi_fs_subtracted2[k].data[index_low:index_high] -= source_subtracted[k].data
    return tdi_fs_subtracted2


# sum the found sources
if mbhbs_removed:
    found_sources = np.load(SAVEPATH+'found_sources_not_anticorrelated_Sangria_12m_no_mbhb_SNR9_seed'+str(seed)+'.pkl', allow_pickle = True)
else:
    found_sources = np.load(SAVEPATH+'found_sources_not_anticorrelated_original_Sangria_'+str(weeks)+'w_mbhb_SNR9_seed'+str(seed)+'.pkl', allow_pickle = True)
found_sources_flat = np.concatenate(found_sources)
# found_sources_flat = np.load(SAVEPATH+'found_sources_Sangria_6m_mbhb_even3_seed1_flat.pkl', allow_pickle = True)
tdi_fs_sum_found = deepcopy(tdi_fs)
for k in ["X", "Y", "Z"]:
    tdi_fs_sum_found[k].data = np.zeros(len(tdi_fs_sum_found[k].data), np.complex128)
for i in range(len(found_sources_flat)):
    # for j in range(len(found_sources_to_subtract[i])):
    freqs = None    
    if GB.T < 365.26*24*3600/4:
            freqs = np.linspace(found_sources_flat[i]['Frequency'], found_sources_flat[i]['Frequency']+0.00001, 16)
    Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_flat[i], oversample=4, freqs=freqs)
    source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
    index_low = np.searchsorted(tdi_fs_sum_found["X"].f, Xs_subtracted.f[0])
    index_high = index_low+len(Xs_subtracted)
    for k in ["X", "Y", "Z"]:
        tdi_fs_sum_found[k].data[index_low:index_high] += source_subtracted[k].data

tdi_fs_subtracted = deepcopy(tdi_fs)
for k in ["X", "Y", "Z"]:
    tdi_fs_subtracted[k].data -= tdi_fs_sum_found[k].data

# plt.figure()
# plt.semilogx(tdi_fs['X'].f, (tdi_fs['X'].data), label = 'original')
# plt.semilogx(tdi_fs_sum_found['X'].f, (tdi_fs_sum_found['X'].data), label = 'sum found')
# plt.semilogx(tdi_fs_subtracted['X'].f, (tdi_fs_subtracted['X'].data), label = 'subtracted')
# plt.legend()
# plt.show()

# plt.figure()
# plt.loglog(tdi_fs['X'].f, np.abs(tdi_fs['X'].data), label = 'original')
# plt.loglog(tdi_fs_sum_found['X'].f, np.abs(tdi_fs_sum_found['X'].data), label = 'sum found')
# plt.loglog(tdi_fs_subtracted['X'].f, np.abs(tdi_fs_subtracted['X'].data),   label = 'subtracted')
# plt.legend()
# plt.savefig(SAVEPATH+'tdi_f_subtracted/tdi_fs_subtracted_found_GB_'+str(weeks)+'w_mbhb_SNR9_seed'+str(seed)+'.png')
# plt.show()

pickle.dump(tdi_fs_subtracted, open(MBHBPATH+'Sangria_tdi_fs_'+str(weeks)+'w_residual.pkl', "wb"))
