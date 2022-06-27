#%%
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Ellipse
#from getdist import plots, MCSamples
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
# from ldc.common.tools import compute_tdi_snr

from fastkde import fastKDE
from sklearn.metrics import mean_squared_error
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
from sources import *

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
SAVEPATH = grandparent+"/LDC/pictures/LDC1-3_v2"

# sangria_fn = DATAPATH + "/dgb-tdi.h5"
sangria_fn = DATAPATH + "/LDC1-3_VGB_v2.hdf5"
# sangria_fn = DATAPATH + "/LDC1-4_GB_v2.hdf5"
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

noise_model = "MRDv1"
Nmodel = get_noise_model(noise_model, np.logspace(-5, -1, 100))

pGB = {}
ind = 0
found_sources = []
target_sources = []
first_start = time.time()
np.random.seed(40) #40
number_of_signals = 1
signals_per_subtraction = 1

f = 0.0115248
def frequency_derivative(f, Mc):
    G = 6.674*10**(-11)
    c = 3*10**8
    Mc = Mc * 2*10**30
    Mc_s = Mc*G/c**3
    return 96/(5*np.pi*Mc_s**2)*(np.pi*Mc_s*f)**(11/3)
def frequency_derivative_Neil(f):
    return 8*10**8*f**(11/3)
def frequency_derivative2(f, Mc_s):
    return 96/5*np.pi**(8/3)*Mc_s**(5/3)*(f)**(11/3)
print('frequency derivative', frequency_derivative(f,0.1),frequency_derivative(f,2),' at f=', f)
chandrasekhar_limit = 1.4
M_chirp_upper_boundary = (chandrasekhar_limit**2)**(3/5)/(2*chandrasekhar_limit)**(1/5)

start_frequency = 0.0005
end_frequency = 0.02
number_of_windows = 0
current_frequency = deepcopy(start_frequency)
while current_frequency < end_frequency:
    current_frequency += 300*current_frequency * 10**3 / 10**9
    number_of_windows += 1

padding = 0.5e-6

save_name = 'LDC1-3'
indexes = np.argsort(cat['Frequency'])
cat_sorted = cat[indexes]

# LDC1-3 ##########################################
target_frequencies = cat_sorted['Frequency']
frequencies = []
window_length = 10**-6 # Hz
for i in range(len(target_frequencies)):
    window_shift = ((np.random.random(1)-0.5)*window_length*0.5)[0]
    frequencies.append([target_frequencies[i]-window_length/2+window_shift,target_frequencies[i]+window_length/2+window_shift])
# frequencies = [frequencies[2]]
frequencies_search = frequencies
do_subtract = False

do_search = False
if do_search:
    MLP = MLP_search(tdi_fs, Tobs, signals_per_window = 1, strategy = 'DE')
    start = time.time()
    pool = mp.Pool(mp.cpu_count())
    found_sources_mp = pool.starmap(MLP.search, frequencies_search)
    pool.close()
    pool.join()
    print('time to search ', number_of_windows, 'windows: ', time.time()-start)
    np.save(SAVEPATH+'/found_sources' +save_name+'.npy', found_sources_mp)
    
do_print = True
if do_print:
    found_sources_mp = np.load(SAVEPATH+'/found_sources' +save_name+'.npy', allow_pickle = True)
    found_sources_mp_best = []
    found_sources_mp_all = []
    frequencies_search = []
    for i in range(len(found_sources_mp)):
        found_sources_mp_best.append(found_sources_mp[i][0])
        found_sources_in_window = []
        for j in range(len(found_sources_mp[i][1])):
            found_sources_in_window.append(found_sources_mp[i][1][j][0][0])
        found_sources_mp_all.append(found_sources_in_window)
        frequencies_search.append(found_sources_mp[i][4])
        
    found_sources_in = []
    for i in range(len(found_sources_mp)):
        found_sources_in.append([])
        for j in range(len(found_sources_mp[i][3])):
            found_sources_in[i].append(found_sources_mp[i][3][j])

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
        pGB_injected_window = []
        pGB_stacked = {}
        for parameter in parameters:
            pGB_stacked[parameter] = cat_sorted[parameter][index_low:index_high][indexesA]
        for i in range(len(cat_sorted['Amplitude'][index_low:index_high])):
            pGBs = {}
            for parameter in parameters:
                pGBs[parameter] = pGB_stacked[parameter][i]
            pGB_injected_window.append(pGBs)
        pGB_injected.append(pGB_injected_window)

# LDC1-3 ####################
start_training_size = 1000
evalutation_times = []
training_times = []
posterior_calculation_input = []
for i in range(len(found_sources_in)):
    # if i != 5:
    #     continue
    for j in range(len(found_sources_in[i])):
        posterior_calculation_input.append((tdi_fs, Tobs, frequencies_search[i], found_sources_in[i][j], pGB_injected[i][j]))
        chain_save_name = SAVEPATH+'/Chain/frequency'+str(int(np.round(frequencies_search[i][0]*10**9)))+'nHz'+save_name+'evaluation_time_seed40.csv'
        mcmc_samples, evalutation_time, training_time = compute_posterior(tdi_fs, Tobs, frequencies_search[i], found_sources_in[i][j], pGB_injected[i][j],
                                            start_training_size, dt, noise_model, parameters, number_of_signals, GB, intrinsic_parameters, 
                                            chain_save_name, save_chain= True)
        evalutation_times.append(evalutation_time)
        training_times.append(training_time)

print('mean evaluation time', np.mean(evalutation_times), 'with ', start_training_size, 'training samples')
print('mean training time', np.mean(training_times), 'with ', start_training_size, 'training samples')

# print('time to search ', number_of_windows, 'windows: ', time.time()-start)
# start = time.time()
# pool = mp.Pool(mp.cpu_count())
# pool = mp.Pool(16)
# mcmc_samples = pool.starmap(compute_posterior,posterior_calculation_input[:16])
# pool.close()
# pool.join()
# print('time to search ', number_of_windows, 'windows: ', time.time()-start)

