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
weeks = 3
Tobs = float(weeks*7*24*3600)
# Tobs = float(2628000)
SNR_threshold = 9
HM = False
mbhbs_removed = bool(int(sys.argv[3]))

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
        
    # Tobs = float(int(np.array(fid['obs/config/t_max']))/reduction)
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
        # for i, k in enumerate(["X", "Y", "Z"]):
        #     # td[k] = td_mbhb[k]
        #     td[k] -= wave[k] 

# Build timeseries and frequencyseries object for X,Y,Z
t_max_index = np.searchsorted(td['t'], Tobs)
tdi_ts = dict([(k, TimeSeries(td[k][:t_max_index], dt=dt, t0=td.t[0])) for k in ["X", "Y", "Z"]])
# tdi_ts = dict([(k, TimeSeries(td[k][:int(len(td[k][:])/reduction)], dt=dt, t0=td.t[0])) for k in ["X", "Y", "Z"]])
tdi_fs = xr.Dataset(dict([(k, tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
GB = fastGB.FastGB(delta_t=dt, T=Tobs)  # in seconds

pGBadded20 = {}
pGBadded20['Amplitude'] = 5e-21
pGBadded20['EclipticLatitude'] = 0.2
pGBadded20['EclipticLongitude'] = 1.5
pGBadded20['Frequency'] = 0.00031
pGBadded20['FrequencyDerivative'] = 5*1e-20
pGBadded20['Inclination'] = 1.2
pGBadded20['InitialPhase'] = 3
pGBadded20['Polarization'] = 2

# add signal
add_signal = False
if add_signal:
    for pGBadding in [pGBadded20]: # faint
    # for pGBadding in [pGBadded11]:              # overlap single signal
        cat = np.hstack((cat,cat[0]))
        for parameter in parameters:
            cat[-1][parameter] = pGBadding[parameter]
        Xs_added, Ys_added, Zs_added = GB.get_fd_tdixyz(template=pGBadding, oversample=4)
        source_added = dict({"X": Xs_added, "Y": Ys_added, "Z": Zs_added})
        index_low = np.searchsorted(tdi_fs["X"].f, Xs_added.f[0])
        index_high = index_low+len(Xs_added)
        # tdi_fs['X'] = tdi_fs['X'] #+ Xs_added
        for k in ["X", "Y", "Z"]:
            tdi_fs[k].data[index_low:index_high] = tdi_fs[k].data[index_low:index_high] + source_added[k].data
    tdi_ts = xr.Dataset(dict([(k, tdi_fs[k].ts.ifft(dt=dt)) for k, n in [["X", 1], ["Y", 2], ["Z", 3]]]))

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
        print('start search', lower_frequency, upper_frequency)
        start_search = time.time()
        initial_guess = []
        if len(self.found_sources_previous) > 0:
            search1 = Search(tdi_fs_search,self.Tobs, lower_frequency, upper_frequency)
            if do_subtract:
                padding_of_initial_guess_range = 0
            else:
                padding_of_initial_guess_range = search1.padding
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
            
            ### sort the initial guesses such that the highest SNR guess comes first
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
        loglikelihood_ratio_threshold = 50
        f_transfer = 19.1*10**-3
        # if lower_frequency > f_transfer:
        #     loglikelihood_ratio_threshold = 200
        # if lower_frequency > 10*10**-3:
        #     SNR_threshold = 15
        # if lower_frequency > 15*10**-3:
        #     SNR_threshold = 20
        # if lower_frequency > 20*10**-3:
        #     SNR_threshold = 100
        if lower_frequency > 10**-2:
            self.signals_per_window = 2
        # current_loglikelihood_ratio = 1000
        ind = 0
        while current_SNR > SNR_threshold and ind < self.signals_per_window:
        # while current_loglikelihood_ratio > loglikelihood_ratio_threshold and ind < self.signals_per_window:
            ind += 1
            search1 = Search(tdi_fs_search,self.Tobs, lower_frequency, upper_frequency)
            # search1.update_noise()
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
                print('which signal per window', ind-1,'and repetition:', i)
                if i == 0:
                    current_SNR = deepcopy(new_SNR)
                    maxpGBsearch = deepcopy(maxpGBsearch_new)
                if new_SNR >= current_SNR:
                    current_SNR = deepcopy(new_SNR)
                    maxpGBsearch = deepcopy(maxpGBsearch_new)
                print('current SNR', np.round(current_SNR,3))
                found_sources_all[-1] = maxpGBsearch_new
                if current_SNR < SNR_threshold-2:
                    break

            if current_SNR < SNR_threshold:
                break

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

            # current_loglikelihood_ratio = search1.loglikelihood(maxpGBsearch[0])
            # print('current loglikelihood ratio', current_loglikelihood_ratio)
            current_SNR = search1.SNR(maxpGBsearch[0])
            print('current SNR ratio', current_SNR)



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
                    freqs = None
                    if GB.T < 365.26*24*3600/4:
                        freqs = np.linspace(found_sources_out[i]['Frequency'], found_sources_out[i]['Frequency']+0.00001, 16)
                    Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_out[i], oversample=4, freqs=freqs)
                    source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
                    index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
                    index_high = index_low+len(Xs_subtracted)
                    for k in ["X", "Y", "Z"]:
                        tdi_fs_subtracted[k].data[index_low:index_high] = tdi_fs_subtracted[k].data[index_low:index_high] - source_subtracted[k].data

                search_out_subtracted = Search(tdi_fs_subtracted,self.Tobs, lower_frequency, upper_frequency)
                # search_out_subtracted.update_noise()
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
                freqs = None
                if GB.T < 365.26*24*3600/4:
                    freqs = np.linspace(found_sources[i]['Frequency'], found_sources[i]['Frequency']+0.00001, 16)
                Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources[i], oversample=4, freqs=freqs)
                source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
                index_low = np.searchsorted(tdi_fs_search["X"].f, Xs_subtracted.f[0])
                index_high = index_low+len(Xs_subtracted)
                for k in ["X", "Y", "Z"]:
                    tdi_fs_search[k].data[index_low:index_high] = tdi_fs_search[k].data[index_low:index_high] - source_subtracted[k].data
        print('search time', time.time()-start_search, 'frequency', lower_frequency, upper_frequency)
        print('found_sources_in',found_sources_in)
        found_sources_mp = [[found_sources, found_sources_all, number_of_evaluations_all, found_sources_in, [lower_frequency, upper_frequency], time.time()-start_search]]
        fn = SAVEPATH+'found_signals/found_sources'+ str(int(np.round(lower_frequency*10**9)))+'nHz_to'+ str(int(np.round(upper_frequency*10**9)))+'nHz_' +save_name+'.pkl'
        # pickle.dump(found_sources_mp, open(fn, "wb"))
        return found_sources, found_sources_all, number_of_evaluations_all, found_sources_in, [lower_frequency, upper_frequency], time.time()-start_search

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
            freqs = None
            if GB.T < 365.26*24*3600/4:
                freqs = np.linspace(found_sources_mp_subtract[i]['Frequency'], found_sources_mp_subtract[i]['Frequency']+0.00001, 16)
            Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_mp_subtract[i], oversample=4, freqs=freqs)
            source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
            index_low = np.searchsorted(tdi_fs_subtracted2["X"].f, Xs_subtracted.f[0])
            index_high = index_low+len(Xs_subtracted)
            for k in ["X", "Y", "Z"]:
                tdi_fs_subtracted2[k].data[index_low:index_high] -= source_subtracted[k].data
    return tdi_fs_subtracted2



#sum the found sources
# if mbhbs_removed:
#     found_sources = np.load(SAVEPATH+'found_sources_not_anticorrelated_Sangria_12m_no_mbhb_SNR9_seed'+str(seed)+'.pkl', allow_pickle = True)
# else:
#     found_sources = np.load(SAVEPATH+'found_sources_not_anticorrelated_original_Sangria_12m_mbhb_SNR9_seed'+str(seed)+'.pkl', allow_pickle = True)
# found_sources_flat = np.concatenate(found_sources)
# # found_sources_flat = np.load(SAVEPATH+'found_sources_Sangria_6m_mbhb_even3_seed1_flat.pkl', allow_pickle = True)
# tdi_fs_sum_found = deepcopy(tdi_fs)
# for k in ["X", "Y", "Z"]:
#     tdi_fs_sum_found[k].data = np.zeros(len(tdi_fs_sum_found[k].data), np.complex128)
# for i in range(len(found_sources_flat)):
#     # for j in range(len(found_sources_to_subtract[i])):
#     Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_flat[i], oversample=4)
#     source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
#     index_low = np.searchsorted(tdi_fs_sum_found["X"].f, Xs_subtracted.f[0])
#     index_high = index_low+len(Xs_subtracted)
#     for k in ["X", "Y", "Z"]:
#         tdi_fs_sum_found[k].data[index_low:index_high] += source_subtracted[k].data
# if mbhbs_removed:
#     pickle.dump(tdi_fs_sum_found, open(SAVEPATH+'tdi_fs_sum_found_12m_no_mbhb_SNR9_seed'+str(seed)+'.pkl', "wb"))
# else:
#     pickle.dump(tdi_fs_sum_found, open(SAVEPATH+'tdi_fs_sum_found_12m_mbhb_SNR9_seed'+str(seed)+'.pkl', "wb"))
        

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
#         cat = np.rec.fromarrays(params, names=list(names))
#     # Sangria
#     else:
#         names_dgb = fid["sky/dgb/cat"].dtype.names # Sangria
#         params_dgb = [np.array(fid["sky/dgb/cat"][k]).squeeze() for k in names_dgb]
#         names_igb = fid["sky/igb/cat"].dtype.names # Sangria
#         params_igb = [np.array(fid["sky/igb/cat"][k]).squeeze() for k in names_igb]
#         names_vgb = fid["sky/vgb/cat"].dtype.names # Sangria
#         params_vgb = [np.array(fid["sky/vgb/cat"][k]).squeeze() for k in names_vgb]

#         cat = np.rec.fromarrays(params_dgb, names=list(names_dgb))
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
if Radler:
    search_range = [0.0003, 0.0319]
# search_range = [0.0001, 0.11]
# search_range = [0.0019935, 0.0020135]
# search_range = [0.0029935, 0.0030135]
# window_length = 1*10**-7 # Hz

frequencies = create_frequency_windows(search_range, Tobs)

frequencies_even = frequencies[::2]
frequencies_odd = frequencies[1::2]

f_bins_per_segment = []
for i in range(len(frequencies)):
    index_low = np.searchsorted(tdi_fs["X"].f, frequencies[i][0])
    index_high = np.searchsorted(tdi_fs["X"].f, frequencies[i][1])
    f_bins_per_segment.append(index_high-index_low)
# counts = np.zeros(len(pGB_injected))
# for i in range(len(pGB_injected)):
#     counts[i] = len(pGB_injected[i])

# frequencies_search = np.asarray(frequencies)
# figure = plt.figure()
# plt.loglog(frequencies_search[:,1],counts, '.')
# plt.xlabel('Frequency [Hz]')
# plt.ylabel('Number of signals')
# plt.show()


# figure = plt.figure()
# plt.loglog(frequencies_search[:,0],frequencies_search[:,1]-frequencies_search[:,0],  linewidth= 4, label= '$B$')
# plt.loglog(frequencies_search[:,0],frequency_derivative(frequencies_search[:,0],2)*Tobs, label= '$B_{F}$')
# plt.loglog(frequencies_search[:,0],frequencies_search[:,0]*3* 10**-4, label= '$3 \cdot B_{O}$')
# plt.loglog(frequencies_search[:,0],np.ones(len(frequencies_search[:,0]))*4*32*10**-9*2, label= '$2 \cdot B_{C}$')
# plt.xlabel('Frequency [Hz]')
# plt.ylabel('Frequency window witdh [Hz]')
# plt.xlim(search_range[0],0.1)
# plt.ylim(bottom=(frequencies_search[0,1]-frequencies_search[0,0])/10**1)
# plt.legend()
# plt.show()
# plt.savefig(SAVEPATH+'bandwidth.png')


if mbhbs_removed:
    save_name = 'Sangria_1w_no_mbhb_SNR'+str(SNR_threshold)+'_odd'
else:
    save_name = 'Sangria_'+str(weeks)+'w_mbhb_SNR'+str(SNR_threshold)+'_even'
frequencies_search = frequencies_even
frequencies_search_full = deepcopy(frequencies_search)
batch_index = int(sys.argv[1])
# batch_index = int(150)
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
batch_size = int(10)
start_index = int(batch_size*batch_index)
print('batch',batch_index, start_index, batch_size)
frequencies_search = frequencies_search[start_index:start_index+batch_size]
# frequencies_search = [frequencies_search[1359],frequencies_search[1370],frequencies_search[1419],frequencies_search[1430],frequencies_search[1659],frequencies_search[1670]]

# print(i, frequencies_search[0])
while frequencies_search[-1][1] + (frequencies_search[-1][1] - frequencies_search[-1][0])/2 > f_Nyquist:
    frequencies_search = frequencies_search[:-1]
# frequencies_search = frequencies_search[70:80]
# frequencies_search = frequencies_search[25:]

search_range = [frequencies_search[0][0],frequencies_search[-1][1]]
# search_range = [1619472*10**-8,2689639*10**-8]
print('search range '+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))))

# save_name_previous = 'found_sources_Radler_12m_even10_first'
# found_sources_mp_subtract = np.load(SAVEPATH+save_name_previous+'.npy', allow_pickle = True)
# frequencies_search_reduced = []
# for i in range(len(found_sources_mp_subtract)):
#     if i in range(start_index,start_index+batch_size):
#         if len(found_sources_mp_subtract[i][0]) > 3:
#             frequencies_search_reduced.append(frequencies_search_full[i])
# frequencies_search = frequencies_search_reduced

# frequencies_search = frequencies_even
do_subtract = True
if do_subtract:
    start = time.time()
    # save_name_previous = 'found_sourcesRadler_half_odd_dynamic_noise'
    # Sangria
    if mbhbs_removed:
        save_name_previous = 'found_sources_original_Sangria_6m_no_mbhb_SNR9_odd_seed'+str(seed)+'_flat'
    else:
        save_name_previous = 'found_sources_original_Sangria_'+str(weeks)+'w_mbhb_SNR9_odd_seed'+str(seed)+'_flat'
    # save_name_previous = 'found_sources_Sangria_12m_even3'
    # save_name_previous = 'found_sources_Radler_12m_odd'
    # save_name_previous = 'found_sources_Radler_half_odd_dynamic_noise'
    # save_name_previous = 'found_sources_Sangria_1_odd_dynamic_noise'
    # save_name_previous = 'found_sourcesSangria_half_odd'
    # save_name_previous = 'found_sourcesSangria_1_full'
    # save_name_previous = 'found_signals_1_even/found_sources1630803to2001429Sangria_1_even'
    # save_name_previous = 'found_sourcesSangria_odd'
    # save_name_previous = 'found_sourcesRadler_1_even3'
    # save_name_previous = 'found_sourcesRadler_1_odd'
    # save_name_previous = 'found_sourcesSangria_1_odd_dynamic_noise'
    # save_name_previous = 'found_sources'+save_name
    found_sources_mp_subtract = np.load(SAVEPATH+save_name_previous+'.pkl', allow_pickle = True)

    found_sources_flat = []
    for i in range(len(found_sources_mp_subtract)):
        # for j in range(len(found_sources_mp_subtract[i][3])):
        #     found_sources_flat.append(found_sources_mp_subtract[i][3][j])
        # for j in range(len(found_sources_mp_subtract[i])):
        #     found_sources_flat.append(found_sources_mp_subtract[i][j])
        found_sources_flat.append(found_sources_mp_subtract[i])
    found_sources_flat = np.asarray(found_sources_flat)
    found_sources_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_flat]) for attribute in found_sources_flat[0].keys()}
    found_sources_flat_df = pd.DataFrame(found_sources_flat_array)
    found_sources_out_flat_df = found_sources_flat_df.sort_values('Frequency')
    for i in range(len(frequencies_search_full)):
        found_sources_out_flat_df = found_sources_out_flat_df[(found_sources_out_flat_df['Frequency']< frequencies_search_full[i][0]) | (found_sources_out_flat_df['Frequency']> frequencies_search_full[i][1])]
    found_sources_out_flat_df = found_sources_out_flat_df.sort_values('Frequency')
    found_sources_out_flat = found_sources_out_flat_df.to_dict(orient='records')
    tdi_fs_subtracted = tdi_subtraction(tdi_fs,found_sources_out_flat, frequencies_search_full)

    print('subtraction time', time.time()-start)
    plot_subtraction = False
    if plot_subtraction:
        i = 6
        # lower_frequency = frequencies_search_full[start_index+i][0]
        # upper_frequency = frequencies_search_full[start_index+i][1]
        lower_frequency = frequencies_search_full[start_index+i][1]
        upper_frequency = frequencies_search_full[start_index+i+1][0]
        found_sources_neighbor = found_sources_flat[(found_sources_flat_df['Frequency']> frequencies_search_full[start_index+i-1][0]) & (found_sources_flat_df['Frequency']< frequencies_search_full[start_index+i+1][1])]
        search1 = Search(tdi_fs,Tobs, lower_frequency, upper_frequency)
        search1.plot(second_data= tdi_fs_subtracted, found_sources_in=found_sources_neighbor)
        # search1.plot()
        
    tdi_fs = deepcopy(tdi_fs_subtracted)

do_not_search_unchanged_even_windows = True
if do_not_search_unchanged_even_windows:
    frequencies_search_reduced = []
    frequencies_search_skipped = []

    if mbhbs_removed:
        save_name_previous = 'found_sources_original_Sangria_'+str(weeks)+'w_no_mbhb_SNR9_even3_seed'+str(seed)+'_flat'
    else:
        save_name_previous = 'found_sources_original_Sangria_'+str(weeks)+'w_mbhb_SNR9_even3_seed'+str(seed)+'_flat'
    found_sources_mp_previous = np.load(SAVEPATH+save_name_previous+'.pkl', allow_pickle = True)
    found_sources_flat = np.asarray(found_sources_mp_previous)
    found_sources_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_flat]) for attribute in found_sources_flat[0].keys()}
    found_sources_flat_df = pd.DataFrame(found_sources_flat_array)
    found_sources_in_skipped = []
    found_sources_in_not_skipped = []
    for i in range(len(frequencies_search)):
        found_sources_out_lower = []
        found_sources_out_upper = []
        try:
            found_sources_out_lower = found_sources_out_flat_df[(found_sources_out_flat_df['Frequency']> frequencies_search[i-1][1]) & (found_sources_out_flat_df['Frequency']< frequencies_search[i][0])]
        except:
            pass
        try:
            found_sources_out_upper = found_sources_out_flat_df[(found_sources_out_flat_df['Frequency']> frequencies_search[i][1]) & (found_sources_out_flat_df['Frequency']< frequencies_search[i+1][0])]
        except:
            pass
        if not(len(found_sources_out_lower) == 0 and len(found_sources_out_upper) == 0):
            frequencies_search_reduced.append(frequencies_search[i])
            found_sources_in_flat_df = found_sources_flat_df[(found_sources_flat_df['Frequency']> frequencies_search[i][0]) & (found_sources_flat_df['Frequency']< frequencies_search[i][1])]
            # found_sources_in_not_skipped.append(found_sources_in_flat_df.to_dict(orient='records'))
        else:
            found_sources_in_flat_df = found_sources_flat_df[(found_sources_flat_df['Frequency']> frequencies_search[i][0]) & (found_sources_flat_df['Frequency']< frequencies_search[i][1])]
            if len(found_sources_in_flat_df) > 2:
                frequencies_search_reduced.append(frequencies_search[i])
                # found_sources_in_not_skipped.append(found_sources_in_flat_df.to_dict(orient='records'))
            else:
                # print('no search', frequencies_search[i])
                frequencies_search_skipped.append(frequencies_search[i])
                found_sources_in_skipped.append(found_sources_in_flat_df.to_dict(orient='records'))
    found_sources_in_skipped = np.concatenate(found_sources_in_skipped)
    # found_sources_in_not_skipped = np.concatenate(found_sources_in_not_skipped)
    # pickle.dump(found_sources_in_skipped, open(SAVEPATH+save_name_previous+'_skipped.pkl', "wb"))
    frequencies_search = frequencies_search_reduced
    # frequencies_search = frequencies_search_skipped

found_sources_sorted = []
use_initial_guess = True
if use_initial_guess:
    # save_name_found_sources_previous = 'found_sources397769to400619LDC1-4_4mHz_half_year_even10'
    # save_name_found_sources_previous = 'found_sources397919to400770LDC1-4_4mHz_half_year_odd'
    # save_name_found_sources_previous = 'found_sources2537595to3305084LDC1-4_4mHz_half_year_even'
    # save_name_found_sources_previous = 'found_sourcesLDC1-4_half_even10'
    if mbhbs_removed:
        save_name_found_sources_previous = 'found_sources_not_anticorrelated_original_Sangria_6m_no_mbhb_SNR9_seed'+str(seed)
    else:
        save_name_found_sources_previous = 'found_sources_not_anticorrelated_original_Sangria_'+str(weeks-1)+'w_mbhb_SNR9_seed'+str(seed)
    # save_name_found_sources_previous = 'found_sourcesLDC1-4_half_odd'
    # found_sources_mp_subtract = np.load(SAVEPATH+save_name_found_sources_previous+'.npy', allow_pickle = True)

    # found_sources_flat = []
    # for i in range(len(found_sources_mp_subtract)):
    #     for j in range(len(found_sources_mp_subtract[i][3])):
    #         found_sources_flat.append(found_sources_mp_subtract[i][3][j])
    # found_sources_flat = np.asarray(found_sources_flat)
    # found_sources_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_flat]) for attribute in found_sources_flat[0].keys()}
    # found_sources_flat_df = pd.DataFrame(found_sources_flat_array)
    # found_sources_flat_df = found_sources_flat_df.sort_values('Frequency')
    # found_sources_sorted = found_sources_flat_df.to_dict(orient='records')

    found_sources_loaded = []
    found_sources_loaded.append(np.load(SAVEPATH+save_name_found_sources_previous+'.pkl', allow_pickle = True))
    found_sources_previous = []
    for i in range(len(found_sources_loaded)):
        for j in range(len(found_sources_loaded[i])):
            for k in range(len(found_sources_loaded[i][j])):
                found_sources_previous.append(found_sources_loaded[i][j][k])

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

# index =0
# search1 = Search(tdi_fs,Tobs, frequencies_search[index][0], frequencies_search[index][1])
# search1.plot()
# print(search1.SNR([pGB_injected[index][0]]))
# search1.update_noise()
# print(search1.SNR([found_sources_mp_loaded[index][0][2]]))
# search1.plot(found_sources_in=found_sources_mp_loaded[index][0], pGB_injected=pGB_injected[index])
# search1.SNR(pGB_injected[index])

do_search = True
if do_search:
    MLP = MLP_search(tdi_fs, Tobs, signals_per_window = 10, found_sources_previous = found_sources_sorted, strategy = 'DE')
    start = time.time()

    # cpu_cores = 16
    # pool = mp.Pool(cpu_cores)
    # pool.starmap(MLP.search, frequencies_search)
    # # found_sources_mp = pool.starmap(MLP.search, frequencies_search)
    # pool.close()
    # pool.join()

    found_sources_mp = []
    for i in range(len(frequencies_search)):
        found_sources_mp.append(MLP.search(*frequencies_search[i]))

    print('time to search ', len(frequencies_search), 'windows: ', time.time()-start)
    if mbhbs_removed:
        directory = SAVEPATH+'found_signals_original_'+save_name+'_seed'+str(seed)
    else:
        directory = SAVEPATH+'found_signals_original_'+save_name+'_seed'+str(seed)
    fn = directory+'/found_sources_batch_index_'+str(batch_index)+'_'+ str(int(np.round(search_range[0]*10**9)))+'nHz_to'+ str(int(np.round(search_range[1]*10**9)))+'nHz_' +save_name+'.pkl'
    if not os.path.exists(directory):
        os.makedirs(directory)
    pickle.dump(found_sources_mp, open(fn, "wb"))
    
    # np.save(SAVEPATH+'found_signals/found_sources'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', found_sources_mp, allow_pickle=True)

final_optimization = False
if final_optimization:
    # found_sources_mp = np.load(SAVEPATH+'found_sources'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'noise_matrix.npy', allow_pickle = True)
    # found_sources_mp = np.load(SAVEPATH+'found_sourcesLDC1-4_2_years_full.npy', allow_pickle = True)
    # found_sources_mp = np.load(SAVEPATH+'found_sources387812to408573LDC1-4_2year_oddnoise_matrix.npy', allow_pickle = True)
    # found_sources_mp = np.load(SAVEPATH+'found_sourcesLDC1-4_2_odd_optimized.npy', allow_pickle = True)
    # found_sources_mp = np.load(SAVEPATH+'found_sourcesLDC1-4_2_even_optimized.npy', allow_pickle = True)
    # found_sources_mp = np.load(SAVEPATH+'found_sourcesLDC1-4_2_optimized.npy', allow_pickle = True)
    found_sources_mp_loaded = pickle.load(open(SAVEPATH+'found_signals/found_sources'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.pkl', 'rb'))
    # found_sources_mp = np.load(SAVEPATH+'found_sources'+save_name+'.npy', allow_pickle = True)
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
    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB_found, oversample=4)
    Xs_injected, Ys_injected, Zs_injected = GB.get_fd_tdixyz(template=pGB_injected, oversample=4)
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
    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB_found, oversample=4)
    Xs_injected, Ys_injected, Zs_injected = GB.get_fd_tdixyz(template=pGB_injected, oversample=4)
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
            Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_out[i][j], oversample=4)
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
                Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=pGB_injected[i][n], oversample=4)
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
                Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_in[i][n], oversample=4)
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
#             Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_in[i][j], oversample=4)
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
#             Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_in[i][j], oversample=4)
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
# f, psd_x_noisy = scipy.signal.welch(tdi_ts["X"], fs=1.0/dt, nperseg=nperseg)
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