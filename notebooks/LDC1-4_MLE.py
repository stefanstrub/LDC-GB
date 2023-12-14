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


# from bbhx.waveformbuild import BBHWaveformFD
cluster = True
from ldc.lisa.noise import get_noise_model
from ldc.common.series import TimeSeries#, TDI
import ldc.waveform.fastGB as fastGB
if cluster:
    from ldc.common.series import window ### manual install of  ldc
else:
    from ldc.common.tools import window ### pip install of ldc

from sources import *

gpu = False
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

seed = int(sys.argv[2])
np.random.seed(seed) #40

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
dataset = 'Sangria'
# dataset = 'Spritz'
VGB = False
add_gaps = False
zero_gaps = False
zero_glitches = False
fill_gaps = False
tdi2 = False
if dataset == 'Radler':
    DATAPATH = grandparent+"/LDC/Radler/data"
    SAVEPATH = grandparent+"/LDC/pictures/LDC1-4/"
elif dataset == 'Sangria':
    DATAPATH = grandparent+"/LDC/Sangria/data"
    SAVEPATH = grandparent+"/LDC/pictures/Sangria/"
    MBHBPATH = grandparent+"/LDC/MBHB/"
elif dataset == 'Spritz':
    DATAPATH = grandparent+"/LDC/Spritz/data"
    SAVEPATH = grandparent+"/LDC/Spritz/"

if dataset == 'Radler':
    data_fn = DATAPATH + "/LDC1-4_GB_v2.hdf5"
    # data_fn = DATAPATH + "/LDC1-3_VGB_v2_FD_noiseless.hdf5"
    if VGB:
        data_fn = DATAPATH + "/LDC1-3_VGB_v2.hdf5"
elif dataset == 'Sangria':
    data_fn = DATAPATH + "/LDC2_sangria_training_v2.h5"
elif dataset == 'Spritz':
    data_fn = DATAPATH + "/LDC2_spritz_vgb_training_v2.h5"
fid = h5py.File(data_fn)

reduction = 1
HM = False
mbhbs_removed = bool(int(sys.argv[3]))

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
    
    mbhb = fid["sky/mbhb/cat"]
    cat_mbhb = []
    for index in range(len(mbhb)):
        pMBHB = dict(zip(mbhb.dtype.names, mbhb[index]))
        for i, name in enumerate(mbhb.dtype.names):
            pMBHB[name] = mbhb[index][name][0]
        cat_mbhb.append(pMBHB)

    td_mbhb = fid["sky/mbhb/tdi"][()]
    # cat_mbhb = fid["sky/mbhb/cat"]
    td_mbhb  = np.rec.fromarrays(list(td_mbhb .T), names=["t", "X", "Y", "Z"])
    td_mbhb  = td_mbhb ['t']

    # td_injected = deepcopy(td)
    Tobs = float(int(np.array(fid['obs/config/t_max']))/reduction)
    if mbhbs_removed:
        for k in ["X", "Y", "Z"]:
            td[k] = td[k] - td_mbhb[k]
            # td_injected[k] -= td_injected[k]
    else:
        td_original = deepcopy(td)
        if reduction == 2:
            # wave = pickle.load(open(MBHBPATH+dataset+"_mbhbh_found_6months.pkl", "rb"))
            wave = pickle.load(open(MBHBPATH+dataset+'_mbhbh_found_6months_seed'+str(seed)+'.pkl', "rb"))
            # wave = pickle.load(open(MBHBPATH+dataset+"_mbhbh_found_6months.pkl", "rb"))
        else:
            wave = pickle.load(open(MBHBPATH+dataset+"_mbhbh_found_12months_seed12.pkl", "rb"))
        if HM:
            wave = pickle.load(open(MBHBPATH+dataset+"_mbhbh_HM_found.pkl", "rb"))
        for i, k in enumerate(["X", "Y", "Z"]):
            # td[k] = td_mbhb[k]
            td[k] -= wave[k] 

    # td_mbhbs_subtracted = deepcopy(td)

    # wave = pickle.load(open(MBHBPATH+dataset+"_mbhbh_injected.pkl", "rb"))
    # if HM:
    #     wave = pickle.load(open(MBHBPATH+dataset+"_mbhbh_HM_injected.pkl", "rb"))
    # for i, k in enumerate(["X", "Y", "Z"]):
    #     # td[k] = td_mbhb[k]
    #     # td[k] += wave[k]
    #     td_injected[k] += wave[k]


elif dataset == 'Spritz':
    names = fid["sky/cat"].dtype.names
    params = [np.array(fid["sky/cat"][k]).squeeze() for k in names]
    cat = np.rec.fromarrays(params, names=list(names))
    indexes = np.argsort(cat['Frequency'])
    cat = cat[indexes]
    # print(cat)
    td = fid["obs/tdi"][()]
    td = np.rec.fromarrays(list(td.T), names=["t", "X", "Y", "Z"])
    td = td['t']
    dt = td["t"][1]-td["t"][0]
    Tobs = float(td['t'][-1]/reduction)

freq_full = np.fft.rfftfreq(int(len(td['t'])), d=dt)

# subtract found MBHBs
# wave_gen = BBHWaveformFD(amp_phase_kwargs=dict(run_phenomd=not(HM)), use_gpu=gpu)

# modes = [(2,2)]
# if HM:
#     modes = [(2,2), (3,3), (4,4), (2,1), (3,2), (4,3)]
# waveform_kwargs = dict( modes=modes, direct=False, fill=True, squeeze=True, length=1024,  shift_t_limits=False) 

# td_found = deepcopy(td)
# for k in ["X", "Y", "Z"]:
#     td_found[k] -= td_found[k]
# if dataset == 'Sangria':
#     if reduction == 2:
#         folderpath = MBHBPATH + 'found_signals_Sangria_HM_false_6months/'
#     else:
#         folderpath = MBHBPATH + 'found_signals_Sangria_HM_false_12months/'
#     onlyfiles = [f for f in os.listdir(folderpath) if os.path.isfile(os.path.join(folderpath, f))]
#     found_sources_mp_even_unsorted = []
#     for s_index2 in range(len(onlyfiles)):
#         signal = pickle.load(open(folderpath+'/'+onlyfiles[s_index2], 'rb'))
#         wave = wave_gen(*signal, freqs=freq_full,
#                                 **waveform_kwargs)[0]
#         found_signal_bbhx_tdi_fd = TDI(dict([(k,FrequencySeries(wave[i], fs=freq_full)) for i,k in enumerate(["A", "E", "T"])]))
#         found_signal_bbhx_tdi_fd.AET2XYZ()
#         found_signal_bbhx_tdi_td = TDI(dict([(k,found_signal_bbhx_tdi_fd[k].ts.ifft()) for k in ["X", "Y", "Z"]]))
#         for i, k in enumerate(["X", "Y", "Z"]):
#             td[k] -= found_signal_bbhx_tdi_td[k]
#             td_found[k] += found_signal_bbhx_tdi_td[k]
# pickle.dump(td_found, open(MBHBPATH+dataset+"_mbhbh_found.pkl", "wb"))


# plt.figure()
# # plt.plot(td_original['t'], td_original['X'])
# plt.plot(td['t'], td['X'])
# plt.plot(td_mbhb['t'], td_mbhb['X'])
# # plt.plot(td_injected['t'], td_injected['X'])
# plt.plot(td_original['t'], wave['X'])
# # plt.plot(td_original['t'], td_original['X']-td_mbhb['X']-td['X'])
# # plt.plot(td_original['t'], td_original['X']-td_injected['X']-td['X'])
# plt.plot(td_original['t'], td_mbhb['X']-wave['X'])
# # plt.plot(td_original['t'], td_injected['X']-wave['X'])
# # plt.plot(td_original['t'], td_mbhb['X']-td_injected['X'])
# plt.show()

# # Build timeseries and frequencyseries object for X,Y,Z
tdi_ts = dict([(k, TimeSeries(td[k][:int(len(td[k][:])/reduction)], dt=dt, t0=td.t[0])) for k in ["X", "Y", "Z"]])
# td_mbhbs_subtracted = dict([(k, TimeSeries(td_mbhbs_subtracted[k][:int(len(td[k][:])/reduction)], dt=dt, t0=td.t[0])) for k in ["X", "Y", "Z"]])
# # tdi_ts_o = deepcopy(tdi_ts)
# tdi_fs_o = xr.Dataset(dict([(k, tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
# tdi_fs_mbhbs_subtracted = xr.Dataset(dict([(k, td_mbhbs_subtracted[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))


# plt.figure()
# plt.plot(tdi_fs_o['f'], np.abs(tdi_fs_o['X'].values))
# plt.plot(tdi_fs_mbhbs_subtracted['f'], np.abs(tdi_fs_mbhbs_subtracted['X'].values))
# plt.show()


GB = fastGB.FastGB(delta_t=dt, T=Tobs)  # in seconds


if add_gaps:
    DATAPATH_spritz = grandparent+"/LDC/Spritz/data"
    data_fn_spritz = DATAPATH_spritz + "/LDC2_spritz_vgb_training_v2.h5"
    fid_spritz = h5py.File(data_fn_spritz)

    td_obs = fid_spritz["obs/tdi"][()]
    td_obs = np.rec.fromarrays(list(td_obs.T), names=["t", "X", "Y", "Z"])
    td_obs = td_obs['t']
    td_clean = fid_spritz["clean/tdi"][()]
    td_clean = np.rec.fromarrays(list(td_clean.T), names=["t", "X", "Y", "Z"])
    td_clean = td_clean['t']
    td_galaxy = fid_spritz["gal/tdi"][()]
    td_galaxy = np.rec.fromarrays(list(td_galaxy.T), names=["t", "X", "Y", "Z"])
    td_galaxy = td_galaxy['t']

    tdi = fid_spritz["obs/tdi"][()].squeeze()
    tdi_nf = fid_spritz["noisefree/tdi"][()].squeeze()
    dt_spritz = tdi['t'][1]-tdi['t'][0]

    # Build timeseries and frequencyseries object for X,Y,Z
    tdi_ts = dict([(k, TimeSeries(td[k][:int(len(td[k][:])/reduction)], dt=dt, t0=td.t[0])) for k in ["X", "Y", "Z"]])
    tdi_ts_obs = dict([(k, TimeSeries(td_obs[k][:int(len(td_obs[k][:])/reduction)], dt=dt_spritz, t0=td_obs.t[0])) for k in ["X", "Y", "Z"]])
    tdi_ts_clean = dict([(k, TimeSeries(td_clean[k][:int(len(td_clean[k][:])/reduction)], dt=dt_spritz, t0=td_clean.t[0])) for k in ["X", "Y", "Z"]])
    tdi_ts_galaxy = dict([(k, TimeSeries(td_galaxy[k][:int(len(td_galaxy[k][:])/reduction)], dt=dt_spritz, t0=td_galaxy.t[0])) for k in ["X", "Y", "Z"]])

    tdi_ts_glitches = deepcopy(tdi_ts_obs)
    for k in ["X", "Y", "Z"]:
        tdi_ts_glitches[k].values = tdi_ts_obs[k].values - tdi_ts_clean[k].values - tdi_ts_galaxy[k].values

    ## add gaps to tdi
    tdi_ts_with_glitches = deepcopy(tdi_ts)
    if dataset != 'Spritz':
        for k in ["X", "Y", "Z"]:
            tdi_ts_with_glitches[k].values = tdi_ts[k].values + tdi_ts_glitches[k].values - tdi_ts_glitches[k].values
    if dataset == 'Spritz':
        for k in ["X", "Y", "Z"]:
            # tdi_ts_with_glitches[k].values = tdi_ts_clean[k].values + tdi_ts_glitches[k].values - tdi_ts_glitches[k].values + tdi_ts_galaxy[k].values
            tdi_ts_with_glitches[k].values = tdi_ts_clean[k].values + tdi_ts_galaxy[k].values ### no gaps
            # tdi_ts_with_glitches[k].values = tdi_ts[k].values - tdi_ts_glitches[k].values


    tdi_ts = deepcopy(tdi_ts_with_glitches)

if zero_gaps:
    gaps = {}
    for k in ["X", "Y", "Z"]:
        gap = np.isnan(tdi_ts_with_glitches[k])
        tdi_ts_with_glitches[k][gap] = 0
        gaps[k] = tdi_ts_with_glitches[k] == 0
        # gaps = np.isnan(tdi_ts_with_glitches[k])
        tdi_ts_with_glitches[k][gaps[k]] = 0

if zero_glitches:
    for k in ["X", "Y", "Z"]:
        mad = scipy.stats.median_abs_deviation(tdi_ts_with_glitches[k])
        peaks, properties = scipy.signal.find_peaks(np.abs(tdi_ts_with_glitches[k]), height=10*mad, threshold=None, distance=1)
        # Turning glitches into gaps
        for pk in peaks:
            tdi_ts_with_glitches[k][pk-10:pk+10] = 0.0

    for k in ["X", "Y", "Z"]:
        tdi_ts_with_glitches[k][:300] = 0



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
    tdi_ts = deepcopy(tdi_ts_with_glitches)


tdi_fs = xr.Dataset(dict([(k, tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
# tdi_fs = TDI(dict([(k, tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
# pickle.dump(tdi_fs, open(MBHBPATH+dataset+"_tdi_fs_residual.pkl", "wb"))

# noise_model = "SciRDv1"
# Nmodel = get_noise_model(noise_model, np.logspace(-5, -1, 100))

pGB = {}
ind = 0
found_sources = []
target_sources = []
first_start = time.time()
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
            search1 = Search(tdi_fs_search,self.Tobs, lower_frequency, upper_frequency, tdi2=tdi2)
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
        SNR_threshold = 10
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
        if lower_frequency > 5*10**-3:
            self.signals_per_window = 3
        if lower_frequency > 10**-2:
            self.signals_per_window = 1
        # current_loglikelihood_ratio = 1000
        ind = 0
        while current_SNR > SNR_threshold and ind < self.signals_per_window:
        # while current_loglikelihood_ratio > loglikelihood_ratio_threshold and ind < self.signals_per_window:
            ind += 1
            search1 = Search(tdi_fs_search,self.Tobs, lower_frequency, upper_frequency, tdi2=tdi2, recombination=1)
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
            if lower_frequency > 10**-2:
                search_repetitions = 6
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
                    for i in range(10):
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
                    Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_out[i], oversample=4)
                    source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
                    index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
                    index_high = index_low+len(Xs_subtracted)
                    for k in ["X", "Y", "Z"]:
                        tdi_fs_subtracted[k].data[index_low:index_high] = tdi_fs_subtracted[k].data[index_low:index_high] - source_subtracted[k].data

                search_out_subtracted = Search(tdi_fs_subtracted,self.Tobs, lower_frequency, upper_frequency, tdi2=tdi2)
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
                Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources[i], oversample=4)
                source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
                index_low = np.searchsorted(tdi_fs_search["X"].f, Xs_subtracted.f[0])
                index_high = index_low+len(Xs_subtracted)
                for k in ["X", "Y", "Z"]:
                    tdi_fs_search[k][index_low:index_high] = tdi_fs_search[k][index_low:index_high] - source_subtracted[k]
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
            search_out_subtracted = Search(tdi_fs_subtracted,self.Tobs, lower_frequency, upper_frequency, tdi2=tdi2)

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

def tdi_subtraction(tdi_fs,found_sources_mp_subtract, frequencies_search=None):

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
            Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_mp_subtract[i], oversample=4)
            source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
            index_low = np.searchsorted(tdi_fs_subtracted2["X"].f, Xs_subtracted.f[0])
            index_high = index_low+len(Xs_subtracted)
            for k in ["X", "Y", "Z"]:
                tdi_fs_subtracted2[k].data[index_low:index_high] -= source_subtracted[k].data
    return tdi_fs_subtracted2


#sum the found sources
# if mbhbs_removed:
#     found_sources = np.load(SAVEPATH+'found_sources_not_anticorrelated_Sangria_6m_no_mbhb_seed'+str(seed)+'.pkl', allow_pickle = True)
# else:
#     found_sources = np.load(SAVEPATH+'found_sources_not_anticorrelated_Sangria_6m_mbhb_seed'+str(seed)+'.pkl', allow_pickle = True)
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
# # if mbhbs_removed:
# #     pickle.dump(tdi_fs_sum_found, open(SAVEPATH+'tdi_fs_sum_found_6m_no_mbhb_seed'+str(seed)+'.pkl', "wb"))
# # else:
# #     pickle.dump(tdi_fs_sum_found, open(SAVEPATH+'tdi_fs_sum_found_6m_mbhb_seed'+str(seed)+'.pkl', "wb"))

# tdi_fs_subtracted = deepcopy(tdi_fs)
# for k in ["X", "Y", "Z"]:
#     tdi_fs_subtracted[k].data -= tdi_fs_sum_found[k].data

# plt.figure()
# plt.semilogx(tdi_fs['X'].f, (tdi_fs['X'].data))
# plt.semilogx(tdi_fs_sum_found['X'].f, (tdi_fs_sum_found['X'].data))
# plt.semilogx(tdi_fs_subtracted['X'].f, (tdi_fs_subtracted['X'].data))
# plt.show()

# plt.figure()
# plt.loglog(tdi_fs['X'].f, np.abs(tdi_fs['X'].data))
# plt.loglog(tdi_fs_sum_found['X'].f, np.abs(tdi_fs_sum_found['X'].data))
# plt.loglog(tdi_fs_subtracted['X'].f, np.abs(tdi_fs_subtracted['X'].data))
# plt.show()

# pickle.dump(tdi_fs_subtracted, open(MBHBPATH+dataset+"_tdi_fs_residual.pkl", "wb"))

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
# search_range = [0.0003, 0.03]
if dataset == 'Radler':
    search_range = [0.0003, 0.0319]
# search_range = [0.0001, 0.11]

# search_range = [0.0019935, 0.0020135]
# search_range = [0.0029935, 0.0030135]
# window_length = 1*10**-7 # Hz

frequencies = create_frequency_windows(search_range, Tobs)

frequencies_even = frequencies[::2]
frequencies_odd = frequencies[1::2]

# counts = np.zeros(len(pGB_injected))
# for i in range(len(pGB_injected)):
#     counts[i] = len(pGB_injected[i])

frequencies_search = np.asarray(frequencies)
# figure = plt.figure()
# plt.loglog(frequencies_search[:,1],counts, '.')
# plt.xlabel('Frequency [Hz]')
# plt.ylabel('Number of signals')
# plt.show()

# figure = plt.figure()
# plt.loglog(frequencies_search[:,0],(frequencies_search[:,1]-frequencies_search[:,0]),  linewidth= 4, label= '$B_{segment}$')
# plt.loglog(frequencies_search[:,0],(frequencies_search[:,1]-frequencies_search[:,0])/2,  linewidth= 4, label= '$B_{max}$')
# plt.loglog(frequencies_search[:,0],frequency_derivative(frequencies_search[:,0],M_chirp_upper_boundary)*Tobs,  linewidth= 4, label= '$B_{max}$')
# # plt.loglog(frequencies_search[:,0],frequency_derivative(frequencies_search[:,0],M_chirp_upper_boundary)*Tobs+frequencies_search[:,0]*2* 10**-4+4/31536000*2, label= '$B_{F}$')
# plt.loglog(frequencies_search[:,0],frequencies_search[:,0]*2* 10**-4, label= '$2 \cdot B_{O}$')
# plt.loglog(frequencies_search[:,0],np.ones(len(frequencies_search[:,0]))*4/31536000*2, label= '$2 \cdot B_{C}$')
# plt.xlabel(r'$f$ (Hz)')
# plt.ylabel(r'window width (Hz)')
# plt.xlim(10**-4,0.1)
# plt.ylim(bottom=(frequencies_search[0,1]-frequencies_search[0,0])/10**1)
# plt.legend()
# plt.show()
# plt.savefig(SAVEPATH+'bandwidth.png')


save_name = 'Sangria_12m_mbhb_even'
# save_name = 'Radler_24m_filled_anticorrelation'
# save_name = 'Spritz'
# save_name = dataset + '_12m_eventest'
# save_name = dataset + '_VGB_gaps'
# for i in range(65):
frequencies_search = frequencies_even
frequencies_search_full = deepcopy(frequencies_search)
batch_index = int(sys.argv[1]) 
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.008545)
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


# frequencies_search = frequencies_even

if VGB:
    frequencies_search = []
    for i in range(len(cat)):
        start_index = np.searchsorted(np.asarray(frequencies)[:,0], cat[i]['Frequency'])-1
        frequencies_search.append(frequencies[start_index])
# frequencies_search = frequencies_search[30:]
### redo the search for the anticorrelated signals
# frequencies_search = []
# # found_sources_anticorrelated_flat = pickle.load(open(SAVEPATH+'found_sources_anticorrelated_Sangria_12m.pkl', 'rb'))
# found_sources_anticorrelated_flat = pickle.load(open(SAVEPATH+'found_sources_anticorrelated_Radler_24m.pkl', 'rb'))
# for i in range(int(len(found_sources_anticorrelated_flat)/2)):
#     start_index = np.searchsorted(np.asarray(frequencies_odd)[:,0], found_sources_anticorrelated_flat[i*2]['Frequency'])-1
#     frequencies_search.append(frequencies_odd[start_index])
# frequencies_search = [frequencies_search[1]]

# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.003977)

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


subtract_all = False
do_subtract = True
if do_subtract:
    start = time.time()
    # save_name_previous = 'found_sourcesRadler_half_odd_dynamic_noise'
    # Sangria
    # save_name_previous = 'found_sources_Sangria_6m_mbhb_odd_seed'+str(seed)+'_flat'
    if mbhbs_removed:
        save_name_previous = 'found_sources_Sangria_12m_no_mbhb_odd_seed'+str(seed)+'_flat'
    else:
        save_name_previous = 'found_sources_Sangria_12m_mbhb_odd_seed'+str(seed)+'_flat'
    # save_name_previous = 'found_sources_Sangria_6m_mbhb_flat'
    # save_name_previous = 'found_sources_not_anticorrelated_Sangria_12m'
    # save_name_previous = 'found_sources_Radler_24m_even'
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
    # found_sources_mp_subtract = pickle.load(open(SAVEPATH+'found_sources_not_anticorrelated_Sangria_12m.pkl', 'rb'))
    # found_sources_mp_subtract = pickle.load(open(SAVEPATH+'found_sources_not_anticorrelated_Radler_24m.pkl', 'rb'))

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
    if not(subtract_all):
        for i in range(len(frequencies_search_full)):
            found_sources_out_flat_df = found_sources_out_flat_df[(found_sources_out_flat_df['Frequency']< frequencies_search_full[i][0]) | (found_sources_out_flat_df['Frequency']> frequencies_search_full[i][1])]
        found_sources_out_flat_df = found_sources_out_flat_df.sort_values('Frequency')
    found_sources_out_flat = found_sources_out_flat_df.to_dict(orient='records')
    tdi_fs_subtracted = tdi_subtraction(tdi_fs,found_sources_out_flat, frequencies_search_full)

    print('subtraction time', time.time()-start)
    plot_subtraction = False
    if plot_subtraction:
        # i = start_index
        i = 4
        # lower_frequency = frequencies_search_full[i][0]
        # upper_frequency = frequencies_search_full[i][1]
        lower_frequency = frequencies_search[i][0]
        upper_frequency = frequencies_search[i][1]
        # lower_frequency = frequencies_search[i][1]
        # upper_frequency = frequencies_search[i+1][0]
        # lower_frequency = frequencies_search_full[i][1]
        # upper_frequency = frequencies_search_full[i+1][0]
        found_sources_neighbor = found_sources_flat[(found_sources_flat_df['Frequency']> frequencies_search[i-1][0]) & (found_sources_flat_df['Frequency']< frequencies_search[i+1][1])]
        search1 = Search(tdi_fs,Tobs, lower_frequency, upper_frequency, tdi2=tdi2)
        search1.plot(second_data= tdi_fs_subtracted, found_sources_in=found_sources_neighbor)
        # search1.plot()
        
    tdi_fs = deepcopy(tdi_fs_subtracted)

do_not_search_unchanged_even_windows = True
if do_not_search_unchanged_even_windows:
    frequencies_search_reduced = []

    # save_name_previous = 'found_sources_Sangria_12m_even3'
    if mbhbs_removed:
        save_name_previous = 'found_sources_Sangria_12m_no_mbhb_even3_seed'+str(seed)+'_flat'
    else:
        save_name_previous = 'found_sources_Sangria_12m_mbhb_even3_seed'+str(seed)+'_flat'
    found_sources_mp_previous = np.load(SAVEPATH+save_name_previous+'.pkl', allow_pickle = True)
    # found_sources_flat = []
    # for j in range(len(found_sources_mp_previous)):
    #     for k in range(len(found_sources_mp_previous[j][3])):
    #         found_sources_flat.append(found_sources_mp_previous[j][3][k])
    # found_sources_flat = np.asarray(found_sources_flat)
    found_sources_flat = np.asarray(found_sources_mp_previous)
    found_sources_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_flat]) for attribute in found_sources_flat[0].keys()}
    found_sources_flat_df = pd.DataFrame(found_sources_flat_array)
    found_sources_in_skipped = []
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
        else:
            found_sources_in_flat_df = found_sources_flat_df[(found_sources_flat_df['Frequency']> frequencies_search[i][0]) & (found_sources_flat_df['Frequency']< frequencies_search[i][1])]
            if len(found_sources_in_flat_df) > 2:
                frequencies_search_reduced.append(frequencies_search[i])
            else:
                # print('no search', frequencies_search[i])
                found_sources_in_skipped.append(found_sources_in_flat_df.to_dict(orient='records'))
    found_sources_in_skipped = np.concatenate(found_sources_in_skipped)
    # pickle.dump(found_sources_in_skipped, open(SAVEPATH+save_name_previous+'_skipped.pkl', "wb"))
    frequencies_search = frequencies_search_reduced

found_sources_sorted = []
use_initial_guess = True
if use_initial_guess:
    # save_name_found_sources_previous = 'found_sources397769to400619LDC1-4_4mHz_half_year_even10'
    # save_name_found_sources_previous = 'found_sources397919to400770LDC1-4_4mHz_half_year_odd'
    # save_name_found_sources_previous = 'found_sources2537595to3305084LDC1-4_4mHz_half_year_even'
    # save_name_found_sources_previous = 'found_sourcesLDC1-4_half_even10'
    # save_name_found_sources_previous = 'found_sources_Radler_12m'
    # save_name_found_sources_previous = 'found_sources_not_anticorrelated_Sangria_6m_mbhb'
    if mbhbs_removed:
        save_name_found_sources_previous = 'found_sources_not_anticorrelated_Sangria_6m_no_mbhb_seed'+str(seed)
    else:
        save_name_found_sources_previous = 'found_sources_not_anticorrelated_Sangria_6m_mbhb_seed'+str(seed)
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

# index =35
# pGB_injected = {}
# for parameter in parameters:
#     pGB_injected[parameter] = cat[index][parameter]
# search1 = Search(tdi_fs,Tobs, frequencies_search[index][0], frequencies_search[index][1], tdi2=tdi2)
# search1.plot(pGB_injected=[pGB_injected])
# print(search1.SNR([pGB_injected]))
# print(search1.SNR([pGB_injected[index][0]]))
# search1.update_noise()
# print(search1.SNR([found_sources_mp_loaded[index][0][2]]))
# search1.plot(found_sources_in=found_sources_mp_loaded[index][0], pGB_injected=pGB_injected[index])
# search1.SNR(pGB_injected[index])
# frequencies_search = [frequencies_search[35]]
do_search = True
if do_search:
    MLP = MLP_search(tdi_fs, Tobs, signals_per_window = 10, found_sources_previous = found_sources_sorted, strategy = 'DE')
    start = time.time()

    # cpu_cores = 16
    # pool = mp.Pool(cpu_cores)
    # # pool.starmap(MLP.search, frequencies_search)
    # found_sources_mp = pool.starmap(MLP.search, frequencies_search)
    # pool.close()
    # pool.join()

    found_sources_mp = []
    for i in range(len(frequencies_search)):
        found_sources_mp.append(MLP.search(*frequencies_search[i]))

    print('time to search ', len(frequencies_search), 'windows: ', time.time()-start)
    if mbhbs_removed:
        directory = SAVEPATH+'found_signals_12m_no_mbhb_seed'+str(seed)
    else:
        directory = SAVEPATH+'found_signals_12m_mbhb_seed'+str(seed)
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
    found_sources_mp_loaded = pickle.load(open(SAVEPATH+'found_signals/found_sources_batch_'+batch_index+'_'+str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.pkl', 'rb'))
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
