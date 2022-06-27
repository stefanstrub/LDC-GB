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
sys.path.append('/cluster/home/sstrub/Repositories/LDC/lib/lib64/python3.8/site-packages/ldc-0.1-py3.8-linux-x86_64.egg')

from ldc.lisa.noise import get_noise_model
from ldc.common.series import TimeSeries, window
import ldc.waveform.fastGB as fastGB
from ldc.common.tools import compute_tdi_snr
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
SAVEPATH = grandparent+"/LDC/pictures"

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
reduction = 4
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
np.random.seed(42) #40
number_of_signals = 1
signals_per_subtraction = 1

f = 0.0115248
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

SNR_threshold = 10

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

        initial_guess = []
        if len(self.found_sources_previous) > 0:
            if do_subtract:
                padding_of_initial_guess_range = 0
            else:
                padding_of_initial_guess_range = (upper_frequency - lower_frequency)/2
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
            
            search1 = Search(tdi_fs_search,self.Tobs, lower_frequency, upper_frequency)
            ### sort the initial guesses such that the highest loglikelihhod guess comes first
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
        ind = 0
        while current_SNR > SNR_threshold and ind < self.signals_per_window:
            ind += 1
            
            search1 = Search(tdi_fs_search,self.Tobs, lower_frequency, upper_frequency)
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
            # print('SNRm', np.round(search1.SNRm([search1.pGB]),3))
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

                print('SNRm of found signal', np.round(search1.SNRm(maxpGBsearch_new[0]),3))
                print('which signal per window', ind,'and repetition:', i)
                found_sources_all.append(maxpGBsearch_new)
                number_of_evaluations_all.append(number_of_evaluations)
                new_SNR = search1.SNR(maxpGBsearch_new[0])
                if i == 0:
                    current_SNR = deepcopy(new_SNR)
                if new_SNR >= current_SNR:
                    current_SNR = deepcopy(new_SNR)
                if current_SNR < SNR_threshold:
                    break
                
                try:
                    for j in range(len(maxpGBsearch_new[0])):
                        A_optimized = search1.calculate_Amplitude([maxpGBsearch_new[0][j]])
                        maxpGBsearch_new[0][j]['Amplitude'] *= A_optimized.values
                    print('in range', maxpGBsearch_new[0][0]['Frequency'] > lower_frequency and maxpGBsearch_new[0][0]['Frequency'] < upper_frequency)
                    # new_SNR = search1.SNR(maxpGBsearch_new[0])
                    if i == 0:
                        maxpGBsearch = deepcopy(maxpGBsearch_new)
                    if new_SNR >= current_SNR:
                        maxpGBsearch = deepcopy(maxpGBsearch_new)
                    found_sources_all[-1] = maxpGBsearch_new
                except:
                    pass
                print('current SNR', current_SNR)

            if current_SNR < SNR_threshold:
                break
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
            if len(found_sources_in) > 1:
                tdi_fs_subtracted = deepcopy(self.tdi_fs)
                for i in range(len(found_sources_out)):
                    Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_out[i], oversample=4, simulator="synthlisa")
                    source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
                    index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
                    index_high = index_low+len(Xs_subtracted)
                    for k in ["X", "Y", "Z"]:
                        tdi_fs_subtracted[k].data[index_low:index_high] = tdi_fs_subtracted[k].data[index_low:index_high] - source_subtracted[k].data

                search_out_subtracted = Search(tdi_fs_subtracted,self.Tobs, lower_frequency, upper_frequency)

                total_boundaries = deepcopy(search1.boundaries)
                amplitudes = []
                for i in range(len(found_sources_in)):
                    amplitudes.append(found_sources_in[i]['Amplitude'])
                total_boundaries['Amplitude'] = [np.min(amplitudes),np.max(amplitudes)]
                amplitudes_length = np.log10(total_boundaries['Amplitude'][1]) - np.log10(total_boundaries['Amplitude'][0])
                total_boundaries['Amplitude'] = [np.log10(total_boundaries['Amplitude'][0]) - amplitudes_length/5, np.log10(total_boundaries['Amplitude'][1]) + amplitudes_length/5]

                start = time.time()
                found_sources_in = search_out_subtracted.optimize([found_sources_in], boundaries= total_boundaries)
                print(time.time()-start)

            found_sources = found_sources_in + found_sources_out

            #subtract the found sources from original
            tdi_fs_search = deepcopy(self.tdi_fs)
            for i in range(len(found_sources)):
                Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources[i], oversample=4, simulator="synthlisa")
                source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
                index_low = np.searchsorted(tdi_fs_search["X"].f, Xs_subtracted.f[0])
                index_high = index_low+len(Xs_subtracted)
                for k in ["X", "Y", "Z"]:
                    tdi_fs_search[k].data[index_low:index_high] = tdi_fs_search[k].data[index_low:index_high] - source_subtracted[k].data
        return found_sources, found_sources_all, number_of_evaluations_all, found_sources_in, [lower_frequency, upper_frequency]

def tdi_subtraction(tdi_fs,found_sources_mp_subtract, frequencies_search):

    found_sources_mp_best = []
    for i in range(len(found_sources_mp_subtract)):
        found_sources_mp_best.append(found_sources_mp_subtract[i][0])

    frequencies_search = np.asarray(frequencies_search)
    found_sources_to_subtract = []
    for i in range(len(found_sources_mp_best)):
        found_sources_to_subtract.append([])
        for j in range(len(found_sources_mp_best[i])):        
            # find closest frequency window
            frequency_window_index = np.searchsorted(frequencies_search[:,0], found_sources_mp_best[i][j]['Frequency'])-1
            if frequency_window_index < 0:
                found_sources_to_subtract[i].append(found_sources_mp_best[i][j])
            elif found_sources_mp_best[i][j]['Frequency'] > frequencies_search[frequency_window_index][1]:
                found_sources_to_subtract[i].append(found_sources_mp_best[i][j])

    #subtract the found sources from original
    tdi_fs_subtracted2 = deepcopy(tdi_fs)
    for i in range(len(found_sources_to_subtract)):
        for j in range(len(found_sources_to_subtract[i])):
            Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_to_subtract[i][j], oversample=4, simulator="synthlisa")
            source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
            index_low = np.searchsorted(tdi_fs_subtracted2["X"].f, Xs_subtracted.f[0])
            index_high = index_low+len(Xs_subtracted)
            for k in ["X", "Y", "Z"]:
                tdi_fs_subtracted2[k].data[index_low:index_high] -= source_subtracted[k].data
    return tdi_fs_subtracted2

padding = 0.5e-6

save_name = 'LDC1-4_4mHz_half_year_odd'
indexes = np.argsort(cat['Frequency'])
cat_sorted = cat[indexes]

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

# frequencies = frequencies[:32]
frequencies_even = frequencies[::2]
frequencies_odd = frequencies[1::2]

##### plot number of signals per frequency window
# frequencies_search = frequencies[::10]
# pGB_injected = []
# for j in range(len(frequencies_search)):
#     padding = (frequencies_search[j][1] - frequencies_search[j][0])/2 *0
#     index_low = np.searchsorted(cat_sorted['Frequency'], frequencies_search[j][0]-padding)
#     index_high = np.searchsorted(cat_sorted['Frequency'], frequencies_search[j][1]+padding)
#     if cat_sorted['Frequency'][index_high] < frequencies_search[j][1]:
#         index_high -= 1
#     indexesA = np.argsort(-cat_sorted[index_low:index_high]['Amplitude'])
#     pGB_injected_window = []
#     pGB_stacked = {}
#     for parameter in parameters:
#         pGB_stacked[parameter] = cat_sorted[parameter][index_low:index_high][indexesA]
#     for i in range(len(cat_sorted['Amplitude'][index_low:index_high])):
#         pGBs = {}
#         for parameter in parameters:
#             pGBs[parameter] = pGB_stacked[parameter][i]
#         pGB_injected_window.append(pGBs)
#     pGB_injected.append(pGB_injected_window)

# counts = np.zeros(len(pGB_injected))
# for i in range(len(pGB_injected)):
#     counts[i] = len(pGB_injected[i])

# frequencies_search = np.asarray(frequencies_search)
# figure = plt.figure()
# plt.loglog(frequencies_search[:,1],counts, '.')
# plt.xlabel('Frequency [Hz]')
# plt.ylabel('Number of signals')
# plt.show()
# figure = plt.figure()
# plt.loglog(frequencies_search[:,1],frequencies_search[:,1]-frequencies_search[:,0],  linewidth= 4, label= 'Frequency window width')
# plt.loglog(frequencies_search[:,1],np.ones(len(frequencies_search[:,1]))*4*32*10**-9*2, label= 'LISA rotation')
# plt.loglog(frequencies_search[:,1],frequencies_search[:,1]*3* 10**-4, label= 'Doppler modulation')
# plt.loglog(frequencies_search[:,1],frequency_derivative(frequencies_search[:,1],2)*Tobs, label= '$\dot{f}_{max} \cdot T_{obs}$')
# plt.xlabel('Frequency [Hz]')
# plt.ylabel('Frequency window witdh [Hz]')
# plt.ylim(bottom=(frequencies_search[0,1]-frequencies_search[0,0])/10**1)
# plt.legend()
# plt.show()

frequencies_search = frequencies_odd
# batch_index = int(sys.argv[1])
batch_index = 65
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.003977)
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], cat_sorted[-2]['Frequency'])-1
# start_index = np.searchsorted(np.asarray(frequencies_search)[:,0], 0.0004)-1
batch_size = 64
start_index = batch_size*batch_index
print('batch',batch_index, start_index)
frequencies_search = frequencies_search[start_index:start_index+10]
### highest + padding has to be less than f Nyqist
while frequencies_search[-1][1] + (frequencies_search[-1][1] - frequencies_search[-1][0])/2 > f_Nyquist:
    frequencies_search = frequencies_search[:-1]
# frequencies_search = frequencies_search[70:80]
# frequencies_search = frequencies_search[25:]

# target_frequencies = []
# index_low = np.searchsorted(cat_sorted['Frequency'], search_range[0])
# for i in range(10):
#     target_frequencies.append(cat_sorted[-10+i-1]['Frequency'])
#     # target_frequencies.append(cat_sorted[index_low+i*100]['Frequency'])
# # target_frequencies = cat_sorted[-17:-1]['Frequency']
# frequencies_search = []
# for i in range(len(target_frequencies)):
#     current_frequency = target_frequencies[i]
#     f_smear = current_frequency *3* 10**-4
#     f_deviation = frequency_derivative(current_frequency,M_chirp_upper_boundary)*Tobs
#     print(current_frequency,frequency_derivative(current_frequency,M_chirp_upper_boundary))
#     window_length = f_smear + f_deviation
#     window_length += 4*32*10**-9*2
#     window_shift = ((np.random.random(1)-0.5)*window_length*0.5)[0]
#     frequencies_search.append([target_frequencies[i]-window_length/2+window_shift,target_frequencies[i]+window_length/2+window_shift])

search_range = [frequencies_search[0][0],frequencies_search[-1][1]]
# search_range = [1619472*10**-8,2689639*10**-8]
print('search range '+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))))

do_subtract = False
if do_subtract:
    # save_name_previous = 'found_sources397769to400619LDC1-4_4mHz_half_year_even3'
    # save_name_previous = 'found_sources397919to400770LDC1-4_4mHz_half_year_odd'
    save_name_previous = 'found_sources397956to401074LDC1-4_4mHz_2_year_initial_half_even3'
    # save_name_previous = 'found_sources397793to400909LDC1-4_4mHz_2_year_initial_half_odd_SNR10'
    # save_name_previous = 'LDC1-4 odd'
    save_name_previous = 'found_sourcesLDC1-4_half_even3'
    found_sources_mp_subtract = np.load(SAVEPATH+'/'+save_name_previous+'.npy', allow_pickle = True)
    tdi_fs_subtracted = tdi_subtraction(tdi_fs,found_sources_mp_subtract, frequencies_search, GB)
    tdi_fs = deepcopy(tdi_fs_subtracted)
# length = 100
# for i in range(length):
#     print(i,found_sources_mp[-i][0])
#     print(i,found_sources_mp_even[-i][0])
# found_sources_mp = found_sources_mp_even[-100:]
# frequencies_search = frequencies_even[-100:]

found_sources_sorted = []
use_initial_guess = False
if use_initial_guess:
    found_sources_loaded = []
    # save_name_found_sources_previous = 'found_sources397769to400619LDC1-4_4mHz_half_year_even10'
    save_name_found_sources_previous = 'found_sources397919to400770LDC1-4_4mHz_half_year_odd'
    found_sources_loaded.append(np.load(SAVEPATH+'/'+save_name_found_sources_previous+'.npy', allow_pickle = True))

    found_sources_previous = []
    for i in range(len(found_sources_loaded)):
        for j in range(len(found_sources_loaded[i])):
            for k in range(len(found_sources_loaded[i][j][0])):
                found_sources_previous.append(found_sources_loaded[i][j][0][k])

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

# search1 = Search(tdi_fs,Tobs, frequencies_search[0][0], frequencies_search[0][1])
# pGBadded['Frequency'] = frequencies_search[0][0]
# start = time.time()
# for i in range(10**3):
#     search1.loglikelihood_SNR([pGBadded])
# print('full time', time.time()-start)
# start = time.time()
# for i in range(10**3):
#     Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBadded, oversample=4, simulator="synthlisa")
# print('full time', time.time()-start)

# MLP = MLP_search(tdi_fs, Tobs, signals_per_window = 1, neighbor_subtracted = do_subtract, recombination=0.75, found_sources_previous = found_sources_sorted, strategy = 'DE')
# found_sources_mp = MLP.search(frequencies_search[4][0], frequencies_search[4][1], dt, noise_model, parameters, number_of_signals, GB, intrinsic_parameters)
# found_sources_mp = [found_sources_mp]
# frequencies_search = [frequencies_search[7]]
search_input = []
for i in range(len(frequencies_search)):
    search_input.append((frequencies_search[i][0],frequencies_search[i][1], tdi_fs, Tobs, 1, do_subtract, dt, noise_model, parameters, number_of_signals, GB, intrinsic_parameters))
do_search = True
if do_search:
    # found_sources_mp = MLP_start_search(frequencies_search[0][0],frequencies_search[0][1], tdi_fs, Tobs, 1, do_subtract, dt, noise_model, parameters, number_of_signals, GB, intrinsic_parameters)
    
    # MLP = MLP_search(tdi_fs, Tobs, signals_per_window = 1, neighbor_subtracted = do_subtract,recombination=0.75, found_sources_previous = found_sources_sorted, strategy = 'DE')
    start = time.time()
    pool = mp.Pool(mp.cpu_count())
    pool = mp.Pool(16)
    found_sources_mp = pool.starmap(MLP_start_search, search_input)
    pool.close()
    pool.join()
    print('time to search ', number_of_windows, 'windows: ', time.time()-start)
    np.save(SAVEPATH+'/found_sources'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', found_sources_mp)
    
do_print = False
if do_print:
    found_sources_mp = np.load(SAVEPATH+'/found_sources'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'.npy', allow_pickle = True)
    found_sources_mp_best = []
    found_sources_mp_all = []
    for i in range(len(found_sources_mp)):
        found_sources_mp_best.append(found_sources_mp[i][0])
        found_sources_in_window = []
        for j in range(len(found_sources_mp[i][1])):
            found_sources_in_window.append(found_sources_mp[i][1][j][0][0])
        found_sources_mp_all.append(found_sources_in_window)

    found_sources_in = []
    found_sources_out = []
    for i in range(len(found_sources_mp_best)):
        found_sources_in.append([])
        found_sources_out.append([])
        for j in range(len(found_sources_mp_best[i])):
            if found_sources_mp_best[i][j]['Frequency'] > frequencies_search[i][0] and found_sources_mp_best[i][j]['Frequency'] < frequencies_search[i][1]:
                found_sources_in[i].append(found_sources_mp_best[i][j])
            else:
                found_sources_out[i].append(found_sources_mp_best[i][j])

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
            pGB_injected[i][j]['SNR'] = float(search1.SNRm([pGB_injected[i][j]])[2].values)
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

    #check SNR
    for i in range(len(found_sources_in)):
        tdi_fs_subtracted_out = deepcopy(tdi_fs)
        if i != 1:
            continue
        for j in range(len( found_sources_out[i])):
            Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_mp_best[i][j], oversample=4, simulator="synthlisa")
            source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
            index_low = np.searchsorted(tdi_fs_subtracted_out["X"].f, Xs_subtracted.f[0])
            index_high = index_low+len(Xs_subtracted)
            for k in ["X", "Y", "Z"]:
                tdi_fs_subtracted_out[k].data[index_low:index_high] = tdi_fs_subtracted_out[k].data[index_low:index_high] - source_subtracted[k].data
        print('frequency range', frequencies_search[i][0],frequencies_search[i][1])
        search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i][1])
        tdi_fs_subtracted = deepcopy(tdi_fs_subtracted_out)
        for j in range(len( pGB_injected[i][:10])):
            #subtract the found sources from original
            for n in range(len( pGB_injected[i][:j])):
                Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=pGB_injected[i][n], oversample=4, simulator="synthlisa")
                source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
                index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
                index_high = index_low+len(Xs_subtracted)
                for k in ["X", "Y", "Z"]:
                    tdi_fs_subtracted[k].data[index_low:index_high] = tdi_fs_subtracted[k].data[index_low:index_high] - source_subtracted[k].data
            search_subtracted = Search(tdi_fs_subtracted,Tobs, frequencies_search[i][0], frequencies_search[i][1])
            print('true subtracted',search_subtracted.SNRm([pGB_injected[i][j]])[2].values, 'original data', search1.SNRm([pGB_injected[i][j]])[2].values)
        tdi_fs_subtracted = deepcopy(tdi_fs_subtracted)
        for j in range(len(found_sources_mp_best[i])):
            search_subtracted = Search(tdi_fs_subtracted,Tobs, frequencies_search[i][0], frequencies_search[i][1])
            print('found subtracted',search_subtracted.SNRm([found_sources_mp_best[i][j]])[2].values, 'original data', search1.SNRm([found_sources_mp_best[i][j]])[2].values, found_sources_mp_best[i][j]['Frequency'])
            #subtract the found sources from original
            Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_mp_best[i][j], oversample=4, simulator="synthlisa")
            source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
            index_low = np.searchsorted(tdi_fs_subtracted["X"].f, Xs_subtracted.f[0])
            index_high = index_low+len(Xs_subtracted)
            for k in ["X", "Y", "Z"]:
                tdi_fs_subtracted[k].data[index_low:index_high] = tdi_fs_subtracted[k].data[index_low:index_high] - source_subtracted[k].data


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
        search1 = Search(tdi_fs_subtracted,Tobs, lower_frequency, upper_frequency)
        if len(pGB_injected[i]) > 0:
            search1.plot(found_sources_in=found_sources_mp_best[i], pGB_injected=pGB_injected[i], saving_label =SAVEPATH+'/strain added'+ str(int(np.round(lower_frequency*10**8))) +save_name+'.png') 
            # search1.plot(found_sources_in=found_sources_in[i], pGB_injected=pGB_injected[i][:10], saving_label =SAVEPATH+'/strain added'+ str(int(np.round(lower_frequency*10**8))) +save_name+'in.png') 

#     #subtract the found sources from original
#     tdi_fs_subtracted = deepcopy(tdi_fs)
#     for i in range(len(found_sources_in)):
#         for j in range(len(found_sources_in[i])):
#             Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_in[i][j], oversample=4, simulator="synthlisa")
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
#             Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_in[i][j], oversample=4, simulator="synthlisa")
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



f_line = np.logspace(-4,-1, num=20)

indexes = np.argsort(cat['Frequency'])
pGB_injected = []
cat_sorted = cat[indexes]
frequencies_plot = []
for i in range(len(f_line)-1):
    frequencies_plot.append([f_line[i],f_line[i+1]])
cat_plot = []
for j in range(len(frequencies_plot)):
    index_low = np.searchsorted(cat_sorted['Frequency'], frequencies_plot[j][0])
    index_high = np.searchsorted(cat_sorted['Frequency'], frequencies_plot[j][1])

    try:
        cat_plot.append(cat_sorted[index_low:index_low+100])
    except:
        cat_plot.append(cat_sorted[index_low:])

def frequency_derivative_thyson(f):
    return 8*10**-7*f**(11/3)
def frequency_derivative_thyson_lower(f):
    return -5*10**-6*f**(13/3)

fig = plt.figure()
parameter_x = 'Frequency'
parameter_y = 'FrequencyDerivative'
# for i in range(len(cat_plot)):
#     for j in range(len(cat_plot[i])):
#         plt.scatter(cat_plot[i][parameter_x][j],cat_plot[i][parameter_y][j])
plt.plot(f_line, frequency_derivative(f_line,0.1))
plt.plot(f_line, frequency_derivative(f_line,M_chirp_upper_boundary))
plt.plot(f_line, frequency_derivative_thyson(f_line), color='red')
# plt.plot(f_line, frequency_derivative(f_line,100))
# plt.plot(f_line, frequency_derivative_Neil(f_line))
plt.hlines(0.01/Tobs**2, xmin=f_line[0], xmax=f_line[-1], linestyles='--')
plt.hlines(0.01/(Tobs/2)**2, xmin=f_line[0], xmax=f_line[-1], linestyles='--')
plt.hlines(0.01/(Tobs*2)**2, xmin=f_line[0], xmax=f_line[-1], linestyles='--')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$\log $f $[$Hz$]$')
plt.ylabel('$\log  \dot{f} [s^{-2}]$')
# plt.savefig(SAVEPATH+'/found_sources_'+save_name+'f-fd.png')


# %%
