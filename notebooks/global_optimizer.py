import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import time
from copy import deepcopy
import multiprocessing as mp
import pandas as pd
import os
import h5py
import pickle

from ldc.lisa.noise import get_noise_model
from ldc.common.series import TimeSeries, window
import ldc.waveform.fastGB as fastGB

from sources2 import *

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
            search_out_subtracted = Search(tdi_fs_subtracted,self.Tobs, lower_frequency, upper_frequency, dt=dt)

            total_boundaries = deepcopy(search_out_subtracted.boundaries)
            start = time.time()
            start_loglikelihood = search_out_subtracted.loglikelihood(found_sources_in)
            found_sources_in_new = search_out_subtracted.optimize([found_sources_in], boundaries= total_boundaries)
            optimized_loglikelihood = search_out_subtracted.loglikelihood(found_sources_in_new)
            if optimized_loglikelihood > start_loglikelihood:
                found_sources_in = found_sources_in_new
            print('global optimization time', np.round(time.time()-start), 'initial loglikelihood', np.round(start_loglikelihood,5), 'optimized_loglikelihood', np.round(optimized_loglikelihood,5), 'difference loglikelihood', np.round(optimized_loglikelihood-start_loglikelihood,5), 'frequency', lower_frequency )

            found_sources = found_sources_in + found_sources_out

        return found_sources_in

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

if Radler:
    sangria_fn = DATAPATH + "/LDC1-4_GB_v2.hdf5"
    # sangria_fn = DATAPATH + "/LDC1-3_VGB_v2_FD_noiseless.hdf5"
    # sangria_fn = DATAPATH + "/LDC1-3_VGB_v2.hdf5"
else:
    sangria_fn = DATAPATH + "/LDC2_sangria_training_v2.h5"
fid = h5py.File(sangria_fn)

reduction = 1

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
tdi_ts = dict([(k, TimeSeries(td[k][:int(len(td[k][:])/reduction)], dt=dt)) for k in ["X", "Y", "Z"]])
tdi_fs = xr.Dataset(dict([(k, tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
GB = fastGB.FastGB(delta_t=dt, T=Tobs)  # in seconds

noise_model = "SciRDv1"
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

chandrasekhar_limit = 1.4
M_chirp_upper_boundary = (chandrasekhar_limit**2)**(3/5)/(2*chandrasekhar_limit)**(1/5)


frequencies = []
frequencies_even = []
frequencies_odd = []
# search_range = [0.00398, 0.0041]
# search_range = [0.0039885, 0.0040205]
# search_range = [0.0039935, 0.0039965]
f_Nyquist = 1/dt/2
search_range = [0.0003, f_Nyquist]
# search_range = [0.0001, 0.11]
# search_range = [0.0019935, 0.0020135]
# search_range = [0.0029935, 0.0030135]
# window_length = 1*10**-7 # Hz

    
frequencies = create_frequency_windows(search_range, Tobs)

# frequencies_half_shifted = []
# for i in range(len(frequencies)-1):
#     frequencies_half_shifted.append([(frequencies[i][1]-frequencies[i][0])/2 +frequencies[i][0],(frequencies[i+1][1]-frequencies[i+1][0])/2 +frequencies[i+1][0]])
# frequencies = frequencies_half_shifted #### if shifted
frequencies_even = frequencies[::2]
frequencies_odd = frequencies[1::2]

frequencies_search = frequencies_even
frequencies_search_full = deepcopy(frequencies_search)


while frequencies_search[-1][1] + (frequencies_search[-1][1] - frequencies_search[-1][0])/2 > f_Nyquist:
    frequencies_search = frequencies_search[:-1]
search_range = [frequencies_search[0][0],frequencies_search[-1][1]]

save_name = 'Sangria_1year_dynamic_noise_opt_odd'

# found_sources_mp = np.load(SAVEPATH+'found_sources_'+save_name+'.npy', allow_pickle = True)
# found_sources_flat = []
# for i in range(len(found_sources_mp)):
#     for j in range(len(found_sources_mp[i])):
#         found_sources_flat.append(found_sources_mp[i][j])
# found_sources_flat = np.asarray(found_sources_flat)

found_sources_flat = np.load(SAVEPATH+'found_sources'+save_name+'_flat.npy', allow_pickle = True)
found_sources_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_flat]) for attribute in found_sources_flat[0].keys()}
found_sources_flat_df = pd.DataFrame(found_sources_flat_array)
found_sources_flat_df = found_sources_flat_df.sort_values('Frequency')
input = []
for i in range(len(frequencies_search)):
    input_dict = found_sources_flat_df[(found_sources_flat_df['Frequency'] > frequencies_search[i][0]) & (found_sources_flat_df['Frequency'] < frequencies_search[i][1])].to_dict(orient='records')
    input.append([frequencies_search[i][0],frequencies_search[i][1],input_dict])



def tdi_subtraction(tdi_fs,found_sources_mp_subtract, frequencies_search):

    #subtract the found sources from original
    tdi_fs_subtracted2 = deepcopy(tdi_fs)
    for i in range(len(found_sources_mp_subtract)):
        # for j in range(len(found_sources_to_subtract[i])):
            Xs_subtracted, Ys_subtracted, Zs_subtracted = GB.get_fd_tdixyz(template=found_sources_mp_subtract[i], oversample=4, simulator="synthlisa")
            source_subtracted = dict({"X": Xs_subtracted, "Y": Ys_subtracted, "Z": Zs_subtracted})
            index_low = np.searchsorted(tdi_fs_subtracted2["X"].f, Xs_subtracted.f[0])
            index_high = index_low+len(Xs_subtracted)
            for k in ["X", "Y", "Z"]:
                tdi_fs_subtracted2[k].data[index_low:index_high] -= source_subtracted[k].data
    return tdi_fs_subtracted2

do_subtract = True
if do_subtract:
    start = time.time()
    found_sources_out_flat = deepcopy(found_sources_flat_df)
    for i in range(len(frequencies_search_full)):
        found_sources_out_flat = found_sources_out_flat[(found_sources_flat_df['Frequency']< frequencies_search_full[i][0]) | (found_sources_flat_df['Frequency']> frequencies_search_full[i][1])]
    found_sources_out_flat = found_sources_out_flat.sort_values('Frequency')
    found_sources_out_flat = found_sources_out_flat.to_dict(orient='records')
    tdi_fs_subtracted = tdi_subtraction(tdi_fs,found_sources_out_flat, frequencies_search_full)

    print('subtraction time', time.time()-start)
    plot_subtraction = True
    if plot_subtraction:
        i = 3000
        lower_frequency = frequencies_search[i][0]
        upper_frequency = frequencies_search[i][1]
        search1 = Search(tdi_fs,Tobs, lower_frequency, upper_frequency)
        # source = [{'Amplitude': 4.500916389929765e-20, 'EclipticLatitude': 0.8528320149942861, 'EclipticLongitude': -0.9418744765040503, 'Frequency': frequencies_search[i][0]+(frequencies_search[i][1]-frequencies_search[i][0])/2, 'FrequencyDerivative': 2.1688352300259018e-22, 'Inclination': 1.343872907043714, 'InitialPhase': 3.583816929574315, 'Polarization': 2.69557290704741}]
        search1.plot(second_data= tdi_fs_subtracted)
        # search1.plot(found_sources_in=source)
        # search1.plot(second_data= tdi_fs_subtracted, found_sources_in=found_sources_mp_o[start_index][0])
        
    tdi_fs = deepcopy(tdi_fs_subtracted)


optimizer = Global_optimizer(tdi_fs, Tobs)

start = time.time()
cpu_cores = 10
pool = mp.Pool(cpu_cores)
found_sources_mp = pool.starmap(optimizer.optimize, input)
pool.close()
pool.join()
print('time to optimize', len(frequencies_search), 'windows: ', time.time()-start)

fn = SAVEPATH+'optimized/found_sources'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'_opt1_even.pkl'
pickle.dump(found_sources_mp, open(fn, "wb"))

# np.save(SAVEPATH+'optimized/found_sources'+ str(int(np.round(search_range[0]*10**8)))+'to'+ str(int(np.round(search_range[1]*10**8))) +save_name+'_opt1_odd.npy', found_sources_mp)
found_sources_mp = pickle.load(open(fn, 'rb'))

found_sources_new_flat = []
for i in range(len(found_sources_mp)):
    for j in range(len(found_sources_mp[i])):
        found_sources_new_flat.append(found_sources_mp[i][j])
found_sources_new_flat = np.asarray(found_sources_new_flat)
found_sources_new_flat_array = {attribute: np.asarray([x[attribute] for x in found_sources_new_flat]) for attribute in found_sources_new_flat[0].keys()}
found_sources_new_flat_df = pd.DataFrame(found_sources_new_flat_array)
found_sources_new_flat_df = found_sources_new_flat_df.sort_values('Frequency')

found_sources_out_flat = pd.DataFrame(found_sources_out_flat)
for i in range(len(found_sources_mp)):
    found_sources_out_flat = found_sources_out_flat[(found_sources_out_flat['Frequency']< frequencies_search[i][0]) | (found_sources_out_flat['Frequency']> frequencies_search[i][1])]
found_sources_combined_flat_df = found_sources_out_flat.append(found_sources_new_flat_df, ignore_index=True)
found_sources_combined_flat_df = found_sources_combined_flat_df.sort_values('Frequency')
found_sources_combined_flat_df = found_sources_combined_flat_df.to_dict(orient='records')

fn = SAVEPATH+'found_sources_'+ save_name+'_opt.npy'
np.save(fn,found_sources_combined_flat_df)