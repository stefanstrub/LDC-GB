import matplotlib.pyplot as plt
import scipy
import numpy as np
import xarray as xr
from astropy import units as u
import pandas as pd
import time
from copy import deepcopy

import ldc.io.hdf5 as hdfio
from ldc.lisa.noise import get_noise_model
from ldc.lisa.orbits import Orbits
from ldc.lisa.projection import ProjectedStrain
from ldc.common.series import TimeSeries, FrequencySeries, window
import ldc.waveform.fastGB as fastGB
from ldc.common.tools import compute_tdi_snr
from ldc.waveform.waveform import HpHc

DATAPATH = "/home/stefan/LDC/Sangria/data"
sangria_fn = DATAPATH+"/dgb-tdi.h5"
# sangria_fn = DATAPATH+"/LDC2_sangria_blind_v1.h5"
# sangria_fn = DATAPATH+"/LDC2_sangria_gdb-tdi_v1_v3U3MxS.h5"
# sangria_fn = DATAPATH+"/LDC2_sangria_idb-tdi_v1_DgtGV85.h5"
# sangria_fn = DATAPATH+"/LDC2_sangria_mbhb-tdi_v1_MN5aIPz.h5"
# sangria_fn = DATAPATH+"/LDC2_sangria_training_v1.h5"
# tdi_ts, tdi_descr = hdfio.load_array(sangria_fn, name="obs/tdi")
sangria_fn = DATAPATH+"/LDC2_sangria_vgb-tdi_v1_sgsEVXb.h5"
tdi_ts, tdi_descr = hdfio.load_array(sangria_fn)
sangria_fn_training = DATAPATH+"/LDC2_sangria_training_v1.h5"
dt = int(1/(tdi_descr["sampling_frequency"]))

# Build timeseries and frequencyseries object for X,Y,Z
# tdi_ts = xr.Dataset(dict([(k,TimeSeries(tdi_ts[k], dt=dt)) for k in ["X", "Y", "Z"]]))
tdi_ts = xr.Dataset(dict([(k,TimeSeries(tdi_ts[k][:,1], dt=dt)) for k in ["X", "Y", "Z"]]))
tdi_fs = xr.Dataset(dict([(k,tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))

# tdi_ts_training, tdi_descr_training = hdfio.load_array(sangria_fn_training, name="obs/tdi")
# tdi_ts_training = xr.Dataset(dict([(k,TimeSeries(tdi_ts_training[k], dt=dt)) for k in ["X", "Y", "Z"]]))
# tdi_fs_training = xr.Dataset(dict([(k,tdi_ts_training[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))


noise_model = "MRDv1"
Nmodel = get_noise_model(noise_model, np.logspace(-5, -1, 100))
Npsd = Nmodel.psd()
# plt.figure(figsize=(12,6))
# plt.subplot(131)
# f, psdX =  scipy.signal.welch(tdi_ts["X"], fs=1.0/dt, window='hanning', nperseg=256*256)
# plt.loglog(f, np.sqrt(psdX), label="TDI X")
# plt.loglog(Nmodel.freq, np.sqrt(Npsd), label=noise_model, alpha=2)
# plt.legend()
# plt.xlabel("freq [Hz]")
# plt.ylabel("PSD")
# plt.axis([1e-5, None, 4e-22, 2e-19])
# plt.subplot(132)
# f, psdX =  scipy.signal.welch(tdi_ts["Y"], fs=1.0/dt, window='hanning', nperseg=256*256)
# plt.loglog(f, np.sqrt(psdX), label="TDI Y")
# plt.loglog(Nmodel.freq, np.sqrt(Npsd), label=noise_model, alpha=2)
# plt.legend()
# plt.xlabel("freq [Hz]")
# plt.ylabel("PSD")
# plt.axis([1e-5, None, 4e-22, 2e-19])
# plt.subplot(133)
# f, psdX =  scipy.signal.welch(tdi_ts["Z"], fs=1.0/dt, window='hanning', nperseg=256*256)
# plt.loglog(f, np.sqrt(psdX), label="TDI Z")
# plt.loglog(Nmodel.freq, np.sqrt(Npsd), label=noise_model, alpha=2)
# plt.legend()
# plt.xlabel("freq [Hz]")
# plt.ylabel("PSD")
# plt.axis([1e-5, None, 4e-22, 2e-19])


vgb, units = hdfio.load_array(sangria_fn_training, name="sky/vgb/cat")
start = time.time()
GB = fastGB.FastGB(delta_t=dt, T=float(tdi_ts["X"].t[-1])) # in seconds
print(time.time()- start)
start = time.time()
pGB = dict(zip(vgb.dtype.names, vgb[8])) # we take the source #8
print(time.time()- start)
#modify pGB
# pGB['InitialPhase'] *= 1.01
# pGB['Amplitude'] *= 1.01
# pGB['Frequency'] *= 1.05
# pGB['FrequencyDerivative'] *= 1.1
# pGB['Polarization'] += 0.01
# pGB['EclipticLatitude'] += 0.03
# pGB['EclipticLongitude'] += 0.01
# pGB['Inclination'] += 0.01
start = time.time()
Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB, oversample=4, simulator='synthlisa')
print(time.time()- start)
fmin, fmax = float(Xs.f[0]) , float(Xs.f[-1]+Xs.attrs['df'])
source = dict({"X":Xs, "Y":Ys, "Z":Zs})
Xs_td, Ys_td, Zs_td = GB.get_td_tdixyz(template=pGB, simulator='synthlisa')


plt.figure(figsize=(12,3))
# plt.plot(Xs_td.t, Xs_td2, label="TDI X")
# plt.plot(tdi_ts_training['X'].t, tdi_ts_training['X'], label="data")
plt.plot(tdi_ts['X'].t/86400, tdi_ts['X'], label="Verification Binaries")
plt.plot(Xs_td.t/86400, Xs_td, label="Binary")
# plt.plot(Xp.t, amplitude_envelope, label="envelope")
# plt.plot(Xs_td.t[::100], amplitude_envelope2[::100], label="envelope 2")
plt.ylabel('X TDI strain')
plt.xlabel('time [days]')
# plt.ylim(-10**-20,10**-20)
plt.legend()
plt.show()

# plt.figure(figsize=(12,6))
# plt.subplot(121)
# plt.title("real part")
# plt.plot(tdi_fs["X"].f, tdi_fs["X"].real, label="TDI X")
# plt.plot(Xs.f, (tdi_fs["X"][Xs.kmin:Xs.kmin+len(Xs)]-Xs.values).real, label="TDI X - fast "+pGB["Name"])
# plt.axis([pGB["Frequency"]-6e-7, pGB["Frequency"]+6e-7, -3e-17, 5e-17])
# plt.legend(loc="lower right")
# plt.xlabel("freq [Hz]")
# plt.subplot(122)
# plt.title("imaginary part")
# plt.plot(tdi_fs["X"].f, tdi_fs["X"].imag, label="TDI X")
# plt.plot(Xs.f, (tdi_fs["X"].isel(f=slice(Xs.kmin, Xs.kmin+len(Xs)))-Xs.values).imag, label="TDI X - fast "+pGB["Name"])
# plt.axis([pGB["Frequency"]-6e-7, pGB["Frequency"]+6e-7, -3e-17, 5e-17])
# plt.legend(loc="lower left")
# plt.xlabel("freq [Hz]")

# vgb, units = hdfio.load_array(sangria_fn_training, name="sky/vgb/cat")
# GB = fastGB.FastGB(delta_t=dt, T=float(tdi_ts["X"].t[-1])) # in seconds
# noise_model = "MRDv1"
# f_noise = np.logspace(-5, -1, 100)
# Nmodel = get_noise_model(noise_model, f_noise)
# SNR2 = np.zeros((len(vgb), 2)) # snr square
# for j,s in enumerate(vgb):
#     pGB = dict(zip(vgb.dtype.names, s))
#     Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB, oversample=4, simulator='synthlisa')
#     fmin, fmax = float(Xs.f[0]) , float(Xs.f[-1]+Xs.attrs['df'])
#     source = dict({"X":Xs, "Y":Ys, "Z":Zs})
#     SNR2[j,1] = compute_tdi_snr(source, Nmodel, data=tdi_fs, fmin=fmin, fmax=fmax)["tot2"]
#     SNR2[j,0] = compute_tdi_snr(source, Nmodel)["tot2"] 

f_noise = np.logspace(-5, -1, 100)
Nmodel = get_noise_model(noise_model, f_noise)
freq = np.array(source["X"].sel(f=slice(fmin, fmax)).f)

Sn = Nmodel.psd(freq=freq, option='X')

# plt.figure(figsize=(12,6))
# plt.semilogx(freq, Sn)
# plt.semilogx(source['X'].f,source['X']/np.sqrt(Nmodel.psd(freq=freq, option='X')))
# plt.show()
######################################################
#%%
#find GB parameters

def likelihood(data, simulation, Sn):
    diff = data - simulation
    p = float(np.mean(np.abs(diff)**2 / Sn/ 10**4))
    return np.exp(-p / 2.0)
print(likelihood(tdi_fs["X"].isel(f=slice(Xs.kmin, Xs.kmin+len(Xs))),Xs.values,Sn))
# Number of histogram bins.
n_bin = 50
# Total number of proposed samples.
number_of_samples = 10 ** 4
number_of_parameters = 3

# Make the first random sample. ------------------------------------
pGBs = deepcopy(pGB)
pGBs['Amplitude'] *= 1.101
pGBs['EclipticLatitude'] = (np.random.random()-0.5) * np.pi 
pGBs['EclipticLongitude'] = np.random.random() * np.pi 
pGBs['InitialPhase'] = np.random.random() * np.pi 
pGBs['Frequency'] *= 1.0001
pGBs['FrequencyDerivative'] *= 1.01
pGBs['Polarization'] = np.random.random() * np.pi 
pGBs['Inclination'] = np.random.random()* np.pi 
Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs, oversample=4, simulator='synthlisa')

dataX = tdi_fs["X"].isel(f=slice(Xs.kmin, Xs.kmin+len(Xs)))
dataY = tdi_fs["Y"].isel(f=slice(Ys.kmin, Ys.kmin+len(Ys)))
dataZ = tdi_fs["Z"].isel(f=slice(Zs.kmin, Zs.kmin+len(Zs)))
# dataX_training = tdi_fs_training["X"].isel(f=slice(Xs.kmin, Xs.kmin+len(Xs)))
# dataY_training = tdi_fs_training["Y"].isel(f=slice(Ys.kmin, Ys.kmin+len(Ys)))
# dataZ_training = tdi_fs_training["Z"].isel(f=slice(Zs.kmin, Zs.kmin+len(Zs)))
fmin, fmax = float(Xs.f[0]) , float(Xs.f[-1]+Xs.attrs['df'])
freq = np.array(Xs.sel(f=slice(fmin, fmax)).f)
Sn = Nmodel.psd(freq=freq, option='X')

# Evaluate posterior for the first sample.
# p1 = likelihood(dataX, Xs.values, Sn)
diff = np.abs(dataX - Xs.values)**2 + np.abs(dataY - Ys.values)**2 + np.abs(dataZ - Zs.values)**2

p = float(np.sum(diff / Sn)*Xs.attrs['df'])
p1 = np.exp(-p / 2.0)
samples = xr.Dataset(dict([(name,xr.DataArray(np.zeros(number_of_samples), dims=('number_of_sample'), coords={"number_of_sample": range(number_of_samples)},
                         )) for name, titles in pGBs.items()]))

for name, titles in pGBs.items():
    if name != 'Name':
        samples[name][0] = pGBs[name]
samples['Likelihood'] = samples['Name']
samples = samples.drop(['Name'])
samples['Likelihood'][0] = p1

# Run Metropolis-Hastings sampler. ---------------------------------
plt.figure()
plt.subplot(231)
# plt.plot(dataX_training.f*1000,dataX_training.values, label='data')
plt.plot(dataX.f*1000,dataX.values, label='binary')
# plt.plot(Xs.f, Xs.values, label='start')
plt.subplot(232)
# plt.plot(dataY_training.f*1000,dataY_training.values, label='data')
plt.plot(dataY.f*1000,dataY.values, label='binary')
# plt.plot(Ys.f, Ys.values, label='start')
plt.subplot(233)
# plt.plot(dataZ_training.f*1000,dataZ_training.values, label='data')
plt.plot(dataZ.f*1000,dataZ.values, label='binary')
# plt.plot(Zs.f, Zs.values, label='start')
plt.subplot(234)
# plt.plot(dataX_training.f*1000,dataX_training.values.imag, label='data')
plt.plot(dataX.f*1000,dataX.values.imag, label='binary')
# plt.plot(Xs.f, Xs.values.imag, label='start')
plt.subplot(235)
# plt.plot(dataY_training.f*1000,dataY_training.values.imag, label='data')
plt.plot(dataY.f*1000,dataY.values.imag, label='binary')
# plt.plot(Ys.f, Ys.values.imag, label='start')
plt.subplot(236)
# plt.plot(dataZ_training.f*1000,dataZ_training.values.imag, label='data')
plt.plot(dataZ.f*1000,dataZ.values.imag, label='binary')
# plt.plot(Zs.f, Zs.values.imag, label='start')

print(p1)
start = time.time()
for i in range(1, number_of_samples):

    # Normal distributed proposal.
    std = np.array(np.eye(8)*[1*10**-22,1/10, 1/10,1*10**-7,10**-19,1/10,1/10,1/10])
    previous_samples =  [samples['Amplitude'][i-1], samples['EclipticLatitude'][i-1], samples['EclipticLongitude'][i-1], samples['Frequency'][0],samples['FrequencyDerivative'][i-1], samples['Inclination'][i-1], samples['InitialPhase'][i-1], samples['Polarization'][i-1]]
    pGBs['Amplitude'], pGBs['EclipticLatitude'], pGBs['EclipticLongitude'], pGBs['Frequency'],pGB['FrequencyDerivative'], pGBs['Inclination'], pGBs['InitialPhase'], pGBs['Polarization'] = np.random.multivariate_normal(previous_samples, std**2)
    if i > 1000 and samples['Frequency'][i-1] == samples['Frequency'][i-401]:
        pGBs['Frequency'] = samples['Frequency'][i-1]
    if pGBs['Frequency'] < 10**-4 or pGBs['Frequency'] > 10**-2:
        pGBs['Frequency'] = 5*10**-4
    if pGBs['Amplitude'] < 10**-23 or pGBs['Amplitude'] > 10**-5:
        pGBs['Amplitude'] = 5*10**-23
    if pGBs['InitialPhase'] < 0 or pGBs['InitialPhase'] > 2*np.pi:
        pGBs['InitialPhase'] = 0
    # print(pGBs['Amplitude'])
    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs, oversample=4, simulator='synthlisa')
    dataX = tdi_fs["X"].isel(f=slice(Xs.kmin, Xs.kmin+len(Xs)))
    dataY = tdi_fs["Y"].isel(f=slice(Ys.kmin, Ys.kmin+len(Ys)))
    dataZ = tdi_fs["Z"].isel(f=slice(Zs.kmin, Zs.kmin+len(Zs)))

    fmin, fmax = float(Xs.f[0]) , float(Xs.f[-1]+Xs.attrs['df'])
    freq = np.array(Xs.sel(f=slice(fmin, fmax)).f)
    Sn = Nmodel.psd(freq=freq, option='X')
    # plt.plot(dataX.f,dataX.values.real)
    # plt.plot(Xs.f, Xs.values.real)

    # Evaluate posterior.
    diff = np.abs(dataX - Xs.values)**2 + np.abs(dataY - Ys.values)**2 + np.abs(dataZ - Zs.values)**2
    p = float(np.sum(diff / Sn)*Xs.attrs['df'])
    p_test = np.exp(-p / 2.0)
    T_inv =  1
    # print(p_test/p1, pGBs['Frequency'])
    # Apply Metropolis rule.
    # print(p_test.values,p.values,(p_test / p) ** T_inv)
    if p1 == 0:
        p1 = p_test
        for name, titles in pGBs.items():
            if name != 'Name':
                samples[name][i] = pGBs[name]
        samples['Likelihood'][i] = p1
    elif (p_test / p1) ** T_inv > np.random.rand():  # L^i/L^j
        p1 = p_test
        for name, titles in pGBs.items():
            if name != 'Name':
                samples[name][i] = pGBs[name]
        samples['Likelihood'][i] = p1
    
    else:
        for name in samples.data_vars:
            samples[name][i] = samples[name][i-1]
        
print(time.time()- start)
# plot
n_max = np.argmax(samples['Likelihood'].values)
pGBmax = deepcopy(pGB)
for name, titles in pGBmax.items():
    if name != 'Name':
        pGBmax[name] = samples[name][n_max].values
Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBmax, oversample=4, simulator='synthlisa')


plt.subplot(231)
# plt.plot(Xs.f, Xs.values, label='optimized')
plt.xlabel('f [mHz]')
plt.ylabel('X-TDI real [1/Hz]')
plt.legend()
plt.subplot(232)
# plt.plot(Ys.f, Ys.values, label='optimized')
plt.xlabel('f [mHz]')
plt.ylabel('Y-TDI real [1/Hz]')
plt.legend()
plt.subplot(233)
# plt.plot(Zs.f, Zs.values, label='optimized')
plt.xlabel('f [mHz]')
plt.ylabel('Z-TDI real [1/Hz]')
plt.legend()
plt.subplot(234)
# plt.plot(Xs.f, Xs.values.imag, label='optimized')
plt.xlabel('f [mHz]')
plt.ylabel('X-TDI imag [1/Hz]')
plt.legend()
plt.subplot(235)
# plt.plot(Ys.f, Ys.values.imag, label='optimized')
plt.xlabel('f [mHz]')
plt.ylabel('Y-TDI imag [1/Hz]')
plt.legend()
plt.subplot(236)
# plt.plot(Zs.f, Zs.values.imag, label='optimized')
plt.xlabel('f [mHz]')
plt.ylabel('Z-TDI imag [1/Hz]')
plt.legend()

n_bin = 50
plt.figure(figsize=(20,30))
i = 0
number_of_parameters = 0
for name in samples.data_vars:
    number_of_parameters += 1
for name in samples.data_vars:
    i += 1
    if name != 'Likelihood':
        plt.suptitle("sampled posterior")
        plt.subplot(3,number_of_parameters,i)
        plt.axvline(x=pGB[name], color='r')
        plt.plot(samples[name],samples['Likelihood'], '.')

        plt.ylabel('Likelihood')
        plt.subplot(3,number_of_parameters,i+number_of_parameters)
        plt.axvline(x=pGB[name], color='r')
        n, bins, patches = plt.hist(samples[name], n_bin, density=True, facecolor='k', alpha=0.5)

        plt.subplot(3,number_of_parameters,i+2*number_of_parameters)
        if name == 'Frequency':
            plt.axvline(x=pGB[name]*1000, color='r')
            plt.plot(samples[name]*1000, range(number_of_samples), 'k')
        else:
            plt.axvline(x=pGB[name], color='r')
            plt.plot(samples[name], range(number_of_samples), 'k')
        plt.xlabel(name)
plt.show()

plt.figure()
plt.suptitle("sampled posterior")
plt.subplot(231)
plt.axvline(x=pGB['Amplitude'], color='r')
plt.plot(samples[:,0],likelihood_values, '.')
plt.xlabel('Amplitude')
plt.ylabel('Likelihood')
plt.subplot(232)
plt.axvline(x=pGB['Frequency'], color='r')
plt.plot(samples[:,1],likelihood_values, '.')
plt.xlabel('Frequeny')
plt.ylabel('Likelihood')
plt.figure()
plt.suptitle("sampled posterior")
plt.subplot(231)
plt.axvline(x=pGB['Amplitude'], color='r')
n, bins, patches = plt.hist(samples[:, 0], n_bin, density=True, facecolor='k', alpha=0.5)
plt.xlabel('Amplitude')

plt.subplot(232)
plt.axvline(x=pGB['Frequency'], color='r')
n, bins, patches = plt.hist(samples[:, 1], n_bin, density=True, facecolor='k', alpha=0.5)
plt.xlabel('Frequency')
# plt.subplot(233)
# plt.axvline(x=1.6, color='r')
# n, bins, patches = plt.hist(samples[:, 1], n_bin, density=True, facecolor='k', alpha=0.5)
# plt.xlabel('phi')

plt.subplot(234)
plt.axvline(x=pGB['Amplitude'], color='r')
plt.plot(samples[:, 0], range(number_of_samples), 'k')
plt.xlabel('Amplitude')
plt.subplot(235)
plt.axvline(x=pGB['Frequency'], color='r')
plt.plot(samples[:, 1], range(number_of_samples), 'k')
plt.xlabel('Frequency [Hz]')
plt.ylabel('sample number')
# plt.subplot(236)
# plt.axvline(x=1.6, color='r')
# plt.plot(samples[:, 2], range(number_of_samples), 'k')
# plt.xlabel('phi')
# plt.ylabel('sample number')
plt.show()

plt.figure()
s = template(Xp.t, samples[-1, 0], samples[-1, 1], samples[-1, 2])
plt.plot(Xp.t, Xp, label="TDI X")
plt.plot(Xp.t, s, 'g', linewidth=1, alpha=0.9, label='Found signal')
plt.xlabel('t [s]')
plt.legend()
plt.xlim(4*10**6, 4*10**6+10000)
plt.title('time series') 
plt.show()
plt.figure()
s = template(Xp.t, samples[-1, 0], samples[-1, 1], samples[-1, 2])
plt.plot(Xp.t, Xp, label="TDI X")
plt.plot(Xp.t, s, 'g', linewidth=1, alpha=0.9, label='Found signal')
plt.xlabel('t [s]')
plt.legend()
plt.xlim(8*10**6, 8*10**6+10000)
plt.title('time series')
plt.show()
plt.figure()
s = template(Xp.t, samples[-1, 0], samples[-1, 1], samples[-1, 2])
plt.plot(Xp.t, Xp, label="TDI X")
plt.plot(Xp.t, s, 'g', linewidth=1, alpha=0.9, label='Found signal')
plt.xlabel('t [s]')
plt.legend()
plt.xlim(9*10**6, 9*10**6+10000)
plt.title('time series')
plt.show()

plt.figure()
s = template( Xp.t, samples[-1, 0], samples[-1, 1], 2.3934989)
plt.plot(Xp.t[:1000], Xp[:1000], label="TDI X")
plt.plot(Xp.t[:1000], s[:1000], 'g', linewidth=1, alpha=0.9, label='Found signal')
plt.xlabel('t [s]')
plt.legend()
plt.title('time series')
plt.show()

plt.figure()
plt.specgram(Xp,Fs=1/(3*dt), NFFT=int(len(Xp.t)/3))
plt.xlabel('Time')
plt.ylabel('Frequency')
plt.ylim(0.0004,0.0005)

plt.figure()
plt.specgram(tdi_ts["X"][window:-window],Fs=1/(dt), NFFT=int(len(tdi_ts.t)))
plt.xlabel('Time')
plt.ylabel('Frequency')
plt.ylim(0.0004,0.0005)

print(samples[-1, :])

#%%
tdi_tsx_invmod = tdi_ts["X"][1:]/amplitude_envelope*amplitude_target
arr = Xs_td.shift(t=2)
subtracted = tdi_ts["X"]-arr
plt.figure(figsize=(24,6))
plt.subplot(131)
# plt.plot(tdi_ts["X"].t[cut:cut+length], tdi_ts["X"][cut:cut+length], label="TDI X")
plt.plot(tdi_ts["X"].t[:-1][cut:cut+length], Xs_td[cut:cut+length], label="TDI X fastGB")
plt.plot(arr.t[cut:cut+length], arr[cut:cut+length], label="TDI X fastGB")
# plt.xlim(cut,cut+length)
# plt.ylim(-10**-22,10**-22)
plt.plot(Xs_td.t[cut:cut+length], amplitude_envelope[cut:cut+length], label="Envelope")
plt.plot(arr.t[cut:cut+length], (arr/amplitude_envelope*amplitude_target/5000)[cut:cut+length], label="TDI X")
# plt.plot(tdi_ts["X"].t[cut+1:cut+length+1], tdi_tsx_invmod[cut:cut+length], label="TDI X inv")
# plt.plot(subtracted.t, subtracted, label="Subtracted")
plt.legend()

window = int(len(arr)/10)
plt.figure()
plt.specgram((arr/amplitude_envelope*amplitude_target/5000)[window:-window],Fs=1/(dt), NFFT=int(len(tdi_ts.t)/1))
plt.xlabel('Time')
plt.ylabel('Frequency')
plt.ylim(0.0004,0.0005)
plt.figure()
plt.specgram(arr[window:-window],Fs=1/(dt), NFFT=int(len(tdi_ts.t)/1))
plt.xlabel('Time')
plt.ylabel('Frequency')
plt.ylim(0.0004,0.0005)
plt.figure()
plt.specgram(tdi_tsx_invmod[window:-window],Fs=1/(dt), NFFT=int(len(tdi_ts.t)))
plt.xlabel('Time')
plt.ylabel('Frequency')
plt.ylim(0.0004,0.0005)
plt.figure()
plt.specgram(data_filtered[window:-window],Fs=1/(dt), NFFT=int(len(tdi_ts.t)))
plt.xlabel('Time')
plt.ylabel('Frequency')
plt.ylim(0.0004,0.0005)

plt.figure()
f, psdX =  scipy.signal.welch(tdi_ts["X"][window:-window], fs=1.0/dt, window='hanning', nperseg=256*256)
plt.loglog(f, np.sqrt(psdX), label="TDI X")
f, psdX =  scipy.signal.welch(tdi_tsx_invmod[window:-window], fs=1.0/dt, window='hanning', nperseg=256*256)
plt.loglog(f, np.sqrt(psdX), label="TDI X inv")
f, psdX =  scipy.signal.welch(arr[window:-window], fs=1.0/dt, window='hanning', nperseg=256*256)
plt.loglog(f, np.sqrt(psdX), label="TDI X")
plt.legend()
plt.xlabel("freq [Hz]")
plt.ylabel("PSD")
plt.xlim(0.0004,0.0005)
# plt.axis([1e-5, None, 4e-22, 2e-19])
#%%
plt.subplot(132)
plt.plot(tdi_ts["Y"].t, tdi_ts["Y"], label="TDI Y")
plt.plot(Ys_td.t, Ys_td, label="TDI X")
plt.plot(Ys_td.t, tdi_ts["Y"]-Ys_td, label="TDI Y")
plt.subplot(133)
plt.plot(tdi_ts["Z"].t, tdi_ts["Z"], label="TDI Z")
plt.plot(Zs_td.t, Zs_td, label="TDI Z")
plt.plot(Zs_td.t, tdi_ts["Z"]-Zs_td, label="TDI Z")

plt.figure(figsize=(12,6))
f, psdX =  scipy.signal.welch(tdi_ts["X"], fs=1.0/dt, window='hanning', nperseg=256*256)
fs, psdXs =  scipy.signal.welch(Xs_td.values, fs=1.0/dt, window='hanning', nperseg=256*256)
plt.loglog(f, np.sqrt(psdX), label="TDI X")
plt.loglog(f, np.sqrt(psdXs), label="TDI X")
plt.legend()
plt.xlabel("freq [Hz]")
plt.ylabel("PSD")
plt.axis([1e-5, None, 4e-22, 2e-19])
plt.show()


data = np.array([vgb["Name"], vgb["Frequency"], SNR2[:,0], SNR2[:,1]]).T
df = pd.DataFrame(data, columns=["Name", "f0", "fastGB", "TDI"])
print(df)
#######################################################################
def semi_fast_tdi(config, pMBHB, t_max, dt):
    hphc = HpHc.type("MBHB-%d"%s_index, "MBHB", "IMRPhenomD")
    hphc.set_param(pMBHB)
    orbits = Orbits.type(config)
    P = ProjectedStrain(orbits)    
    yArm = P.arm_response(0, t_max, dt, [hphc], tt_order=1)
    X = P.compute_tdi_x(np.arange(0, t_max, dt))
    return TimeSeries(X, dt=dt)

default_units = {'EclipticLatitude':'rad','EclipticLongitude':'rad',
         'PolarAngleOfSpin1':'rad','PolarAngleOfSpin2':'rad',
         'Spin1': '1','Spin2':'1',
         'Mass1':'Msun','Mass2':'Msun',
         'CoalescenceTime': 's','PhaseAtCoalescence':'rad',
         'InitialPolarAngleL':'rad','InitialAzimuthalAngleL':'rad',
         'Cadence': 's','Redshift': '1','Distance': 'Gpc',
         'ObservationDuration':'s'}

mbhb, units = hdfio.load_array(sangria_fn_training, name="sky/mbhb/cat")
print(units)
if not units:
    units = default_units
config = hdfio.load_config(sangria_fn, name="obs/config")
print(config)
s_index = 0
pMBHB = dict(zip(mbhb.dtype.names, mbhb[s_index]))
# units = default_units
for k,v in pMBHB.items():
    print(k)
    pMBHB[k] *= u.Unit(units[k])

t_max = float(tdi_ts["X"].t[-1]+tdi_ts["X"].attrs["dt"])
Xs = semi_fast_tdi(config, pMBHB, t_max, dt)

plt.figure(figsize=(12,6))
plt.plot(tdi_ts["X"].t, tdi_ts["X"], label="TDI X")
plt.plot(Xs.t, (tdi_ts["X"]-Xs), label="TDI X - fast %d"%s_index)
plt.axis([pMBHB["CoalescenceTime"]-600, pMBHB["CoalescenceTime"]+600, None, None])
plt.legend(loc="lower right")
plt.xlabel("time [s]")

config["t_max"].unit

def mldc_fast_tdi(pMBHB, t_max, dt):
    from GenerateFD_SignalTDIs import ComputeMBHBXYZ_FD
    from LISAhdf5 import ParsUnits
    hphc = HpHc.type("MBHB-%d"%s_index, "MBHB", "IMRPhenomD")
    hphc.set_param(pMBHB)
    pMBHB["Cadence"] = dt
    pMBHB["ObservationDuration"] = t_max/2
    pu = ParsUnits(pars_i=pMBHB, units_i=hphc.info())
    fr, Xs, Ys, Zs = ComputeMBHBXYZ_FD(pu)
    return (FrequencySeries(Xs, df=1/t_max), 
            FrequencySeries(Ys, df=1/t_max),
            FrequencySeries(Zs, df=1/t_max))

mbhb = hdfio.load_array(sangria_fn, name="sky/mbhb/cat", full_output=False)
config = hdfio.load_config(sangria_fn, name="obs/config")
noise_model = "MRDv1"
Nmodel = get_noise_model(noise_model, np.logspace(-5, -1, 100))
SNR2 = np.zeros((len(mbhb), 2)) 
for j,s in enumerate(mbhb):
    pMBHB = dict(zip(mbhb.dtype.names, s))
    Xs,Ys,Zs = mldc_fast_tdi(pMBHB, t_max, dt)
    source = dict({"X":Xs, "Y":Ys, "Z":Zs})
    SNR2[j,1] = compute_tdi_snr(source, Nmodel, data=tdi_fs)["tot2"]
    SNR2[j,0] = compute_tdi_snr(source, Nmodel)["tot2"]

import pandas as pd
data = np.array([range(len(mbhb)), mbhb["CoalescenceTime"], SNR2[:,0], SNR2[:,1]]).T
pd.DataFrame(data, columns=["index", "tc", "fast", "TDI"])
# %%
