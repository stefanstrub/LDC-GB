#%%
import matplotlib.pyplot as plt
import scipy
import numpy as np
import xarray as xr
from astropy import units as u
import pandas as pd

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

# filtered_fft = tdi_fs.copy()
# peak_freq = 0.0004726
# filtered_fft["X"][tdi_fs.f < peak_freq * 0.99] = 0+0j
# filtered_fft["X"][tdi_fs.f > peak_freq * 1.01] = 0+0j
# data_filtered = filtered_fft["X"].ts.ifft()
# data_filtered2 = TimeSeries(np.fft.irfft(filtered_fft["X"]),dt=dt)

plt.figure
plt.semilogx(tdi_fs["X"].f, tdi_fs["X"].real, label="TDI X")
# plt.loglog(filtered_fft["X"].f, filtered_fft["X"].real, label="TDI X")
plt.legend()
plt.xlabel("freq [Hz]")
plt.ylabel("PSD")
# plt.axis([1e-5, None, 4e-22, 2e-14])
plt.show()

template = lambda t, a, f, phi : a * np.sin(f* 2 * np.pi * t + phi)
sine_func = lambda t, a, f, phi: a * np.sin(f* 2 * np.pi * t + phi)

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
GB = fastGB.FastGB(delta_t=dt, T=float(tdi_ts["X"].t[-1])) # in seconds
pGB = dict(zip(vgb.dtype.names, vgb[8])) # we take the source #8
#modify pGB
# pGB['InitialPhase'] *= 1.001
# pGB['Amplitude'] *= 1.001
# pGB['Frequency'] *= 1.01
# pGB['FrequencyDerivative'] *= 1.1
# pGB['Polarization'] += 0.01
# pGB['EclipticLatitude'] += 0.01
# pGB['EclipticLongitude'] += 0.01
# pGB['Inclination'] += 0.01
Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB, oversample=4, simulator='synthlisa')
Xs_td, Ys_td, Zs_td = GB.get_td_tdixyz(template=pGB, simulator='synthlisa')

analytic_signal = scipy.signal.hilbert(Xs_td)
amplitude_envelope = np.abs(analytic_signal)
frequency_target = pGB['Frequency']
amplitude_target = pGB['Amplitude']
phi_target = pGB['InitialPhase']

#%%
Xp = Xs_td/amplitude_envelope*amplitude_target
number_segments = 24
Xp_segments = []
params = np.zeros((number_segments,3))
for i in range(number_segments):
    Xp_segments.append(Xp[int(i/number_segments*len(Xp.t)):int((1+i)/number_segments*len(Xp.t))])
    params[i], params_covariance = scipy.optimize.curve_fit(sine_func, Xp_segments[-1].t, Xp_segments[-1]/amplitude_target, p0=[1,frequency_target*1.1,0], sigma=0.000001*np.ones(len(Xp_segments[-1].t)))
    # print((abs(params)- [vgb[8][1], vgb[8][4], vgb[8][6]]))
    print(params[i])
    # print((abs(params[i])- [vgb[8][1], vgb[8][4], vgb[8][6]])/[vgb[8][1], vgb[8][4], vgb[8][6]])
    plt.figure()
    s = template( Xp_segments[-1].t, params[i][0], params[i][1], params[i][2])
    plt.plot(Xp_segments[-1][:10000].t, Xp_segments[-1][:10000]/amplitude_target, label="TDI X")
    plt.plot(Xp_segments[-1][:10000].t, s[:10000], 'g', linewidth=1, alpha=0.9, label='Found signal')
    plt.xlabel('t [s]')
    plt.legend()
    plt.title('time series')
    plt.show()
    frequency_correction = params[:,1]- pGB['Frequency']


#%%
#modify pGB
# pGB['InitialPhase'] *= 1.01
# pGB['Amplitude'] *= 1.01
# pGB['Frequency'] *= 1.01
# pGB['FrequencyDerivative'] *= 1.1
# pGB['Polarization'] *= 1.001
# pGB['EclipticLatitude'] += 0.05
# pGB['EclipticLongitude'] += 0.01
# pGB['Inclination'] *= 1.001
Xs_td2, Ys_td2, Zs_td2 = GB.get_td_tdixyz(template=pGB, simulator='synthlisa')
analytic_signal2 = scipy.signal.hilbert(Xs_td2)
amplitude_envelope2 = np.abs(analytic_signal2)

pGB2['EclipticLatitude'] -= 0.05
#modify pGB
# pGB['InitialPhase'] *= 1.01
# pGB['Amplitude'] *= 1.01
pGB2['Frequency'] *= 1.1
# pGB['FrequencyDerivative'] *= 1.1
# pGB['Polarization'] *= 1.001
# pGB['EclipticLatitude'] += 0.05
# pGB['EclipticLongitude'] += 0.01
# pGB['Inclination'] *= 1.001
Xs_td2, Ys_td2, Zs_td2 = GB.get_td_tdixyz(template=pGB2, simulator='synthlisa')
analytic_signal2 = scipy.signal.hilbert(Xs_td2)
amplitude_envelope3 = np.abs(analytic_signal2)


# plt.figure(figsize=(24,6))
# plt.subplot(131)
# # plt.plot(Xs_td.t, Xs_td2, label="TDI X")
# plt.plot(Xs_td.t, Xs_td, label="TDI X")
# # plt.plot(Xp.t, amplitude_envelope, label="envelope")
# # plt.plot(Xs_td.t[::100], amplitude_envelope2[::100], label="envelope 2")
# plt.legend()
# plt.show()


# plt.figure(figsize=(24,6))
# plt.subplot(132)
# plt.plot(Ys_td.t, Ys_td, label="TDI Y")
# plt.subplot(133)
# plt.plot(Zs_td.t, Zs_td, label="TDI Z")
# plt.show()

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

#%%

vgb, units = hdfio.load_array(sangria_fn_training, name="sky/vgb/cat")
GB = fastGB.FastGB(delta_t=dt, T=float(tdi_ts["X"].t[-1])) # in seconds
noise_model = "MRDv1"
Nmodel = get_noise_model(noise_model, np.logspace(-5, -1, 100))
SNR2 = np.zeros((len(vgb), 2)) # snr square
Xs_td.values = np.zeros(len(Xs_td.values))
Ys_td.values = np.zeros(len(Ys_td.values))
Zs_td.values = np.zeros(len(Zs_td.values))
for j,s in enumerate(vgb):
    pGB2 = dict(zip(vgb.dtype.names, s))
    Xs2, Ys2, Zs2 = GB.get_fd_tdixyz(template=pGB2, oversample=4, simulator='synthlisa')
    fmin, fmax = float(Xs2.f[0]) , float(Xs2.f[-1]+Xs2.attrs['df'])
    source = dict({"X":Xs2, "Y":Ys2, "Z":Zs2})
    SNR2[j,1] = compute_tdi_snr(source, Nmodel, data=tdi_fs, fmin=fmin, fmax=fmax)["tot2"]
    SNR2[j,0] = compute_tdi_snr(source, Nmodel)["tot2"] 

    Xs_tdvb, Ys_tdvb, Zs_tdvb = GB.get_td_tdixyz(template=pGB2, simulator='synthlisa')
    Xs_td += Xs_tdvb
    Ys_td += Ys_tdvb
    Zs_td += Zs_tdvb

plt.figure(figsize=(12,6))
plt.subplot(121)
plt.title("real part")
plt.plot(tdi_fs["X"].f, tdi_fs["X"].real, label="TDI X")
plt.plot(Xs.f, (tdi_fs["X"][Xs.kmin:Xs.kmin+len(Xs)]-Xs.values).real, label="TDI X - fast "+pGB["Name"])
plt.axis([pGB["Frequency"]-6e-7, pGB["Frequency"]+6e-7, -3e-17, 5e-17])
plt.legend(loc="lower right")
plt.xlabel("freq [Hz]")
plt.subplot(122)
plt.title("imaginary part")
plt.plot(tdi_fs["X"].f, tdi_fs["X"].imag, label="TDI X")
plt.plot(Xs.f, (tdi_fs["X"].isel(f=slice(Xs.kmin, Xs.kmin+len(Xs)))-Xs.values).imag, label="TDI X - fast "+pGB["Name"])
plt.axis([pGB["Frequency"]-6e-7, pGB["Frequency"]+6e-7, -3e-17, 5e-17])
plt.legend(loc="lower left")
plt.xlabel("freq [Hz]")
#%%
Xp = Xs_td/amplitude_envelope*amplitude_target
Xp_segments2 = []
params2 = np.zeros((number_segments,3))
for i in range(number_segments):
    Xp_segments2.append(Xp[int(i/number_segments*len(Xp.t)):int((1+i)/number_segments*len(Xp.t))])
    params2[i], params_covariance = scipy.optimize.curve_fit(sine_func, Xp_segments2[-1].t, Xp_segments2[-1], p0=samples[0])
frequency_calculated = params2[:,1] - frequency_correction
print(frequency_calculated)
#%%
tdi_tsc = tdi_ts.copy(deep=True)
tdi_fsc = tdi_fs.copy(deep=True)
tdi_tsc["X"] = Xs_td#/amplitude_envelope*amplitude_target
# Xp = tdi_ts["X"]/amplitude_envelope*amplitude_target
# Xp = Xp[int(0.1*len(Xp.t)):-int(0.1*len(Xp.t))]
# Xp = Xp[int(0.499*len(Xp.t)):-int(0.499*len(Xp.t))]
# Xp = Xp_segments[2]
# Xp = Xp[:int(1/1000*len(Xp.t))]
# Xp = data_filtered2/amplitude_envelope*amplitude_target
plt.figure(figsize=(12,6))
plt.plot(tdi_tsc["X"].t/(24*3600), tdi_tsc["X"], label="VGB")
plt.plot(tdi_tsc["X"][:-1].t/(24*3600), tdi_tsc["X"][:-1]/amplitude_envelope*amplitude_target/5000, label="VGB-inv")
plt.plot(tdi_tsc["X"][:-1].t/(24*3600), amplitude_envelope, label="Envelope")
plt.plot(tdi_tsc["X"][:-1].t/(24*3600), amplitude_envelope2, label="Envelope_lat+0.05")
plt.plot(tdi_tsc["X"][:-1].t/(24*3600), amplitude_envelope3, label="Envelope_f*1.1")
# plt.plot(tdi_tsc8.t/(24*3600), tdi_tsc8/amplitude_envelope[::3]*amplitude_target/5000, label="VGB8-inv")
plt.xlabel("Time [days]")
plt.ylabel("X-TDI strain")
plt.legend()
# plt.savefig("pictures/VGB88-inv.png")
plt.show()

tdi_tsc_inv = tdi_tsc.copy(deep=True)
tdi_tsc_inv["X"][:-1].values = tdi_tsc["X"][:-1].values/amplitude_envelope*amplitude_target/5000
tdi_fsc_inv = xr.Dataset(dict([(k,tdi_tsc_inv[k][:-1].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
tdi_fsc = xr.Dataset(dict([(k,tdi_tsc[k][:-1].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))

plt.figure(figsize=(12,6))
plt.subplot(121)
plt.title("real part")
plt.plot(Xs.f, (tdi_fsc["X"][Xs.kmin:Xs.kmin+len(Xs)]-Xs.values).real, label="TDI X - fast "+pGB["Name"])
plt.plot(Xs.f, (tdi_fsc_inv["X"][Xs.kmin:Xs.kmin+len(Xs)]).real, label="TDI X - fast "+pGB["Name"])
plt.axis([pGB["Frequency"]-6e-7, pGB["Frequency"]+6e-7, -3e-17, 5e-17])
plt.legend(loc="lower right")
plt.xlabel("freq [Hz]")
plt.subplot(122)
plt.title("imaginary part")
plt.plot(tdi_fsc.f, (tdi_fsc_inv.isel(f=slice(Xs.kmin, Xs.kmin+len(Xs)))-Xs.values).imag, label="TDI X - fast "+pGB["Name"])
plt.axis([pGB["Frequency"]-6e-7, pGB["Frequency"]+6e-7, -3e-17, 5e-17])
plt.legend(loc="lower left")
plt.xlabel("freq [Hz]")
plt.savefig("pictures/fs.png")

plt.figure(figsize=(12,6))
plt.plot(tdi_tsc_inv.t[cut:cut+length]/(24*3600), tdi_tsc_inv[cut:cut+length], label="VGB8")
plt.xlabel("Time [days]")
plt.ylabel("X-TDI strain")
plt.legend()
plt.savefig("pictures/VGB88-inv.png")
plt.show()
######################################################
#%%
#find GB parameters

def prior_data(data, s2, noise_level):

    p = np.sum((data - s2) ** 2) / (len(data) * noise_level ** 2)
    # print(p)
    return np.exp(-p / 2.0)


# Number of histogram bins.
n_bin = 50
# Total number of proposed samples.
number_of_samples = 10 ** 2
number_of_parameters = 3

# Make the first random sample. ------------------------------------

samples = np.zeros((number_of_samples, number_of_parameters))

# Random time shift between 10 and 30.
samples[0, 0] = 10**-20 * np.random.rand()
samples[0, 1] = 10**-3 * np.random.rand()
samples[0, 2] = 2 * np.pi * np.random.rand()
# samples[0, 0] = amplitude_target
# samples[0, 1] = frequency_target
# samples[0, 2] = 2

params, params_covariance = scipy.optimize.curve_fit(sine_func, Xp.t, Xp, p0=samples[0])
# print(params, vgb[8][1], vgb[8][4], vgb[6])
# Make test time series by shifting s2.
s = template(Xp.t, samples[0, 0], samples[0, 1], samples[0, 2])

noise_level = 1*10**-17
# Evaluate posterior for the first sample.
p = prior_data(s, Xp, noise_level)

# Run Metropolis-Hastings sampler. ---------------------------------

for i in range(1, number_of_samples):

    # Normal distributed proposal.
    std = np.array(np.eye(3)*[5*10**-22,5*10**-5,0.1])
    a_test, f_test, phi_test = np.random.multivariate_normal(samples[i-1,:], std**2)
    f_test = 10**-3 * np.random.rand()
    a_test = 10**-20 * np.random.rand()
    phi_test = 2 * np.pi * np.random.rand()
    if a_test < 10**-21:
        a_test = 10**-21
    s = template( Xp.t,a_test, f_test, phi_test)
    
    # Evaluate posterior.
    p_test = prior_data(s, Xp, noise_level)
    T_inv =  10**4
    # Apply Metropolis rule.
    # print(p_test.values,p.values,(p_test / p) ** T_inv)
    if (p_test / p) ** T_inv > np.random.rand():  # L^i/L^j
        p = p_test
        samples[i, 0] = a_test
        samples[i, 1] = f_test
        samples[i, 2] = phi_test
    else:
        samples[i] = samples[i - 1]

# plot
n_bin = 50
plt.figure()
plt.suptitle("sampled posterior")
plt.subplot(231)
plt.axvline(x=amplitude_target, color='r')
n, bins, patches = plt.hist(samples[:, 0], n_bin, density=True, facecolor='k', alpha=0.5)
plt.xlabel('amplitude')


plt.subplot(232)
plt.axvline(x=frequency_target, color='r')
n, bins, patches = plt.hist(samples[:, 1], n_bin, density=True, facecolor='k', alpha=0.7)
plt.xlabel('frequency')

plt.subplot(233)
plt.axvline(x=1.6, color='r')
n, bins, patches = plt.hist(samples[:, 2], n_bin, density=True, facecolor='k', alpha=0.5)
plt.xlabel('phi')

plt.subplot(234)
plt.axvline(x=amplitude_target, color='r')
plt.plot(samples[:, 0], range(number_of_samples), 'k')
plt.xlabel('amplitude')
plt.subplot(235)
plt.axvline(x=frequency_target, color='r')
plt.plot(samples[:, 1], range(number_of_samples), 'k')
plt.xlabel('frequency [Hz]')
plt.ylabel('sample number')
plt.subplot(236)
plt.axvline(x=1.6, color='r')
plt.plot(samples[:, 2], range(number_of_samples), 'k')
plt.xlabel('phi')
plt.ylabel('sample number')

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
