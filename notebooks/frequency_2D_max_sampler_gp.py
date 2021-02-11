#%%
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

import torch
import gpytorch
from sklearn.metrics import mean_squared_error

from botorch.models.gpytorch import GPyTorchModel
from botorch.acquisition.monte_carlo import qExpectedImprovement, qUpperConfidenceBound
from botorch.optim import optimize_acqf



# use a GPU if available
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
dtype = torch.float

DATAPATH = "/home/stefan/LDC/Sangria/data"
sangria_fn = DATAPATH+"/dgb-tdi.h5"
# sangria_fn = DATAPATH+"/LDC2_sangria_blind_v1.h5"
# sangria_fn = DATAPATH+"/LDC2_sangria_gdb-tdi_v1_v3U3MxS.h5"
# sangria_fn = DATAPATH+"/LDC2_sangria_idb-tdi_v1_DgtGV85.h5"
# sangria_fn = DATAPATH+"/LDC2_sangria_mbhb-tdi_v1_MN5aIPz.h5"
sangria_fn = DATAPATH+"/LDC2_sangria_training_v1.h5"
tdi_ts, tdi_descr = hdfio.load_array(sangria_fn, name="obs/tdi")
# sangria_fn = DATAPATH+"/LDC2_sangria_vgb-tdi_v1_sgsEVXb.h5"
# tdi_ts, tdi_descr = hdfio.load_array(sangria_fn)
sangria_fn_training = DATAPATH+"/LDC2_sangria_training_v1.h5"
dt = int(1/(tdi_descr["sampling_frequency"]))

# Build timeseries and frequencyseries object for X,Y,Z
tdi_ts = xr.Dataset(dict([(k,TimeSeries(tdi_ts[k], dt=dt)) for k in ["X", "Y", "Z"]]))
# tdi_ts = xr.Dataset(dict([(k,TimeSeries(tdi_ts[k][:,1], dt=dt)) for k in ["X", "Y", "Z"]]))
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
GB = fastGB.FastGB(delta_t=dt, T=float(tdi_ts["X"].t[-1])) # in seconds
pGB = dict(zip(vgb.dtype.names, vgb[8])) # we take the source #8
#modify pGB
# pGB['InitialPhase'] = np.random.rand()*np.pi*2
# pGB['Amplitude'] *= 1.01
# pGB['Frequency'] *= 1.05
# pGB['FrequencyDerivative'] *= 1.1
# pGB['Polarization'] = np.random.rand()*np.pi*2
# pGB['EclipticLatitude'] += 0.5
# pGB['EclipticLongitude'] += 0.5
# pGB['Inclination'] += 0.1
Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB, oversample=4, simulator='synthlisa')
fmin, fmax = float(Xs.f[0]) , float(Xs.f[-1]+Xs.attrs['df'])
source = dict({"X":Xs, "Y":Ys, "Z":Zs})
start = time.time()
Xs_td, Ys_td, Zs_td = GB.get_td_tdixyz(template=pGB, simulator='synthlisa')


# plt.figure(figsize=(12,3))
# # plt.plot(Xs_td.t, Xs_td2, label="TDI X")
# # plt.plot(tdi_ts_training['X'].t, tdi_ts_training['X'], label="data")
# plt.plot(tdi_ts['X'].t/86400, tdi_ts['X'], label="Verification Binaries")
# plt.plot(Xs_td.t/86400, Xs_td, label="Binary")
# # plt.plot(Xp.t, amplitude_envelope, label="envelope")
# # plt.plot(Xs_td.t[::100], amplitude_envelope2[::100], label="envelope 2")
# plt.ylabel('X TDI strain')
# plt.xlabel('time [days]')
# # plt.ylim(-10**-20,10**-20)
# plt.legend()
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
def loglikelihood(pGBs):
    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs, oversample=4, simulator='synthlisa')
    index_low = np.searchsorted(Xs.f, dataX.f[0])
    Xs = Xs[index_low:index_low+len(dataX)]
    Ys = Ys[index_low:index_low+len(dataY)]
    Zs = Zs[index_low:index_low+len(dataZ)]
    diff = np.abs(dataX - Xs.values)**2 + np.abs(dataY - Ys.values)**2 + np.abs(dataZ - Zs.values)**2
    # p1 = -float(np.sum(diff / Sn)*Xs.attrs['df'])/2.0
    p1 = -float(np.sum(diff / (Sn+noise))/len(diff))/2.0
    # p1 = np.exp(p1)
    return p1

# Number of histogram bins.
n_bin = 50
# Total number of proposed samples.
number_of_samples = 8*10 **1
cutoff_ratio = 1000

parameters = ['Amplitude','EclipticLatitude','EclipticLongitude','Frequency','FrequencyDerivative','Inclination','InitialPhase','Polarization']
parametersfd = ['Amplitude','EclipticLatitude','EclipticLongitude','Frequency','Inclination','InitialPhase','Polarization']
boundaries = {'Amplitude': [10**-22.0, 5*10**-21.0],'EclipticLatitude': [-1.0, 1.0],
'EclipticLongitude': [0.0, 2.0*np.pi],'Frequency': [0.0004725, 0.0004727],'FrequencyDerivative': [10**-20.0, 10**-18.0],
'Inclination': [-1.0, 1.0],'InitialPhase': [0.0, 2.0*np.pi],'Polarization': [0.0, 1.0*np.pi]}

# previous_max = [0.2090, 0.1000, 0.8469, 0.5276, 0.7168, 0.9667, 0.0970, 0.0000]
# previous_max = [0.2090, 0.2667, 0.7333, 0.5276, 0.7168, 0.9667, 0.0970, 0.0000]
# previous_max = [0.2090, 0.3333, 0.7333, 0.5276, 0.7168, 0.8667, 0.0970, 0.0000]
# previous_max = [0.2090, 0.3667, 0.6667, 0.5276, 0.7168, 0.8667, 0.0970, 0.0000]
# previous_max = [0.4000, 0.3667, 0.6667, 0.5276, 0.7168, 0.8667, 0.0970, 0.0000]
# previous_max = [0.4000, 0.3667, 0.6667, 0.5276, 0.7667, 0.8667, 0.0970, 0.0000]
# previous_max = [0.2333, 0.3667, 0.6667, 0.5276, 0.9667, 0.8667, 0.0970, 0.0000]
previous_max = [0.333,0.54,0.134,0.7345,0.6456,0.2645,0.8216,0.000]
i = 0
for parameter in parameters:
    if parameter in ['FrequencyDerivative']:
        i -= 1
    elif parameter in ['EclipticLatitude']:
        pGB[parameter] = np.arcsin((previous_max[i]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
    elif parameter in ['Inclination']:
        pGB[parameter] = np.arccos((previous_max[i]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
    else:
        pGB[parameter] = (previous_max[i]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
        print(parameter, pGB[parameter],previous_max[i])
    i += 1

# print(pGB)
# pGB['EclipticLatitude'] = 0.0
# pGB = {'Amplitude': 4.150747034e-21, 'EclipticLatitude': 0.0, 'EclipticLongitude': 1.6707963267948967, 'Frequency': 0.000472612, 'FrequencyDerivative': 1.440030107143638e-19, 'Inclination': 1.4468779499032993, 'InitialPhase': 1.7, 'Name': 'CD-30o11223', 'Polarization': 1.7}
# print(pGB)
# pGB = {'Amplitude': 1.1240340491284739e-21, 'EclipticLatitude': 0.48579612338022954, 'EclipticLongitude': 5.3215032928260895, 'Frequency': 0.0004726055124556648, 'FrequencyDerivative': 7.422865032688086e-19, 'Inclination': 1.12226944780215, 'InitialPhase': 2.391197244408699, 'Name': 'CD-30o11223', 'Polarization': 0.6092901760377935}
# Make the first random sample. ------------------------------------
pGBs = deepcopy(pGB)
# pGBs['Amplitude'] *= 1.0
# pGBs['EclipticLatitude'] = (np.random.random()-0.5) * np.pi 
# pGBs['EclipticLongitude'] = np.random.random() * np.pi 
# pGBs['InitialPhase'] = np.random.random() * np.pi 
# pGBs['Frequency'] *= 1.001
# pGBs['FrequencyDerivative'] *= 1.01
# pGBs['Polarization'] = np.random.random() * np.pi 
# pGBs['Inclination'] = np.random.random()* np.pi 
Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBs, oversample=4, simulator='synthlisa')
psd_signal = np.abs(Xs.values)**2 + np.abs(Ys.values)**2 + np.abs(Zs.values)**2
highSNR = psd_signal > np.max(psd_signal)/cutoff_ratio
dataX = tdi_fs["X"].isel(f=slice(Xs.kmin, Xs.kmin+len(Xs)))[highSNR]
dataY = tdi_fs["Y"].isel(f=slice(Ys.kmin, Ys.kmin+len(Ys)))[highSNR]
dataZ = tdi_fs["Z"].isel(f=slice(Zs.kmin, Zs.kmin+len(Zs)))[highSNR]
spd_data = np.abs(dataX)**2 + np.abs(dataY)**2 + np.abs(dataZ)**2
noise = (np.mean(spd_data[:2])+np.mean(spd_data[-3:])).values/2
Xs, Ys, Zs = Xs[highSNR], Ys[highSNR], Zs[highSNR]
fmin, fmax = float(Xs.f[0]) , float(Xs.f[-1]+Xs.attrs['df'])
freq = np.array(Xs.sel(f=slice(fmin, fmax)).f)
Sn = Nmodel.psd(freq=freq, option='X')
diff = np.abs(dataX - Xs.values)**2 + np.abs(dataY - Ys.values)**2 + np.abs(dataZ - Zs.values)**2
p1 = -float(np.sum(diff / (Sn+noise))/len(diff))/2.0
# p1 = np.exp(p1)

frequency_lower_boundary = dataX.f[0].values+(dataX.f[-1].values - dataX.f[0].values)*4/10
frequency_upper_boundary = dataX.f[-1].values-(dataX.f[-1].values - dataX.f[0].values)*3/10
# [0.0004725, 0.0004727]
# boundaries_small = deepcopy(boundaries)
# part_ratio = 10
# for parameter in parameters:
#     if parameter in ['EclipticLongitude','Frequency']:
#         boundaries_small[parameter] = [pGB[parameter]-(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio,pGB[parameter]+(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio]
#     if parameter == 'EclipticLatitude':
#         boundaries_small[parameter] = [np.sin(pGB[parameter])-(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio,np.sin(pGB[parameter])+(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio]
#     elif parameter == 'Inclination':
#         boundaries_small[parameter] = [np.cos(pGB[parameter])-(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio,np.cos(pGB[parameter])+(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio]
#     else:
#         boundaries_small[parameter] = [pGB[parameter]-(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio,pGB[parameter]+(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio]
# boundaries = boundaries_small



samples = xr.Dataset(dict([(name,xr.DataArray(np.zeros(number_of_samples), dims=('number_of_sample'), coords={"number_of_sample": range(number_of_samples)},
                         )) for name, titles in pGBs.items()]))
pGBs01 = {}
for parameter in parameters:
    if parameter in ['EclipticLatitude']:
        samples[parameter][0] = (np.sin(pGBs[parameter])-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0])
    elif parameter in ['Inclination']:
        samples[parameter][0] = (np.cos(pGBs[parameter])-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0])
    else:
        samples[parameter][0] = (pGBs[parameter]-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0])
    pGBs01[parameter] = samples[parameter][0]
samples['Likelihood'] = samples['Name']
samples = samples.drop(['Name'])
samples['Likelihood'][0] = p1


# Run Metropolis-Hastings sampler. ---------------------------------
plt.figure()
ax1=plt.subplot(231)
# plt.plot(dataX_training.f*1000,dataX_training.values, label='data')
ax1.plot(dataX.f*1000,dataX.values.real, label='binary')
ax1.plot(Xs.f*1000, Xs.values.real, label='start')
ax2=plt.subplot(232)
# plt.plot(dataY_training.f*1000,dataY_training.values, label='data')
ax2.plot(dataY.f*1000,dataY.values.real, label='binary')
ax2.plot(Ys.f*1000, Ys.values.real, label='start')
ax3=plt.subplot(233)
# plt.plot(dataZ_training.f*1000,dataZ_training.values, label='data')
ax3.plot(dataZ.f*1000,dataZ.values.real, label='binary')
ax3.plot(Zs.f*1000, Zs.values.real, label='start')
ax4=plt.subplot(234)
# plt.plot(dataX_training.f*1000,dataX_training.values.imag, label='data')
ax4.plot(dataX.f*1000,dataX.values.imag, label='binary')
ax4.plot(Xs.f*1000, Xs.values.imag, label='start')
ax5=plt.subplot(235)
# plt.plot(dataY_training.f*1000,dataY_training.values.imag, label='data')
ax5.plot(dataY.f*1000,dataY.values.imag, label='binary')
ax5.plot(Ys.f*1000, Ys.values.imag, label='start')
ax6=plt.subplot(236)
# plt.plot(dataZ_training.f*1000,dataZ_training.values.imag, label='data')
ax6.plot(dataZ.f*1000,dataZ.values.imag, label='binary')
ax6.plot(Zs.f*1000, Zs.values.imag, label='start')
print('p1',p1)

def sampler(number_of_samples,parameters,pGB,boundaries,p1, uniform=False, MCMC=False, only=False, onlyparameter='Frequency', twoD=False, secondparameter='Amplitude'):
    samples = xr.Dataset(dict([(name,xr.DataArray(np.zeros(number_of_samples), dims=('number_of_sample'), coords={"number_of_sample": range(number_of_samples)},
                         )) for name, titles in pGB.items()]))
    samples = {}
    pGBs01 = {}
    for parameter in parameters:
        samples[parameter] = []
        if parameter in ['EclipticLatitude']:
            samples[parameter].append((np.sin(pGB[parameter])-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0]))
        elif parameter in ['Inclination']:
            samples[parameter].append((np.cos(pGB[parameter])-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0]))
        else:
            samples[parameter].append((pGB[parameter]-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0]))
        pGBs01[parameter] = samples[parameter][0]
    samples['Likelihood'] = []
    samples['Likelihood'].append(p1)
    start = time.time()
    j = 0
    number_of_sampels_sqrt = np.sqrt(number_of_samples)
    for i in range(1, number_of_samples):
        if only:
            parameter = onlyparameter
            if uniform:
                pGBs01[parameter] = i/number_of_samples
            else:
                pGBs01[parameter] = np.random.rand()
            if parameter != 'InitialPhase':
                pGBs01['InitialPhase'] = samples['InitialPhase'][0]
            if parameter != 'Polarization':
                pGBs01['Polarization'] = samples['Polarization'][0]
        elif twoD:
            parameter = onlyparameter
            parameter2 = secondparameter
            if uniform:
                if i % number_of_sampels_sqrt == 0:
                    j += 1
                pGBs01[parameter] = ((i-1)%number_of_sampels_sqrt)/number_of_sampels_sqrt
                pGBs01[parameter2] = j/number_of_sampels_sqrt
            else: 
                pGBs01[parameter] = np.random.rand()
                pGBs01[parameter2] = np.random.rand()
            # if parameter != 'InitialPhase':
            #     pGBs01['InitialPhase'] = samples['InitialPhase'][0]
        else:
            for parameter in parameters:
                pGBs01[parameter] = np.random.rand()
                # if parameter in ['Amplitude']:#,'FrequencyDerivative','Amplitude','EclipticLongitude']:
                #     pGBs01[parameter] = i/number_of_samples
                # elif parameter in ['FrequencyDerivative']:
                #     pass
                # if parameter in ['Inclination']:
                #     pGBs01[parameter] = np.random.rand()
                # elif parameter in ['Inclination']:
                #     pGBs01[parameter] = np.random.rand()
                # else:
                #     pGBs01[parameter] = np.random.rand()
        for parameter in parameters:
            # if parameter in ['FrequencyDerivative']:
            #     pass
            if parameter in ['EclipticLatitude']:
                pGBs[parameter] = np.arcsin((pGBs01[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
            elif parameter in ['Inclination']:
                pGBs[parameter] = np.arccos((pGBs01[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
            else:
                pGBs[parameter] = (pGBs01[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
        p_test = loglikelihood(pGBs)
        p1 = p_test
        s = 1/p1
        for parameter in parameters:
            samples[parameter].append(pGBs01[parameter])
        samples['Likelihood'].append(p1)
    print('sampler time',time.time()- start)
    for parameter in parameters:
        samples[parameter] = np.asarray(samples[parameter])
    samples['Likelihood'] = np.asarray(samples['Likelihood'])
    return samples

samples = sampler(number_of_samples,parameters,pGBs,boundaries,p1,uniform=False, MCMC=False)
print('samples',len(samples['Amplitude']))

#%%

class ExactGPModel(gpytorch.models.ExactGP, GPyTorchModel):
    _num_outputs = 1  # to inform GPyTorchModel API
    def __init__(self, train_x, train_y, likelihood, kernel):
        super(ExactGPModel, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ConstantMean()
        kernelA = gpytorch.kernels.RBFKernel(active_dims=(0),lengthscale_constraint=gpytorch.constraints.GreaterThan(0.2))
        kernelLat = gpytorch.kernels.PeriodicKernel(active_dims=torch.tensor([1]))
        kernelLong = gpytorch.kernels.PeriodicKernel(active_dims=torch.tensor([2]))
        kernelF = gpytorch.kernels.RBFKernel(active_dims=(3),lengthscale_constraint=gpytorch.constraints.Interval(0.06,0.1))
        # kernelFD = gpytorch.kernels.RBFKernel(active_dims=(4))
        kernelI = gpytorch.kernels.PolynomialKernel(power= 2,active_dims=(4))
        # kernelIP = gpytorch.kernels.PeriodicKernel(active_dims=torch.tensor([6]))
        kernelIP = gpytorch.kernels.ScaleKernel(gpytorch.kernels.CosineKernel(active_dims=torch.tensor([5]),period_length_constraint=gpytorch.constraints.Interval(0.499,0.501)), outputscale_constraint=gpytorch.constraints.Interval(0.9,1.1))
        # kernelP = gpytorch.kernels.PeriodicKernel(active_dims=torch.tensor([7]))
        # kernelP = gpytorch.kernels.CosineKernel(active_dims=torch.tensor([6]))*gpytorch.kernels.CosineKernel(active_dims=torch.tensor([6]))
        kernelP = gpytorch.kernels.ScaleKernel(gpytorch.kernels.CosineKernel(active_dims=torch.tensor([6]),period_length_constraint=gpytorch.constraints.Interval(0.249,0.251)), outputscale_constraint=gpytorch.constraints.Interval(0.9,1.1))
        # kernel = kernelA + kernelLat + kernelLong + kernelF + kernelFD + kernelI + kernelIP + kernelP
        # kernel = kernelA * kernelLat * kernelLong * kernelF * kernelI * kernelIP * kernelP
        self.covar_module = gpytorch.kernels.ScaleKernel(gpytorch.kernels.RBFKernel(ard_num_dims=8))
        self.to(train_x)  # make sure we're on the right device/dtype
  
    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)


train_x = np.zeros((len(samples['Amplitude'])-1,len(parameters)))
i = 0
for name in parametersfd:
    train_x[:,i] = samples[name][1:]
    i +=1
nu = np.mean(samples['Likelihood'][1:])
sigma = np.std(samples['Likelihood'][1:])
train_y = (samples['Likelihood'][1:]-nu)/sigma
train_x = torch.from_numpy(train_x).float()
train_y = torch.from_numpy(train_y).float()

number_of_test_samples = 20
test_x = {}
test_y = {}
for parameter in parametersfd:
    test_samples = sampler(number_of_test_samples,parameters,pGB,boundaries,p1, uniform= True, only=True, onlyparameter=parameter)
    test_x[parameter] = np.zeros((number_of_test_samples,len(parameters)))
    i = 0
    for name in parametersfd:
        test_x[parameter][:,i] = test_samples[name]
        i +=1
    test_y[parameter] = test_samples['Likelihood']
    test_x[parameter] = torch.from_numpy(test_x[parameter]).float()
    test_x[parameter] = test_x[parameter]
    test_y[parameter] = torch.from_numpy(test_y[parameter]).float()

# test_samples = sampler(number_of_test_samples,parameters,pGB,boundaries,p1)
# parameter = 'random'
# test_x[parameter] = np.zeros((number_of_test_samples,len(parameters)))
# i = 0
# for name in parametersfd:
#     test_x[parameter][:,i] = test_samples[name]
#     i +=1
# test_y[parameter] = test_samples['Likelihood']
# test_x[parameter] = torch.from_numpy(test_x[parameter]).float()
# test_x[parameter] = test_x[parameter]
# test_y[parameter] = torch.from_numpy(test_y[parameter]).float()

def planemaxsearch(maxpGB,parameterstocheck, parameter2, resolution):
    for parameter in parameterstocheck:
        # if parameter != parameter2:
        train_samples = sampler(resolution,parameters,maxpGB,boundaries,p1, uniform= False, twoD = True, onlyparameter=parameter, secondparameter=parameter2)
        train_x = np.zeros((resolution,len(parameters)))
        i = 0
        for name in parametersfd:
            train_x[:,i] = train_samples[name]
            i +=1
        train_y = train_samples['Likelihood']
        train_x = torch.from_numpy(train_x).float()
        train_y = torch.from_numpy(train_y).float()
    for parameter in parameterstocheck:
        # if parameter != parameter2:
        test_samples = sampler(resolution,parameters,maxpGB,boundaries,p1, uniform= True, twoD = True, onlyparameter=parameter, secondparameter=parameter2)
        test_x[parameter+parameter2] = np.zeros((resolution,len(parameters)))
        i = 0
        for name in parametersfd:
            test_x[parameter+parameter2][:,i] = test_samples[name]
            i +=1
        test_y[parameter+parameter2] = test_samples['Likelihood']
        test_x[parameter+parameter2] = torch.from_numpy(test_x[parameter+parameter2]).float()
        test_y[parameter+parameter2] = torch.from_numpy(test_y[parameter+parameter2]).float()

    kernel = gpytorch.kernels.RBFKernel(ard_num_dims=2)
    
    nu = np.mean(train_y.numpy())
    sigma = np.std(train_y.numpy())
    train_y = (train_y-nu)/sigma
    model, likelihood = traingpmodel(train_x,train_y, kernel)
    model.eval()
    likelihood.eval()
    observed_pred = {}
    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        observed_pred[parameter+parameter2] = likelihood(model(test_x[parameter+parameter2]))
        observed_pred_print = (observed_pred[parameter+parameter2].mean*sigma)+nu
    print('sqrt(MSE) ',parameter+parameter2, np.sqrt(mean_squared_error(test_y[parameter+parameter2].numpy(),observed_pred_print.numpy())))

    flatsamples = np.zeros((len(parameterstocheck),resolution))
    flatsamplesparameters = []
    i = 0
    for parameter in parameterstocheck:
        if parameter != parameter2:
            flatsamples[i,:] = test_y[parameter+parameter2].numpy()
            flatsamplesparameters.append(test_x[parameter+parameter2].numpy())
            i+=1
    maxindx = np.unravel_index(flatsamples.argmax(), flatsamples.shape)
    max_parameters = flatsamplesparameters[maxindx[0]][maxindx[1]]
    max_loglike = flatsamples.max()
    return max_parameters, max_loglike


def scaletooriginal(previous_max,boundaries):
    i = 0
    maxpGB = deepcopy(pGB)
    for parameter in parameters:
        if parameter in ['FrequencyDerivative']:
            i -= 1
        elif parameter in ['EclipticLatitude']:
            maxpGB[parameter] = np.arcsin((previous_max[i]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
        elif parameter in ['Inclination']:
            maxpGB[parameter] = np.arccos((previous_max[i]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
        else:
            maxpGB[parameter] = (previous_max[i]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
        i += 1
    return maxpGB
def plotplanes(parameterstocheck,parameter2):
    fig, ax = plt.subplots(2, 4,figsize=(15,15))
    plt.suptitle("loglikelihood true")
    i = 0    
    for parameter in parameterstocheck:
        if parameter != parameter2:
            j = 0
            if i > 3:
                j = 1
            with torch.no_grad():
                ax[j,i%4].axvline(x=previous_max[parametersfd.index(parameter)], color='k')
                ax[j,i%4].axhline(y=previous_max[parametersfd.index(parameter2)], color='k')
                im = ax[j,i%4].scatter(test_x[parameter+parameter2].numpy()[:,parametersfd.index(parameter)],test_x[parameter+parameter2].numpy()[:,parametersfd.index(parameter2)],c=test_y[parameter+parameter2].numpy()[:])
                ax[j,i%4].set_xlabel(parameter)
                ax[j,i%4].set_ylabel(parameter2)
                fig.colorbar(im, ax=ax[j,i%4])
        i += 1
    plt.show()





train_x = train_x
train_y = train_y
def traingpmodel(train_x,train_y,kernel):
    bo_iterations = 2
    for bo_iter in range(bo_iterations):
        # initialize likelihood and model
        likelihood = gpytorch.likelihoods.GaussianLikelihood()
        model = ExactGPModel(train_x, train_y, likelihood, kernel)

        training_iter = 50

        hypers = {
            'likelihood.noise': torch.tensor(0.0001),
        }
        model.initialize(**hypers)

        # Find optimal model hyperparameters
        model.train()
        likelihood.train()

        # Use the adam optimizer
        optimizer = torch.optim.Adam([
            {"params": model.mean_module.parameters()},
            {"params": model.covar_module.parameters()},
        ], lr=0.1)  # Includes GaussianLikelihood parameters
        # optimizer = torch.optim.Adam(model.parameters(), lr=0.1)

        # "Loss" for GPs - the marginal log likelihood
        mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)

        for i in range(training_iter):
            torch.cuda.empty_cache()
            # Zero gradients from previous iteration
            optimizer.zero_grad()
            # Output from model
            output = model(train_x)
            # Calc loss and backprop gradients
            loss = -mll(output, train_y)
            loss.backward()
            print('Iter %d/%d - Loss: %.3f     noise: %.3f' % (
            i + 1, training_iter, loss.item(),
            model.likelihood.noise.item()
            ))
            optimizer.step()

        best_value = train_y.max()
        print('best value',bo_iter, (best_value*sigma)+nu, )
        best_value = best_value
        if bo_iter < bo_iterations-1:
            qEI = qExpectedImprovement(model=model, best_f=best_value)
            qUCB = qUpperConfidenceBound(model=model, beta=1)
            candidates = 200
            new_point_analytic, _ = optimize_acqf(
                acq_function=qUCB,
                bounds=torch.tensor([[0.0] * 8, [1.0] * 8]),
                q=candidates,
                num_restarts=1,
                raw_samples=100,
                options={},
            )
            train_x = torch.cat((train_x,new_point_analytic),0)
            new_point_analytic = new_point_analytic.cpu().numpy()
            for candi_n in range(candidates):
                para_n = 0
                for parameter in parameters:
                    if parameter in ['EclipticLatitude']:
                        pGBs[parameter] = np.arcsin((new_point_analytic[candi_n][para_n]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
                    elif parameter in ['Inclination']:
                        pGBs[parameter] = np.arccos((new_point_analytic[candi_n][para_n]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
                    else:
                        pGBs[parameter] = (new_point_analytic[candi_n][para_n]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
                    para_n += 1
                loglike = loglikelihood(pGBs)
                # print('new point',i, new_point_analytic[candi_n], loglike)

                loglike = (loglike-nu)/sigma
                loglike = torch.tensor([loglike]).float()
                train_y = torch.cat((train_y,loglike),0)
    return model, likelihood

maxpGB = deepcopy(pGB)
for i in range(40):
    resolution = 21**2
    parameter2 = 'Polarization'
    if i > 3:
        parameter2 = 'EclipticLatitude'
        resolution = 19**2
    if i > 7:
        parameter2 = 'EclipticLatitude'
        resolution = 23**2
    if i > 11:
        parameter2 = 'Frequency'
        resolution = 30**2
    parameter1 = parametersfd[np.random.randint(0,6)]
    parameter2 = parametersfd[np.random.randint(0,6)]
    parameter2 = 'InitialPhase'
    while parameter2 == parameter1:
        parameter2 = parametersfd[np.random.randint(0,6)]
    parametersreduced = [parameter1]
    previous_max, max_loglike = planemaxsearch(maxpGB,parametersreduced,parameter2,resolution)
    maxpGB = scaletooriginal(previous_max,boundaries)
    # plotplanes(parametersreduced,parameter2)
    
    print(i,max_loglike,maxpGB)
    print(previous_max)


# Get into evaluation (predictive posterior) mode
model.eval()
likelihood.eval()
observed_pred = {}


for parameter in parametersfd:
    # Make predictions by feeding model through likelihood
    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        observed_pred[parameter] = likelihood(model(test_x[parameter]))
        observed_pred_print = (observed_pred[parameter]*sigma)+nu
    print('sqrt(MSE) ',parameter, np.sqrt(mean_squared_error(test_y[parameter].numpy(),observed_pred_print.mean.cpu().numpy())))
    if parameter != parameter2:
        with torch.no_grad(), gpytorch.settings.fast_pred_var():
            observed_pred[parameter+parameter2] = likelihood(model(test_x[parameter+parameter2]))
            observed_pred_print = (observed_pred[parameter+parameter2]*sigma)+nu
        print('sqrt(MSE) ',parameter+parameter2, np.sqrt(mean_squared_error(test_y[parameter+parameter2].numpy(),observed_pred_print.mean.cpu().numpy())))

parameter = 'random'
with torch.no_grad(), gpytorch.settings.fast_pred_var():
    observed_pred[parameter] = likelihood(model(test_x[parameter]))
observed_pred_print = (observed_pred[parameter]*sigma)+nu
print('sqrt(MSE) ','random', np.sqrt(mean_squared_error(test_y[parameter].numpy(),observed_pred_print.mean.cpu().numpy())))


for parameter in parameters:
    if parameter in ['EclipticLatitude']:
        samples[parameter] = np.arcsin((samples[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
        test_samples[parameter] = np.arcsin((test_samples[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
    elif parameter in ['Inclination']:
        samples[parameter] = np.arccos((samples[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
        test_samples[parameter] = np.arccos((test_samples[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
    else:
        samples[parameter] = (samples[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
        test_samples[parameter] = (test_samples[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]


# plt.plot(Xs.f, Xs, label='optimized')
ax6.set_xlabel('f [mHz]')
ax6.set_ylabel('X-TDI real [1/Hz]')
ax1.legend()
# ax1.plot(Ys.f, Ys, label='optimized')
ax2.set_xlabel('f [mHz]')
ax2.set_ylabel('Y-TDI real [1/Hz]')
ax2.legend()
# ax1.plot(Zs.f, Zs, label='optimized')
ax3.set_xlabel('f [mHz]')
ax3.set_ylabel('Z-TDI real [1/Hz]')
ax3.legend()
# ax1.plot(Xs.f, Xs.imag, label='optimized')
ax4.set_xlabel('f [mHz]')
ax4.set_label('X-TDI imag [1/Hz]')
ax4.legend()
# ax1.plot(Ys.f, Ys.imag, label='optimized')
ax5.set_xlabel('f [mHz]')
ax5.set_ylabel('Y-TDI imag [1/Hz]')
ax5.legend()
# plt.plot(Zs.f, Zs.imag, label='optimized')
ax6.set_xlabel('f [mHz]')
ax6.set_ylabel('Z-TDI imag [1/Hz]')
ax6.legend()

train_x = train_x.cpu()
train_y = train_y.cpu()
# train_y = (train_y*sigma)+nu

pGB01 = {}
for parameter in parameters:
    if parameter in ['EclipticLatitude']:
        pGB01[parameter] = ((np.sin(pGB[parameter])-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0]))
    elif parameter in ['Inclination']:
        pGB01[parameter] = ((np.cos(pGB[parameter])-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0]))
    else:
        pGB01[parameter] = ((pGB[parameter]-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0]))


fig, ax = plt.subplots(2, 4,figsize=(15,15))
plt.suptitle("loglikelihood")
i = 0    
mean = {}
for parameter in parametersfd:
    j = 0
    if i > 3:
        j = 1
    with torch.no_grad():
        # Get upper and lower confidence bounds
        lower, upper = observed_pred[parameter].confidence_region()
        mean2 = observed_pred[parameter].mean
        mean2 = mean2.cpu()
        lower = lower.cpu()
        upper = upper.cpu()
        test_x[parameter] = test_x[parameter].cpu()
        # mean[parameter] = (mean2*sigma)+nu
        mean[parameter] = mean2
        # lower = (lower*sigma)+nu
        # upper = (upper*sigma)+nu
        # Plot training data as black stars
        ax[j,i%4].axvline(x=pGB01[parameter], color='k',label='True')
        # ax[j,i%4].plot(train_x.numpy()[:,i], train_y.numpy(), 'k*')
        ax[j,i%4].plot(test_x[parameter].numpy()[1:,i], test_y[parameter].numpy()[1:], 'g.',label='True')
        ax[j,i%4].plot(train_x.numpy()[:,i], train_y.numpy()[:], 'r.',label='Train')
        # Plot predictive means as blue line
        ax[j,i%4].plot(test_x[parameter].numpy()[1:,i], mean[parameter].numpy()[1:], 'b.',label='Mean')
        # Shade between the lower and upper confidence bounds
        ax[j,i%4].fill_between(test_x[parameter].numpy()[1:,i], lower.numpy()[1:], upper.numpy()[1:], alpha=0.5,label='Confidence')
        ax[j,i%4].legend()
        ax[j,i%4].set_xlabel(parameter)
    i += 1
fig, ax = plt.subplots(2, 4,figsize=(15,15))
plt.suptitle("loglikelihood true")
i = 0    
for parameter in parametersfd:
    if parameter != parameter2:
        j = 0
        if i > 3:
            j = 1
        with torch.no_grad():
            ax[j,i%4].axvline(x=pGB01[parameter], color='k')
            ax[j,i%4].axhline(y=pGB01[parameter2], color='k')
            test_x[parameter+parameter2] = test_x[parameter+parameter2].cpu()
            im = ax[j,i%4].scatter(test_x[parameter+parameter2].numpy()[:,i],test_x[parameter+parameter2].numpy()[:,parametersfd.index(parameter2)],c=test_y[parameter+parameter2].numpy()[:])
            ax[j,i%4].set_xlabel(parameter)
            ax[j,i%4].set_ylabel(parameter2)
            fig.colorbar(im, ax=ax[j,i%4])
    else:
        ax[j,i%4].plot(test_x[parameter].numpy()[1:,i], test_y[parameter].numpy()[1:], 'g.')
        ax[j,i%4].set_xlabel(parameter)
        ax[j,i%4].set_ylabel('loglikelihood')
        ax[j,i%4].legend(['True'])
    i += 1
    
# fig, ax = plt.subplots(2, 4,figsize=(15,15))
# plt.suptitle("loglikelihood predicted mean")
# i = 0    
# for parameter in parametersfd:
#     if parameter != parameter2:
#         j = 0
#         if i > 3:
#             j = 1
#         with torch.no_grad():
#             # Get upper and lower confidence bounds
#             mean2 = observed_pred[parameter+parameter2].mean
#             mean2 = mean2.cpu()
#             test_x[parameter+parameter2] = test_x[parameter+parameter2].cpu()
#             mean[parameter+parameter2] = (mean2*sigma)+nu
#             ax[j,i%4].axvline(x=pGB01[parameter], color='k')
#             ax[j,i%4].axhline(y=pGB01[parameter2], color='k')
#             im = ax[j,i%4].scatter(test_x[parameter+parameter2].numpy()[:,i],test_x[parameter+parameter2].numpy()[:,parametersfd.index(parameter2)],c=mean[parameter+parameter2][:])
#             ax[j,i%4].set_xlabel(parameter)
#             ax[j,i%4].set_ylabel(parameter2)
#             fig.colorbar(im, ax=ax[j,i%4])
#     else:
#         with torch.no_grad():
#             lower, upper = observed_pred[parameter].confidence_region()
#             lower = lower.cpu()
#             upper = upper.cpu()
#             lower = (lower*sigma)+nu
#             upper = (upper*sigma)+nu
#             ax[j,i%4].plot(test_x[parameter].numpy()[1:,i], test_y[parameter].numpy()[1:], 'g.')
#             ax[j,i%4].plot(test_x[parameter].numpy()[1:,i], mean[parameter].numpy()[1:], 'b.')
#             ax[j,i%4].fill_between(test_x[parameter].numpy()[1:,i], lower.numpy()[1:], upper.numpy()[1:], alpha=0.5)
#             ax[j,i%4].set_xlabel(parameter)
#             ax[j,i%4].set_ylabel('loglikelihood')
#             ax[j,i%4].legend(['True','Mean', 'Confidence'])
#     i += 1
plt.show()

prediction = observed_pred['Frequency'].mean.cpu().numpy()
n_bin = 50
i = 0    
fig, axes = plt.subplots(2, 4,figsize=(15,15))
plt.suptitle("sampled posterior")
for parameter in parameters:
    j = 0
    if i > 3:
        j = 1
    axes[j,i%4].axvline(x=pGB[parameter], color='r')
    axes[j,i%4].plot(samples[parameter],samples['Likelihood'], '.',label='train')
    axes[j,i%4].plot(test_samples[parameter],test_samples['Likelihood'],'.',label='true')
    axes[j,i%4].errorbar(test_samples[parameter][1:],prediction[1:],prediction[1:]-lower.numpy()[1:], fmt='.',label='prediction')
    # plt.fill_between(test_samples[parameter][1:],prediction[1:]-p_std[1:],prediction[1:]+p_std[1:], color='g', alpha= 0.4)
    if i == 0:
        axes[j,i%4].set_ylabel('log-Likelihood')
    # plt.subplot(3,number_of_parameters,i+number_of_parameters)
    # plt.axvline(x=pGB[parameter], color='r')
    # n, bins, patches = plt.hist(samples[parameter], n_bin, density=True, facecolor='k', alpha=0.5)

    # plt.subplot(3,number_of_parameters,i+2*number_of_parameters)
    # if parameter == 'Frequency':
    #     plt.axvline(x=pGB[parameter]*1000, color='r')
    #     plt.plot(samples[parameter]*1000, range(number_of_samples), 'k')
    # else:
    #     plt.axvline(x=pGB[parameter], color='r')
    #     plt.plot(samples[parameter], range(number_of_samples), 'k')
    axes[j,i%4].set_xlabel(parameter)
    i += 1
plt.legend()

name = 'Frequency'
plt.figure(figsize=(10,8))
plt.suptitle("sampled posterior")
plt.subplot(1,1,1)
plt.axvline(x=pGB[name], color='r')
plt.plot(samples[name],samples['Likelihood'], '.',label='train',zorder=1)
plt.plot(test_samples[name],test_samples['Likelihood'],'.',label='true')
plt.plot(test_samples[name][1:],prediction[1:],label='prediction')
plt.fill_between(test_samples[name][1:],lower.numpy()[1:],upper.numpy()[1:], color='g', alpha= 0.4)
plt.ylabel('log-Likelihood')
plt.xlabel(name)
plt.legend()
# plt.figure()
# plt.scatter(samples['Amplitude'],samples['Frequency'],c=samples['Likelihood'])
# plt.xlim(min(samples['Amplitude']),max(samples['Amplitude']))
# plt.ylim(min(samples['Frequency']),max(samples['Frequency']))
# plt.colorbar()
# plt.figure()
# # plt.title('sqrt(MSE)',np.round(np.sqrt(mean_squared_error(test_y,prediction)),2))
# plt.scatter(samples['EclipticLatitude'],samples['EclipticLongitude'],c=samples['Likelihood'])
# plt.xlim(min(samples['EclipticLatitude']),max(samples['EclipticLatitude']))
# plt.ylim(min(samples['EclipticLongitude']),max(samples['EclipticLongitude']))
# plt.colorbar()
plt.show()