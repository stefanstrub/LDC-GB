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
Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB, oversample=4, simulator='synthlisa')
fmin, fmax = float(Xs.f[0]) , float(Xs.f[-1]+Xs.attrs['df'])
source = dict({"X":Xs, "Y":Ys, "Z":Zs})
start = time.time()
Xs_td, Ys_td, Zs_td = GB.get_td_tdixyz(template=pGB, simulator='synthlisa')

print('ftransform',time.time()- start)
print(len(Xs))

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
    p = -float(np.sum(diff / Sn)*Xs.attrs['df'])/2.0
    return p

# Number of histogram bins.
n_bin = 50
# Total number of proposed samples.
number_of_samples = 1*10 **3
cutoff_ratio = 1000

parameters = ['Amplitude','EclipticLatitude','EclipticLongitude','Frequency','FrequencyDerivative','Inclination','InitialPhase','Polarization']
boundaries = {'Amplitude': [10**-22.0, 5*10**-21.0],'EclipticLatitude': [-1.0, 1.0],
'EclipticLongitude': [0.0, 2.0*np.pi],'Frequency': [0.0004725, 0.0004727],'FrequencyDerivative': [10**-20.0, 10**-18.0],
'Inclination': [-1.0, 1.0],'InitialPhase': [0.0, 2.0*np.pi],'Polarization': [0.0, 2.0*np.pi]}

# boundaries_small = deepcopy(boundaries)
# part_ratio = 5
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
Xs, Ys, Zs = Xs[highSNR], Ys[highSNR], Zs[highSNR]
fmin, fmax = float(Xs.f[0]) , float(Xs.f[-1]+Xs.attrs['df'])
freq = np.array(Xs.sel(f=slice(fmin, fmax)).f)
Sn = Nmodel.psd(freq=freq, option='X')
diff = np.abs(dataX - Xs.values)**2 + np.abs(dataY - Ys.values)**2 + np.abs(dataZ - Zs.values)**2
p1 = -float(np.sum(diff / Sn)*Xs.attrs['df'])/2.0
p1 = p1

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
def sampler(number_of_samples,parameters,pGB,boundaries,p1, uniform=False, MCMC=False, only=False, onlyparameter='Frequency'):
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
    
    for i in range(1, number_of_samples):
        if only:
            parameter = onlyparameter
            if uniform:
                pGBs01[parameter] = i/number_of_samples
            pGBs[parameter] = (pGBs01[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
            # pGBs[parameter] = np.arccos((pGBs01[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
        else:
            for parameter in parameters:
                if parameter in ['Frequency']:#,'FrequencyDerivative','Amplitude','EclipticLongitude']:
                    pGBs01[parameter] = np.random.rand()
                elif parameter in ['FrequencyDerivative']:
                    pass
                elif parameter in ['EclipticLatitude']:
                    pGBs01[parameter] = np.random.rand()
                elif parameter in ['Inclination']:
                    pGBs01[parameter] = np.random.rand()
                else:
                    pGBs01[parameter] = np.random.rand()
        for parameter in parameters:
            if parameter in ['FrequencyDerivative']:
                pass
            if parameter in ['EclipticLatitude']:
                pGBs[parameter] = np.arcsin((pGBs01[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
            elif parameter in ['Inclination']:
                pGBs[parameter] = np.arccos((pGBs01[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
            else:
                pGBs[parameter] = (pGBs01[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
        p_test = loglikelihood(pGBs)
        Tinv = 1
        if i > number_of_samples/5:
            Tinv =  i/10
        if not(MCMC):
            Tinv = 0 
        # print(p_test/p1, pGBs['Frequency'])
        # Apply Metropolis rule.
        # print(p_test.values,p.values,(p_test / p) ** T_inv)
        if (p_test / p1) ** (-Tinv) > np.random.rand():  # L^i/L^j
        # elif True:  # L^i/L^j
            p1 = p_test
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
# plot
n_max = np.argmax(samples['Likelihood'][1:])+1
pGBmax = deepcopy(pGB)
for name, titles in pGBmax.items():
    if name != 'Name':
        pGBmax[name] = deepcopy(samples[name][n_max])
Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGBmax, oversample=4, simulator='synthlisa')


class ExactGPModel(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood):
        super(ExactGPModel, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ConstantMean()
        kernelA = gpytorch.kernels.RBFKernel(active_dims=(0))
        kernelLat = gpytorch.kernels.PeriodicKernel(active_dims=torch.tensor([1]))
        kernelLong = gpytorch.kernels.PeriodicKernel(active_dims=torch.tensor([2]))
        kernelF = gpytorch.kernels.RBFKernel(active_dims=(3))
        kernelFD = gpytorch.kernels.RBFKernel(active_dims=(4))
        kernelI = gpytorch.kernels.PeriodicKernel(active_dims=(5))
        kernelIP = gpytorch.kernels.PeriodicKernel(active_dims=torch.tensor([6]))
        kernelP = gpytorch.kernels.PeriodicKernel(active_dims=torch.tensor([7]))
        kernel = kernelA + kernelLat + kernelLong + kernelF + kernelFD + kernelI + kernelIP + kernelP
        kernel = kernelA * kernelLat * kernelLong * kernelF * kernelFD * kernelI * kernelIP * kernelP
        self.covar_module = gpytorch.kernels.ScaleKernel(kernel)

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)


train_x = np.zeros((len(samples['Amplitude'])-1,len(parameters)))
i = 0
for name in parameters:
    train_x[:,i] = samples[name][1:]
    i +=1
nu = np.mean(samples['Likelihood'][1:])
sigma = np.std(samples['Likelihood'][1:])
train_y = (samples['Likelihood'][1:]-nu)/sigma

number_of_test_samples = 300
test_x = {}
test_y = {}
parameter = 'Frequency'
test_samples = sampler(number_of_test_samples,parameters,pGB,boundaries,p1, uniform= True, only=True, onlyparameter=parameter)
test_x[parameter] = np.zeros((number_of_test_samples,len(parameters)))
i = 0
for name in parameters:
    test_x[parameter][:,i] = test_samples[name]
    i +=1
test_y[parameter] = test_samples['Likelihood']
test_x[parameter] = torch.from_numpy(test_x[parameter]).float()
test_x[parameter] = test_x[parameter]
test_y[parameter] = torch.from_numpy(test_y[parameter]).float()

train_x = torch.from_numpy(train_x).float()
train_y = torch.from_numpy(train_y).float()
# initialize likelihood and model
likelihood = gpytorch.likelihoods.GaussianLikelihood()
model = ExactGPModel(train_x, train_y, likelihood)

training_iter = 50


# Find optimal model hyperparameters
model.train()
likelihood.train()


hypers = {
    'likelihood.noise': torch.tensor(0.0001),
    # 'covar_module.base_kernel.lengthscale': torch.tensor(0.08)
    # 'covar_module.base_kernel.kernels.1.period_length': torch.tensor(0.5)
    # list(covar_module.base_kernel.kernels)[1].period_length
}
model.initialize(**hypers)
# Polarization
model.covar_module.base_kernel.kernels[1].period_length = torch.tensor([[0.5]])
list(model.covar_module.base_kernel.kernels[1].parameters())[1].requires_grad=False
# InitialPhase
model.covar_module.base_kernel.kernels[0].kernels[1].period_length = torch.tensor([[1.0]])
list(model.covar_module.base_kernel.kernels[0].kernels[1].parameters())[1].requires_grad=False
# Frequency
model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[1].lengthscale = torch.tensor([[0.066]])
list(model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[1].parameters())[0].requires_grad=False
# Longitude
model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[1].period_length = torch.tensor([[1.0]])
list(model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[1].parameters())[1].requires_grad=False
# Latitude
model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[1].period_length = torch.tensor([[2.0]])
list(model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[1].parameters())[1].requires_grad=False
# Amplitude
model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].lengthscale = torch.tensor([[0.4]])
list(model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[1].parameters())[0].requires_grad=False


# Use the adam optimizer
optimizer = torch.optim.Adam([
    {"params": model.mean_module.parameters()},
    {"params": model.covar_module.parameters()},
], lr=0.1)  # Includes GaussianLikelihood parameters
# optimizer = torch.optim.Adam(model.parameters(), lr=0.1)

# "Loss" for GPs - the marginal log likelihood
mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)

for i in range(training_iter):
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
    print('Iter %d/%d - Loss: %.3f   Al: %.3f Lal: %.3f Lap: %.3f  Lol: %.3f Lop: %.3f  Fl: %.3f  FDl: %.3f  Il: %.3f  IPl: %.3f IPp: %.3f  Pl: %.3f Pp: %.3f   noise: %.3f' % (
        i + 1, training_iter, loss.item(),
        model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].lengthscale.item(),
        model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[1].lengthscale.item(),
        model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[1].period_length.item(),
        model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[1].lengthscale.item(),
        model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[1].period_length.item(),
        model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[1].lengthscale.item(),
        model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[1].lengthscale.item(),
        model.covar_module.base_kernel.kernels[0].kernels[0].kernels[1].lengthscale.item(),
        model.covar_module.base_kernel.kernels[0].kernels[1].lengthscale.item(),
        model.covar_module.base_kernel.kernels[0].kernels[1].period_length.item(),
        model.covar_module.base_kernel.kernels[1].lengthscale.item(),
        model.covar_module.base_kernel.kernels[1].period_length.item(),
        model.likelihood.noise.item()
    ))
    optimizer.step()

# Get into evaluation (predictive posterior) mode
model.eval()
likelihood.eval()

observed_pred = {}
parameter = 'Frequency'
# Make predictions by feeding model through likelihood
with torch.no_grad(), gpytorch.settings.fast_pred_var():
    observed_pred[parameter] = likelihood(model(test_x[parameter]))
print('sqrt(MSE) ',parameter, np.sqrt(mean_squared_error(test_y[parameter].cpu().numpy(),observed_pred[parameter].mean.cpu().numpy())))
parameter = 'random'
with torch.no_grad(), gpytorch.settings.fast_pred_var():
    observed_pred[parameter] = likelihood(model(test_x[parameter]))
print('sqrt(MSE) ','random', np.sqrt(mean_squared_error(test_y[parameter].cpu().numpy(),observed_pred[parameter].mean.cpu().numpy())))

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


with torch.no_grad():
    # Initialize plot
    f, ax = plt.subplots(1, 1, figsize=(4, 3))

    # Get upper and lower confidence bounds
    lower, upper = observed_pred.confidence_region()
    # Plot training data as black stars
    ax.plot(train_x.numpy()[:,3], train_y.numpy(), 'k*')
    ax.plot(test_x.numpy()[1:,3], test_y.numpy()[1:], 'g.')
    # Plot predictive means as blue line
    ax.plot(test_x.numpy()[1:,3], observed_pred.mean.numpy()[1:], 'b.')
    # Shade between the lower and upper confidence bounds
    ax.fill_between(test_x.numpy()[1:,3], lower.numpy()[1:], upper.numpy()[1:], alpha=0.5)
    # ax.set_ylim([-3, 3])
    ax.legend(['Observed Data','True', 'Mean', 'Confidence'])
    

prediction = observed_pred.mean.numpy()
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
plt.plot(samples[name],samples['Likelihood'], '.',label='train',zorder=3)
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
plt.figure()
# plt.title('sqrt(MSE)',np.round(np.sqrt(mean_squared_error(test_y,prediction)),2))
plt.scatter(samples['EclipticLatitude'],samples['EclipticLongitude'],c=samples['Likelihood'])
plt.xlim(min(samples['EclipticLatitude']),max(samples['EclipticLatitude']))
plt.ylim(min(samples['EclipticLongitude']),max(samples['EclipticLongitude']))
plt.colorbar()
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
plt.plot(samples[:,1],np.log(likelihood_values), '.')
plt.xlabel('Frequeny')
plt.ylabel('log-Likelihood')
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

## Constant for Cost function
THRESHOLD = 0.5
W1 = 1
W2 = 20
W3 = 100
W4 = 0.04


def cost_function(true, predicted):
    """
        true: true values in 1D numpy array
        predicted: predicted values in 1D numpy array

        return: float
    """
    cost = (true - predicted)**2

    # true above threshold (case 1)
    mask = true > THRESHOLD
    mask_w1 = np.logical_and(predicted>=true,mask)
    mask_w2 = np.logical_and(np.logical_and(predicted<true,predicted >=THRESHOLD),mask)
    mask_w3 = np.logical_and(predicted<THRESHOLD,mask)

    cost[mask_w1] = cost[mask_w1]*W1
    cost[mask_w2] = cost[mask_w2]*W2
    cost[mask_w3] = cost[mask_w3]*W3

    # true value below threshold (case 2)
    mask = true <= THRESHOLD
    mask_w1 = np.logical_and(predicted>true,mask)
    mask_w2 = np.logical_and(predicted<=true,mask)

    cost[mask_w1] = cost[mask_w1]*W1
    cost[mask_w2] = cost[mask_w2]*W2

    reward = W4*np.logical_and(predicted < THRESHOLD,true<THRESHOLD)
    if reward is None:
        reward = 0
    return np.mean(cost) - np.mean(reward)

"""
Fill in the methods of the Model. Please do not change the given methods for the checker script to work.
You can add new methods, and make changes. The checker script performs:


    M = Model()
    M.fit_model(train_x,train_y)
    prediction = M.predict(test_x)

It uses predictions to compare to the ground truth using the cost_function above.
"""


if __name__ == "__main__":
    main()
