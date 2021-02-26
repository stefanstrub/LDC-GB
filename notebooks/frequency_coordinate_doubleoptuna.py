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

import optuna


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

noise_model = "MRDv1"
Nmodel = get_noise_model(noise_model, np.logspace(-5, -1, 100))
Npsd = Nmodel.psd()

vgb, units = hdfio.load_array(sangria_fn_training, name="sky/vgb/cat")
GB = fastGB.FastGB(delta_t=dt, T=float(tdi_ts["X"].t[-1])) # in seconds
pGB = dict(zip(vgb.dtype.names, vgb[8])) # we take the source #8

Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB, oversample=4, simulator='synthlisa')
fmin, fmax = float(Xs.f[0]) , float(Xs.f[-1]+Xs.attrs['df'])
source = dict({"X":Xs, "Y":Ys, "Z":Zs})
start = time.time()
Xs_td, Ys_td, Zs_td = GB.get_td_tdixyz(template=pGB, simulator='synthlisa')


f_noise = np.logspace(-5, -1, 100)
Nmodel = get_noise_model(noise_model, f_noise)
freq = np.array(source["X"].sel(f=slice(fmin, fmax)).f)

Sn = Nmodel.psd(freq=freq, option='X')

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
    p1 = float(np.sum(diff / (Sn+noise))/len(diff))/2.0
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
'EclipticLongitude': [0.0, 2.0*np.pi],'Frequency': [pGB['Frequency']-0.0000001,pGB['Frequency']+0.0000001],'FrequencyDerivative': [10**-20.0, 10**-16.0],
'Inclination': [-1.0, 1.0],'InitialPhase': [0.0, 2.0*np.pi],'Polarization': [0.0, 1.0*np.pi]}

previous_max = [0.2090, 0.1000, 0.8469, 0.5276, 0.7168, 0.9667, 0.0970, 0.0000]
# previous_max = [0.2090, 0.2667, 0.7333, 0.5276, 0.7168, 0.9667, 0.0970, 0.0000]
# previous_max = [0.2090, 0.3333, 0.7333, 0.5276, 0.7168, 0.8667, 0.0970, 0.0000]
# previous_max = [0.2090, 0.3667, 0.6667, 0.5276, 0.7168, 0.8667, 0.0970, 0.0000]
# previous_max = [0.4000, 0.3667, 0.6667, 0.5276, 0.7168, 0.8667, 0.0970, 0.0000]
# previous_max = [0.4000, 0.3667, 0.6667, 0.5276, 0.7667, 0.8667, 0.0970, 0.0000]
# previous_max = [0.2333, 0.3667, 0.6667, 0.5276, 0.9667, 0.8667, 0.0970, 0.0000]
# previous_max = [0.333,0.54,0.134,0.7345,0.6456,0.2645,0.8216,0.000]
# previous_max = [0.2090, 0.1000, 0.8469, 0.3276, 0.7168, 0.9667, 0.0970, 0.0000]
previous_max = [0.45090, 0.5600, 0.123469, 0.87276, 0.2341168, 0.56667, 0.5689970, 0.0000]
previous_max = [0.45090, 0.5600, 0.123469, 0.3276, 0.2341168, 0.56667, 0.9689970, 0.0000]
# previous_max = [0.86436456, 0.3156825,  0.6350386,  0.55715334, 0.5604474,  0.7789818, 0.03608589, 0.0]
# previous_max =[0.80888367, 0.35581076, 0.62365836, 0.5551591,  0.5607991,  0.76172084,0.03608589, 0.0]
# previous_max = np.random.rand(8)
i = 0
pGBs = deepcopy(pGB)
for parameter in parameters:
    if parameter in ['FrequencyDerivative']:
        i -= 1
    elif parameter in ['EclipticLatitude']:
        pGBs[parameter] = np.arcsin((previous_max[i]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
    elif parameter in ['Inclination']:
        pGBs[parameter] = np.arccos((previous_max[i]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
    else:
        pGBs[parameter] = (previous_max[i]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
        print(parameter, pGBs[parameter],previous_max[i])
    i += 1

# print(pGB)
# pGB['EclipticLatitude'] = 0.0
# pGB = {'Amplitude': 4.150747034e-21, 'EclipticLatitude': 0.0, 'EclipticLongitude': 1.6707963267948967, 'Frequency': 0.000472612, 'FrequencyDerivative': 1.440030107143638e-19, 'Inclination': 1.4468779499032993, 'InitialPhase': 1.7, 'Name': 'CD-30o11223', 'Polarization': 1.7}
# print(pGB)
# pGB = {'Amplitude': 1.1240340491284739e-21, 'EclipticLatitude': 0.48579612338022954, 'EclipticLongitude': 5.3215032928260895, 'Frequency': 0.0004726055124556648, 'FrequencyDerivative': 7.422865032688086e-19, 'Inclination': 1.12226944780215, 'InitialPhase': 2.391197244408699, 'Name': 'CD-30o11223', 'Polarization': 0.6092901760377935}
# Make the first random sample. ------------------------------------
# pGBs = deepcopy(pGB)
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
lowerindex = np.where(highSNR)[0][0]-10
higherindex = np.where(highSNR)[0][-1]+10
dataX = tdi_fs["X"].isel(f=slice(Xs.kmin, Xs.kmin+len(Xs)))[lowerindex:higherindex]
dataY = tdi_fs["Y"].isel(f=slice(Ys.kmin, Ys.kmin+len(Ys)))[lowerindex:higherindex]
dataZ = tdi_fs["Z"].isel(f=slice(Zs.kmin, Zs.kmin+len(Zs)))[lowerindex:higherindex]
spd_data = np.abs(dataX)**2 + np.abs(dataY)**2 + np.abs(dataZ)**2
noise = (np.mean(spd_data[:5])+np.mean(spd_data[-5:])).values/2
Xs, Ys, Zs = Xs[lowerindex:higherindex], Ys[lowerindex:higherindex], Zs[lowerindex:higherindex]
fmin, fmax = float(Xs.f[0]) , float(Xs.f[-1]+Xs.attrs['df'])
freq = np.array(Xs.sel(f=slice(fmin, fmax)).f)
Sn = Nmodel.psd(freq=freq, option='X')
diff = np.abs(dataX - Xs.values)**2 + np.abs(dataY - Ys.values)**2 + np.abs(dataZ - Zs.values)**2
p1 = float(np.sum(diff / (Sn+noise))/len(diff))/2.0
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

print('p start',p1, 'p true', loglikelihood(pGB))

def sampler(number_of_samples,parameters,pGB,boundaries,p1, uniform=False, MCMC=False, only=False, onlyparameter='Frequency', twoD=False, secondparameter='Amplitude', calculate_loglikelihood=True):
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
    p1 = loglikelihood(pGB)
    samples['Likelihood'].append(p1)
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
                pGBs01[parameter] = ((i-1)%number_of_sampels_sqrt+1)/(number_of_sampels_sqrt+2)
                pGBs01[parameter2] = (j+1)/(number_of_sampels_sqrt+2)
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
        # pGBs01array = np.zeros(len(parametersfd))
        # i = 0
        # for name in parametersfd:
        #     pGBs01array[i] = pGBs01[name]
        #     i +=1
        # pGBscaled = scaletooriginal(pGBs01array ,boundaries)
        # print(loglikelihood(pGBs), loglikelihood(pGBscaled))
        if calculate_loglikelihood:
            p1 = loglikelihood(pGBs)
            samples['Likelihood'].append(p1)
        for parameter in parameters:
            samples[parameter].append(pGBs01[parameter])
    for parameter in parameters:
        samples[parameter] = np.asarray(samples[parameter])
    samples['Likelihood'] = np.asarray(samples['Likelihood'])
    return samples


#%%

class ExactGPModel(gpytorch.models.ExactGP, GPyTorchModel):
    _num_outputs = 1  # to inform GPyTorchModel API
    def __init__(self, train_x, train_y, likelihood, kernel):
        super(ExactGPModel, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ConstantMean()
        kernelA = gpytorch.kernels.RBFKernel(active_dims=(0),lengthscale_constraint=gpytorch.constraints.GreaterThan(0.2))
        kernelLat = gpytorch.kernels.PeriodicKernel(active_dims=torch.tensor([1]),lengthscale_constraint=gpytorch.constraints.GreaterThan(0.8))
        kernelLong = gpytorch.kernels.PeriodicKernel(active_dims=torch.tensor([2]),lengthscale_constraint=gpytorch.constraints.GreaterThan(0.8))
        kernelF = gpytorch.kernels.RBFKernel(active_dims=(3),lengthscale_constraint=gpytorch.constraints.Interval(0.06,0.1))
        # kernelFD = gpytorch.kernels.RBFKernel(active_dims=(4))
        kernelI = gpytorch.kernels.PolynomialKernel(power= 2,active_dims=(4))
        # kernelIP = gpytorch.kernels.PeriodicKernel(active_dims=torch.tensor([6]))
        kernelIP = gpytorch.kernels.ScaleKernel(gpytorch.kernels.CosineKernel(active_dims=torch.tensor([5]),period_length_constraint=gpytorch.constraints.Interval(0.499,0.501)), outputscale_constraint=gpytorch.constraints.Interval(0.9,1.1))
        # kernelP = gpytorch.kernels.PeriodicKernel(active_dims=torch.tensor([7]))
        # kernelP = gpytorch.kernels.CosineKernel(active_dims=torch.tensor([6]))*gpytorch.kernels.CosineKernel(active_dims=torch.tensor([6]))
        kernelP = gpytorch.kernels.ScaleKernel(gpytorch.kernels.CosineKernel(active_dims=torch.tensor([6]),period_length_constraint=gpytorch.constraints.Interval(0.249,0.251)), outputscale_constraint=gpytorch.constraints.Interval(0.9,1.1))
        # kernel = kernelA + kernelLat + kernelLong + kernelF + kernelFD + kernelI + kernelIP + kernelP
        kernel = kernelA * kernelLat * kernelLong * kernelF * kernelI * kernelIP * kernelP
        # self.covar_module = gpytorch.kernels.ScaleKernel(gpytorch.kernels.RBFKernel(ard_num_dims=8))
        self.covar_module = gpytorch.kernels.ScaleKernel(gpytorch.kernels.RBFKernel(ard_num_dims=8))
        self.to(train_x)  # make sure we're on the right device/dtype
  
    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)


number_of_test_samples = 20
test_x = {}
test_y = {}
# for parameter in parametersfd:
#     test_samples = sampler(number_of_test_samples,parameters,pGBs,boundaries,p1, uniform= True, only=True, onlyparameter=parameter)
#     test_x[parameter] = np.zeros((number_of_test_samples,len(parameters)))
#     i = 0
#     for name in parametersfd:
#         test_x[parameter][:,i] = test_samples[name]
#         i +=1
#     test_y[parameter] = test_samples['Likelihood']
#     test_x[parameter] = torch.from_numpy(test_x[parameter]).float()
#     test_x[parameter] = test_x[parameter]
#     test_y[parameter] = torch.from_numpy(test_y[parameter]).float()

def planeAdam(minpGB,parameterstocheck, parameter2, resolution, boundaries):
    for parameter in parameterstocheck:
        # if parameter != parameter2:
        resolution_reduced = int(20**2)
        resolution_reduced = resolution
        # if 'Frequency' in [parameter,parameter2]:
        #     resolution_reduced = int(15**2)
        train_samples = sampler(resolution_reduced,parameters,minpGB,boundaries,p1, uniform= False, twoD = True, onlyparameter=parameter, secondparameter=parameter2)
        train_x = np.zeros((resolution_reduced,len(parameters)))
        i = 0
        for name in parametersfd:
            train_x[:,i] = train_samples[name]
            i +=1
        train_y = train_samples['Likelihood']
        train_x = torch.from_numpy(train_x).float()
        train_y = torch.from_numpy(train_y).float()
    for parameter in parameterstocheck:
        # if parameter != parameter2:
        test_samples = sampler(resolution,parameters,minpGB,boundaries,p1, uniform= True, twoD = True, onlyparameter=parameter, secondparameter=parameter2, calculate_loglikelihood=False)
        test_x[parameter+parameter2] = np.zeros((resolution,len(parameters)))
        i = 0
        for name in parametersfd:
            test_x[parameter+parameter2][:,i] = test_samples[name]
            i +=1
        # test_y[parameter+parameter2] = test_samples['Likelihood']
        test_x[parameter+parameter2] = torch.from_numpy(test_x[parameter+parameter2]).float()
        test_x[parameter+parameter2] = train_x
        test_y[parameter+parameter2] = train_y
        # test_y[parameter+parameter2] = torch.from_numpy(test_y[parameter+parameter2]).float()
    
    # kernel = gpytorch.kernels.RBFKernel(ard_num_dims=2)
    
    # nu = np.mean(train_y.numpy())
    # sigma = np.std(train_y.numpy())
    # train_y = (train_y-nu)/sigma
    # model, likelihood = traingpmodel(train_x,train_y, kernel, sigma, nu)
    # model.eval()
    # likelihood.eval()
    # observed_pred = {}
    # observed_pred_mean = {}
    # with torch.no_grad(), gpytorch.settings.fast_pred_var():
    #     observed_pred[parameter+parameter2] = likelihood(model(test_x[parameter+parameter2]))
    #     observed_pred_mean[parameter+parameter2] = (observed_pred[parameter+parameter2].mean*sigma)+nu
    # print('sqrt(MSE) ',parameter+parameter2, np.sqrt(mean_squared_error(test_y[parameter+parameter2].numpy(),observed_pred_mean[parameter+parameter2].numpy())))

    flatsamples = np.zeros((len(parameterstocheck),resolution_reduced))
    flatsamplesparameters = []
    i = 0
    for parameter in parameterstocheck:
        if parameter != parameter2:
            flatsamples[i,:] = train_y.numpy()
            # flatsamples[i,:] = observed_pred_mean[parameter+parameter2].numpy()
            flatsamplesparameters.append(train_x.numpy())
            i+=1
    minindx = np.unravel_index(flatsamples.argmin(), flatsamples.shape)
    min_parameters = flatsamplesparameters[minindx[0]][minindx[1]]
    min_loglike = flatsamples.min()
    print(min_loglike,loglikelihood(scaletooriginal(min_parameters,boundaries)))
    return min_parameters, min_loglike, test_y

def objective(trial, changeableparameters):
    for parameter in changeableparameters:
        parametervalue = trial.suggest_uniform(parameter,boundaries[parameter][0],boundaries[parameter][1])
        minpGB2[parameter] = parametervalue
        if parameter in ['EclipticLatitude']:
            minpGB2[parameter] = np.arcsin(parametervalue)
        elif parameter in ['Inclination']:
            minpGB2[parameter] = np.arccos(parametervalue)
    p = loglikelihood(minpGB2)
    return p

def model(optimize_parameters):
    i = 0
    for parameter in changeableparameters:
        parametervalue = optimize_parameters[i]
        minpGB2[parameter] = (parametervalue*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
        if parameter in ['EclipticLatitude']:
            minpGB2[parameter] = troch.arcsin(parametervalue)
        elif parameter in ['Inclination']:
            minpGB2[parameter] = torch.arccos(parametervalue)
    p = loglikelihood(minpGB2)
    i += 1
    return p


def scaletooriginal(previous_min,boundaries):
    i = 0
    minpGB = deepcopy(pGBs)
    for parameter in parametersfd:
        if parameter in ['EclipticLatitude']:
            minpGB[parameter] = np.arcsin((previous_min[parametersfd.index(parameter)]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
        elif parameter in ['Inclination']:
            minpGB[parameter] = np.arccos((previous_min[parametersfd.index(parameter)]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
        else:
            minpGB[parameter] = (previous_min[parametersfd.index(parameter)]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
    return minpGB

def plotplanes(parameterstocheck,parameter2, plot_x, plot_y):
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
                im = ax[j,i%4].scatter(plot_x[parameter+parameter2].numpy()[:,parametersfd.index(parameter)],plot_x[parameter+parameter2].numpy()[:,parametersfd.index(parameter2)],c=plot_y[parameter+parameter2].numpy()[:])
                ax[j,i%4].set_xlabel(parameter)
                ax[j,i%4].set_ylabel(parameter2)
                fig.colorbar(im, ax=ax[j,i%4])
        i += 1
    plt.show()





def traingpmodel(train_x,train_y,kernel, sigma, nu):
    bo_iterations = 1
    for bo_iter in range(bo_iterations):
        # initialize likelihood and model
        likelihood = gpytorch.likelihoods.GaussianLikelihood()
        model = ExactGPModel(train_x, train_y, likelihood, kernel)

        training_iter = 50

        hypers = {
            'likelihood.noise': torch.tensor(0.0001),
        }
        model.initialize(**hypers)
        # Polarization
        # model.covar_module.base_kernel.kernels[1].outputscale = torch.tensor([[1.0]])
        # list(model.covar_module.base_kernel.kernels[1].parameters())[0].requires_grad=False
        # model.covar_module.base_kernel.kernels[1].base_kernel.period_length = torch.tensor([[0.25]])
        # list(model.covar_module.base_kernel.kernels[1].base_kernel.parameters())[0].requires_grad=False
        # # InitialPhase
        # # model.covar_module.base_kernel.kernels[0].kernels[1].kernels[0].period_length = torch.tensor([[1.0]])
        # # list(model.covar_module.base_kernel.kernels[0].kernels[1].kernels[0].parameters())[0].requires_grad=False
        # # model.covar_module.base_kernel.kernels[0].kernels[1].kernels[1].period_length = torch.tensor([[1.0]])
        # # list(model.covar_module.base_kernel.kernels[0].kernels[1].kernels[1].parameters())[0].requires_grad=False
        # # Frequency
        # model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[1].lengthscale = torch.tensor([[0.066]])
        # list(model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[1].parameters())[0].requires_grad=False
        # # # Longitude
        # model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[1].period_length = torch.tensor([[1.5]])
        # # list(model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[1].parameters())[1].requires_grad=False
        # # Latitude
        # model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[1].period_length = torch.tensor([[2.0]])
        # list(model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[1].parameters())[1].requires_grad=False
        # # Amplitude
        # model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].lengthscale = torch.tensor([[0.4]])
        # list(model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].parameters())[0].requires_grad=False

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
            print('Iter %d/%d - Loss: %.3f      noise: %.3f' % (
            i + 1, training_iter, loss.item(),
            model.likelihood.noise.item(),
            ))
            # print('Iter %d/%d - Loss: %.3f   Al: %.3f  Lal: %.3f Lap: %.3f  Lol: %.3f Lop: %.3f  Fl: %.3f  Il: %.3f  IPp: %.3f Po: %.3f  Pp: %.3f   out: %.3f' % (
            # i + 1, training_iter, loss.item(),
            # model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].lengthscale.item(),
            # model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[1].lengthscale.item(),
            # model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[0].kernels[1].period_length.item(),
            # model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[1].lengthscale.item(),
            # model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[0].kernels[1].period_length.item(),
            # model.covar_module.base_kernel.kernels[0].kernels[0].kernels[0].kernels[1].lengthscale.item(),
            # model.covar_module.base_kernel.kernels[0].kernels[0].kernels[1].offset.item(),
            # model.covar_module.base_kernel.kernels[0].kernels[1].base_kernel.period_length.item(),
            # model.covar_module.base_kernel.kernels[1].outputscale.item(),
            # model.covar_module.base_kernel.kernels[1].base_kernel.period_length.item(),
            # # model.mean_module.constant.item(),
            # model.covar_module.outputscale.item()
            # ))
            optimizer.step()

        best_value = train_y.min()
        print('best value',bo_iter, (best_value*sigma)+nu, train_x[train_y.argmin()])
        best_value = best_value
        if bo_iter < bo_iterations-1:
            qEI = qExpectedImprovement(model=model, best_f=best_value)
            qUCB = qUpperConfidenceBound(model=model, beta=1)
            candidates = 10
            new_point_analytic, _ = optimize_acqf(
                acq_function=qEI,
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

minpGB = deepcopy(pGBs)
minpGB2 = deepcopy(pGBs)
notreduced = False
notreduced2 = False
boundaries_reduced = deepcopy(boundaries)
previous_best = p1
def outobjective(outtrial):
    for parameter in parametersfd:
        parametervalue = outtrial.suggest_uniform(parameter,boundaries[parameter][0],boundaries[parameter][1])
        minpGB[parameter] = parametervalue
        if parameter in ['EclipticLatitude']:
            minpGB[parameter] = np.arcsin(parametervalue)
        elif parameter in ['Inclination']:
            minpGB[parameter] = np.arccos(parametervalue)
    for i in range(50):
        resolution = 100
        parameter1 = parametersfd[i%7]
        parameter2 = parametersfd[np.random.randint(0,6)]
        parameter3 = parametersfd[np.random.randint(0,6)]
        # parameter2 = 'InitialPhase'
        # parameter1 = 'Inclination'
        while parameter2 == parameter1:
            parameter2 = parametersfd[np.random.randint(0,6)]
        while parameter3 == parameter1 or parameter3 == parameter2:
            parameter3 = parametersfd[np.random.randint(0,6)]
        # if parameter1 == 'Frequency':
        #     parameter2 = 'Polarization'
        parametersreduced = [parameter1]
        changeableparameters = [parameter1,parameter2]
        params = np.zeros(len(changeableparameters))

        # optuna.logging.set_verbosity(optuna.logging.WARNING)
        # start = time.time()
        # study = optuna.create_study(sampler=optuna.samplers.RandomSampler())
        # study.optimize(objective(), n_trials=50)
        # print('optuna time', time.time()-start)
        # if study.best_value < previous_best:
        #     previous_best = study.best_value
        #     for parameter in changeableparameters:
        #         minpGB[parameter] = study.best_params[parameter]
        #         if parameter in ['EclipticLatitude']:
        #             minpGB[parameter] = np.arcsin(study.best_params[parameter])
        #         elif parameter in ['Inclination']:
        #             minpGB[parameter] = np.arccos(study.best_params[parameter])
        # start = time.time()
        # print(i, previous_best, loglikelihood(minpGB), minpGB)
        # minpGB2 = deepcopy(minpGB)
        start = time.time()
        previous_min, min_loglike, observed_pred = planeAdam(minpGB,parametersreduced,parameter2,resolution, boundaries)
        minpGB3 = scaletooriginal(previous_min,boundaries)
        for parameter in changeableparameters:
            minpGB[parameter] = minpGB3[parameter]
        # print('random time',time.time()-start)
        # print(i,min_loglike, minpGB)
        # print(previous_min)
        real_loglike = loglikelihood(minpGB)
    return min_loglike

start = time.time()
outstudy = optuna.create_study(sampler=optuna.samplers.RandomSampler(), pruner= optuna.pruners.HyperbandPruner())
outstudy.optimize(outobjective, n_trials=20)
print('optuna time', time.time()-start)


# plotplanes(parametersreduced,parameter2, test_x, test_y)
# if parameter1 == 'Frequency':
# plotplanes(parametersreduced,parameter2, test_x, observed_pred)
# if real_loglike > -0.9 and notreduced:
#     notreduced = False
#     ratio = 0.2
#     for parameter in parameters:
#         length = boundaries[parameter][1]-boundaries[parameter][0]
#         if parameter == 'EclipticLatitude':
#             boundaries_reduced[parameter] = [np.sin(minpGB[parameter])-length*ratio/2, np.sin(minpGB[parameter])+length*ratio/2] 
#         elif parameter == 'Inclination':
#             boundaries_reduced[parameter] = [np.cos(minpGB[parameter])-length*ratio/2, np.cos(minpGB[parameter])+length*ratio/2] 
#         else:
#             boundaries_reduced[parameter] = [minpGB[parameter]-length*ratio/2, minpGB[parameter]+length*ratio/2] 
#         if boundaries_reduced[parameter][0] <  boundaries[parameter][0]:
#             boundaries_reduced[parameter][0] =  boundaries[parameter][0]
#         if boundaries_reduced[parameter][1] >  boundaries[parameter][1]:
#             boundaries_reduced[parameter][1] =  boundaries[parameter][1]
#     boundaries = boundaries_reduced
# if real_loglike > -0.8 and notreduced2:
#     notreduced2 = False
#     ratio = 0.5
#     for parameter in parameters:
#         length = boundaries[parameter][1]-boundaries[parameter][0]
#         if parameter == 'EclipticLatitude':
#             boundaries_reduced[parameter] = [np.sin(minpGB[parameter])-length*ratio/2, np.sin(minpGB[parameter])+length*ratio/2] 
#         elif parameter == 'Inclination':
#             boundaries_reduced[parameter] = [np.cos(minpGB[parameter])-length*ratio/2, np.cos(minpGB[parameter])+length*ratio/2] 
#         else:
#             boundaries_reduced[parameter] = [minpGB[parameter]-length*ratio/2, minpGB[parameter]+length*ratio/2] 
#         if boundaries_reduced[parameter][0] <  boundaries[parameter][0]:
#             boundaries_reduced[parameter][0] =  boundaries[parameter][0]
#         if boundaries_reduced[parameter][1] >  boundaries[parameter][1]:
#             boundaries_reduced[parameter][1] =  boundaries[parameter][1]
#     boundaries = boundaries_reduced

ratio = 0.2
for parameter in parameters:
    length = boundaries[parameter][1]-boundaries[parameter][0]
    if parameter == 'EclipticLatitude':
        boundaries_reduced[parameter] = [np.sin(minpGB[parameter])-length*ratio/2, np.sin(minpGB[parameter])+length*ratio/2] 
    elif parameter == 'Inclination':
        boundaries_reduced[parameter] = [np.cos(minpGB[parameter])-length*ratio/2, np.cos(minpGB[parameter])+length*ratio/2] 
    else:
        boundaries_reduced[parameter] = [minpGB[parameter]-length*ratio/2, minpGB[parameter]+length*ratio/2] 
    if boundaries_reduced[parameter][0] <  boundaries[parameter][0]:
        boundaries_reduced[parameter][0] =  boundaries[parameter][0]
    if boundaries_reduced[parameter][1] >  boundaries[parameter][1]:
        boundaries_reduced[parameter][1] =  boundaries[parameter][1]
boundaries = boundaries_reduced
resolution = 800
# if parameter != parameter2:
resolution_reduced = int(20**2)
resolution_reduced = resolution
# if 'Frequency' in [parameter,parameter2]:
#     resolution_reduced = int(15**2)
train_samples = sampler(resolution_reduced,parameters,minpGB,boundaries,p1, uniform= False, twoD = False)
train_x = np.zeros((resolution_reduced,len(parameters)))
i = 0
for name in parametersfd:
    train_x[:,i] = train_samples[name]
    i +=1
train_y = train_samples['Likelihood']
train_x = torch.from_numpy(train_x).float()
train_y = torch.from_numpy(train_y).float()
# if parameter != parameter2:
resolution = 100
test_samples = sampler(resolution,parameters,minpGB,boundaries,p1, uniform= False, twoD = False, onlyparameter=parameter, secondparameter=parameter2, calculate_loglikelihood=True)
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
model, likelihood = traingpmodel(train_x,train_y, kernel, sigma, nu)
model.eval()
likelihood.eval()
observed_pred = {}
observed_pred_mean = {}
with torch.no_grad(), gpytorch.settings.fast_pred_var():
    observed_pred[parameter+parameter2] = likelihood(model(test_x[parameter+parameter2]))
    observed_pred_mean[parameter+parameter2] = (observed_pred[parameter+parameter2].mean*sigma)+nu
print('sqrt(MSE) ',parameter+parameter2, np.sqrt(mean_squared_error(test_y[parameter+parameter2].numpy(),observed_pred_mean[parameter+parameter2].numpy())))

resolution = 10**5
test_samples = sampler(resolution,parameters,minpGB,boundaries,p1, uniform= False, twoD = False, onlyparameter=parameter, secondparameter=parameter2, calculate_loglikelihood=False)
test_x = np.zeros((resolution,len(parameters)))
i = 0
for name in parametersfd:
    test_x[:,i] = test_samples[name]
    i +=1
test_x = torch.from_numpy(test_x).float()
with torch.no_grad(), gpytorch.settings.fast_pred_var():
    observed_pred = likelihood(model(test_x))
    observed_pred_mean = (observed_pred.mean*sigma)+nu


flatsamples = np.zeros((1,resolution))
flatsamplesparameters = []
i = 0
# flatsamples[i,:] = train_y.numpy()
flatsamples[i,:] = observed_pred_mean.numpy()
flatsamplesparameters.append(test_x.numpy())
i+=1
minindx = np.unravel_index(flatsamples.argmin(), flatsamples.shape)
min_parameters = flatsamplesparameters[minindx[0]][minindx[1]]
min_loglike = flatsamples.min()
minpGB = scaletooriginal(min_parameters,boundaries)
print('pred',min_loglike,'true',loglikelihood(scaletooriginal(min_parameters,boundaries)), minpGB)


Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB, oversample=4, simulator='synthlisa')
# Run Metropolis-Hastings sampler. ---------------------------------
plt.figure()
ax1=plt.subplot(231)
# plt.plot(dataX_training.f*1000,dataX_training.values, label='data')
ax1.plot(dataX.f*1000,dataX.values.real, label='data')
ax1.plot(Xs.f*1000, Xs.values.real, label='binary')
ax2=plt.subplot(232)
# plt.plot(dataY_training.f*1000,dataY_training.values, label='data')
ax2.plot(dataY.f*1000,dataY.values.real, label='data')
ax2.plot(Ys.f*1000, Ys.values.real, label='binary')
ax3=plt.subplot(233)
# plt.plot(dataZ_training.f*1000,dataZ_training.values, label='data')
ax3.plot(dataZ.f*1000,dataZ.values.real, label='data')
ax3.plot(Zs.f*1000, Zs.values.real, label='binary')
ax4=plt.subplot(234)
# plt.plot(dataX_training.f*1000,dataX_training.values.imag, label='data')
ax4.plot(dataX.f*1000,dataX.values.imag, label='data')
ax4.plot(Xs.f*1000, Xs.values.imag, label='binary')
ax5=plt.subplot(235)
# plt.plot(dataY_training.f*1000,dataY_training.values.imag, label='data')
ax5.plot(dataY.f*1000,dataY.values.imag, label='data')
ax5.plot(Ys.f*1000, Ys.values.imag, label='binary')
ax6=plt.subplot(236)
# plt.plot(dataZ_training.f*1000,dataZ_training.values.imag, label='data')
ax6.plot(dataZ.f*1000,dataZ.values.imag, label='data')
ax6.plot(Zs.f*1000, Zs.values.imag, label='binary')
Xs, Ys, Zs = GB.get_fd_tdixyz(template=minpGB, oversample=4, simulator='synthlisa')
ax1.plot(Xs.f*1000, Xs, label='solution')
ax1.set_xlabel('f [mHz]')
ax1.set_ylabel('X-TDI real [1/Hz]')
ax1.legend()
ax2.plot(Ys.f*1000, Ys, label='solution')
ax2.set_xlabel('f [mHz]')
ax2.set_ylabel('Y-TDI real [1/Hz]')
ax2.legend()
ax3.plot(Zs.f*1000, Zs, label='solution')
ax3.set_xlabel('f [mHz]')
ax3.set_ylabel('Z-TDI real [1/Hz]')
ax3.legend()
ax4.plot(Xs.f*1000, Xs.imag, label='solution')
ax4.set_xlabel('f [mHz]')
ax4.set_label('X-TDI imag [1/Hz]')
ax4.legend()
ax5.plot(Ys.f*1000, Ys.imag, label='solution')
ax5.set_xlabel('f [mHz]')
ax5.set_ylabel('Y-TDI imag [1/Hz]')
ax5.legend()
ax6.plot(Zs.f*1000, Zs.imag, label='solution')
ax6.set_xlabel('f [mHz]')
ax6.set_ylabel('Z-TDI imag [1/Hz]')
ax6.legend()
plt.show()

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