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

def semi_fast_tdi(config, pMBHB, t_min, t_max, dt):
    hphc = HpHc.type("MBHB-%d"%s_index, "MBHB", "IMRPhenomD")
    hphc.set_param(pMBHB)
    orbits = Orbits.type(config)
    P = ProjectedStrain(orbits)    
    yArm = P.arm_response(t_min, t_max, dt, [hphc], tt_order=1)
    X = P.compute_tdi_x(np.arange(t_min, t_max, dt))
    return TimeSeries(X, dt=dt)

def AziPolAngleL2PsiIncl(bet, lam, theL, phiL):
    """
    Convert Polar and Azimuthal angles of zS (typically orbital angular momentum L)
    to polarisation and inclination (see doc)
    @param bet is the ecliptic latitude of the source in sky [rad]
    @param lam is the ecliptic longitude of the source in sky [rad]
    @param theL is the polar angle of zS [rad]
    @param phiL is the azimuthal angle of zS [rad]
    @return polarisation and inclination
    """
    #inc = np.arccos( np.cos(theL)*np.sin(bet) + np.cos(bet)*np.sin(theL)*np.cos(lam - phiL) )
    #up_psi = np.cos(theL)*np.cos(bet) - np.sin(bet)*np.sin(theL)*np.cos(lam - phiL)
    #down_psi = np.sin(theL)*np.sin(lam - phiL)
    #psi = np.arctan2(up_psi, down_psi)

    inc = np.arccos( - np.cos(theL)*np.sin(bet) - np.cos(bet)*np.sin(theL)*np.cos(lam - phiL) )
    down_psi = np.sin(theL)*np.sin(lam - phiL)
    up_psi = -np.sin(bet)*np.sin(theL)*np.cos(lam - phiL) + np.cos(theL)*np.cos(bet)
    psi = np.arctan2(up_psi, down_psi)

    return psi, inc
    
def scaletooriginal(psample,boundaries):
    pMBHBnew = deepcopy(pMBHBs)
    for parameter in parametersboundary:
        if parameter in parameters:
            pMBHBnew[parameter] = (psample[parametersboundary.index(parameter)]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
        elif parameter in ['ChirpMass']:
            ChirpMass = (psample[parametersboundary.index(parameter)]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
        elif parameter in ['MassRatio']:
            MassRatio = (psample[parametersboundary.index(parameter)]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
        elif parameter in ['sinEclipticLatitude']:
            pMBHBnew['EclipticLatitude'] = np.arcsin((psample[parametersboundary.index(parameter)]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
        else:
            print('parameter missed',parameter)
    pMBHBnew['Mass1'], pMBHBnew['Mass2'] = funcm1m2ofMchirpq(ChirpMass,MassRatio)
    pMBHBnew['Distance'] = 10**pMBHBnew['Distance']
    # v = hubbleconstant*(10**3*pMBHBnew['Distance'])
    # pMBHBnew['Redshift'] = np.sqrt((1+v/speedoflight)/(1-v/speedoflight)) - 1
    return pMBHBnew

hubbleconstant = 70000
speedoflight = 299792458
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


default_units = {'EclipticLatitude':'rad','EclipticLongitude':'rad',
         'PolarAngleOfSpin1':'rad','PolarAngleOfSpin2':'rad',
         'Spin1': '1','Spin2':'1',
         'Mass1':'Msun','Mass2':'Msun',
         'CoalescenceTime': 's','PhaseAtCoalescence':'rad',
         'InitialPolarAngleL':'rad','InitialAzimuthalAngleL':'rad',
         'Cadence': 's','Redshift': '1','Distance': 'Gpc',
         'ObservationDuration':'s'}

mbhb, units = hdfio.load_array(sangria_fn, name="sky/mbhb/cat")
print(units)
if not units:
    units = default_units
config = hdfio.load_config(sangria_fn, name="obs/config")
print(config)
secondsperyear = 60*60*24*365.25
s_index = 0
pMBHB = dict(zip(mbhb.dtype.names, mbhb[s_index]))
dt = 5 # waveform sampling
def funcMchirpofm1m2(m1, m2):
    return pow(m1*m2, 3./5) / pow(m1+m2, 1./5)
def funcm1m2ofMchirpq(Mc, q):
    m2 = (1+q)**(1/5)*Mc/q**(3/5)
    return m2*q, m2
Mc = funcMchirpofm1m2(pMBHB['Mass1'], pMBHB['Mass2'])
q = pMBHB['Mass1']/pMBHB['Mass2']
m1, m2 = funcm1m2ofMchirpq(Mc,q)
######################################################
#%%
#find GB parameters
def loglikelihood(pMBHBs):
    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pMBHBs, oversample=4, simulator='synthlisa')
    index_low = np.searchsorted(Xs.f, dataX.f[0])
    Xs = Xs[index_low:index_low+len(dataX)]
    Ys = Ys[index_low:index_low+len(dataY)]
    Zs = Zs[index_low:index_low+len(dataZ)]
    diff = np.abs(dataX - Xs.values)**2 + np.abs(dataY - Ys.values)**2 + np.abs(dataZ - Zs.values)**2
    # p1 = -float(np.sum(diff / Sn)*Xs.attrs['df'])/2.0
    p1 = float(np.sum(diff / (Sn+noise))/len(diff))/2.0
    # p1 = np.exp(p1)
    return p1

# Total number of proposed samples.
number_of_samples = 8*10 **1
cutoff_ratio = 1000

parameters = ['EclipticLatitude','EclipticLongitude',
         'PolarAngleOfSpin1','PolarAngleOfSpin2',
         'Spin1','Spin2',
         'Mass1','Mass2',
         'CoalescenceTime','PhaseAtCoalescence',
         'InitialPolarAngleL','InitialAzimuthalAngleL',
         'Redshift','Distance',
         'ObservationDuration']
boundaries = {'ChirpMass': [Mc*(1.0 - 0.01), Mc*(1.0 + 0.01)], ### broadening just lead to longer run time
'MassRatio': [1.0, 10.0],
'CoalescenceTime': [pMBHB['CoalescenceTime'] - 100., pMBHB['CoalescenceTime'] + 100.], ### could be reduced to this range by computing the slide
'Spin1': [-0.99, 0.99],
'Spin2': [-0.99, 0.99],
'PolarAngleOfSpin1': [0, np.pi],
'PolarAngleOfSpin1': [0, np.pi],
'Distance': [np.log10(pMBHB['Distance'])*(1.0), np.log10(pMBHB['Distance'])*(1.0)],
'InitialPolarAngleL': [0.0, np.pi],
'sinEclipticLatitude': [-1.0, 1.0],
'EclipticLongitude': [0.0, 2.0*np.pi],
'InitialAzimuthalAngleL': [0.0, 2.0*np.pi],
'PhaseAtCoalescence': [0.0, 2.0*np.pi]}
parametersboundary = ['ChirpMass','MassRatio','CoalescenceTime','Spin1','Spin2','PolarAngleOfSpin1','PolarAngleOfSpin1','Distance','InitialPolarAngleL',
'sinEclipticLatitude','EclipticLongitude','InitialAzimuthalAngleL','PhaseAtCoalescence']
psample = np.random.rand(13)
pMBHBs = deepcopy(pMBHB)
pMBHBs['PolarAngleOfSpin1'] = 0
# pMBHBs = scaletooriginal(psample,boundaries)

t_max = pMBHB["CoalescenceTime"]+1000#*1.0001#60*60*24*365 # time of observation = 1yr
t_min = pMBHB["CoalescenceTime"]-2000#*0.9998
shift = t_min#int(pMBHB["CoalescenceTime"]*0.9998)
t_max -= shift
t_min -= shift 
# t_min = 0
coalescencetime = pMBHB['CoalescenceTime']
pMBHB['CoalescenceTime'] -= shift
pMBHBs['CoalescenceTime'] -= shift
initial_position = 2*np.pi*(((shift)/secondsperyear)%1)
config = {"initial_position": initial_position, "initial_rotation": 0, 
          "nominal_arm_length": 2500000000, "orbit_type": 'analytic'}
lisa_orbits = Orbits.type(config)
pMBHB["ObservationDuration"] = t_max
pMBHB["Cadence"] = dt
pMBHBs["ObservationDuration"] = t_max
pMBHBs["Cadence"] = dt

trangeshift = np.arange(t_min, t_max, dt)
start = time.time()
tdi_X = semi_fast_tdi(config, pMBHB, t_min, t_max, dt)
print(time.time()- start)
start = time.time()
Xs = semi_fast_tdi(config, pMBHBs, t_min, t_max, dt)
print(time.time()- start)

index_low = np.searchsorted(tdi_ts.t,  t_min+shift)

plt.figure()
plt.plot(tdi_ts.t[index_low:index_low+len(tdi_X)], tdi_ts['X'][index_low:index_low+len(tdi_X)], label="data")
plt.plot(trangeshift+shift, tdi_X, label="binary", alpha=0.5)
plt.plot(trangeshift+shift, Xs, label="sample", alpha=0.5)
# plt.plot(trange-t_min, Xs, label="strain to TDI", alpha=0.5)
plt.xlabel("Time [s]")
# plt.axis([coalescencetime-1000, coalescencetime+600, None, None])
plt.ylabel("TDI X")
plt.legend()
tdi_Xfd = tdi_X.ts.fft(win=window)
Xsfd = Xs.ts.fft(win=window)
tdi_fdX = tdi_ts['X'][index_low:index_low+len(tdi_X)].ts.fft(win=window)
plt.figure()
plt.semilogx(tdi_fdX.f, tdi_fdX)
plt.semilogx(tdi_Xfd.f, tdi_Xfd, label="strain to TDI", alpha=0.5)
plt.semilogx(Xsfd.f, Xsfd, label="sample", alpha=0.5)
plt.xlabel("f[Hz]")
# plt.axis([coalescencetime-1000, coalescencetime+600, None, None])
plt.ylabel("TDI X")
plt.legend()
plt.show()

# psd_signal = np.abs(Xs.values)**2 + np.abs(Ys.values)**2 + np.abs(Zs.values)**2
# highSNR = psd_signal > np.max(psd_signal)/cutoff_ratio
# lowerindex = np.where(highSNR)[0][0]-10
# higherindex = np.where(highSNR)[0][-1]+10
# dataX = tdi_fs["X"].isel(f=slice(Xs.kmin, Xs.kmin+len(Xs)))[lowerindex:higherindex]
# dataY = tdi_fs["Y"].isel(f=slice(Ys.kmin, Ys.kmin+len(Ys)))[lowerindex:higherindex]
# dataZ = tdi_fs["Z"].isel(f=slice(Zs.kmin, Zs.kmin+len(Zs)))[lowerindex:higherindex]
# spd_data = np.abs(dataX)**2 + np.abs(dataY)**2 + np.abs(dataZ)**2
# noise = (np.mean(spd_data[:5])+np.mean(spd_data[-5:])).values/2
# Xs, Ys, Zs = Xs[lowerindex:higherindex], Ys[lowerindex:higherindex], Zs[lowerindex:higherindex]
fmin, fmax = float(Xsfd.f[0]) , float(Xsfd.f[-1]+Xsfd.attrs['df'])
freq = np.array(Xsfd.sel(f=slice(fmin, fmax)).f)
Sn = Nmodel.psd(freq=freq, option='X')
diff = np.abs(tdi_fdX - Xsfd.values)**2
p1 = float(np.sum(diff / (Sn))/len(diff))/2.0
# p1 = np.exp(p1)

frequency_lower_boundary = dataX.f[0].values+(dataX.f[-1].values - dataX.f[0].values)*4/10
frequency_upper_boundary = dataX.f[-1].values-(dataX.f[-1].values - dataX.f[0].values)*3/10
# [0.0004725, 0.0004727]
# boundaries_small = deepcopy(boundaries)
# part_ratio = 10
# for parameter in parameters:
#     if parameter in ['EclipticLongitude','Frequency']:
#         boundaries_small[parameter] = [pMBHB[parameter]-(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio,pMBHB[parameter]+(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio]
#     if parameter == 'EclipticLatitude':
#         boundaries_small[parameter] = [np.sin(pMBHB[parameter])-(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio,np.sin(pMBHB[parameter])+(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio]
#     elif parameter == 'Inclination':
#         boundaries_small[parameter] = [np.cos(pMBHB[parameter])-(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio,np.cos(pMBHB[parameter])+(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio]
#     else:
#         boundaries_small[parameter] = [pMBHB[parameter]-(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio,pMBHB[parameter]+(boundaries[parameter][1]-boundaries[parameter][0])/part_ratio]
# boundaries = boundaries_small



samples = xr.Dataset(dict([(name,xr.DataArray(np.zeros(number_of_samples), dims=('number_of_sample'), coords={"number_of_sample": range(number_of_samples)},
                         )) for name, titles in pMBHBs.items()]))
pMBHBs01 = {}
for parameter in parameters:
    if parameter in ['EclipticLatitude']:
        samples[parameter][0] = (np.sin(pMBHBs[parameter])-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0])
    elif parameter in ['Inclination']:
        samples[parameter][0] = (np.cos(pMBHBs[parameter])-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0])
    else:
        samples[parameter][0] = (pMBHBs[parameter]-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0])
    pMBHBs01[parameter] = samples[parameter][0]
samples['Likelihood'] = samples['Name']
samples = samples.drop(['Name'])
samples['Likelihood'][0] = p1

print('p start',p1, 'p true', loglikelihood(pMBHB))

def sampler(number_of_samples,parameters,pMBHB,boundaries,p1, uniform=False, MCMC=False, only=False, onlyparameter='Frequency', twoD=False, secondparameter='Amplitude', calculate_loglikelihood=True):
    samples = xr.Dataset(dict([(name,xr.DataArray(np.zeros(number_of_samples), dims=('number_of_sample'), coords={"number_of_sample": range(number_of_samples)},
                         )) for name, titles in pMBHB.items()]))
    samples = {}
    pMBHBs01 = {}
    for parameter in parameters:
        samples[parameter] = []
        if parameter in ['EclipticLatitude']:
            samples[parameter].append((np.sin(pMBHB[parameter])-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0]))
        elif parameter in ['Inclination']:
            samples[parameter].append((np.cos(pMBHB[parameter])-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0]))
        else:
            samples[parameter].append((pMBHB[parameter]-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0]))
        pMBHBs01[parameter] = samples[parameter][0]
    samples['Likelihood'] = []
    p1 = loglikelihood(pMBHB)
    samples['Likelihood'].append(p1)
    j = 0
    number_of_sampels_sqrt = np.sqrt(number_of_samples)
    for i in range(1, number_of_samples):
        if only:
            parameter = onlyparameter
            if uniform:
                pMBHBs01[parameter] = i/number_of_samples
            else:
                pMBHBs01[parameter] = np.random.rand()
            if parameter != 'InitialPhase':
                pMBHBs01['InitialPhase'] = samples['InitialPhase'][0]
            if parameter != 'Polarization':
                pMBHBs01['Polarization'] = samples['Polarization'][0]
        elif twoD:
            parameter = onlyparameter
            parameter2 = secondparameter
            if uniform:
                if i % number_of_sampels_sqrt == 0:
                    j += 1
                pMBHBs01[parameter] = ((i-1)%number_of_sampels_sqrt+1)/(number_of_sampels_sqrt+2)
                pMBHBs01[parameter2] = (j+1)/(number_of_sampels_sqrt+2)
            else: 
                pMBHBs01[parameter] = np.random.rand()
                pMBHBs01[parameter2] = np.random.rand()
            # if parameter != 'InitialPhase':
            #     pMBHBs01['InitialPhase'] = samples['InitialPhase'][0]
        else:
            for parameter in parameters:
                pMBHBs01[parameter] = np.random.rand()
                # if parameter in ['Amplitude']:#,'FrequencyDerivative','Amplitude','EclipticLongitude']:
                #     pMBHBs01[parameter] = i/number_of_samples
                # elif parameter in ['FrequencyDerivative']:
                #     pass
                # if parameter in ['Inclination']:
                #     pMBHBs01[parameter] = np.random.rand()
                # elif parameter in ['Inclination']:
                #     pMBHBs01[parameter] = np.random.rand()
                # else:
                #     pMBHBs01[parameter] = np.random.rand()
        for parameter in parameters:
            # if parameter in ['FrequencyDerivative']:
            #     pass
            if parameter in ['EclipticLatitude']:
                pMBHBs[parameter] = np.arcsin((pMBHBs01[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
            elif parameter in ['Inclination']:
                pMBHBs[parameter] = np.arccos((pMBHBs01[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
            else:
                pMBHBs[parameter] = (pMBHBs01[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
        # pMBHBs01array = np.zeros(len(parametersfd))
        # i = 0
        # for name in parametersfd:
        #     pMBHBs01array[i] = pMBHBs01[name]
        #     i +=1
        # pMBHBscaled = scaletooriginal(pMBHBs01array ,boundaries)
        # print(loglikelihood(pMBHBs), loglikelihood(pMBHBscaled))
        if calculate_loglikelihood:
            p1 = loglikelihood(pMBHBs)
            samples['Likelihood'].append(p1)
        for parameter in parameters:
            samples[parameter].append(pMBHBs01[parameter])
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
#     test_samples = sampler(number_of_test_samples,parameters,pMBHBs,boundaries,p1, uniform= True, only=True, onlyparameter=parameter)
#     test_x[parameter] = np.zeros((number_of_test_samples,len(parameters)))
#     i = 0
#     for name in parametersfd:
#         test_x[parameter][:,i] = test_samples[name]
#         i +=1
#     test_y[parameter] = test_samples['Likelihood']
#     test_x[parameter] = torch.from_numpy(test_x[parameter]).float()
#     test_x[parameter] = test_x[parameter]
#     test_y[parameter] = torch.from_numpy(test_y[parameter]).float()

def objective(trial):
    for parameter in changeableparameters:
        parametervalue = trial.suggest_uniform(parameter,boundaries[parameter][0],boundaries[parameter][1])
        minpMBHB2[parameter] = parametervalue
        if parameter in ['EclipticLatitude']:
            minpMBHB2[parameter] = np.arcsin(parametervalue)
        elif parameter in ['Inclination']:
            minpMBHB2[parameter] = np.arccos(parametervalue)
    p = loglikelihood(minpMBHB2)
    return p


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
                        pMBHBs[parameter] = np.arcsin((new_point_analytic[candi_n][para_n]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
                    elif parameter in ['Inclination']:
                        pMBHBs[parameter] = np.arccos((new_point_analytic[candi_n][para_n]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
                    else:
                        pMBHBs[parameter] = (new_point_analytic[candi_n][para_n]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
                    para_n += 1
                loglike = loglikelihood(pMBHBs)
                # print('new point',i, new_point_analytic[candi_n], loglike)

                loglike = (loglike-nu)/sigma
                loglike = torch.tensor([loglike]).float()
                train_y = torch.cat((train_y,loglike),0)
    return model, likelihood

minpMBHB = deepcopy(pMBHBs)
minpMBHB2 = deepcopy(pMBHBs)
notreduced = False
notreduced2 = False
boundaries_reduced = deepcopy(boundaries)
best_params = deepcopy(pMBHBs)
best_value = p1
for n in range(3):
    no_improvement_counter = 0
    for parameter in parametersfd:
        parametervalue = np.random.uniform(boundaries[parameter][0],boundaries[parameter][1])
        minpMBHB[parameter] = parametervalue
        if parameter in ['EclipticLatitude']:
            minpMBHB[parameter] = np.arcsin(parametervalue)
        elif parameter in ['Inclination']:
            minpMBHB[parameter] = np.arccos(parametervalue)
    previous_best = p1
    for i in range(100):
        resolution = 300
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

        optuna.logging.set_verbosity(optuna.logging.WARNING)
        start = time.time()
        study = optuna.create_study(sampler=optuna.samplers.RandomSampler())
        study.optimize(objective, n_trials=50)
        print('optuna time', time.time()-start)
        if study.best_value < previous_best:
            no_improvement_counter = 0
            previous_best = study.best_value
            for parameter in changeableparameters:
                minpMBHB[parameter] = study.best_params[parameter]
                if parameter in ['EclipticLatitude']:
                    minpMBHB[parameter] = np.arcsin(study.best_params[parameter])
                elif parameter in ['Inclination']:
                    minpMBHB[parameter] = np.arccos(study.best_params[parameter])
        else:
            no_improvement_counter += 1
        if no_improvement_counter > 16:
            print('no improvement')
            if previous_best < best_value:
                best_value = previous_best
                best_params = deepcopy(minpMBHB)
            break 
        start = time.time()
        print(i, previous_best, loglikelihood(minpMBHB), minpMBHB)
        minpMBHB2 = deepcopy(minpMBHB)

minpMBHB = best_params
ratio = 0.2
for parameter in parameters:
    length = boundaries[parameter][1]-boundaries[parameter][0]
    if parameter == 'EclipticLatitude':
        boundaries_reduced[parameter] = [np.sin(minpMBHB[parameter])-length*ratio/2, np.sin(minpMBHB[parameter])+length*ratio/2] 
    elif parameter == 'Inclination':
        boundaries_reduced[parameter] = [np.cos(minpMBHB[parameter])-length*ratio/2, np.cos(minpMBHB[parameter])+length*ratio/2] 
    else:
        boundaries_reduced[parameter] = [minpMBHB[parameter]-length*ratio/2, minpMBHB[parameter]+length*ratio/2] 
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
train_samples = sampler(resolution_reduced,parameters,minpMBHB,boundaries,p1, uniform= False, twoD = False, onlyparameter=parameter, secondparameter=parameter2)
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
test_samples = sampler(resolution,parameters,minpMBHB,boundaries,p1, uniform= False, twoD = False, onlyparameter=parameter, secondparameter=parameter2, calculate_loglikelihood=True)
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
test_samples = sampler(resolution,parameters,minpMBHB,boundaries,p1, uniform= False, twoD = False, onlyparameter=parameter, secondparameter=parameter2, calculate_loglikelihood=False)
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
minpMBHB = scaletooriginal(min_parameters,boundaries)
print('pred',min_loglike,'true',loglikelihood(scaletooriginal(min_parameters,boundaries)), minpMBHB)


Xs, Ys, Zs = GB.get_fd_tdixyz(template=pMBHB, oversample=4, simulator='synthlisa')
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
Xs, Ys, Zs = GB.get_fd_tdixyz(template=minpMBHB, oversample=4, simulator='synthlisa')
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

pMBHB01 = {}
for parameter in parameters:
    if parameter in ['EclipticLatitude']:
        pMBHB01[parameter] = ((np.sin(pMBHB[parameter])-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0]))
    elif parameter in ['Inclination']:
        pMBHB01[parameter] = ((np.cos(pMBHB[parameter])-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0]))
    else:
        pMBHB01[parameter] = ((pMBHB[parameter]-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0]))


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
        ax[j,i%4].axvline(x=pMBHB01[parameter], color='k',label='True')
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
            ax[j,i%4].axvline(x=pMBHB01[parameter], color='k')
            ax[j,i%4].axhline(y=pMBHB01[parameter2], color='k')
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
#             ax[j,i%4].axvline(x=pMBHB01[parameter], color='k')
#             ax[j,i%4].axhline(y=pMBHB01[parameter2], color='k')
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
    axes[j,i%4].axvline(x=pMBHB[parameter], color='r')
    axes[j,i%4].plot(samples[parameter],samples['Likelihood'], '.',label='train')
    axes[j,i%4].plot(test_samples[parameter],test_samples['Likelihood'],'.',label='true')
    axes[j,i%4].errorbar(test_samples[parameter][1:],prediction[1:],prediction[1:]-lower.numpy()[1:], fmt='.',label='prediction')
    # plt.fill_between(test_samples[parameter][1:],prediction[1:]-p_std[1:],prediction[1:]+p_std[1:], color='g', alpha= 0.4)
    if i == 0:
        axes[j,i%4].set_ylabel('log-Likelihood')
    # plt.subplot(3,number_of_parameters,i+number_of_parameters)
    # plt.axvline(x=pMBHB[parameter], color='r')
    # n, bins, patches = plt.hist(samples[parameter], n_bin, density=True, facecolor='k', alpha=0.5)

    # plt.subplot(3,number_of_parameters,i+2*number_of_parameters)
    # if parameter == 'Frequency':
    #     plt.axvline(x=pMBHB[parameter]*1000, color='r')
    #     plt.plot(samples[parameter]*1000, range(number_of_samples), 'k')
    # else:
    #     plt.axvline(x=pMBHB[parameter], color='r')
    #     plt.plot(samples[parameter], range(number_of_samples), 'k')
    axes[j,i%4].set_xlabel(parameter)
    i += 1
plt.legend()

name = 'Frequency'
plt.figure(figsize=(10,8))
plt.suptitle("sampled posterior")
plt.subplot(1,1,1)
plt.axvline(x=pMBHB[name], color='r')
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