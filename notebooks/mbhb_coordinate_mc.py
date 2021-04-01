#%%
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rcParams
import scipy
import numpy as np
import xarray as xr
from astropy import units as u
import pandas as pd
import time
from copy import deepcopy
import Cosmology

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
fig_width_pt = 464.0  # Get this from LaTeX using \showthe\columnwidth
golden_mean = (np.sqrt(5.0) - 1.0) / 2.0  # Aesthetic ratio
ratio = golden_mean
inches_per_pt = 1.0 / 72.27  # Convert pt to inches
fig_width = fig_width_pt * inches_per_pt*1.3  # width in inches
fig_height = fig_width * ratio  # height in inches
fig_size = [fig_width, fig_height]
fig_size_squared = [fig_width, fig_width]
rcParams.update({"figure.figsize": fig_size})

np.random.seed(40)

def semi_fast_tdi(config, pMBHBi, t_min, t_max, dt):
    hphc = HpHc.type("MBHB-%d"%s_index, "MBHB", "IMRPhenomD")
    hphc.set_param(pMBHBi)
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
    
def scaletooriginal(psample,pMBHBs,boundaries):
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
    # pMBHBnew['Distance'] = 10**pMBHBnew['Distance']
    # v = hubbleconstant*(10**3*pMBHBnew['Distance'])
    # pMBHBnew['Redshift'] = np.sqrt((1+v/speedoflight)/(1-v/speedoflight)) - 1
    return pMBHBnew    

def scaleto01(pMBHBs,boundaries):
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
sangria_fn = DATAPATH+"/dMBHB-tdi.h5"
# sangria_fn = DATAPATH+"/LDC2_sangria_blind_v1.h5"
# sangria_fn = DATAPATH+"/LDC2_sangria_gdb-tdi_v1_v3U3MxS.h5"
# sangria_fn = DATAPATH+"/LDC2_sangria_idb-tdi_v1_DgtGV85.h5"
# sangria_fn = DATAPATH+"/LDC2_sangria_mbhb-tdi_v1_MN5aIPz.h5"
sangria_fn = DATAPATH+"/LDC2_sangria_training_v1.h5"
tdi_ts, tdi_descr = hdfio.load_array(sangria_fn, name="obs/tdi")
# sangria_fn = DATAPATH+"/LDC2_sangria_vMBHB-tdi_v1_sgsEVXb.h5"
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
#find MBHB parameters
def loglikelihood(pMBHBs):
    Xs = semi_fast_tdi(config, pMBHBs, t_min, t_max, dt)
    Xsfd = Xs.ts.fft(win=window)
    index_low = np.searchsorted(Xsfd.f, dataX.f[0])
    Xsfd = Xsfd[index_low:index_low+len(dataX)]
    diff = np.abs(dataX - Xsfd.values)**2
    # p1 = -float(np.sum(diff / Sn)*Xsfd.attrs['df'])/2.0
    p1 = float(np.sum(diff / Sn)/len(diff))/2.0
    # p1 = np.exp(p1)
    return -p1

def changeparameterstoinput(pMBHB01, boundaries):
    pMBHBnew = deepcopy(pMBHBs)
    for parameter in parametersboundary:
        if parameter in parameters:
            pMBHBnew[parameter] = (pMBHB01[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
        elif parameter in ['ChirpMass']:
            ChirpMass = (pMBHB01[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
        elif parameter in ['MassRatio']:
            MassRatio = (pMBHB01[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
        elif parameter in ['sinEclipticLatitude']:
            pMBHBnew['EclipticLatitude'] = np.arcsin((pMBHB01[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
        else:
            print('parameter missed',parameter)
    pMBHBnew['Mass1'], pMBHBnew['Mass2'] = funcm1m2ofMchirpq(ChirpMass,MassRatio)
    # pMBHBnew['Distance'] = 10**pMBHB01['Distance']
    # v = hubbleconstant*(10**3*pMBHBnew['Distance'])
    # pMBHBnew['Redshift'] = np.sqrt((1+v/speedoflight)/(1-v/speedoflight)) - 1
    return pMBHBnew    
def inputto01(pMBHB, boundaries):
    pMBHB01 = {}
    Mc = funcMchirpofm1m2(pMBHB['Mass1'], pMBHB['Mass2'])
    q = pMBHB['Mass1']/pMBHB['Mass2']
    for parameter in parametersboundary:
        if parameter in parameters:
            pMBHB01[parameter] = (pMBHB[parameter]-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0])
        elif parameter in ['ChirpMass']:
            pMBHB01[parameter] = (Mc-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0])
        elif parameter in ['MassRatio']:
            pMBHB01[parameter] = (q-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0])
        elif parameter in ['sinEclipticLatitude']:
            pMBHB01['sinEclipticLatitude'] = (np.sin(pMBHB['EclipticLatitude'])-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0])
        else:
            print('parameter missed',parameter)
    # parameter = 'Distance'
    # pMBHB01[parameter] = (np.log10(pMBHB[parameter])-boundaries[parameter][0])/(boundaries[parameter][1]-boundaries[parameter][0])
    # v = hubbleconstant*(10**3*pMBHBnew['Distance'])
    # pMBHBnew['Redshift'] = np.sqrt((1+v/speedoflight)/(1-v/speedoflight)) - 1
    return pMBHB01  

# Total number of proposed samples.
number_of_samples = 8*10 **1
cutoff_ratio = 1000

parameters = ['EclipticLatitude','EclipticLongitude',
         'PolarAngleOfSpin1','PolarAngleOfSpin2',
         'Spin1','Spin2',
         'Mass1','Mass2',
         'CoalescenceTime','PhaseAtCoalescence',
         'InitialPolarAngleL','InitialAzimuthalAngleL',
         'Distance']

pMBHBs = deepcopy(pMBHB)
# pMBHBs['PolarAngleOfSpin1'] = 0
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

boundaries = {'ChirpMass': [Mc*(1.0 - 0.01), Mc*(1.0 + 0.01)], ### broadening just lead to longer run time
'MassRatio': [1.0, 10.0],
'CoalescenceTime': [pMBHB['CoalescenceTime'] - 50., pMBHB['CoalescenceTime'] + 50.], ### could be reduced to this range by computing the slide
'Spin1': [-0.99, 0.99],
'Spin2': [-0.99, 0.99],
'PolarAngleOfSpin1': [0, np.pi],
'PolarAngleOfSpin2': [0, np.pi],
# 'Distance': [np.log10(pMBHB['Distance'])*(0.9), np.log10(pMBHB['Distance'])*(1.1)],
'sinEclipticLatitude': [-1.0, 1.0],
'EclipticLongitude': [0.0, 2.0*np.pi],
'InitialAzimuthalAngleL': [0.0, 2.0*np.pi],
'InitialPolarAngleL': [0.0, np.pi],
'PhaseAtCoalescence': [0.0, 1.0*np.pi]}
parametersboundary = ['ChirpMass','MassRatio','CoalescenceTime','Spin1','Spin2','PolarAngleOfSpin1','PolarAngleOfSpin2','InitialPolarAngleL',
'sinEclipticLatitude','EclipticLongitude','InitialAzimuthalAngleL','PhaseAtCoalescence']
psample = np.random.rand(13)
samples1 = {}
for parameter in parametersboundary:
    samples1[parameter] = psample[parametersboundary.index(parameter)]
pMBHBs = changeparameterstoinput(samples1, boundaries)
pMBHBs['Distance'] = pMBHB['Distance']

# redshift = Cosmology.zofDl(pMBHB['Distance'],w=0, tolerance=0.00001)
Distance = Cosmology.DL(pMBHB['Redshift'],w=0)[0]
trangeshift = np.arange(t_min, t_max, dt)
start = time.time()
tdi_X = semi_fast_tdi(config, pMBHB, t_min, t_max, dt)
print('pMBHB time',time.time()- start)
start = time.time()
Xs = semi_fast_tdi(config, pMBHBs, t_min, t_max, dt)
print('pMBHBs time',time.time()- start)

index_low = np.searchsorted(tdi_ts.t,  t_min+shift)

secondsperday = 24*3600
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
plt.figure(figsize=fig_size)
plt.plot(tdi_ts.t[index_low:index_low+len(tdi_X)]/secondsperday, tdi_ts['X'][index_low:index_low+len(tdi_X)], label="Data")
plt.plot((trangeshift+shift)/secondsperday, tdi_X, label="MBHB", alpha=1, color= colors[2])
plt.plot((trangeshift+shift)/secondsperday, Xs, label="sample", alpha=0.5)
# plt.plot(trangeshift-t_min, Xs, label="strain to TDI", alpha=0.5)
plt.xlabel("Time (days)")
# plt.axis([coalescencetime-1000, coalescencetime+600, None, None])
plt.locator_params(nbins=4)
plt.ylabel("TDI X strain")
plt.tight_layout()
plt.legend()
tdi_Xfd = tdi_X.ts.fft(win=window)
Xsfd = Xs.ts.fft(win=window)
tdi_fdX = tdi_ts['X'][index_low:index_low+len(tdi_X)].ts.fft(win=window)
plt.figure()
plt.semilogx(tdi_fdX.f, tdi_fdX, label="Data")
plt.semilogx(tdi_Xfd.f, tdi_Xfd, label="MBHB", alpha=1, color=colors[2])
plt.semilogx(Xsfd.f, Xsfd, label="sample", alpha=0.5)
plt.xlabel("f (Hz)")
plt.ylabel("TDI X real (1/Hz)")
plt.legend()
plt.tight_layout()
plt.show()

abovethreshold = tdi_fdX.f > 1.5e-2
higherindex = np.where(abovethreshold)[0][1]
dataX = tdi_fdX[1:higherindex]
Xsfd = Xsfd[1:higherindex]
fmin, fmax = float(Xsfd.f[0]) , float(Xsfd.f[-1]+Xsfd.attrs['df'])
freq = np.array(Xsfd.sel(f=slice(fmin, fmax)).f)
Sn = Nmodel.psd(freq=freq, option='X')
diff = np.abs(dataX - Xsfd.values)**2
p1 = -float(np.sum(diff / (Sn))/len(diff))/2.0


plt.figure()
plt.semilogx(tdi_fdX.f, tdi_fdX, label="Data")
plt.semilogx(dataX.f, dataX, label="dataX", alpha=1, color=colors[2])
plt.semilogx(Xsfd.f, Xsfd, label="sample2", alpha=0.5)
plt.xlabel("f (Hz)")
plt.ylabel("TDI X real (1/Hz)")
plt.legend()
plt.tight_layout()
plt.show()

samples1 = inputto01(pMBHBs,boundaries)
samples = {}
for parameter in parametersboundary:
    samples[parameter] = []
    samples[parameter].append(samples1[parameter])
samples['Likelihood'] = []
samples['Likelihood'].append(p1)

print('p start',p1, 'p true', loglikelihood(pMBHB))

def sampler(number_of_samples,parameters,pMBHB,boundaries,p1, uniform=False, MCMC=False, only=False, onlyparameter='Frequency', twoD=False, secondparameter='Amplitude', calculate_loglikelihood=True):
    samples1 = inputto01(pMBHB,boundaries)
    samples = {}
    for parameter in parametersboundary:
        samples[parameter] = []
        samples[parameter].append(samples1[parameter])
    p1 = loglikelihood(pMBHB)
    samples['Likelihood'] = []
    samples['Likelihood'].append(p1)
    j = 0
    pMBHBs01 = samples1
    number_of_sampels_sqrt = np.sqrt(number_of_samples)
    for i in range(1, number_of_samples):
        if only:
            parameter = onlyparameter
            if uniform:
                pMBHBs01[parameter] = i/number_of_samples
            else:
                pMBHBs01[parameter] = np.random.rand()
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
            for parameter in parametersboundary:
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
        pMBHBs = changeparameterstoinput(pMBHBs01, boundaries)
        # pMBHBs01array = np.zeros(len(parametersfd))
        # i = 0
        # for name in parametersfd:
        #     pMBHBs01array[i] = pMBHBs01[name]
        #     i +=1
        # pMBHBscaled = scaletooriginal(pMBHBs01array ,boundaries)
        # print(loglikelihood(pMBHBs), loglikelihood(pMBHBscaled))
        pMBHBs['Distance'] = pMBHB['Distance']
        if calculate_loglikelihood:
            p1 = loglikelihood(pMBHBs)
            samples['Likelihood'].append(p1)
        for parameter in parametersboundary:
            samples[parameter].append(pMBHBs01[parameter])
    for parameter in parametersboundary:
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
        # kernel = kernelA * kernelLat * kernelLong * kernelF * kernelI * kernelIP * kernelP
        # self.covar_module = gpytorch.kernels.ScaleKernel(gpytorch.kernels.RBFKernel(ard_num_dims=8))
        self.covar_module = gpytorch.kernels.ScaleKernel(kernel)
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
        if parameter in ['ChirpMass']:
            ChirpMass = parametervalue
            MassRatio = maxpMBHB2['Mass1'] / maxpMBHB2['Mass2']
        elif parameter in ['MassRatio']:
            MassRatio = parametervalue
            try:
                ChirpMass
            except:
                ChirpMass = funcMchirpofm1m2(maxpMBHB2['Mass1'], maxpMBHB2['Mass2'])
        elif parameter in ['sinEclipticLatitude']:
            maxpMBHB2['EclipticLatitude'] = np.arcsin(parametervalue)
        else:
            maxpMBHB2[parameter] = parametervalue
    try:
        maxpMBHB2['Mass1'], maxpMBHB2['Mass2'] = funcm1m2ofMchirpq(ChirpMass,MassRatio)
    except:
        pass
    # maxpMBHB2['Distance'] = 10**parametervlue
    maxpMBHB2['Distance'] = pMBHB['Distance']
    p = loglikelihood(maxpMBHB2)
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

        training_iter = 30

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
        ], lr=0.01)  # Includes GaussianLikelihood parameters
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

maxpMBHB = deepcopy(pMBHBs)
notreduced = False
notreduced2 = False
best_params = deepcopy(pMBHBs)
best_value = p1
parameters_recorded = []
best_run = 0
for n in range(2):
    parameters_recorded1 = {}
    no_improvement_counter = 0
    psample = {}
    for parameter in parametersboundary:
        psample[parameter] = np.random.uniform()
    
    maxpMBHB = changeparameterstoinput(psample,boundaries)
    for parameter in parameters:
        parameters_recorded1[parameter] = []
        parameters_recorded1[parameter].append(maxpMBHB[parameter])
    parameters_recorded1['Loglikelihood'] = []
    if n == 0:
        maxpMBHB = deepcopy(pMBHBs)
    maxpMBHB['Distance'] = pMBHB['Distance']
    previous_best = loglikelihood(maxpMBHB)
    maxpMBHB2 = deepcopy(maxpMBHB)
    for i in range(100):
        resolution = 300
        parameter1 = parametersboundary[i%len(parametersboundary)]
        parameter2 = parametersboundary[np.random.randint(0,len(parametersboundary)-1)]
        parameter3 = parametersboundary[np.random.randint(0,len(parametersboundary)-1)]
        # parameter2 = 'InitialPhase'
        # parameter1 = 'Inclination'
        while parameter2 == parameter1:
            parameter2 = parametersboundary[np.random.randint(0,len(parametersboundary)-1)]
        while parameter3 == parameter1 or parameter3 == parameter2:
            parameter3 = parametersboundary[np.random.randint(0,len(parametersboundary)-1)]
        # if parameter1 == 'Frequency':
        #     parameter2 = 'Polarization'
        parametersreduced = [parameter1]
        changeableparameters = [parameter1,parameter2]
        params = np.zeros(len(changeableparameters))

        optuna.logging.set_verbosity(optuna.logging.WARNING)
        start = time.time()
        study = optuna.create_study(sampler=optuna.samplers.RandomSampler(),direction='maximize')
        study.optimize(objective, n_trials=50)
        print('optuna time', time.time()-start)
        if study.best_value > previous_best:
            no_improvement_counter = 0
            previous_best = study.best_value
            for parameter in changeableparameters:
                parametervalue = study.best_params[parameter]
                if parameter in ['ChirpMass']:
                    ChirpMass = parametervalue
                    MassRatio = maxpMBHB['Mass1'] / maxpMBHB['Mass2']
                elif parameter in ['MassRatio']:
                    MassRatio = parametervalue
                    try:
                        ChirpMass
                    except:
                        ChirpMass = funcMchirpofm1m2(maxpMBHB['Mass1'], maxpMBHB['Mass2'])
                elif parameter in ['sinEclipticLatitude']:
                    maxpMBHB['EclipticLatitude'] = np.arcsin(parametervalue)
                else:
                    maxpMBHB[parameter] = parametervalue
            try:
                maxpMBHB['Mass1'], maxpMBHB['Mass2'] = funcm1m2ofMchirpq(ChirpMass,MassRatio)
            except:
                pass
            maxpMBHB['Distance'] = pMBHB['Distance']
        else:
            no_improvement_counter += 1
        if no_improvement_counter > 16:
            past_mean = 0
            sum_count = 0
            for l in range(n):
                s = 0
                try:
                    past_mean += parameters_recorded[l]["Loglikelihood"][i]
                    sum_count += 1
                except:
                    pass
            try:
                past_mean = past_mean / sum_count
            except:
                no_improvement_counter = 0
            if previous_best > past_mean:
                no_improvement_counter = 0
            # elif sum_count == 0:
            #     no_improvement_counter = 0
            else:
                print("no improvement")
                break

        print(i, previous_best, loglikelihood(maxpMBHB), maxpMBHB)
        parameters_recorded1['Loglikelihood'].append(loglikelihood(maxpMBHB))
        for parameter in parameters:
            parameters_recorded1[parameter].append(maxpMBHB[parameter])

        maxpMBHB2 = deepcopy(maxpMBHB)
        x = 1
    if previous_best > best_value:
        best_run = n
        best_value = previous_best
        best_params = deepcopy(maxpMBHB)
    parameters_recorded.append(parameters_recorded1)

maxpMBHB = best_params
ratio = 0.1
boundaries_reduced = deepcopy(boundaries)
Mc = funcMchirpofm1m2(pMBHB['Mass1'], pMBHB['Mass2'])
q = pMBHB['Mass1']/pMBHB['Mass2']
for parameter in parametersboundary:
    length = boundaries[parameter][1]-boundaries[parameter][0]
    print(parameter)
    if parameter in parameters:
        boundaries_reduced[parameter] = [maxpMBHB[parameter]-length*ratio/2, maxpMBHB[parameter]+length*ratio/2] 
    elif parameter in ['ChirpMass']:
        boundaries_reduced[parameter] = [Mc-length*ratio/2, Mc+length*ratio/2] 
    elif parameter in ['MassRatio']:
        boundaries_reduced[parameter] = [q-length*ratio/2, q+length*ratio/2]
    if parameter in ['Frequency']:
        boundaries_reduced[parameter] = [maxpMBHB[parameter]-length*ratio/2, maxpMBHB[parameter]+length*ratio/2] 
    elif parameter in ['EclipticLongitude']:
        print('yes')
        boundaries_reduced[parameter] = [maxpMBHB[parameter]-length*ratio/8, maxpMBHB[parameter]+length*ratio/8] 
    elif parameter in ['sinEclipticLatitude']:
        boundaries_reduced[parameter] = [np.sin(maxpMBHB['EclipticLatitude'])-length*ratio/2, np.sin(maxpMBHB['EclipticLatitude'])+length*ratio/2] 
    if boundaries_reduced[parameter][0] < boundaries[parameter][0]:
        boundaries_reduced[parameter][0] = boundaries[parameter][0]
    if boundaries_reduced[parameter][1] > boundaries[parameter][1]:
        boundaries_reduced[parameter][1] = boundaries[parameter][1]

resolution = 1000
# if parameter != parameter2:
resolution_reduced = int(20**2)
resolution_reduced = resolution
# if 'Frequency' in [parameter,parameter2]:
#     resolution_reduced = int(15**2)
start = time.time()
parameter = 'Frequency'
train_samples = sampler(resolution_reduced,parameters,maxpMBHB,boundaries_reduced,p1, uniform= False, twoD = False, onlyparameter=parameter, secondparameter=parameter2)
print(time.time()-start)
train_x = np.zeros((resolution_reduced,len(parametersboundary)))
i = 0
for name in parametersboundary:
    train_x[:,i] = train_samples[name]
    i +=1
train_y = train_samples['Likelihood']
train_x = torch.from_numpy(train_x).float()
train_y = torch.from_numpy(train_y).float()
# if parameter != parameter2:
resolution = 500
test_samples = sampler(resolution,parameters,maxpMBHB,boundaries_reduced,p1, uniform= False, twoD = False, onlyparameter=parameter, secondparameter=parameter2, calculate_loglikelihood=True)
test_x[parameter+parameter2] = np.zeros((resolution,len(parametersboundary)))
i = 0
for name in parametersboundary:
    test_x[parameter+parameter2][:,i] = test_samples[name]
    i +=1
test_y[parameter+parameter2] = test_samples['Likelihood']
test_x[parameter+parameter2] = torch.from_numpy(test_x[parameter+parameter2]).float()
test_y[parameter+parameter2] = torch.from_numpy(test_y[parameter+parameter2]).float()

kernel = gpytorch.kernels.RBFKernel(ard_num_dims=12)

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
test_samples = sampler(resolution,parameters,maxpMBHB,boundaries_reduced,p1, uniform= False, twoD = False, onlyparameter=parameter, secondparameter=parameter2, calculate_loglikelihood=False)
test_x = np.zeros((resolution,len(parametersboundary)))
i = 0
for name in parametersboundary:
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
maxindx = np.unravel_index(flatsamples.argmax(), flatsamples.shape)
max_parameters = flatsamplesparameters[maxindx[0]][maxindx[1]]
max_loglike = flatsamples.max()
maxpMBHBpredicted = scaletooriginal(max_parameters,pMBHB,boundaries_reduced)
if loglikelihood(maxpMBHBpredicted) > loglikelihood(maxpMBHB):
    maxpMBHB = maxpMBHBpredicted 
print('pred',max_loglike,'true',loglikelihood(scaletooriginal(max_parameters,pMBHB,boundaries_reduced)),'max', loglikelihood(maxpMBHB), maxpMBHB)
parameters_recorded[best_run]['Loglikelihood'].append(loglikelihood(maxpMBHB))
for parameter in parameters:
    parameters_recorded[best_run][parameter].append(maxpMBHB[parameter])

fig, ax = plt.subplots(2, 4,figsize=np.asarray(fig_size)*2)
# plt.suptitle("Intermediate Parameters")
for n in range(len(parameters_recorded)):
    i = 0    
    for parameter in parameters+['Loglikelihood']:
        j = 0
        if i > 3:
            j = 1
        if parameter in ['Frequency']:
            ax[j,i%4].plot(np.asarray(parameters_recorded[n][parameter])*10**3,np.arange(0,len(parameters_recorded[n][parameter])))
        else:
            ax[j,i%4].plot(parameters_recorded[n][parameter],np.arange(0,len(parameters_recorded[n][parameter])))
        i += 1
i = 0
for parameter in parameters+['Log-likelihood']:
    j = 0
    if i > 3:
        j = 1
    if parameter in ['Log-likelihood']:
        ax[j,i%4].axvline(x=loglikelihood(pMBHB), color='k',label='True')
    elif parameter in ['Frequency']:
        ax[j,i%4].axvline(x=pMBHB[parameter]*10**3, color='k',label='True')
    else:
        ax[j,i%4].axvline(x=pMBHB[parameter], color='k',label='True')
    if parameter in ['Amplitude']:
        ax[j,i%4].set_xlabel(parameter, ha='right', va='top')
    elif parameter in ['Frequency']:
        ax[j,i%4].xaxis.set_major_locator(plt.MaxNLocator(3))
        ax[j,i%4].set_xlabel(parameter + ' (mHz)')
    else:
        ax[j,i%4].xaxis.set_major_locator(plt.MaxNLocator(3))
        ax[j,i%4].set_xlabel(parameter)
    # if parameter in ['Log-likelihood']:
    #     ax[j,i%4].set_xlabel('$$')
    ax[j,i%1].set_ylabel('Step')
    ax[j,i%4].set_ylim(0,100)
    # if parameter in ['Log-likelihood']:
    #     pass
    # elif parameter == 'EclipticLatitude':
    #     ax[j,i%4].set_xlim(np.arcsin(boundaries[parameter][0]),np.arcsin(boundaries[parameter][1]))
    # elif parameter == 'Inclination':
    #     ax[j,i%4].set_xlim(np.arccos(boundaries[parameter][1]),np.arccos(boundaries[parameter][0]))
    # elif parameter in ['Frequency']:
    #     ax[j,i%4].set_xlim(boundaries[parameter][0]*10**3,boundaries[parameter][1]*10**3)
    # else:
    #     ax[j,i%4].set_xlim(boundaries[parameter][0],boundaries[parameter][1])
    i += 1
ax[0,0].legend()
plt.tight_layout()
plt.show()

Xs = semi_fast_tdi(config, pMBHB, t_min, t_max, dt)
Xsfd = Xs.ts.fft(win=window)
index_low = np.searchsorted(Xsfd.f, dataX.f[0])
Xsfd = Xsfd[index_low:index_low+len(dataX)]
plt.figure(figsize=np.asarray(fig_size)*2)
ax1=plt.subplot(231)
# plt.plot(dataX_training.f*1000,dataX_training.values, label='data')
ax1.plot(dataX.f*1000,dataX.values.real, label='Data')
ax1.plot(Xsfd.f*1000, Xsfd.values.real, label='MBHB')
ax4=plt.subplot(234)
# plt.plot(dataX_training.f*1000,dataX_training.values.imag, label='data')
ax4.plot(dataX.f*1000,dataX.values.imag, label='Data')
ax4.plot(Xsfd.f*1000, Xsfd.values.imag, label='True Response')
Xs = semi_fast_tdi(config, maxpMBHB, t_min, t_max, dt)
Xsfd = Xs.ts.fft(win=window)
index_low = np.searchsorted(Xsfd.f, dataX.f[0])
Xsfd = Xsfd[index_low:index_low+len(dataX)]
ax1.plot(Xsfd.f*1000, Xsfd, label='Found MBHB',alpha=0.5)
ax1.xaxis.set_major_locator(plt.MaxNLocator(3))
ax1.set_xlabel('f (mHz)')
ax1.set_ylabel('X-TDI real (1/Hz)')
ax1.legend(loc='upper left',numpoints=0.5)
# ax3.legend()
ax4.plot(Xsfd.f*1000, Xsfd.imag, label='Found Response')
ax4.xaxis.set_major_locator(plt.MaxNLocator(3))
ax4.set_xlabel('f (mHz)')
ax4.set_ylabel('X-TDI imag (1/Hz)')
# ax4.legend()
# ax6.legend()
plt.tight_layout()
plt.show()

# Get into evaluation (predictive posterior) mode
model.eval()
likelihood.eval()
observed_pred = {}

number_of_test_samples = 100
test_x = {}
test_y = {}
for parameter in parametersboundary:
    test_samples = sampler(number_of_test_samples,parameters,pMBHB,boundaries_reduced,p1, uniform= True, only=True, onlyparameter=parameter)
    test_x[parameter] = np.zeros((number_of_test_samples,len(parametersboundary)))
    i = 0
    for name in parametersboundary:
        test_x[parameter][:,i] = test_samples[name]
        i +=1
    test_y[parameter] = test_samples['Likelihood']
    test_x[parameter] = torch.from_numpy(test_x[parameter]).float()
    test_x[parameter] = test_x[parameter]
    test_y[parameter] = torch.from_numpy(test_y[parameter]).float()
# number_of_test_samples = 100
# test_x2 = {}
# test_y2 = {}
# for parameter in parametersboundary:
#     test_samples = sampler(number_of_test_samples,parameters,pMBHB,boundaries,p1, uniform= True, only=True, onlyparameter=parameter)
#     test_x2[parameter] = np.zeros((number_of_test_samples,len(parametersboundary)))
#     i = 0
#     for name in parametersboundary:
#         test_x2[parameter][:,i] = test_samples[name]
#         i +=1
#     test_y2[parameter] = test_samples['Likelihood']
#     test_x2[parameter] = torch.from_numpy(test_x2[parameter]).float()
#     test_x2[parameter] = test_x2[parameter]
#     test_y2[parameter] = torch.from_numpy(test_y2[parameter]).float()
# number_of_test_samples2d = 20**2
# parameter2 = 'EclipticLatitude'
# for parameter in parametersboundary:
#     if parameter != parameter2:
#         test_samples = sampler(number_of_test_samples2d,parameters,pMBHB,boundaries_reduced,p1, uniform= True, twoD = True, onlyparameter=parameter, secondparameter=parameter2, calculate_loglikelihood=False)
#         test_x[parameter+parameter2] = np.zeros((number_of_test_samples2d,len(parametersboundary)))
#         i = 0
#         for name in parametersboundary:
#             test_x[parameter+parameter2][:,i] = test_samples[name]
#             i +=1
#         test_x[parameter+parameter2] = torch.from_numpy(test_x[parameter+parameter2]).float()

for parameter in parametersboundary:
    # Make predictions by feeding model through likelihood
    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        observed_pred[parameter] = likelihood(model(test_x[parameter]))
        observed_pred_print = (observed_pred[parameter].mean.numpy()*sigma)+nu
    print('sqrt(MSE) ',parameter, np.sqrt(mean_squared_error(test_y[parameter].numpy(),observed_pred_print)))
    # if parameter != parameter2:
    #     with torch.no_grad(), gpytorch.settings.fast_pred_var():
    #         observed_pred[parameter+parameter2] = likelihood(model(test_x[parameter+parameter2]))

# parameter = 'random'
# with torch.no_grad(), gpytorch.settings.fast_pred_var():
#     observed_pred[parameter] = likelihood(model(test_x[parameter]))
# observed_pred_print = (observed_pred[parameter]*sigma)+nu
# print('sqrt(MSE) ','random', np.sqrt(mean_squared_error(test_y[parameter].numpy(),observed_pred_print.mean.cpu().numpy())))


# for parameter in parameters:
#     if parameter in ['EclipticLatitude']:
#         samples[parameter] = np.arcsin((samples[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
#         test_samples[parameter] = np.arcsin((test_samples[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
#     elif parameter in ['Inclination']:
#         samples[parameter] = np.arccos((samples[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
#         test_samples[parameter] = np.arccos((test_samples[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0])
#     else:
#         samples[parameter] = (samples[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]
#         test_samples[parameter] = (test_samples[parameter]*(boundaries[parameter][1]-boundaries[parameter][0]))+boundaries[parameter][0]


train_x = train_x.cpu()
train_y = train_y.cpu()
# train_y = (train_y*sigma)+nu

fig, ax = plt.subplots(3, 4,figsize=np.asarray(fig_size)*2.5)
# plt.suptitle("loglikelihood")
i = 0    
mean = {}
for parameter in parametersboundary:
    j = 0
    if i > 3:
        j = 1
    if i > 7:
        j = 2
    with torch.no_grad():
        # Get upper and lower confidence bounds
        lower, upper = observed_pred[parameter].confidence_region()
        mean2 = observed_pred[parameter].mean
        mean[parameter] = (mean2*sigma)+nu
        lower = (lower*sigma)+nu
        upper = (upper*sigma)+nu
        # Plot training data as black stars
        if parameter == 'Frequency':
            ax[j,i%4].axvline(x=inputto01(pMBHB, boundaries_reduced)[parameter]*10**3, color='k',label='MBHB')
            ax[j,i%4].plot(np.asarray(test_x2[parameter].numpy()[1:])*10**3, test_y2[parameter].numpy()[1:], 'g',label='True')
            ax[j,i%4].plot(np.asarray(test_x[parameter].numpy()[1:])*10**3, mean[parameter].numpy()[1:], 'b',label='Mean')
            ax[j,i%4].fill_between(np.asarray(test_x[parameter].numpy()[1:])*10**3, lower.numpy()[1:], upper.numpy()[1:], alpha=0.5,label='Confidence')
        else:
            ax[j,i%4].axvline(x=inputto01(pMBHB, boundaries_reduced)[parameter], color='k',label='MBHB')
            ax[j,i%4].plot(test_x2[parameter].numpy()[1:,parametersboundary.index(parameter)], test_y2[parameter].numpy()[1:], 'g',label='True')
            ax[j,i%4].plot(test_x[parameter].numpy()[1:,parametersboundary.index(parameter)], mean[parameter].numpy()[1:], 'b',label='Mean')
            ax[j,i%4].fill_between(test_x[parameter].numpy()[1:,parametersboundary.index(parameter)], lower.numpy()[1:], upper.numpy()[1:], alpha=0.5,label='Confidence')
        ax[0,0].legend()
    if parameter in ['Amplitude']:
        ax[j,i%4].set_xlabel(parameter, ha='right', va='top')
    elif parameter in ['Frequency']:
        ax[j,i%4].set_xlabel(parameter + ' (mHz)')
    else:
        ax[j,i%4].set_xlabel(parameter)
    ax[j,i%4].set_ylabel('Log-likelihood')
    i += 1
plt.tight_layout()
plt.show()
# fig, ax = plt.subplots(2, 4,figsize=(15,15))
# plt.suptitle("loglikelihood true")
# i = 0    
# for parameter in parametersfd:
#     if parameter != parameter2:
#         j = 0
#         if i > 3:
#             j = 1
#         with torch.no_grad():
#             ax[j,i%4].axvline(x=pMBHB01[parameter], color='k')
#             ax[j,i%4].axhline(y=pMBHB01[parameter2], color='k')
#             im = ax[j,i%4].scatter(test_x[parameter+parameter2].numpy()[:,i],test_x[parameter+parameter2].numpy()[:,parametersfd.index(parameter2)],c=test_y[parameter+parameter2].numpy()[:])
#             ax[j,i%4].set_xlabel(parameter)
#             ax[j,i%4].set_ylabel(parameter2)
#             fig.colorbar(im, ax=ax[j,i%4])
#     else:
#         ax[j,i%4].plot(test_x[parameter].numpy()[1:,i], test_y[parameter].numpy()[1:], 'g.')
#         ax[j,i%4].set_xlabel(parameter)
#         ax[j,i%4].set_ylabel('loglikelihood')
#         ax[j,i%4].legend(['True'])
#     i += 1
    
fig, ax = plt.subplots(2, 4,figsize=(15,15))
plt.suptitle("loglikelihood predicted mean")
i = 0    
for parameter in parametersfd:
    if parameter != parameter2:
        j = 0
        if i > 3:
            j = 1
        with torch.no_grad():
            # Get upper and lower confidence bounds
            mean2 = observed_pred[parameter+parameter2].mean
            mean[parameter+parameter2] = (mean2*sigma)+nu
            test_x_rescaled1 = []
            test_x_rescaled2 = []
            for k in range(len(test_x[parameter+parameter2])):
                test_x_rescaled1.append(scaletooriginal(test_x[parameter+parameter2][k].numpy(),boundaries_reduced)[parameter])
                test_x_rescaled2.append(scaletooriginal(test_x[parameter+parameter2][k].numpy(),boundaries_reduced)[parameter2])
            ax[j,i%4].axvline(x=pMBHB[parameter], color='k')
            ax[j,i%4].axhline(y=pMBHB[parameter2], color='k')
            im = ax[j,i%4].scatter(test_x_rescaled1,test_x_rescaled2,c=mean[parameter+parameter2][:])
            ax[j,i%4].set_xlabel(parameter)
            ax[j,i%4].set_ylabel(parameter2)
            fig.colorbar(im, ax=ax[j,i%4])

            if parameter == 'EclipticLatitude':
                ax[j,i%4].set_xlim(np.arcsin(boundaries_reduced[parameter][0]),np.arcsin(boundaries_reduced[parameter][1]))
            elif parameter == 'Inclination':
                ax[j,i%4].set_xlim(np.arccos(boundaries_reduced[parameter][1]),np.arccos(boundaries_reduced[parameter][0]))
            elif parameter in ['Frequency']:
                ax[j,i%4].set_xlim(boundaries_reduced[parameter][0]*10**3,boundaries_reduced[parameter][1]*10**3)
            else:
                ax[j,i%4].set_xlim(boundaries_reduced[parameter][0],boundaries_reduced[parameter][1])
            if parameter2 == 'EclipticLatitude':
                ax[j,i%4].set_ylim(np.arcsin(boundaries_reduced[parameter2][0]),np.arcsin(boundaries_reduced[parameter2][1]))
            elif parameter2 == 'Inclination':
                ax[j,i%4].set_ylim(np.arccos(boundaries_reduced[parameter2][1]),np.arccos(boundaries_reduced[parameter2][0]))
            elif parameter2 in ['Frequency']:
                ax[j,i%4].set_ylim(boundaries_reduced[parameter2][0]*10**3,boundaries_reduced[parameter2][1]*10**3)
            else:
                ax[j,i%4].set_ylim(boundaries_reduced[parameter2][0],boundaries_reduced[parameter2][1])
    else:
        with torch.no_grad():
            lower, upper = observed_pred[parameter].confidence_region()
            lower = lower.cpu()
            upper = upper.cpu()
            lower = (lower*sigma)+nu
            upper = (upper*sigma)+nu
            ax[j,i%4].plot(test_x[parameter].numpy()[1:,i], test_y[parameter].numpy()[1:], 'g.')
            ax[j,i%4].plot(test_x[parameter].numpy()[1:,i], mean[parameter].numpy()[1:], 'b.')
            ax[j,i%4].fill_between(test_x[parameter].numpy()[1:,i], lower.numpy()[1:], upper.numpy()[1:], alpha=0.5)
            ax[j,i%4].set_xlabel(parameter)
            ax[j,i%4].set_ylabel('loglikelihood')
            ax[j,i%4].legend(['True','Mean', 'Confidence'])
    i += 1
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
# plt.xlim(max(samples['Amplitude']),max(samples['Amplitude']))
# plt.ylim(max(samples['Frequency']),max(samples['Frequency']))
# plt.colorbar()
# plt.figure()
# # plt.title('sqrt(MSE)',np.round(np.sqrt(mean_squared_error(test_y,prediction)),2))
# plt.scatter(samples['EclipticLatitude'],samples['EclipticLongitude'],c=samples['Likelihood'])
# plt.xlim(max(samples['EclipticLatitude']),max(samples['EclipticLatitude']))
# plt.ylim(max(samples['EclipticLongitude']),max(samples['EclipticLongitude']))
# plt.colorbar()
plt.show()