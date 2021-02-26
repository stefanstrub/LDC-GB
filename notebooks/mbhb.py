import matplotlib.pyplot as plt
import scipy
import numpy as np
import xarray as xr
from astropy import units as u

import ldc.io.hdf5 as hdfio
from ldc.lisa.noise import get_noise_model
from ldc.lisa.orbits import Orbits
from ldc.lisa.projection import ProjectedStrain
from ldc.common.series import TimeSeries, FrequencySeries, window
import ldc.waveform.fastGB as fastGB
from ldc.common.tools import compute_tdi_snr
from ldc.waveform.waveform import HpHc
import time

DATAPATH = "/home/stefan/LDC/Sangria/data"
sangria_fn = DATAPATH+"/LDC2_sangria_training_v1.h5"
tdi_ts, tdi_descr = hdfio.load_array(sangria_fn, name="obs/tdi")
dt = int(1/(tdi_descr["sampling_frequency"]))

# Build timeseries and frequencyseries object for X,Y,Z
tdi_ts = xr.Dataset(dict([(k,TimeSeries(tdi_ts[k], dt=dt)) for k in ["X", "Y", "Z"]]))
tdi_fs = xr.Dataset(dict([(k,tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))

def semi_fast_tdi(config, pMBHB, t_max, dt):
    hphc = HpHc.type("MBHB-%d"%s_index, "MBHB", "IMRPhenomD")
    hphc.set_param(pMBHB)
    orbits = Orbits.type(config)
    P = ProjectedStrain(orbits)    
    yArm = P.arm_response(0, t_max, dt, [hphc], tt_order=1)
    X = P.compute_tdi_x(np.arange(0, t_max, dt))
    return TimeSeries(X, dt=dt)

def mldc_fast_tdi(pMBHB, t_max, dt):
    from GenerateFD_SignalTDIs import ComputeMBHBXYZ_FD
    from LISAhdf5 import ParsUnits
    hphc = HpHc.type("MBHB-%d"%s_index, "MBHB", "IMRPhenomD")
    hphc.set_param(pMBHB)
    pMBHB["Cadence"] = dt
    pMBHB["ObservationDuration"] = t_max/2
    pu = ParsUnits(pars_i=pMBHB, units_i=hphc.info())
    fr, Xs, Ys, Zs = ComputeMBHBXYZ_FD(pu)
    return Xs
            

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
s_index = 0
pMBHB = dict(zip(mbhb.dtype.names, mbhb[s_index]))
units = default_units
#for k,v in pMBHB.items():
#   print(k)
# #  pMBHB[k] *= u.Unit(units[k])

t_max = float(tdi_ts["X"].t[-1]+tdi_ts["X"].attrs["dt"])

t_max = 10000
t_min = 0
pMBHB['ColaescenceTime'] = 500
pMBHB["ObservationDuration"] = t_max
config = {"initial_position": 0, "initial_rotation": 0, 
          "nominal_arm_length": 2500000000, "orbit_type": 'analytic'}
# pMBHB["CoalescenceTime"] = 100
# t_max = pMBHB["CoalescenceTime"]+100#60*60*24*365 # time of observation = 1yr
start = time.time()
Xs = semi_fast_tdi(config, pMBHB, t_max, dt)
print(time.time()- start)
print(Xs)



noise_model = "MRDv1"
Nmodel = get_noise_model(noise_model, np.logspace(-5, -1, 100))
SNR2 = np.zeros((2,2)) 
start = time.time()
Xs,Ys,Zs = mldc_fast_tdi(pMBHB, t_max, dt)
print(time.time()-start)
source = dict({"X":Xs, "Y":Ys, "Z":Zs})
start = time.time()
SNR2[0,1] = compute_tdi_snr(source, Nmodel, data=tdi_fs)["tot2"]
print(time.time()-start)
SNR2[0,0] = compute_tdi_snr(source, Nmodel)["tot2"]
loglikelihood = SNR2[0,1] - 0.5*SNR2[0,0]
SNR2[1,1] = compute_tdi_snr(source, Nmodel, data=tdi_fs)["tot2"]
print(time.time()-start)
SNR2[1,0] = compute_tdi_snr(source, Nmodel)["tot2"]
loglikelihood = SNR2[1,1] - 0.5*SNR2[1,0]
# Xsfastfd,Ysfastfd,Zsfastfd   = mldc_fast_tdi(pMBHB, t_max, dt)
# Xsfast = Xsfastfd.ts.ifft()
# Xsfd = Xs.ts.fft(win=window)

SNR2 = np.zeros(2) # snr square
fmin, fmax = float(Xs.f[0]) , float(Xs.f[-1]+Xs.attrs['df'])
source = dict({"X":Xs, "Y":Ys, "Z":Zs})
SNR2[1] = compute_tdi_snr(source, Nmodel, data=tdi_fs, fmin=fmin, fmax=fmax)["tot2"]
SNR2[0] = compute_tdi_snr(source, Nmodel)["tot2"] 

# print(Xsfast)
plt.figure(figsize=(12,6))
# plt.semilogx(tdi_fs["X"].f, tdi_fs["X"], label="TDI X")
plt.semilogx(Xsfastfd.f, (Xsfastfd), label=" fast %d"%s_index)
plt.semilogx(Xsfd.f, (Xsfd), label="semi fast %d"%s_index)
# plt.axis([pMBHB["CoalescenceTime"]-1000, pMBHB["CoalescenceTime"]+600, None, None])
# plt.ylim(-3*10**-19,3*10**-19)
plt.legend(loc="lower right")
# plt.xlabel("time [s]")
plt.figure(figsize=(12,6))
plt.plot(tdi_ts["X"].t, tdi_ts["X"], label="TDI X")
plt.plot(Xs.t, (tdi_ts["X"]-Xs), label="TDI X - fast %d"%s_index)
plt.plot(Xs.t, (Xs), label="semi fast %d"%s_index)
plt.plot(Xs.t[:-1], (Xsfast), label=" fast %d"%s_index)
plt.axis([pMBHB["CoalescenceTime"]-1000, pMBHB["CoalescenceTime"]+600, None, None])
plt.ylim(-3*10**-19,3*10**-19)
plt.legend(loc="lower right")
plt.xlabel("time [s]")
plt.show()