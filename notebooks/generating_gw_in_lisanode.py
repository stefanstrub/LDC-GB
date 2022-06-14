import matplotlib.pyplot as plt
import scipy
import numpy as np
import xarray as xr
from astropy import units as u
import pandas as pd

import ldc.io.hdf5 as hdfio
from ldc.waveform.source import load_gb_catalog, load_mbhb_catalog
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

def semi_fast_tdi(config, pMBHB, t_min, t_max, dt):
    hphc = HpHc.type("MBHB-%d"%s_index, "MBHB", "IMRPhenomD")
    hphc.set_param(pMBHB)
    orbits = Orbits.type(config)
    P = ProjectedStrain(orbits)    
    yArm = P.arm_response(t_min, t_max, dt, [hphc], tt_order=1)
    X = P.compute_tdi_x(np.arange(t_min, t_max, dt))
    return TimeSeries(X, dt=dt)

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
t_max = pMBHB["CoalescenceTime"]+2000#*1.0001#60*60*24*365 # time of observation = 1yr
t_min = pMBHB["CoalescenceTime"]-2000#*0.9998
shift = t_min#int(pMBHB["CoalescenceTime"]*0.9998)
t_max -= shift
t_min -= shift 
# t_min = 0
coalescencetime = pMBHB['CoalescenceTime']
pMBHB['CoalescenceTime'] -= shift
initial_position = 2*np.pi*(((shift)/secondsperyear)%1)
config = {"initial_position": initial_position, "initial_rotation": 0, 
          "nominal_arm_length": 2500000000, "orbit_type": 'analytic'}
lisa_orbits = Orbits.type(config)
pMBHB["ObservationDuration"] = t_max
pMBHB["Cadence"] = dt
pMBHB['Redshift'] += 0.0001

trangeshift = np.arange(t_min, t_max, dt)
start = time.time()
tdi_X = semi_fast_tdi(config, pMBHB, t_min, t_max, dt)
print(time.time()- start)

index_low = np.searchsorted(tdi_ts.t,  t_min+shift)

plt.figure()
plt.plot(tdi_ts.t[index_low:index_low+len(tdi_X)], tdi_ts['X'][index_low:index_low+len(tdi_X)])
plt.plot(trangeshift+shift, tdi_X, label="strain to TDI shifted", alpha=0.5)
# plt.plot(trange-t_min, Xs, label="strain to TDI", alpha=0.5)
plt.xlabel("Time [s]")
# plt.axis([coalescencetime-1000, coalescencetime+600, None, None])
plt.ylabel("TDI X")
plt.legend()

tdi_Xfd = tdi_X.ts.fft(win=window)
tdi_fsX = tdi_ts['X'][index_low:index_low+len(tdi_X)].ts.fft(win=window)
plt.figure()
plt.semilogx(tdi_fsX.f, tdi_fsX)
plt.semilogx(tdi_Xfd.f, tdi_Xfd, label="strain to TDI", alpha=0.5)
plt.xlabel("f[Hz]")
# plt.axis([coalescencetime-1000, coalescencetime+600, None, None])
plt.ylabel("TDI X")
plt.legend()
plt.show()
print('s')
