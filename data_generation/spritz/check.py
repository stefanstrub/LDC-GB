import matplotlib.pyplot as plt
import scipy
import numpy as np
import xarray as xr
from astropy import units as u
import pandas as pd

import ldc.io.hdf5 as hdfio
#from ldc.waveform.source import load_gb_catalog, load_mbhb_catalog
#from ldc.lisa.noise import get_noise_model
from ldc.lisa.orbits import Orbits
from ldc.lisa.projection import ProjectedStrain
#from ldc.common.series import TimeSeries, FrequencySeries
#import ldc.waveform.fastGB as fastGB
#from ldc.common.tools import compute_tdi_snr, window
#from ldc.waveform.waveform import HpHc
from ldc.common import constants
CLIGHT = constants.Nature.VELOCITYOFLIGHT_CONSTANT_VACUUM
from ldc.waveform.waveform.hphc import HpHc
import ldc.io.hdf5 as h5io
dt = 5 # waveform sampling
t_max = 60*60*24*7 #365 # time of observation = 1yr
t_min = 0

config = {"initial_position": 0, "initial_rotation": 0, 
          "nominal_arm_length": 2500000000, "orbit_type": 'analytic'}
lisa_orbits = Orbits.type(config)

projector = ProjectedStrain(lisa_orbits)
projector.from_file("run1/vgb-y.h5")
#GWs = HpHc.type("test", "GB", "TD_fdot")
#cat, units = h5io.load_array("run1/vgb-cat.h5")
#GWs.set_param(cat, units=units)
#GWs = GWs.split()
#yArm = projector.arm_response(t_min, t_max, dt, GWs, tt_order=1)

if 1:
    tdi_li, attr = hdfio.load_array("run1/vgb-lisainstrument-noisefree-tdi.h5")
    trange0 = tdi_li[:,0]
    tdi_X0 = tdi_li[:,1]
if 1:
    tdi_ln, attr = hdfio.load_array("run1/vgb-lisanode-noisefree-tdi.h5")
    trange1 = tdi_ln[:,0]
    tdi_X1 = tdi_ln[:,1]

    #istart = 0#np.where(tdi_ln["X"][:,0]>=0)[0][0]
    #tdi_X1 = tdi_ln["X"][istart:,1]
    #trange1 = tdi_ln["X"][istart:,0]
    
    #tdi_X0 = tdi_ln[:,1]/2.816E14

tdi_X0_ = projector.compute_tdi_x(trange0, tdi2=True)
tdi_X1_ = projector.compute_tdi_x(trange1, tdi2=True)

plt.figure(figsize=(15,10))
plt.subplot(211)
#plt.axvline(x=196912.03946846604, alpha=0.5, color='k')
plt.plot(trange0, tdi_X0, label="lisa instrument", alpha=0.5)
plt.plot(trange1, tdi_X1, label="lisa node", alpha=0.5)
plt.plot(trange1, tdi_X1_, label="strain to TDI", alpha=0.5)
plt.xlabel("Time [s]")
plt.ylabel("TDI X")
plt.legend()
plt.axis([5000, 7000, -1e-22, 1e-22])

plt.subplot(212)
#plt.axvline(x=196912.03946846604, alpha=0.5, color='k')
plt.plot(trange0, tdi_X0-tdi_X0_, label="diff with lisa instrument", alpha=0.5)
plt.plot(trange1, tdi_X1-tdi_X1_, label="diff with lisanode", alpha=0.5)
plt.xlabel("Time [s]")
plt.ylabel("TDI X")
plt.legend()
#plt.axis([196000, 198000, None, None])
plt.axis([5000, 7000, -2e-24, 2e-24])

