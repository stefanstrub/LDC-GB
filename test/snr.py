import matplotlib.pyplot as plt
from astropy import units as un
import numpy as np
from ldc.lisa import orbits
from ldc.lisa.projection import ProjectedStrain
from ldc.waveform.waveform import HpHc
from ldc.common.series import TimeSeries, FrequencySeries
from ldc.common import constants, tools
import ldc.waveform.fastGB as fastGB
from ldc.lisa.noise import get_noise_model
import xarray as xr
from ldc.common.tools import compute_tdi_snr
from ldc.common.series import TDI, XYZ2AET


lisa_orbits = orbits.Orbits.type(dict({"nominal_arm_length":2.5e9*un.m,
                                       "initial_rotation":0*un.rad,
                                       "initial_position":0*un.rad,
                                       "orbit_type":"analytic"}))
Tobs = 31536000. * (2/3.)
dt = 5.

vgb = np.rec.fromarrays(['CD-30o11223', 4.15074703e-21, -0.28998647, 3.86029155,
                         0.00047261, 1.44003011e-19, 1.44687795, 1.7, 1.7],
                        names=['Name', 'Amplitude', 'EclipticLatitude',
                               'EclipticLongitude', 'Frequency',
                               'FrequencyDerivative', 'Inclination',
                               'InitialPhase', 'Polarization'])

GB = fastGB.FastGB(delta_t=dt, T=Tobs)
X, Y, Z = GB.get_fd_tdixyz(template=vgb, oversample=4, simulator='synthlisa')
XYZ = TDI([X, Y, Z], ["X", "Y", "Z"])
AET = XYZ2AET(XYZ)

fr = np.logspace(-5, 0, 10000)
model = 'SciRDv1'
noise = get_noise_model(model, fr, wd=1)

snr_xyz = compute_tdi_snr(XYZ, noise, AET=False, full_output=True)
snr_aet = compute_tdi_snr(AET, noise, AET=True, full_output=True)

print(snr_xyz["tot2"])
print(snr_xyz["tot2"])

plt.figure()
plt.plot(X.f, snr_xyz["cumsum"])
plt.plot(X.f, snr_aet["cumsum"])
