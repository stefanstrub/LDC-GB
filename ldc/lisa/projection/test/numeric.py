import numpy as np
from ldc.lisa.orbits import Orbits
from ldc.waveform.waveform import NumericHpHc, HpHc
from ldc.lisa.projection import ProjectedStrain
import matplotlib.pyplot as plt

config = {'nominal_arm_length':2.5e9, #"m", 
          'initial_rotation':0, #'rad', 
          'initial_position':0,#'rad')]
          'orbit_type':'analytic'} 

t_min = 0
t_max = 60*60*24*30*6
dt = 15
trange = np.arange(t_min, t_max, dt)

pGB = dict({'Amplitude': 1.07345e-22,#, "strain"), 
            'EclipticLatitude': 0.312414,#, "radian"),
            'EclipticLongitude': -2.75291,# "radian"), 
            'Frequency': 0.00135962,# "Hz"),
            'FrequencyDerivative': 8.94581279e-19,# "Hz^2"), 
            'Inclination': 0.523599 ,# "radian"), 
            'InitialPhase': 3.0581565,# "radian"), 
            'Polarization': 3.5621656})#,# "radian")})

GB = HpHc.type("my-galactic-binary", "GB", "TD_fdot")
GB.set_param(pGB)

orbits = Orbits.type(config)
P = ProjectedStrain(orbits)
yArm = P.arm_response(t_min, t_max, dt, [GB], tt_order=1)
tdi_X = P.compute_tdi_x(trange)

extended_trange = np.arange(t_min-500, t_max+500, dt/5)
hp, hc = GB.compute_hphc_td(extended_trange)
hphc_num = NumericHpHc(extended_trange, hp, hc, pGB['EclipticLongitude'], pGB['EclipticLatitude'])
yArm = P.arm_response(t_min, t_max, dt, [hphc_num], tt_order=1)
tdi_X_num = P.compute_tdi_x(trange)

plt.figure()
plt.subplot(211)
plt.plot(trange, tdi_X, label='through GB type')
plt.plot(trange, tdi_X_num, label="numeric")
plt.legend()
plt.axis([674000, 676500, None, None])
plt.subplot(212)
plt.plot(trange, tdi_X-tdi_X_num)
plt.axis([674000, 676500, None, None])
