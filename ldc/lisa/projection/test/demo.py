import numpy as np
from ldc.lisa.orbits import Orbits
from ldc.waveform.waveform import HpHc
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

#compare to fastGB
import ldc.waveform.fastGB as FB
GB = FB.FastGB(delta_t=dt, T=t_max) # in seconds
X, Y, Z = GB.get_td_tdixyz(template=pGB, oversample=4)

plt.figure()
plt.plot(trange, tdi_X, label='projected strain to TDI')
plt.plot(X.t, X, label="fastGB")
plt.legend()
plt.axis([674000, 676500, None, None])


from lisainstrument import Instrument
gwfile = 'gb-strain.h5'
P.to_file(gwfile, fmt="sim")
instrument = Instrument(aafilter=None, dt=dt, size=1000, gws=gwfile)
instrument.disable_all_noises()
instrument.disable_dopplers()
instrument.simulate()

