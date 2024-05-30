import numpy as np
import matplotlib.pyplot as plt
from ldc.lisa.orbits import Orbits
from ldc.waveform.waveform import HpHc
from ldc.lisa.projection import ProjectedStrain
import ldc.io.yml as ymlio
import os
import h5py
import ldc.waveform.fastGB as FB
from ldc.common.series import TimeSeries
from ldc.common import constants

config = {"dt":5.0, "initial_position": 0, "initial_rotation": 0, 
          "nominal_arm_length": 2500000000, "orbit_type": 'analytic', 
          "t_max": 60*60*24*30*6, "t_min": 0, "travel_time_order": 1}

cat = np.array([(8837553, 4.688322047828039e-21, -0.7072278874089373,
                 3.5318990119515874,
                 0.01002696164199913, 7.0274293836061735e-15,
                 0.8456506362930373, 0.21987979316696454,
                 2.481152331028798)],
               dtype=[('Name', '<i8'), ('Amplitude', '<f8'),
                      ('EclipticLatitude', '<f8'),
                      ('EclipticLongitude', '<f8'),
                      ('Frequency', '<f8'), ('FrequencyDerivative', '<f8'),
                      ('Inclination', '<f8'), ('InitialPhase', '<f8'),
                      ('Polarization', '<f8')])

GB = HpHc.type("GB", "GB", "TD_fdot")
GB.set_param(cat)
GWs = GB.split()

globals().update(config)
trange = np.arange(t_min, t_max, dt)
orbits = Orbits.type(config)
P = ProjectedStrain(orbits)
yArm = P.arm_response(t_min, t_max, dt, GWs, tt_order=travel_time_order)

simple_tdi_X = P.compute_tdi_x(trange, tdi2=False)
simple_tdi_X = TimeSeries(simple_tdi_X, dt=dt)
simple_tdi_X = simple_tdi_X.ts.fft()

simple_tdi2_X = P.compute_tdi_x(trange, tdi2=True)
simple_tdi2_X = TimeSeries(simple_tdi2_X, dt=dt)
simple_tdi2_X = simple_tdi2_X.ts.fft()

GB = FB.FastGB(delta_t=dt, T=t_max)
pGB = dict(zip(cat.dtype.names, cat[0]))
X, Y, Z = GB.get_fd_tdixyz(template=pGB)
X2, Y2, Z2 = GB.get_fd_tdixyz(template=pGB, tdi2=True)


plt.figure()
i = 1
for x, s, g in zip([X, X2], [simple_tdi_X, simple_tdi2_X], ["1.5", "2.0"]):
    plt.subplot(2,2,i)
    plt.title("real part")
    plt.plot(x.f, x.real, alpha=0.5, label="fastGB")
    plt.plot(s.f, s.real, alpha=0.5, label=f"simple TDI {g}")
    plt.axis([pGB["Frequency"]-1e-6, pGB["Frequency"]+3e-6, None, None])
    plt.legend()
    plt.subplot(2,2,i+1)
    plt.title("imaginary part")
    plt.plot(x.f, x.imag, alpha=0.5, label="fastGB")
    plt.plot(s.f, s.imag, alpha=0.5, label=f"simple TDI {g}")
    plt.axis([pGB["Frequency"]-1e-6, pGB["Frequency"]+3e-6, None, None])
    i = i+2


