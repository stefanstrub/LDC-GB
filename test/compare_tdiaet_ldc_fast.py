"""Compare hp and hc for SBBH between LDC and Fast waveform generation routines.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from ldc.lisa.orbits import Orbits
from ldc.waveform.lisabeta import FastBHB
from ldc.waveform.waveform import HpHc
from ldc.lisa.projection import ProjectedStrain
from ldc.common.series import TimeSeries, FrequencySeries, TDI


if __name__ == "__main__":

    # general params
    approx = "IMRPhenomD"
    source_type = "SBBH"
    tdi2 = False
    dt = 5 # s
    # t_max = 3600*24*365*0.1 # = 0.1yr
    t_max = 3600*24*365*0.5 # = 0.5yr
    t_min = 0

    # source parameters
    params =  dict({
         "Mass1":                 50.,
         "Spin1":                 0.0,
         "Mass2":                 40.,
         "Spin2":                 0.0,
         "EclipticLatitude":      1.7,
         "EclipticLongitude":     1.0471975511965976,
         "Inclination":           1.0471975511965976,
         "InitialFrequency":      1.2e-2,
         "InitialPhase":          0.7,
         "Polarization":          1.2,
         "Redshift":              2.0,
         "Distance":              15974.456786495544,
         'Cadence':               dt,
         'ObservationDuration':   t_max-t_min,
         'PolarAngleOfSpin1': 0.0,
         'PolarAngleOfSpin2': 0.0,
    })

    # define orbit.
    orbits = Orbits.type(
        dict({'orbit_type':'analytic',
              'nominal_arm_length':2.5e9,
              "initial_position": 0,
              "initial_rotation": 0})
    )

    # compare TDI AET for Fast
    fast = FastBHB(source_type, T=t_max, delta_t=dt, approx=approx, orbits=orbits)
    A_fast, E_fast, T_fast = fast.get_td_tdiaet(template=params, tdi2=tdi2)

    # compare TDI AET for LDC
    tvec = np.arange(t_min, t_max, dt)
    SOBHB = HpHc.type("my-sobhb", source_type, approx)
    SOBHB.set_param(params)
    projector = ProjectedStrain(orbits)
    yArm = projector.arm_response(t_min, t_max, dt, SOBHB.split())

    X_ldc = TimeSeries(projector.compute_tdi_x(tvec, tdi2=tdi2), t0=t_min, dt=dt)
    Y_ldc = TimeSeries(projector.compute_tdi_y(tvec, tdi2=tdi2), t0=t_min, dt=dt)
    Z_ldc = TimeSeries(projector.compute_tdi_z(tvec, tdi2=tdi2), t0=t_min, dt=dt)
    tdi = TDI(dict(zip(["X", "Y", "Z"], [X_ldc, Y_ldc, Z_ldc])))
    tdi.XYZ2AET()

    A_ldc, E_ldc, T_ldc = tdi.A, tdi.E, tdi.T

    # plot and compare
    fig, axs = plt.subplots(3, sharex=True)

    axs[0].plot(A_ldc.t, A_ldc, label='ldc')
    axs[0].plot(A_fast.t, A_fast, label='fast')
    axs[0].set_ylabel("A")
    axs[0].legend()

    axs[1].plot(E_ldc.t, E_ldc, label='ldc')
    axs[1].plot(E_fast.t, E_fast, label='fast')
    axs[1].set_ylabel("E")
    axs[1].legend()

    axs[2].plot(T_ldc.t, T_ldc, label='ldc')
    axs[2].plot(T_fast.t, T_fast, label='fast')
    axs[2].set_ylabel("T")
    axs[2].set_title("Time (s)")
    axs[2].legend()

    plt.show()
