import matplotlib.pyplot as plt
from astropy import units as un
import numpy as np
from ldc.waveform.lisabeta import FastBHB
from ldc.lisa import orbits
from ldc.lisa.projection import ProjectedStrain
from ldc.waveform.waveform import HpHc
from ldc.common.series import TimeSeries, FrequencySeries
from ldc.common import tools

# params_3 = dict(zip(['EclipticLatitude', 'EclipticLongitude',
#                      'PolarAngleOfSpin1', 'PolarAngleOfSpin2',
#                      'Spin1', 'Spin2', 'Mass1', 'Mass2', 'CoalescenceTime',
#                      'PhaseAtCoalescence', 'InitialPolarAngleL',
#                      'InitialAzimuthalAngleL', 'Redshift', 'Distance',
#                      'ObservationDuration', 'Cadence'],
#                     [0.42199769, 3.2910523,
#                      2.16808415, 2.47809227,
#                      0.892382, 0.570879, 298995., 287572.,
#                      3409125.34065294, 3.07758541,
#                      1.25648348, 1.65710114,
#                      2.83462, 24363.40064515,
#                      31558149.7635456, 5.]))

lisa_orbits = orbits.Orbits.type(dict({"nominal_arm_length":2.5e9*un.m,
                                       "initial_rotation":0*un.rad,
                                       "initial_position":0*un.rad,
                                       "orbit_type":"analytic"}))
Tobs = 31536000. * (2/3.)
dt = 5.


def AET(X,Y,Z):
    return ((Z - X)/np.sqrt(2.0),
            (X - 2.0*Y + Z)/np.sqrt(6.0),
            (X + Y + Z)/np.sqrt(3.0))


def semi_fast_tdi(orbits, pMBHB, t_max, dt):
    hphc = HpHc.type("MBHB", "MBHB", "IMRPhenomD")
    hphc.set_param(pMBHB)
    P = ProjectedStrain(orbits)    
    yArm = P.arm_response(0, t_max, dt, [hphc], tt_order=0)
    X = P.compute_tdi_x(np.arange(0, t_max, dt))
    Z = P.compute_tdi_z(np.arange(0, t_max, dt))
    Y = P.compute_tdi_y(np.arange(0, t_max, dt))
    return TimeSeries(X, dt=dt), TimeSeries(Y, dt=dt), TimeSeries(Z, dt=dt)

params_2 = dict(zip(['EclipticLatitude','EclipticLongitude',
                     'PolarAngleOfSpin1', 'PolarAngleOfSpin2',
                     'Spin1', 'Spin2', 'Mass1', 'Mass2',
                     'CoalescenceTime', 'PhaseAtCoalescence',
                     'InitialPolarAngleL', 'InitialAzimuthalAngleL',
                     'Redshift', 'Distance', 'ObservationDuration', 'Cadence'],
                    [1.28888303, 2.10716542, 0.62478351, 1.60719115, 0.953675,
                     0.954511, 294021., 260563., 20426222.35914461,
                     3.43517638, 2.06463032, 4.28787516, 3.58897,
                     32.30154329*1000, 31558149.7635456, 5.]))

params_2["ObservationDuration"] = Tobs
params_2["Cadence"] = dt

params_1 = params_2.copy()
#params_1["Mass1"] *= (1+ params_2["Redshift"])
#params_1["Mass2"] *= (1+ params_2["Redshift"])
psi, incl = tools.aziPolAngleL2PsiIncl(params_2["EclipticLatitude"],
                                       params_2["EclipticLongitude"],
                                       params_2['InitialPolarAngleL'],
                                       params_2['InitialAzimuthalAngleL'])
params_1['Polarization'] = psi
params_1['Inclination'] = incl
params_1["Spin1"] = params_2['Spin1']*np.cos(params_2['PolarAngleOfSpin1'])
params_1["Spin2"] = params_2['Spin2']*np.cos(params_2['PolarAngleOfSpin2'])
params_1.pop('PolarAngleOfSpin1')
params_1.pop('PolarAngleOfSpin2')
params_1.pop('InitialPolarAngleL')
params_1.pop('InitialAzimuthalAngleL')
params_1.pop('ObservationDuration')
params_1.pop('Cadence')
params_1.pop('Redshift')

if 1:

    # LDC
    X2, Y2, Z2 = semi_fast_tdi(lisa_orbits, params_2, Tobs, dt)
    A2, E2, T2 = AET(X2,Y2,Z2)

    # lisabeta
    phenomD = FastBHB(approx="IMRPhenomD", T=Tobs, delta_t=dt, orbits=lisa_orbits, bbh_type='mbhb')
    A1, E1, T1 = phenomD.get_td_tdiaet(template=params_2)

    plt.figure()
    plt.subplot(311)
    A1.plot(label="lisabeta")
    A2.plot(alpha=0.5, label='LDC')
    plt.axis([5900+2.042e7, 6500+2.042e7, -2.5e-19, 2.5e-19])
    plt.subplot(312)
    E1.plot(label="lisabeta")
    E2.plot(alpha=0.5, label='LDC')
    plt.axis([5900+2.042e7, 6500+2.042e7, -2.5e-19, 2.5e-19])
    plt.subplot(313)
    T1.plot(label="lisabeta")
    T2.plot(alpha=0.5, label='LDC')
    plt.axis([5900+2.042e7, 6500+2.042e7, -2.5e-19, 2.5e-19])
    plt.legend()
    plt.show()
    
