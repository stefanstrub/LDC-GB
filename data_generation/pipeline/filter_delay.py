import ldc.io.hdf5 as hdfio
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.interpolate import InterpolatedUnivariateSpline as spline
import os


duration = 24*60*60
window = signal.gaussian(3*duration, std=2000)
arr = window

if 0:
    hdfio.save_array("input.h5", arr, name="strain")

# get delay of lisanode @dt=5 where 
# passband_freq: 0.40*(1/5)
# stopband_freq: 0.48*(1/5)
elif 0:#se: 

    out, jk = hdfio.load_array("Response_20200721101931.h5", name="out") #14.61
    out2, jk = hdfio.load_array("Response_20200721103641.h5", name="out") # 14.535
    window2 = signal.gaussian(duration/5, std=2000/15)

    interpolator = spline(out[:,0], out[:,1])
    out_i = interpolator(np.arange(0, duration, 5))

    interpolator = spline(out2[:,0], out2[:,1])
    out_i_2 = interpolator(np.arange(0, duration, 5))
    
    interpolator = spline(out[:,0], out[:,1])
    out_i_shift = interpolator(np.arange(0, duration, 5) - 0.075)

    interpolator = spline(np.arange(0, duration, 1/3), window)
    true_i = interpolator(np.arange(0, duration, 5))
    
    
    plt.figure()
    plt.subplot(211)
    plt.plot(np.arange(0, duration, 1/3), window, label="true 3Hz")
    plt.plot(np.arange(0, duration, 5), window2, label="true 0.2Hz")
    plt.plot(np.arange(0, duration, 5), true_i, label="interpolated 3Hz")
    plt.plot(np.arange(0, duration, 5), out_i, label="out")
    plt.plot(np.arange(0, duration, 5), out_i_shift, label="out shifted")
    plt.plot(np.arange(0, duration, 5), out_i_2, label="out2")
    plt.legend()
    plt.subplot(212)
    plt.plot(np.arange(0, duration, 5), true_i-out_i, label="diff")
    plt.plot(np.arange(0, duration, 5), true_i-out_i_shift, label="diff shited")
    plt.plot(np.arange(0, duration, 5), true_i-out_i_2, label="diff 2")
    plt.legend()


# get filter coeff to have delay=10 
# passband_freq: 0.40*(1/3.45)
# stopband_freq: 0.48*(1/3.45)   
elif 1:
    os.system("rm out.h5")
    os.system("lisanode run filter_graph.py:Response -d 86410 -o out.h5")
    out2, jk = hdfio.load_array("out.h5", name="out") # 14.535
    
    out, jk = hdfio.load_array("Response_20200721101931.h5", name="out") #14.61
    interpolator = spline(out[:,0], out[:,1])
    out_i = interpolator(np.arange(0, duration, 5) - 0.075)
    
    out_i_shift = out2[2:, 1]
    
    interpolator = spline(np.arange(0, duration, 1/3), window)
    true_i = interpolator(np.arange(0, duration, 5))
    
    plt.figure()
    plt.subplot(211)
    plt.plot(np.arange(0, duration, 1/3), window, label="true 3Hz")
    plt.plot(np.arange(0, duration, 5), true_i, label="interpolated 3Hz")
    plt.plot(np.arange(0, duration, 5), out_i, label="out")
    plt.plot(np.arange(0, duration, 5), out_i_shift, label="out shifted")
    plt.legend()
    plt.subplot(212)
    plt.plot(np.arange(0, duration, 5), true_i-out_i, label="diff")
    plt.plot(np.arange(0, duration, 5), true_i-out_i_shift, label="diff shited")
    plt.legend()
    plt.show()

    
if 0: # get delay of lisacode @dt=5
    from RunSimuLC2 import RunSimuLC2
    from lc import ConfigureInstrument
    from ldc.lisa.orbits import Orbits
    from ldc.waveform.waveform import HpHc
    from ldc.lisa.projection import ProjectedStrain
    import os
    tr = np.arange(0, duration, 5)

    pMBHB = dict({'EclipticLatitude': 0.312414, #"radian"),
                  'EclipticLongitude': -2.75291,# "radian"),
                  'CoalescenceTime': 28086000.0,# 's'),
                  'Distance':  9.14450149011798,# 'Gpc'),
                  'InitialAzimuthalAngleL': 3.9,# 'radian'),
                  'InitialPolarAngleL': 2.3535, #'radian'),
                  'Mass1': 132628.202,# "SolarMass"),
                  'Mass2': 30997.2481,# "SolarMass"),
                  'PhaseAtCoalescence':  3.8, #'Radian'),
                  'PolarAngleOfSpin1': 0.0,#'Radian'),
                  'PolarAngleOfSpin2': 0.0,#'Radian'),
                  'Redshift': 1.2687,# 'dimensionless'),
                  'Spin1': 0.9481998052314212, #'MassSquared'),
                  'Spin2': 0.9871324769575264, #'MassSquared'),
                  'Cadence': 5.,
                  'ObservationDuration':duration})#, 's') })
    
    
    GW = HpHc.type("MBHB", "MBHB", "IMRPhenomD")
    pdict = dict(zip(cat.dtype.names, cat[0]))
    GW.set_param(pMBHB)
    GW.hp = 50000*signal.gaussian(duration/5, std=200)#2000/16)
    GW.hc = 0*GW.hp 
    GW.t = tr
    GW.interp_hphc(tr)

    hphcfile = "hphc.hdf5";     os.system("rm %s"%hphcfile)
    GW.to_file(hphcfile) # save hp hc to file

    config = dict({"nominal_arm_length":2.5e9,#meter
                   "initial_rotation":0,      #rad
                   "initial_position":0,      #rad
                   "orbit_type":"analytic"})
    orbits = Orbits.type(config)
    P = ProjectedStrain(orbits)
    yArm = P.arm_response(0, duration, 5, [GW], tt_order=1, precomputed=True)
    simple_tdi_X = P.compute_tdi_x(tr)
    
    ConfigureInstrument(hphcfile, scriptPath="",
                        options=[], TDI="X,Y,Z",
                        duration=duration, timeStep=5, orbits="MLDC_Orbits")
    
    RunSimuLC2(hphcfile, scriptPath="", 
               options=[],
               debug=True, verbose=True,
               path2LISACode="/usr/local/bin/",
               NoNoise=True, NoGW=False, seed=-1)#, source_index=0)

    X = np.loadtxt("TmpLC2_hphc-TDI.txt")
    interpolator = spline(X[:,0], X[:,1])
    shifted = interpolator(np.arange(0, duration, 5)-251.75)

    plt.figure()
    plt.subplot(211)
    plt.plot(tr, simple_tdi_X, label='simple tdi')   
    plt.plot(X[:,0], X[:,1], label='lisacode')
    plt.plot(tr, shifted, label='lisacode shitfed')
    plt.legend()
    plt.subplot(212)
    plt.plot(tr, simple_tdi_X-shifted, label='diff shifted')   
    plt.plot(tr, simple_tdi_X-shifted2, label='diff shifted')   
    plt.legend()
    
