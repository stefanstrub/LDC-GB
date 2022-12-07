import lisagwresponse
import lisaorbits
#from pytdi.michelson import X1, Y1, Z1
from pytdi.michelson import X1_ETA as X1
from pytdi.michelson import Y1_ETA as Y1
from pytdi.michelson import Z1_ETA as Z1

import h5py
from scipy.interpolate import InterpolatedUnivariateSpline
import numpy as np
from pytdi import Data

pGB = dict({"Frequency":0.0109966533684,
            "FrequencyDerivative": 2.95291919174e-14, 
            "EclipticLatitude":-0.150495923831, 
            "EclipticLongitude":4.60685356883, 
            "Amplitude":2.60186425363e-22,
            "Inclination":0.389613740033,
            "Polarization":0.605754423063,
            "InitialPhase":4.54588971721})

tobs = (365*24*60*60)

if 0:
    lisaorbits.KeplerianOrbits(L=2.5e9, dt=86400, size=400, tt_order=0, t0=-100).write('orbits.h5')

dt = 15
size = int(tobs/dt)

if 0:
    GB = lisagwresponse.GalacticBinary(A=pGB["Amplitude"], f=pGB["Frequency"], df=pGB["FrequencyDerivative"],
                                       phi0=pGB["InitialPhase"], iota=pGB["Inclination"], psi=pGB["Polarization"],
                                       orbits='orbits.h5',
                                       gw_beta=pGB["EclipticLatitude"], gw_lambda=pGB["EclipticLongitude"], dt=dt,
                                       size=size)
    GB.write('gw.h5')

with h5py.File('gw.h5', 'r') as fid:

    o = h5py.File("orbits.h5", 'r')
    data = {}
    traveltimes = {}
    for j, n in enumerate(["12", "21", "13", "31", "23", "32"]):
        r,s = int(n[0]), int(n[-1])
        data[f"eta_{r}{s}"] = np.array(fid[f'l_{n}'])
        traveltimes[f"d_{r}{s}"] = InterpolatedUnivariateSpline(o[f"tcb/t"], o[f"tcb/l_{n}"]['tt'])(fid["t"])
        
    fs = 1/dt
    data = Data(data, traveltimes, fs)
    build = X1.build(**data.args_nodoppler)
    data_x = build(data.measurements)
    build = Y1.build(**data.args_nodoppler)
    data_y = build(data.measurements)
    build = Z1.build(**data.args_nodoppler)
    data_z = build(data.measurements)

    tvec = np.array(fid['t'])
    
xyz = np.rec.fromarrays([tvec, data_x, data_y, data_z], names=["t", "X", "Y", "Z"])
np.save("XYZ.npy", xyz)
