from ldc.lisa import orbits
import numpy as np

config = dict({"nominal_arm_length":2.5e9,#meter
               "initial_rotation":0,      #rad
               "initial_position":0,      #rad
               "orbit_type":"analytic"})
    
X = orbits.Orbits.type(config)

t = np.arange(0, 1000, 10)
P = X.compute_position(1, t) # Compute position for 1st spacecraft
