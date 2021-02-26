from ldc.lisa import orbits
import numpy as np
import matplotlib.pyplot as plt

config = dict({"nominal_arm_length":2.5e9,#meter
               "initial_rotation":0,      #rad
               "initial_position":0,      #rad
               "orbit_type":"analytic"})
    
X = orbits.Orbits.type(config)
secondsperyear = 60*60*24*365.25


t = np.arange(0, secondsperyear, 100)
P = X.compute_position(1, t) # Compute position for 1st spacecraft
plt.figure()
plt.plot(t,P[0])
plt.plot(t,P[1])
plt.plot(t,P[2])
config = dict({"nominal_arm_length":2.5e9,#meter
               "initial_rotation":secondsperyear/2,      #rad
               "initial_position":0,      #rad
               "orbit_type":"analytic"})
    
X = orbits.Orbits.type(config)
secondsperyear = 60*60*24*365.25


t = np.arange(secondsperyear/4, secondsperyear, 100)
P = X.compute_position(1, t) # Compute position for 1st spacecraft
plt.plot(t,P[0])
plt.plot(t,P[1])
plt.plot(t,P[2])
plt.show()
print(P)