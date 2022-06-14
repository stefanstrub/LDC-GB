from ldc.lisa import orbits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

config = dict({"nominal_arm_length":2.5e9,#meter
               "initial_rotation":0,      #rad
               "initial_position":0,      #rad
               "orbit_type":"analytic"})
    
X = orbits.Orbits.type(config)
secondsperyear = 60*60*24*365.25


t = np.arange(0, secondsperyear, 100)
P = X.compute_position(1, t) # Compute position for 1st spacecraft
P2 = X.compute_position(2, t) # Compute position for 1st spacecraft
plt.figure()
plt.plot(t,P[0])
plt.plot(t,P2[0])
# plt.plot(t,P[2])
shift = secondsperyear/10
initial_position = 2*np.pi*(((shift)/secondsperyear)%1)
samples = 1
for i in range(samples):
    config = dict({"nominal_arm_length":2.5e9,#meter
                "initial_rotation":0,      #rad
                "initial_position":initial_position,      #rad
                "orbit_type":"analytic"})
        
    X = orbits.Orbits.type(config)

    P = X.compute_position(1, t) # Compute position for 1st spacecraft
    plt.plot(t+shift,P[0], color= cm.gray(i/samples))
# plt.plot(t+shift,P[1])
# plt.plot(t+shift,P[2])
plt.show()
print(P)