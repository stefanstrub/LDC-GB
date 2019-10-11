from LISAhdf5 import ParsUnits
from ldc.lisa import orbits

# Instance of a configuration class
lconfig = [('nominal_arm_length', 2.5e9, "m"),
           ('initial_rotation', 0, 'rad'), 
           ('initial_position', 0, 'rad')]
conf = ParsUnits(name='orbit_type', value='analytic')
for k,v,u in lconfig:
    conf.addPar(k,v,u)
    
X = orbits.Orbits.type(conf)

t = np.arange(0, 1000, 10)
P = X.compute_position(1, t) # Compute position for 1st spacecraft
