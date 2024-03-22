import math
import numpy as np
import lisaconstants as constants

from ldc.lisa.orbits import Orbits
from ldc.lisa.noise import simple_snr

YEAR = constants.SIDEREALYEAR_J2000DAY*24*60*60

def get_default_orbits(orbits=None):
    """ Provide default LISA orbits
    """
    if orbits is not None:
        if isinstance(orbits, dict):
            orbits = Orbits.type(orbits)
    else:
        arm_length = 2.5e9 # m
        init_rotation = 0 # rad
        init_position = 0 # rad
        orbits = Orbits.type(dict({'orbit_type':'analytic',
                                   'nominal_arm_length':2.5e9, # m
                                   "initial_position": 0, 
                                   "initial_rotation": 0, }))
        return orbits

def get_buffersize(T, f0, ampl, oversample):
    """Get array dimension needed to compute TDI.
    """
    Acut = simple_snr(f0,ampl,years=T/YEAR)
    mult = 8
    if((T/YEAR) <= 8.0): mult = 8
    if((T/YEAR) <= 4.0): mult = 4
    if((T/YEAR) <= 2.0): mult = 2
    if((T/YEAR) <= 1.0): mult = 1
    N = 32*mult
    if(f0 > 0.001): N = 64*mult
    if(f0 > 0.01):  N = 256*mult
    if(f0 > 0.03):  N = 512*mult
    if(f0 > 0.1):   N = 1024*mult
    M = int(math.pow(2.0,1 + int(np.log(Acut)/np.log(2.0))))
    if(M > 8192):
        M = 8192
    M = N = max(M,N)
    N *= oversample
    return N
