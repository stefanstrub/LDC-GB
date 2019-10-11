""" Compute LISA s/c dynamics based on analytic orbits or orbits from file """

from abc import ABC, abstractmethod
import numpy as np
from itertools import permutations
from LISAhdf5 import ParsUnits 
import LISAConstants as LC 
from _orbits import pyAnalyticOrbits
from ldc.common import constants

C = constants.Nature
AU_IN_M = C.ASTRONOMICALUNIT_METER
# ASTRONOMICAL_YEAR = C.SIDEREALYEAR_J2000DAY*60*60*24
# POST_NEWTONIAN_CONSTANT = 1.0
# HALF_SCHWARZSCHILD_RADIUS = 1.47664E3
# BARYCENTER_ANGULAR_VELOCITY = 2. * np.pi / ASTRONOMICAL_YEAR

class Orbits(ABC):
    """ This abstract base class is the gateway to the LISA set of orbits functions """
    
    def __init__(self, config):
        """Initialize orbits from a configuration, given as ParsUnits object.
        
        Configuration should include:
        - orbit type in ['analytic', 'file']
        - nominal_arm_length

        Orbits are given in SSB reference frame. 
        >>> X = Orbits.type(config)
        """
        if not isinstance(config, ParsUnits):
            raise "Configuration should be a ParsUnit object"
        
        self.orbit_type = config.get('orbit_type')
        self.reference_frame = 'SSB'
        self.number_of_spacecraft = 3
        self.number_of_arms = 6
        self.arm_length = config.getConvert('nominal_arm_length', LC.convDistance, "m")
        self.eccentricity = self.arm_length/(2*np.sqrt(3)*AU_IN_M)
        self.tt = None

    @abstractmethod
    def compute_travel_time(self, emitter, receiver, receiver_time, order):
        """ Returns travel times for an emitter-receiver pair.
        
        Receiver_time can be a collection.
        """
        pass
    
    @abstractmethod
    def compute_position(self, spacecraft_index, time):
        """ Returns spacecraft positions.
        
        Time can be a collection. 
        """
        pass
    
    @abstractmethod
    def compute_velocity(self, spacecraft_index, time):
        """ Returns spacecraft velocities.

        Time can be a collection. 
        """
        pass

    def get_pairs(self):
        """ Returns arm pairs of indices, in the order used as a convention.
        
        >>> X = Orbits.type(config)
        >>> X.get_pairs()
        [(1, 2), (1, 3), (2, 1), (2, 3), (3, 1), (3, 2)]
        """
        perms = list(permutations(range(1,self.number_of_spacecraft+1), 2))
        return perms

    @classmethod
    def type(cls, config):
        """ Switches to an analytic or a numerical orbits class. """
        if config.get('orbit_type') == "file":
            return OrbitsFromFile(config)
        if config.get('orbit_type') == "analytic":
            return AnalyticOrbits(config)
        raise ValueError("Invalid orbit type.")

class AnalyticOrbits(pyAnalyticOrbits, Orbits):
    """ Computes Keplerian orbits. 
    """
    
    def __init__(self, config):
        pyAnalyticOrbits.__init__(self, config)
        Orbits.__init__(self, config)

class OrbitsFromFile(Orbits):
    """ Computes the orbits from an orbits file """
    
    def __init__(self, config):
        super().__init__(config)
        self.filename = config.get_parameter('filename')

    def compute_position(self, spacecraft_index, time):
        """ Computes spacecraft position from an orbits file """
        return "To be implemented."

    def compute_velocity(self, spacecraft_index, time):
        """ Computes spacecraft velocity from an orbits file """
        return "To be implemented."

    def compute_travel_time(self, emitter, receiver, receiver_time, order):
        """ Computes travel time for a single arm from an orbits file """
        return "To be implemented."




if __name__ == "__main__":
    from LISAhdf5 import ParsUnits
    import doctest

    lconfig = [('nominal_arm_length', 2.5e9, "m"), 
               ('initial_rotation', 0, 'rad'), 
               ('initial_position', 0, 'rad')]
    config = ParsUnits(name='orbit_type', value='analytic')
    for k,v,u in lconfig:
        config.addPar(k,v,u)
    doctest.testmod()
