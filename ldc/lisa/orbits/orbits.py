""" Compute LISA s/c dynamics based on analytic orbits or orbits from file """

from abc import ABC, abstractmethod
import numpy as np
from astropy import units as un
from itertools import permutations
from _orbits import pyAnalyticOrbits
import lisaconstants
import scipy.interpolate
import h5py

AU_IN_M = lisaconstants.ASTRONOMICAL_UNIT

class Orbits(ABC):
    """ This abstract base class is the gateway to the LISA set of orbits functions """

    def __init__(self, config):
        """Initialize orbits from a configuration.

        Configuration should include:
        - orbit type in ['analytic', 'file']
        - nominal_arm_length

        Orbits are given in SSB reference frame.
        >>> X = Orbits.type(config)
        """
        config = self.check_units(config)
        self.orbit_type = config.get('orbit_type')
        self.reference_frame = 'SSB'
        self.number_of_spacecraft = 3
        self.number_of_arms = 6
        self.arm_length = config['nominal_arm_length']
        self.eccentricity = self.arm_length/(2*np.sqrt(3)*AU_IN_M)
        self.tt = None

    def check_units(self, config):
        """ Convert parameters in expected units.
        """
        self.units = dict({'nominal_arm_length':'m'})
        for k,u in self.units.items():
            if isinstance(config[k], un.Quantity):
                config[k] = config[k].to(un.Unit(u))
                config[k] = config[k].value
        return config


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
        """Returns arm pairs of indices (receiver-emitter), in the order used
        as a convention.

        >>> X = Orbits.type(config)
        >>> X.get_pairs()
        [(1, 2), (1, 3), (2, 1), (2, 3), (3, 1), (3, 2)]

        """
        perms = list(permutations(range(1,self.number_of_spacecraft+1), 2))
        return perms

    # TODO PB: base version identification from installed lisaorbits version
    # directly instead of reading an extra argument.
    @classmethod
    def type(cls, config):
        """ Switches to an analytic or a numerical orbits class. 
        
        :param config: configuration for orbit
        :type config: dict with structure     
        
        config = {
            'nominal_arm_length': 2.5e9, 
            'orbit_type':'file', 
            'filename': <hdf5_orbit_filename>,
            'lisaorbit_old_version': <True/False>,
        }
        """
        if config.get('orbit_type') == 'file':
            if config.get('lisaorbit_old_version'):
                print("Selected lisaorbit 1.0+ config...")
                return OrbitsFromFileDeprecated(config)
            else:
                print("Falling back to lisaorbit 2.0+ config...")
                return OrbitsFromFile(config)
        if config.get('orbit_type') == "analytic":
            return AnalyticOrbits(config)
        raise ValueError("Invalid orbit type.")

class AnalyticOrbits(pyAnalyticOrbits, Orbits):
    """ Computes Keplerian orbits.
    """

    def __init__(self, config):
        """ Set instrumental configuration for orbits.
        """
        Orbits.__init__(self, config)
        config = self.check_units(config)
        arm_length_meter = config["nominal_arm_length"]
        irot_rad = config["initial_rotation"]
        ipos_rad = config["initial_position"]
        pyAnalyticOrbits.__init__(self, arm_length_meter, irot_rad, ipos_rad)

    def check_units(self, config):
        """ Convert parameters in expected units.
        """
        self.units = dict({'nominal_arm_length':'m',
                           'initial_rotation':'rad',
                           'initial_position':'rad'})
        for k,u in self.units.items():
            if isinstance(config[k], un.Quantity): 
                config[k] = config[k].to(un.Unit(u))
                config[k] = config[k].value
        return config


class OrbitsFromFile(Orbits):
    """ Computes the orbits from an orbits file from lisaorbits.
    This class is to be used with version 2.0 of lisaorbits since
    the HDF5 file format has changed since version 1.X.X"""

    def __init__(self, config, read_t0=True):
        super().__init__(config)
        self.filename = config['filename']
        self.pos_interp = [None]*self.number_of_spacecraft
        self.vel_interp = [None]*self.number_of_spacecraft
        self.tt_interp = [None]*self.number_of_arms
        with h5py.File(self.filename, 'r') as fid:
            #self.t = fid["tcb"]["t"][:]
            if read_t0:
                self.t0 = fid.attrs["t0"]
            else:
                self.t0 = 0
            self.dt = fid.attrs["dt"]
            self.size = fid.attrs["size"]
            self.t = np.arange(self.t0, self.t0+self.dt*self.size, self.dt)
            
            
    def compute_position(self, spacecraft_index, time):
        """ Computes spacecraft position from an orbits file """
        if self.pos_interp[spacecraft_index-1] is None:
            with h5py.File(self.filename, 'r') as fid:
                self.pos_interp[spacecraft_index-1] = [
                    scipy.interpolate.InterpolatedUnivariateSpline(self.t, #fid["tcb"]["t"],
                                                                   #fid["tcb"]['sc_%d'%spacecraft_index][k])
                                                                   fid["tcb"]["x"][:,spacecraft_index-1, j])
                    for j,k in enumerate(["x", "y", "z"])]
        x = self.pos_interp[spacecraft_index-1][0](time)
        y = self.pos_interp[spacecraft_index-1][1](time)
        z = self.pos_interp[spacecraft_index-1][2](time)
        return np.array([x,y,z])

    def compute_velocity(self, spacecraft_index, time):
        """ Computes spacecraft velocity from an orbits file """
        if self.vel_interp[spacecraft_index-1] is None:
            with h5py.File(self.filename, 'r') as fid:
                self.vel_interp[spacecraft_index-1] = [
                    scipy.interpolate.InterpolatedUnivariateSpline(self.t,#fid["tcb"]["t"],
                                                                   #fid["tcb"]['sc_%d'%spacecraft_index][k])
                                                                   fid["tcb"]["v"][:,spacecraft_index-1, j])
                    
                    for j,k in enumerate(["vx", "vy", "vz"])]
        vx = self.vel_interp[spacecraft_index-1][0](time)
        vy = self.vel_interp[spacecraft_index-1][1](time)
        vz = self.vel_interp[spacecraft_index-1][2](time)
        return np.array([vx,vy,vz])

    def compute_travel_time(self, emitter, receiver, receiver_time, **kwargs):
        """ Computes travel time for a single arm from an orbits file """
        itt = self.get_pairs().index((receiver, emitter))
        LINKS = [12, 23, 31, 13, 32, 21] #receiver/emitter
        idx = LINKS.index(int(f"{receiver}{emitter}"))
        if self.tt_interp[itt] is None:
            with h5py.File(self.filename, 'r') as fid:
                self.tt_interp[itt] = \
                    scipy.interpolate.InterpolatedUnivariateSpline(self.t,#fid["tcb"]["t"],
                                                                   #fid["tcb"]['l_%d%d'%(receiver, emitter)]['tt'])
                                                                   fid["tcb"]['ltt'][:,idx])
                
        tt = self.tt_interp[itt](receiver_time)
        return tt


class OrbitsFromFileDeprecated(Orbits):
    """ Computes the orbits from an orbits file from lisaorbits.
    This class is to be used with version 1.0 of lisaorbits since
    the HDF5 file format has changed since version 1.X.X"""

    def __init__(self, config, read_t0=True):
        super().__init__(config)
        self.filename = config['filename']
        self.pos_interp = [None]*self.number_of_spacecraft
        self.vel_interp = [None]*self.number_of_spacecraft
        self.tt_interp = [None]*self.number_of_arms
        with h5py.File(self.filename, 'r') as fid:
            #self.t = fid["tcb"]["t"][:]
            if read_t0:
                self.t0 = fid.attrs["t0"]
            else:
                self.t0 = 0
            self.dt = fid.attrs["dt"]
            self.size = fid["tcb"]["t"].size
            self.t = np.arange(self.t0, self.t0+self.dt*self.size, self.dt)

    def compute_position(self, spacecraft_index, time):
        """ Computes spacecraft position from an orbits file """
        if self.pos_interp[spacecraft_index-1] is None:
            with h5py.File(self.filename, 'r') as fid:
                self.pos_interp[spacecraft_index-1] = [
                    scipy.interpolate.InterpolatedUnivariateSpline(
                        self.t, fid["tcb"][f'sc_{spacecraft_index}'][k])
                    for j,k in enumerate(["x", "y", "z"])]
        x = self.pos_interp[spacecraft_index-1][0](time)
        y = self.pos_interp[spacecraft_index-1][1](time)
        z = self.pos_interp[spacecraft_index-1][2](time)
        return np.array([x,y,z])

    def compute_velocity(self, spacecraft_index, time):
        """ Computes spacecraft velocity from an orbits file """
        if self.vel_interp[spacecraft_index-1] is None:
            with h5py.File(self.filename, 'r') as fid:
                self.vel_interp[spacecraft_index-1] = [
                    scipy.interpolate.InterpolatedUnivariateSpline(
                        self.t, fid["tcb"][f'sc_{spacecraft_index}'][k])
                    for j,k in enumerate(["vx", "vy", "vz"])]
        vx = self.vel_interp[spacecraft_index-1][0](time)
        vy = self.vel_interp[spacecraft_index-1][1](time)
        vz = self.vel_interp[spacecraft_index-1][2](time)
        return np.array([vx,vy,vz])

    def compute_travel_time(self, emitter, receiver, receiver_time, **kwargs):
        """ Computes travel time for a single arm from an orbits file """
        itt = self.get_pairs().index((receiver, emitter))
        LINKS = [12, 23, 31, 13, 32, 21] #receiver/emitter
        idx = LINKS.index(int(f"{receiver}{emitter}"))
        if self.tt_interp[itt] is None:
            with h5py.File(self.filename, 'r') as fid:
                self.tt_interp[itt] = \
                    scipy.interpolate.InterpolatedUnivariateSpline(
                        self.t, fid["tcb"][f'l_{receiver}{emitter}']['tt'])
        tt = self.tt_interp[itt](receiver_time)
        return tt


if __name__ == "__main__":
    import doctest
    config = dict({"nominal_arm_length":2.5e9*un.m,
                   "initial_rotation":0*un.rad,
                   "initial_position":0*un.rad,
                   "orbit_type":"analytic"})
    doctest.testmod()
