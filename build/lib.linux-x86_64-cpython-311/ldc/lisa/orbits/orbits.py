""" Compute LISA s/c dynamics based on analytic orbits or orbits from file """

from abc import ABC, abstractmethod
from itertools import permutations
import numpy as np
from astropy import units as un
from _orbits import pyAnalyticOrbits
import lisaconstants
import scipy.interpolate
import h5py

#pylint:disable=C0103

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
        if 'nominal_arm_length' in config.keys():
            self.arm_length = config['nominal_arm_length']
        else:
            self.arm_length = 2.5e9
            self.eccentricity = self.arm_length/(2*np.sqrt(3)*AU_IN_M)
        self.tt = None

    def check_units(self, config):
        """ Convert parameters in expected units.
        """
        self.units = dict({'nominal_arm_length':'m'})
        for k,u in self.units.items():
            if k in config.keys() and isinstance(config[k], un.Quantity):
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

    @classmethod
    def type(cls, config):
        """ Switches to an analytic or a numerical orbits class.
        """
        if config.get('orbit_type') == 'file':
            return OrbitsFromFile(config)
        if config.get('orbit_type') in ["analytic", "equalarmlength",
                                        'equalarmlength-orbits']:
            return AnalyticOrbits(config)
        raise ValueError("Invalid orbit type.")

class AnalyticOrbits(pyAnalyticOrbits, Orbits):
    """ Computes equal armlength orbits.
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
    """

    def __init__(self, config, read_t0=True):
        """
        >>> X = Orbits.type({"orbit_type":'file', 'filename':'orbits.h5'})
        """
        super().__init__(config)
        self.filename = config['filename']
        self.pos_interp = [None]*self.number_of_spacecraft
        self.vel_interp = [None]*self.number_of_spacecraft
        self.tt_interp = [None]*self.number_of_arms
        with h5py.File(self.filename, 'r') as fid:
            self.deprecated = not 'size' in fid.attrs
            if read_t0:
                self.t0 = fid.attrs["t0"]
            else:
                self.t0 = 0
            self.dt = fid.attrs["dt"]
            if self.deprecated:
                self.size = fid["tcb"]["t"].size
            else:
                self.size = fid.attrs["size"]
            self.t = np.arange(self.t0, self.t0+self.dt*self.size, self.dt)


    def compute_position(self, spacecraft_index, time):
        """ Computes spacecraft position from an orbits file

        >>> X = Orbits.type({"orbit_type":'file', 'filename':'orbits.h5'})
        >>> X.compute_position(1,np.arange(100)).shape
        (3, 100)
        """
        if self.pos_interp[spacecraft_index-1] is None:
            if self.deprecated:
                with h5py.File(self.filename, 'r') as fid:
                    self.pos_interp[spacecraft_index-1] = [
                        scipy.interpolate.InterpolatedUnivariateSpline(
                            self.t, fid["tcb"][f'sc_{spacecraft_index}'][k])
                        for j,k in enumerate(["x", "y", "z"])]
            else:
                with h5py.File(self.filename, 'r') as fid:
                    self.pos_interp[spacecraft_index-1] = [
                        scipy.interpolate.InterpolatedUnivariateSpline(
                            self.t, fid["tcb"]["x"][:,spacecraft_index-1, j])
                        for j,k in enumerate(["x", "y", "z"])]
        x = self.pos_interp[spacecraft_index-1][0](time)
        y = self.pos_interp[spacecraft_index-1][1](time)
        z = self.pos_interp[spacecraft_index-1][2](time)
        return np.array([x,y,z])

    def compute_velocity(self, spacecraft_index, time):
        """ Computes spacecraft velocity from an orbits file

        >>> X = Orbits.type({"orbit_type":'file', 'filename':'orbits.h5'})
        >>> X.compute_velocity(1,np.arange(100)).shape
        (3, 100)
        """
        if self.vel_interp[spacecraft_index-1] is None:
            if self.deprecated:
                with h5py.File(self.filename, 'r') as fid:
                    self.vel_interp[spacecraft_index-1] = [
                        scipy.interpolate.InterpolatedUnivariateSpline(
                            self.t, fid["tcb"][f'sc_{spacecraft_index}'][k])
                        for j,k in enumerate(["vx", "vy", "vz"])]
            else:
                with h5py.File(self.filename, 'r') as fid:
                    self.vel_interp[spacecraft_index-1] = [
                        scipy.interpolate.InterpolatedUnivariateSpline(
                            self.t, fid["tcb"]["v"][:,spacecraft_index-1, j])
                        for j,k in enumerate(["vx", "vy", "vz"])]
        vx = self.vel_interp[spacecraft_index-1][0](time)
        vy = self.vel_interp[spacecraft_index-1][1](time)
        vz = self.vel_interp[spacecraft_index-1][2](time)
        return np.array([vx,vy,vz])

    def compute_travel_time(self, emitter, receiver, receiver_time, **kwargs):
        """ Computes travel time for a single arm from an orbits file

        >>> X = Orbits.type({"orbit_type":'file', 'filename':'orbits.h5'})
        >>> X.compute_travel_time(1, 2, np.arange(100)).shape
        (100,)
        """
        itt = self.get_pairs().index((receiver, emitter))
        LINKS = [12, 23, 31, 13, 32, 21] #receiver/emitter
        idx = LINKS.index(int(f"{receiver}{emitter}"))
        if self.tt_interp[itt] is None:
            if self.deprecated:
                with h5py.File(self.filename, 'r') as fid:
                    self.tt_interp[itt] = \
                        scipy.interpolate.InterpolatedUnivariateSpline(
                            self.t, fid["tcb"][f'l_{receiver}{emitter}']['tt'])
            else:
                with h5py.File(self.filename, 'r') as fid:
                    self.tt_interp[itt] = \
                        scipy.interpolate.InterpolatedUnivariateSpline(
                            self.t, fid["tcb"]['ltt'][:,idx])
        tt = self.tt_interp[itt](receiver_time)
        return tt


if __name__ == "__main__":
    import doctest
    config = dict({"nominal_arm_length":2.5e9*un.m,
                   "initial_rotation":0*un.rad,
                   "initial_position":0*un.rad,
                   "orbit_type":"analytic"})
    import lisaorbits
    my_orbits = lisaorbits.EqualArmlengthOrbits(dt=100, size=100)
    my_orbits.write(f'orbits.h5')

    doctest.testmod()
