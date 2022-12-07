# distutils: sources = orbits.cc common.cc
# distutils: language = c++
# cython: language_level=3
cdef extern from "orbits.hpp":
    cdef cppclass AnalyticOrbits:
        AnalyticOrbits() except +
        AnalyticOrbits(double, double, double) except +
        void position_x(int, double*, double*, int)
        void position_y(int, double*, double*, int)
        void position_z(int, double*, double*, int) 
        void velocity_x(int, double*, double*, int)
        void velocity_y(int, double*, double*, int)
        void velocity_z(int, double*, double*, int) 
        void get_travel_time(int, int, double*, double*, int, int)
        void get_travel_time(int, int, double*, double*, double*, double*, double*, double*, double*, double*, int, int)


import numpy as np
cimport numpy as np
import lisaconstants as constants
ASTRONOMICAL_YEAR = constants.ASTRONOMICAL_YEAR
BARYCENTER_ANGULAR_VELOCITY = (2.0*np.pi)/(constants.SIDEREALYEAR_J2000DAY*24*60*60)

cdef class pyAnalyticOrbits:
    cdef AnalyticOrbits O
    cdef public double initial_rotation
    cdef public double initial_position
    cdef public double arm_length

    def __init__(self, armlength_meter=2.5e9, irotation_rad=0, ipos_rad=0):
        """ Set instrumental configuration for orbits. 
        """
        
        self.arm_length = armlength_meter
        self.initial_rotation = irotation_rad
        self.initial_position = ipos_rad
        self.O = AnalyticOrbits(self.arm_length, self.initial_position, self.initial_rotation)


    def compute_position(self, spacecraft_index, itime):
        """ Return (x_position, y_position, z_position) in the SSB. 
        
        Spacecraft index goes from 1 to 3. 
        """
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] time = itime.astype("double")
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] x = np.zeros(len(time))
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] y = np.zeros(len(time))
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] z = np.zeros(len(time))
        self.O.position_x(spacecraft_index, &time[0], &x[0], len(time))
        self.O.position_y(spacecraft_index, &time[0], &y[0], len(time))
        self.O.position_z(spacecraft_index, &time[0], &z[0], len(time))
        return np.array([x,y,z])

    def compute_velocity(self, spacecraft_index, itime):#np.ndarray[double, ndim=1, mode="c"] time not None):
        """ Return (x_velocity, y_velocity, z_velocity) in the SSB
        
        Spacecraft index goes from 1 to 3. 
        These are the derivatives w.r.t. time of the positions.
        """
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] time = itime.astype("double")
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] vx = np.zeros(len(time))
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] vy = np.zeros(len(time))
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] vz = np.zeros(len(time))
        self.O.velocity_x(spacecraft_index, &time[0], &vx[0], len(time))
        self.O.velocity_y(spacecraft_index, &time[0], &vy[0], len(time))
        self.O.velocity_z(spacecraft_index, &time[0], &vz[0], len(time))
        return np.array([vx,vy,vz])

    def compute_travel_time(self, emitter, receiver, receiver_itime, order=2,
                            position_receiver=None, position_emitter=None):
        """ Compute travel time between emitter and receiver at time receiver_time 
        
        Spacecraft index goes from 1 to 3. 
        """
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] receiver_time = receiver_itime.astype("double")
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] tt = np.zeros(len(receiver_time))

        cdef np.ndarray[np.double_t, ndim=1, mode="c"] x_emitter
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] y_emitter
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] z_emitter
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] x_receiver
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] y_receiver
        cdef np.ndarray[np.double_t, ndim=1, mode="c"] z_receiver
        if position_emitter is not None:
            x_emitter = np.ascontiguousarray(position_emitter[0,:])
            y_emitter = np.ascontiguousarray(position_emitter[1,:])
            z_emitter = np.ascontiguousarray(position_emitter[2,:])
            x_receiver = np.ascontiguousarray(position_receiver[0,:])
            y_receiver = np.ascontiguousarray(position_receiver[1,:])
            z_receiver = np.ascontiguousarray(position_receiver[2,:])
            self.O.get_travel_time(emitter, receiver, &receiver_time[0],
                                   &x_emitter[0], &y_emitter[0], &z_emitter[0],
                                   &x_receiver[0], &y_receiver[0], &z_receiver[0],
                                   &tt[0], len(tt), order)
        else:
            self.O.get_travel_time(emitter, receiver, &receiver_time[0], &tt[0], len(tt), order)
        return tt

    def alpha(self, itime):
        """ Return BARYCENTER_ANGULAR_VELOCITY * (init_time + time)
        """
        init_time = self.initial_position * ASTRONOMICAL_YEAR / (2 * np.pi)
        return BARYCENTER_ANGULAR_VELOCITY * (init_time + itime)
