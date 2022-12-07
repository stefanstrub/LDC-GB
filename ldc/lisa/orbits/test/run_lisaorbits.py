from ldc.lisa.orbits import Orbits
import lisaorbits
import numpy as np
import matplotlib.pyplot as plt
import time
import lisaconstants

ASTRONOMICAL_YEAR = lisaconstants.ASTRONOMICAL_YEAR


ARM_LENGTH =  2.5E9
config = dict({"nominal_arm_length":ARM_LENGTH,#meter
               "initial_rotation":np.pi/3.,      #rad
               "initial_position":np.pi/5.,     #rad
               "orbit_type":"analytic"})
orbits = Orbits.type(config)

config2 = dict({"nominal_arm_length":ARM_LENGTH,#meter
                "initial_rotation":0,      #rad
                "initial_position":0,      #rad
                "orbit_type":"file",
                "filename":"orbits.h5"})
orbits2 = Orbits.type(config2)

dt = 500
tmax = 60*60*24*365
trange = np.arange(0, tmax, dt)
tinit = config["initial_position"] /(2 * np.pi / ASTRONOMICAL_YEAR)
sc = 2
sc1 = sc
sc2 = 3

if 1: # look at position
    t0 = time.time()
    xyz1 = orbits.compute_position(sc, trange)
    t1 = time.time()
    print(f"took {t1-t0} s")

    t0 = time.time()
    EAO = lisaorbits.EqualArmlengthOrbit(lambda1=config['initial_rotation'],
                                         tinit=-tinit, m_init1=-config['initial_rotation'])
    xyz2 = EAO.compute_spacecraft_position(sc, trange)
    t1 = time.time()
    print(f"took {t1-t0} s")

    t0 = time.time()
    xyz3 = orbits2.compute_position(sc, trange)
    t1 = time.time()
    print(f"took {t1-t0} s")
    
    
if 1: # plot it
    plt.figure()
    for i in range(3):
        plt.subplot(2,3,i+1)
        plt.plot(trange, xyz1[i,:], label="ldc")
        plt.plot(trange, xyz2[:,i], label="lisa orbits")
        #plt.plot(trange, xyz3[i,:], label="ldc from file")
        plt.subplot(2,3,3+i+1)
        plt.plot(trange, xyz1[i,:]-xyz2[:,i], label="diff")
        #plt.plot(trange, xyz1[i,:]-xyz3[i,:], label="diff")

if 0: # look at travel time
    order = 2
    t0 = time.time()
    tt1 = orbits.compute_travel_time(sc1, sc2, trange, order)
    t1 = time.time()
    print(f"took {t1-t0} s")

    t0 = time.time()
    EAO = lisaorbits.EqualArmlengthOrbit(tt_order=2, tt_niter=0)
    tt2 = EAO.compute_light_travel_time(sc1, sc2, trange)
    t1 = time.time()
    print(f"took {t1-t0} s")

    t0 = time.time()
    tt3 = orbits2.compute_travel_time(sc1, sc2, trange)
    t1 = time.time()
    print(f"took {t1-t0} s")

    
if 0: # plot it
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(trange, tt1, label="ldc")
    plt.plot(trange, tt2, label="lisa orbits")
    plt.plot(trange, tt3, label="ldc from file")

    plt.legend()
    plt.subplot(2,1,2)
    plt.plot(trange, tt1-tt2, label="diff")
    plt.plot(trange, tt1-tt3, label="diff")
    plt.legend()

if 0: # look at velocities
    t0 = time.time()
    xyz1 = orbits.compute_velocity(sc, trange)
    t1 = time.time()
    print(f"took {t1-t0} s")

    t0 = time.time()
    EAO = lisaorbits.EqualArmlengthOrbit()
    xyz2 = EAO.compute_spacecraft_velocity(sc, trange)
    t1 = time.time()
    print(f"took {t1-t0} s")

    t0 = time.time()
    xyz3 = orbits2.compute_velocity(sc, trange)
    t1 = time.time()
    print(f"took {t1-t0} s")

    
if 0: # plot it
    plt.figure()
    for i in range(3):
        plt.subplot(2,3,i+1)
        plt.plot(trange, xyz1[i,:], label="ldc")
        plt.plot(trange, xyz2[:,i], label="lisa orbits")
        plt.plot(trange, xyz3[i,:], label="ldc from file")

        plt.subplot(2,3,3+i+1)
        plt.plot(trange, xyz1[i,:]-xyz2[:,i], label="diff")
        plt.plot(trange, xyz1[i,:]-xyz3[i,:], label="diff")

