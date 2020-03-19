#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['axes.formatter.min_exponent'] = 2
from graph import TT_ORDER, ARM_LENGTH
from ldc.lisa.orbits import Orbits

def run_lisanode(case="Reforbits", duration=3600*24*365):
    flags = '-I../nodes  -I../lib -L../lib -lorbits -I../../../common/constants'
    os.system("lisanode run -o %s -f text --flags='%s' graph.py:%s -d %d"%(case, flags, case,duration))

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--pos', action="store_true", help="Plot orbits position")
    parser.add_argument('--vel', action="store_true", help="Plot orbits velocity")
    parser.add_argument('--tt', action="store_true", help="Plot travel times")

    args = parser.parse_args()

    cases = ["RefOrbits", "TestingLDCOrbits", "TestingTravelTimes_old", "TestingTravelTimes"]
    
    for case in cases:
        run_lisanode(case)

    config = dict({"nominal_arm_length":ARM_LENGTH,#meter
                   "initial_rotation":0,      #rad
                   "initial_position":0,      #rad
                   "orbit_type":"analytic"})
    orbits = Orbits.type(config)
        
    LC_links= dict({1:(3,2), 2:(1,3), 3:(2,1), 4:(2,3), 5:(3,1), 6:(1,2)})
    
    if args.tt:
        
        plt.figure() # plot tt
        plt.subplot(2,1,1)
        for link in range(1):
            for case in ["TestingLDCOrbits", "TestingTravelTimes"]:
                nt = np.loadtxt("%s/tt%d.txt"%(case, link+1))
                plt.plot(nt[:,0], nt[:,1], label=case)
            emitter, receiver = LC_links[link+1]
            trange = nt[:,0]
            tt_ldc = orbits.compute_travel_time(emitter, receiver, trange,TT_ORDER)
            plt.plot(nt[:,0], tt_ldc, label="LDC")
        plt.ylabel("Travel time for link 1 [s]")
        plt.legend()
        plt.subplot(2,1,2)
        for link in range(1):
            nts = list()
            for case in ["TestingLDCOrbits", "TestingTravelTimes"]:
                nt = np.loadtxt("%s/tt%d.txt"%(case, link+1))
                nts.append(nt)
                
            plt.plot(nts[0][:,0], 1e9*(nts[1][:,1]-nts[0][:,1]), label="LDC - LISANode")
            emitter, receiver = LC_links[link+1]
            trange = nt[:,0]
            tt_ldc = orbits.compute_travel_time(emitter, receiver, trange,TT_ORDER)
            plt.plot(nts[0][:,0], 1e9*(nts[1][:,1]-tt_ldc), label="LISANode - LDC")
            plt.plot(nts[0][:,0], 1e9*(nts[0][:,1]-tt_ldc), label="LDC - LDC")
        plt.xlabel("Time [s]")
        plt.ylabel("Tt diff for link 1 [ns]")
        plt.legend()
        #plt.savefig("tt_order%s.png"%TT_ORDER)

    if args.pos or args.vel:
        plt.figure() # plot x
        plt.subplot(2,1,1)
        for link in range(1):
            for p, color in zip(["x", "y", "z"], ["b", "orange", "g"]):
                for case in ["TestingTravelTimes_old", "TestingLDCOrbits"]:
                    suffix = "v" if args.vel else ""
                    nt = np.loadtxt("%s/Orbit%s%s%d.txt"%(case, suffix, p, link))
                    label= p if case=="TestingTravelTimes" else ""
                    plt.plot(nt[:,0], nt[:,1], label=label, color=color)
                    #print(nt)
        plt.legend()
        if args.pos:
            plt.ylabel("Positions for sc 1 [m]")
        else:
            plt.ylabel("Velocity for sc 1 [m]")
            
        #plt.figure() # plot x-x
        plt.subplot(2,1,2)
        for link in range(1):
            for p, color in zip(["x", "y", "z"], ["b", "orange", "g"]):
                nts = list()
                for case in ["TestingTravelTimes_old", "TestingLDCOrbits"]:
                    suffix = "v" if args.vel else ""
                    nts.append(np.loadtxt("%s/Orbit%s%s%d.txt"%(case, suffix, p, link)))
                label="LDC-LISANode" if p=="x" else ""
                plt.plot(nts[0][:,0], np.abs(nts[1][:,1]-nts[0][:,1]), label=label, color=color)

        plt.legend()
        plt.xlabel("Time [s]")
        if args.pos:
            plt.ylabel("Pos abs diff for link 1 [m]")
        else:
            plt.ylabel("Vel abs diff for link 1 [m]")
        #plt.savefig("pos.png")


