#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['axes.formatter.min_exponent'] = 2
from graph import TT_ORDER

def run_lisanode(case="Reforbits", duration=3600*24*365):
    flags = '-I../nodes  -I../lib -L../lib -lorbits'
    os.system("lisanode run -o %s --flags='%s' graph.py:%s -d %d"%(case, flags, case,duration))

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--orbits', action="store_true", help="Plot orbits")
    parser.add_argument('--tt', action="store_true", help="Plot travel times")

    args = parser.parse_args()

    cases = ["RefOrbits", "TestingLDCOrbits", "TestingTravelTimes"]
    
    for case in cases:
        run_lisanode(case)

    
    if args.tt:
        plt.figure() # plot tt
        plt.subplot(2,1,1)
        for link in range(1):
            for case in ["RefOrbits", "TestingTravelTimes"]:
                nt = np.loadtxt("%s/tt%d.txt"%(case, link+1))
                plt.plot(nt[:,0], nt[:,1], label=case)
        plt.ylabel("Travel time for link 1 [s]")
        plt.legend()
        plt.subplot(2,1,2)
        for link in range(1):
            nts = list()
            for case in ["RefOrbits", "TestingTravelTimes"]:
                nt = np.loadtxt("%s/tt%d.txt"%(case, link+1))
                nts.append(nt)
                
            plt.plot(nts[0][:,0], 1e9*(nts[1][:,1]-nts[0][:,1]), label="LDC - LISANode")
        plt.xlabel("Time [s]")
        plt.ylabel("Tt diff for link 1 [ns]")
        plt.legend()
        plt.savefig("tt_order%s.png"%TT_ORDER)

    if args.orbits:
        plt.figure() # plot x
        plt.subplot(2,1,1)
        for link in range(1):
            for p, color in zip(["x", "y", "z"], ["b", "orange", "g"]):
                for case in ["TestingTravelTimes", "TestingLDCOrbits"]:
                    nt = np.loadtxt("%s/Orbit%s%d.txt"%(case, p, link))
                    label= p if case=="TestingTravelTimes" else ""
                    plt.plot(nt[:,0], nt[:,1], label=label, color=color)
                    #print(nt)
        plt.legend()
        plt.ylabel("Positions for link 1 [m]")
        
        #plt.figure() # plot x-x
        plt.subplot(2,1,2)
        for link in range(1):
            for p, color in zip(["x", "y", "z"], ["b", "orange", "g"]):
                nts = list()
                for case in ["TestingTravelTimes", "TestingLDCOrbits"]:
                    nts.append(np.loadtxt("%s/Orbit%s%d.txt"%(case, p, link)))
                label="LDC-LISANode" if p=="x" else ""
                plt.plot(nts[0][:,0], np.abs(nts[1][:,1]-nts[0][:,1]), label=label, color=color)

        plt.legend()
        plt.xlabel("Position [m]")
        plt.ylabel("Pos abs diff for link 1 [m]")
        plt.savefig("pos.png")
