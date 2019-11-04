#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import matplotlib.pyplot as plt
import numpy as np

def run_lisanode(case="Reforbits", duration=3600*24*365):
    flags = '-I../nodes  -I../lib -L../lib -L. -lorbits -lhdf5_serial -lhdf5_cpp'
    os.system("lisanode run -o %s --flags='%s' graph.py:%s -d %d"%(case, flags, case,duration))

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--orbits', action="store_true", help="Plot orbits")
    parser.add_argument('--tt', action="store_true", help="Plot travel times")

    args = parser.parse_args()

    for case in ["RefOrbits", "TestingLDCOrbits"]:
        run_lisanode(case)

    
    if args.tt:
        plt.figure() # plot tt
        for link in range(1):
            for case in ["RefOrbits", "TestingLDCOrbits"]:
                nt = np.loadtxt("%s/tt%d.txt"%(case, link+1))
                plt.plot(nt[:,0], nt[:,1], label=case)

    if args.orbits:
        plt.figure() # plot x
        for link in range(1):
            for p, color in zip(["x", "y", "z"], ["b", "orange", "g"]):
                for case in ["RefOrbits", "TestingLDCOrbits"]:
                    nt = np.loadtxt("%s/Orbit%s%d.txt"%(case, p, link))
                    plt.plot(nt[:,0], nt[:,1], label="%s-%s"%(case,p), color=color)
                    #print(nt)
        plt.legend()
        
        plt.figure() # plot x-x
        for link in range(1):
            for p, color in zip(["x", "y", "z"], ["b", "orange", "g"]):
                nts = list()
                for case in ["RefOrbits", "TestingLDCOrbits"]:
                    nts.append(np.loadtxt("%s/Orbit%s%d.txt"%(case, p, link)))
                plt.plot(nts[0][:,0], np.abs(nts[0][:,1]-nts[1][:,1]), label="%s"%(p), color=color)

