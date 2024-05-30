import matplotlib.pyplot as plt
import numpy as np
import xml.dom.minidom
from ldc.lisa.orbits import Orbits
import os

ARM_LENGTH =  2.5E9

def get_param(xml_content, name):
    params = xml_content.getElementsByTagName("Param")
    value =  [p.childNodes[0].data for p in params if p.getAttribute("Name")==name][0]
    value = value.strip()
    return value

def set_param(xml_content, name, value):
    params = xml_content.getElementsByTagName("Param")
    for p in params:
        if p.getAttribute("Name")==name:
            p.childNodes[0].data = value
    
def run_lisacode(case="", tt_order=0, duration=3600*24*365):
    """ case in : 
    - default: lisacode config: orbits and travel time
    """
    filename = "test_orbits.xml"
    basename = case
    xml_content = xml.dom.minidom.parse(filename)
    params = xml_content.getElementsByTagName("Param")
    armlength = get_param(xml_content, "Armlength")
    set_param(xml_content, "OrbitApproximation", "Conventional")
    set_param(xml_content, "OrderDelay", str(tt_order))
    set_param(xml_content, "Armlength", str(ARM_LENGTH))
    xmlfile = "test_orbits_tmp.xml"
    fid = open(xmlfile, "w")
    xml_content.writexml(fid)
    fid.close()
    os.system("LC2Orbits %s %s"%(basename,xmlfile))


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--pos', action="store_true", help="Plot position")
    parser.add_argument('--vel', action="store_true", help="Plot velocity")
    parser.add_argument('--tt', action="store_true", help="Plot travel times")
    parser.add_argument('--tt-order', type=int, default=2,
                        help="Travel time order in [0,1,2]")
    args = parser.parse_args()

    cases = ["testorb"]

    config = dict({"nominal_arm_length":ARM_LENGTH,#meter
                   "initial_rotation":0,      #rad
                   "initial_position":0,      #rad
                   "orbit_type":"analytic"})
    orbits = Orbits.type(config)
    for case in cases:
        run_lisacode(case, tt_order=args.tt_order)

    LC_links= dict({1:(3,2), 2:(1,3), 3:(2,1), 4:(2,3), 5:(3,1), 6:(1,2)})

    if args.tt:
        for case in cases:
            for link in [1]:#, 2, 3, 4, 5, 6]:
                plt.figure() # plot tt
                plt.subplot(2,1,1)
                nt = np.loadtxt("%s-Arm.txt"%(case))
                trange = nt[:,0]
                emitter, receiver = LC_links[link]
                tt_ldc = orbits.compute_travel_time(emitter, receiver, trange,
                                                    args.tt_order) 
                plt.plot(nt[:,0], -1*nt[:,link], label="LISACode")
                plt.plot(nt[:,0], tt_ldc, label="LDC")
                plt.ylabel("Travel time for link %d [s]"%link)
                plt.legend()
                plt.subplot(2,1,2)
                plt.plot(nt[:,0], 1e9*(tt_ldc+nt[:,link]), label="LDC - LISACode")
                plt.xlabel("Time [s]")
                plt.ylabel("Tt diff for link 1 [ns]")
                plt.legend()
                #plt.savefig("tt_order%s.png"%TT_ORDER)

    if args.pos or args.vel:
        for case in cases:
            for sc in [1]:
                plt.figure() # plot tt
                plt.subplot(2,1,1)
                suffix = "Pos" if args.pos else "Vel"
                nt = np.loadtxt("%s-%s.txt"%(case, suffix))
                trange = nt[:,0]
                if args.pos:
                    xyz_ldc = orbits.compute_position(sc, trange)
                elif args.vel:
                    xyz_ldc = orbits.compute_velocity(sc, trange)
                for j,(p,color) in enumerate(zip(["x", "y", "z"], ["b", "orange", "g"])):
                    plt.plot(nt[:,0], nt[:,sc+j], label=p, color=color)
                    plt.plot(nt[:,0], xyz_ldc[j,:])
                plt.ylabel("Positions for sc %d [m]"%sc)
                plt.legend()
                plt.subplot(2,1,2)
                for j,(p,color) in enumerate(zip(["x", "y", "z"], ["b", "orange", "g"])):
                    plt.plot(nt[:,0],xyz_ldc[j,:]-nt[:,sc+j] , label=p, color=color)
                plt.xlabel("Time [s]")
                if args.pos:
                    plt.ylabel("Pos abs diff for sc %d [m]"%sc)
                else:
                    plt.ylabel("Vel abs diff for link 1 [m]")
                plt.legend()
                #plt.savefig("tt_order%s.png"%TT_ORDER)
