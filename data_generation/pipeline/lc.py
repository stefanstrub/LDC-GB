import os
from LISAhdf5 import LISAhdf5,ParsUnits
from RunSimuLC2 import RunSimuLC2
import numpy as np

MOSAlist=["1","1s","2","2s","3","3s"]
distMOSA = {"1":"2s","2":"3s","3":"1s","1s":"3","2s":"1","3s":"2"}
joinMOSA = {"1":"1s","2":"2s","3":"3s","1s":"1","2s":"2","3s":"3"}


# cfgpu = dict({"dt":5, "initial_position":0, "initial_rotation":0,
#               "interp_order":3, "nominal_arm_length":2.5e9, "orbit_type":"analytic",
#               "t_min":0, "t_max":60*60*24})

def ConfigureInstrument(fNh5,\
                        scriptPath="",\
                        options=[],\
                        TDI="X,Y,Z",\
                        duration=125829105.0,\
                        timeStep=15.0,\
                        orbits="MLDC_Orbits",\
                        psdfile="None"):
    """
    Configure the instrument
    @param fNh5 is the filename of the hdf5
    @param scriptPath is the script path
    @param options is the list of options from the command line
    @param TDI is string describing the TDI generator
    @param duration is the duration (sec)
    @param timeStep is the time step (sec)
    @param orbits is the string describing the orbits ("LISACode_Orbits")
    """

    LH = LISAhdf5(fNh5)
    LH.addHistory(scriptPath, options, [fNh5])

    #### Configure duration and time steps
    pS1 = ParsUnits()
    pS1.addPar("Cadence",timeStep,"second")
    pS1.addPar("StepPhysic",timeStep,"second")
    pS1.addPar("Duration",duration,"second")
    pS1.addPar("UpdateShotNoise","Off","string")
    pS1.addPar("InterpolationNoises",7,"Lagrange")

    LH.addSimulation(pS1)


    #### Configure orbits
    pO = ParsUnits()
    pO.addPar("TimeOffset",0.,"second")
    pO.addPar("InitialPosition",0.,"radian")
    pO.addPar("InitialRotation",0.,"radian")
    pO.addPar("Armlength",2.5e9,"meter")
    #pO.addPar("OrbitApproximation","Eccentric","string")
    pO.addPar("OrbitApproximation","Conventional","string")
    pO.addPar("OrderDelay","1","string")
    LH.addLISADataSource("Orbits",orbits,pO,overwrite=True)


    #### Configure particular pieces of the instrument
    for xMOSA in MOSAlist:

        pL = ParsUnits()
        pL.addPar("Power",2.,"watt")
        pL.addPar("Wavelength",1064.e-9,"meter")
        LH.addLISADataSource("laser"+xMOSA,"Laser",pL,overwrite=True)

        pT = ParsUnits()
        pT.addPar("Diameter",0.3,"meter")
        LH.addLISADataSource("tel"+xMOSA,"Telescope",pT,overwrite=True)

        pOB = ParsUnits()
        pOB.addPar("OpticalBenchType","Std2002","string")
        pOB.addPar("FilterType","Weak","None")
        Inputs = {"beam1":"laser"+distMOSA[xMOSA]+",tel"+distMOSA[xMOSA]+",tel"+xMOSA}
        Inputs.update({"beam2":"laser"+xMOSA})
        LH.addLISADataModel("sci"+xMOSA,"OpticalBenchPhasemeter",pOB,InputsData=Inputs,OutputsData={},overwrite=True)


    #### Configure TDI
    pTDI = ParsUnits()
    pTDI.addPar("TDIInterpolationData",31,"Lagrange")
    pTDI.addPar("TDIInterpolationDelay",0,"Linear")
    pTDI.addPar("TDIGenerator",TDI,"string")
    if psdfile == "None":
        LH.addPreProcess(pTDI)
    else:
        d = np.loadtxt(psdfile,skiprows=2)
        LH.addPreProcess(pTDI,PSDdata=d)

def source_to_file(filename, GW):
    from LISAhdf5 import LISAhdf5, ParsUnits
    h5 = LISAhdf5(filename)
    units = GW.source_parameters.copy()
    for k, v in GW.units.items():
        units[k] = v
    pu = ParsUnits(pars_i=GW.source_parameters, units_i=units)
    h5.addSource(GW.source_name, pu,
                 overwrite=True, hphcData=np.vstack([GW.t, GW.hp, GW.hc]).T)
    
        
def run_lisacode(GWs, t_min, t_max, dt_source, dt_obs=None):

    if dt_obs is None:
        dt_obs = dt_source
    
    hphcfile = "hphc.hdf5"
    os.system("rm %s"%hphcfile)
    for GW1 in GWs:
        GW1.compute_hphc_td(np.arange(t_min, t_max, dt_source), set_attr=True)
        source_to_file(hphcfile, GW1)


    output_file = "TmpLC2_hphc-TDI.txt"
    os.system("rm %s"%output_file)
    
    ConfigureInstrument(hphcfile, scriptPath="",
                        options=[], TDI="X,Y,Z",
                        duration=t_max, timeStep=dt_obs, orbits="MLDC_Orbits")
    
        
    RunSimuLC2(hphcfile, scriptPath="", 
               options=[],
               debug=True, verbose=True,
               path2LISACode="/usr/local/bin/",
               NoNoise=True, NoGW=False, seed=-1)#, source_index=0)

    X = np.loadtxt("TmpLC2_hphc-TDI.txt")
    return X

