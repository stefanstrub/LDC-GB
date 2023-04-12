import numpy as np
import matplotlib.pyplot as plt
from ldc.lisa.orbits import Orbits
from ldc.common.series import TDI
from ldc.lisa.projection import ProjectedStrain
from ldc.waveform.waveform import HpHc
from ldc.common.series import TimeSeries
import ldc.io.hdf5 as h5io

param = {'EclipticLatitude': -0.30300442294174235,
         'EclipticLongitude': 1.2925183861048521,
         'PolarAngleOfSpin1': 1.2031361791056812,
         'PolarAngleOfSpin2': 2.097303543065685,
         'Spin1': 0.747377,
         'Spin2': 0.8388,
         'Mass1': 1323277.47932,
         'Mass2': 612485.5060299999,
         'CoalescenceTime': 11526944.921879262,
         'PhaseAtCoalescence': 1.2201968860015653,
         'InitialPolarAngleL': 2.6919824500032945,
         'InitialAzimuthalAngleL': 1.808398497592109,
         'Redshift': 1.73941,
         'Distance': 13470.983558972537,
         'ObservationDuration': 31558149.763545603,
         'Cadence': 3.0}

config = {"initial_position": 0, "initial_rotation": 0, 
          "nominal_arm_length": 2500000000, "orbit_type": 'analytic'}
lisa_orbits = Orbits.type(config)

cat, descr = h5io.load_array('/home/maude/data/LDC/sangria/mbhb-big/mbhb1-cat.h5')

param = dict(zip(cat.dtype.names, cat[()]))

t_min = 0
t_max = 365*24*60*60
dt = 5

Proj = ProjectedStrain(lisa_orbits)
GW = HpHc.type("debug", "MBHB", "IMRPhenomD")
GW.set_param(param)
Proj.from_file('/home/maude/data/LDC/sangria/mbhb-big/mbhb1-y.h5')
#yArm = Proj.arm_response(t_min, t_max, dt, [GW])

if 0:
    tdi_x = TimeSeries(Proj.compute_tdi_x(np.arange(t_min, t_max, dt)), dt=dt)
    tdi_y = TimeSeries(Proj.compute_tdi_y(np.arange(t_min, t_max, dt)), dt=dt)
    tdi_z = TimeSeries(Proj.compute_tdi_z(np.arange(t_min, t_max, dt)), dt=dt)
    tdi = TDI(dict({"X":tdi_x, "Y":tdi_y, "Z":tdi_z}))
    #tdi.XYZ2AET()

if 1:
    tdi_lisanode = TDI.load("/home/maude/data/LDC/sangria/mbhb-big/new/mbhb1-lisanode-noisefree-tdi.h5")
    dt = tdi_lisanode["t"][1]-tdi_lisanode["t"][0]
    #subset = tdi_lisanode.sel(t=tdi.t, method="nearest")
    tdi_x = TimeSeries(Proj.compute_tdi_x(tdi_lisanode.t.values), dt=dt, t0=tdi_lisanode["t"][0])
    tdi_y = TimeSeries(Proj.compute_tdi_y(tdi_lisanode.t.values), dt=dt, t0=tdi_lisanode["t"][0])
    tdi_z = TimeSeries(Proj.compute_tdi_z(tdi_lisanode.t.values), dt=dt, t0=tdi_lisanode["t"][0])
    tdi1 = TDI(dict({"X":tdi_x, "Y":tdi_y, "Z":tdi_z}))

if 1:
    tdi_, descr = h5io.load_array("/home/maude/data/LDC/sangria/1.7/mbhb-tdi-XYZ.h5")
    dt = tdi_["t"][1]-tdi_["t"][0]
    tdi_sangria = TDI(dict({"X":TimeSeries(tdi_["X"], dt=dt),
                            "Y":TimeSeries(tdi_["Y"], dt=dt),
                            "Z":TimeSeries(tdi_["Z"], dt=dt)}))
    tdi_x = TimeSeries(Proj.compute_tdi_x(tdi_sangria.t.values), dt=dt, t0=tdi_sangria["t"][0])
    tdi_y = TimeSeries(Proj.compute_tdi_y(tdi_sangria.t.values), dt=dt, t0=tdi_sangria["t"][0])
    tdi_z = TimeSeries(Proj.compute_tdi_z(tdi_sangria.t.values), dt=dt, t0=tdi_sangria["t"][0])
    tdi2 = TDI(dict({"X":tdi_x, "Y":tdi_y, "Z":tdi_z}))


if 0:
    merger_time = np.arange(param["CoalescenceTime"]-1000, param["CoalescenceTime"]+500, dt)
    merger = tdi_lisanode.sel(t=merger_time, method="nearest")
    ldc_merger = tdi.sel(t=merger_time, method="nearest")

    plt.figure()
    for time_delay in [0.29,0.3,0.31]:#,30]:
        tdi_candidate = merger.assign_coords(t=(merger_time-time_delay))
        s = tdi_candidate["X"].interp(t=merger_time)#, method="nearest")
        plt.plot(merger_time, ldc_merger["X"].values-s.values, label='LDC-LISANode')
    plt.axis([param["CoalescenceTime"]-1000, param["CoalescenceTime"]+500, -4e-20, 4e-20])
    stop

if 0:
    plt.figure()
    plt.plot(tdi.t, tdi["X"], label='X from LDC')
    plt.plot(tdi_lisanode.t, tdi_lisanode["X"], label='X from LISANode')
    plt.plot(tdi.t, tdi["X"].values-tdi_lisanode["X"].values, label='LDC-LISANode', color='grey')
    plt.axis([param["CoalescenceTime"]-1000, param["CoalescenceTime"]+500, -4e-19, 4e-19])
    plt.legend()
    plt.savefig("comp5.png")

if 0:
    plt.figure()
    plt.plot(tdi.t, tdi["X"], label='X from LDC')
    plt.plot(tdi_sangria.t, tdi_sangria["X"], label='X from sangria')
    plt.plot(tdi.t, tdi["X"].values-subset_sangria["X"].values, label='LDC-sangria', color='grey')
    plt.axis([param["CoalescenceTime"]-1000, param["CoalescenceTime"]+500, -4e-19, 4e-19])
    plt.legend()
    plt.savefig("comp2.png")
if 1:
    plt.figure()
    plt.plot(tdi1.t, tdi1["X"].values-tdi_lisanode["X"].values, label='LDC-LISANode')
    plt.plot(tdi2.t, tdi2["X"].values-tdi_sangria["X"].values, label='LDC-sangria')
    plt.axis([param["CoalescenceTime"]-1000, param["CoalescenceTime"]+500, -1e-20, 1e-20])
    plt.legend()
    plt.savefig("comp5-zoom.png")
#/work/SC/lisa/lejeune/sangria2/run1/mbhb1-lisanode-noisefree-tdi.h5
