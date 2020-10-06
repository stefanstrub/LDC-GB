import numpy as np
import matplotlib.pyplot as plt
import glob
from ldc.waveform.waveform.hphc import HpHc
from ldc.waveform.source import SourceMaker
import ldc.io.hdf5 as h5io

cfg = dict({"approximant": "IMRPhenomD",
            "catalogs": '/home/maude/soft/scienceperfs_evalfom/SOA/data/Q3d_complete',
            "coalescence_time":(0.1, 0.95),
            "mass_ratio":( 1, 10),
            "nsource": 200, 
            "seed": 1234, 
            "source_type": "MBHB",
            "spin1":(0.5, 0.99),
            "spin2":(0.5, 0.99)})

cfg["redshifted_mass"] = True
cfg["non_precessing"] = False

source_maker = SourceMaker.type(cfg["source_type"],
                                cfg["approximant"],
                                catalogs=sorted(glob.glob(cfg["catalogs"]+"/*")))
cat, units = source_maker.choose_from_catalog(**cfg)
hphc = HpHc.type('test', cfg["source_type"], cfg["approximant"])

ldc_cat, jk = h5io.load_array("/home/maude/data/LDC/sangria/1.6/sangria_with_units.h5",
                              name="sky/mbhb/cat")


ldc_cat = np.load("/home/maude/data/LDC/sangria/blind-1.6/mbhb.npy") 
indices = []
for i,s in enumerate(ldc_cat):
    loc = np.where(s["CoalescenceTime"]==cat["CoalescenceTime"])
    print(loc)
    indices.append(loc[0][0])
#[159, 29, 87, 0, 153, 156]
#"indices": [41, 125, 176, 150, 190, 88, 67, 10, 73, 131, 26, 97, 24, 63, 45],

cat = cat[np.array(indices)]

if 1:
    plt.figure()
    for i,s in enumerate(cat):
        print(i)
        p = dict(zip(cat.dtype.names, s))
        hphc.set_param(p, units=units)
        frq, amp, phase = hphc._IMRPhenomD_waveform()
        plt.loglog(frq, amp, label="%d"%i)
    plt.legend(loc="upper left")
