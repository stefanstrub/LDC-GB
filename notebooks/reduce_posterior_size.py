
from fastkde import fastKDE
from getdist import plots, MCSamples
import numpy as np
import h5py
import corner
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd
import os

plt.style.use(['seaborn-ticks','seaborn-deep'])

parameters = [
    "Amplitude",
    "EclipticLatitude",
    "EclipticLongitude",
    "Frequency",
    "FrequencyDerivative",
    "Inclination",
    "InitialPhase",
    "Polarization",
]

# get current directory
path = os.getcwd()
 
# parent directory
parent = os.path.dirname(path)
# grandparent directory
grandparent = os.path.dirname(parent)

DATAPATH = "/home/stefan/LDC/Radler/data"
DATAPATH = grandparent+"/LDC/Radler/data"
SAVEPATH = grandparent+"/LDC/pictures/LDC1-3_v2"

# sangria_fn = DATAPATH + "/dgb-tdi.h5"
sangria_fn = DATAPATH + "/LDC1-3_VGB_v2.hdf5"
# sangria_fn = DATAPATH + "/LDC1-4_GB_v2.hdf5"
# sangria_fn = DATAPATH + "/LDC1-3_VGB_v2_FD_noiseless.hdf5"
fid = h5py.File(sangria_fn)
# get the source parameters
names = np.array(fid['H5LISA/GWSources/GalBinaries'])
params = [fid['H5LISA/GWSources/GalBinaries'][k] for k in names]
reduced_names = []
i = 0
for p in params:
    i += 1
    if p.shape:
        reduced_names.append(names[i-1])
params = [np.array(p) for p in params if p.shape]
cat = np.rec.fromarrays(params, names=list(reduced_names))

save_name = 'compare_LDC1-3'
indexes = np.argsort(cat['Frequency'])
cat_sorted = cat[indexes]

pGB_injected = []
pGB_injected_window = []
pGB_stacked = cat_sorted
for i in range(len(pGB_stacked)):
    pGBs = {}
    for parameter in parameters:
        pGBs[parameter] = pGB_stacked[parameter][i]
    pGB_injected_window.append(pGBs)
pGB_injected.append(pGB_injected_window)

### Dictionary of submitted results

VGBs = {'GW135962': {'BC':3, 'MM':0, 'BH':2},  
        'GW125313': {'BC':6, 'MM':9, 'BH':5},
        'GW181324': {'BC':8, 'MM':3, 'BH':1},
        'GW166667': {'BC':1, 'MM':7, 'BH':0},
        'GW194414': {'BC':5, 'MM':1, 'BH':6},
        'GW322061': {'BC':4, 'MM':6, 'BH':7},
        'GW351250': {'BC':0, 'MM':5, 'BH':9}, 
        'GW168350': {'BC':2, 'MM':4, 'BH':8},
        'GW622028': {'BC':7, 'MM':2, 'BH':4},
        'GW261301': {'BC':9, 'MM':8, 'BH':3}}

print('end')