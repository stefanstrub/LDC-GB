
from fastkde import fastKDE
from getdist import plots, MCSamples
import numpy as np
import h5py
import corner
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd
import os
import glob

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
pGB_stacked = cat_sorted
for i in range(len(pGB_stacked)):
    pGBs = {}
    for parameter in parameters:
        pGBs[parameter] = pGB_stacked[parameter][i]
    pGB_injected.append(pGBs)

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

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

data_set_names = ['ETH','BC']
names = ['EclipticLongitude','EclipticLatitude','Frequency','FrequencyDerivative','Inclination','Amplitude']
file_paths = ['/home/stefan/Repositories/ldc1_evaluation_data/submission/ETH_LDC1-3_v2/*.csv','/home/stefan/Repositories/ldc1_evaluation_data/submission/Barcelona/ldc1-3_posteriors_ice-csic-ieec_v2/*.csv']
data = {}
for j in range(2):
    data_set = []
    fls = glob.glob(file_paths[j])
    print (fls)

    for fl in fls:
        df = pd.read_csv(fl)[:10**4]
        df['Inclination'] = np.cos(df['Inclination'].values)
        df['EclipticLatitude'] = np.sin(df['EclipticLatitude'].values)
        df['FrequencyDerivative'] = np.log10(df['FrequencyDerivative'].values)
        df['Amplitude'] = np.log10(df['Amplitude'].values)
        df = df.rename(columns = {'# Frequency':'Frequency'})
        df['Frequency'] *= 1000
        df = df.drop(labels= ['InitialPhase', 'Polarization'], axis=1)
        df = df.reindex(columns = names)
        # df['Frequency'] *= 1000
        data_set.append(df)

    frequencies = np.zeros(len(pGB_injected))
    for i in range(len(pGB_injected)):
        frequencies[i] = pGB_injected[i]['Frequency']
    frequencies = np.sort(frequencies)
    data_sorted = []
    Barc_frquencies = []
    for i in range(len(pGB_injected)):
        Barc_frquencies.append(data_set[i]['Frequency'][0])
    for i in range(len(pGB_injected)):
        frequency = pGB_injected[i]['Frequency']
        idx = find_nearest(Barc_frquencies,frequency*10**3)
        print(idx)
        data_sorted.append(data_set[idx])
    data[data_set_names[j]] = data_sorted


rng = [0.999, 0.999, 0.999, 0.999, (0, 1.1), (-22,-21.4), (0.8, np.pi), (0, 2*np.pi)]

# lbls = [ r'\log A', r'\sin \beta',r'\lambda', 'f ($mHz$)', '\log \dot{f} $ $ ($Hz/s$)', r'\cos \iota', r'\phi', r'\Phi']
lbls = [r'\lambda', r'\sin \beta', r'f ($mHz$)', r'\log \dot{f}$ $ ($Hz/s$)', r'\cos \iota', r'\log A']

#########################################
# plot overlap
for k in range(len(pGB_injected)):
    m = 0
    samples = []
    for j in range(2):
        samples.append(MCSamples(samples=data[data_set_names[j]][k].to_numpy(), names = names, labels = lbls))
        samples[-1].updateSettings({'contours': [0.68, 0.95]})
        m += 1
    g = plots.get_subplot_plotter(subplot_size=0.9)
    g.settings.num_plot_contours = 2
    g.triangle_plot(samples, shaded=True, legend_labels=[])

    tr_s = np.zeros(len(parameters))
    maxvalues = np.zeros(len(parameters))
    i = 0
    pGB = pGB_injected[k]
    for parameter in names:
        if parameter in ['Amplitude','FrequencyDerivative']:
            tr_s[i] = np.log10(pGB[parameter])
        elif parameter in ['Inclination']:
            tr_s[i] = np.cos(pGB[parameter])
        elif parameter in ['EclipticLatitude']:
            tr_s[i] = np.sin(pGB[parameter])
        elif parameter in ['Frequency']:
            tr_s[i] = pGB[parameter]*10**3
        else:
            tr_s[i] = pGB[parameter]
        i += 1
    if tr_s[0] > np.pi:
        tr_s[0] -= 2*np.pi
    #markers vertical
    ndim= 6
    for i in range(ndim):
        for ax in g.subplots[i:,i]:
            ax.axvline(tr_s[i], color='black', lw = 1)
        i += 1
    #markers horizontal
    for i in range(ndim):
        for ax in g.subplots[i,:i]:
            ax.axhline(tr_s[i], color='black', lw = 1)
        i += 1
    g.export('/home/stefan/LDC/LDC/pictures/corner overlap '+save_name+'_'+str(int(np.round(pGB_injected[k]['Frequency']*10**8)))+'.png')
print('end')