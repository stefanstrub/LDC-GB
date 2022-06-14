
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

from scipy import stats
def jensen_shannon_divergence(
    samples, kde=stats.gaussian_kde, decimal=5, base=np.e, **kwargs
):
    """Calculate the JS divergence between two sets of samples

    Parameters
    ----------
    samples: list
        2d list containing the samples drawn from two pdfs
    kde: func
        function to use when calculating the kde of the samples
    decimal: int, float
        number of decimal places to round the JS divergence to
    base: float, optional
        optional base to use for the scipy.stats.entropy function. Default
        np.e
    kwargs: dict
        all kwargs are passed to the kde function
    """
    try:
        kernel = [kde(i, **kwargs) for i in samples]
    except np.linalg.LinAlgError:
        return float("nan")
    x = np.linspace(
        np.min([np.min(i) for i in samples]),
        np.max([np.max(i) for i in samples]),
        100
    )
    a, b = [k(x) for k in kernel]
    a = np.asarray(a)
    b = np.asarray(b)
    a /= a.sum()
    b /= b.sum()
    m = 1. / 2 * (a + b)
    kl_forward = stats.entropy(a, qk=m, base=base)
    kl_backward = stats.entropy(b, qk=m, base=base)
    return np.round(kl_forward / 2. + kl_backward / 2., decimal)

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

#########################################
# plot overlap
lbls = [r'\lambda', r'\sin \beta', r'f ($mHz$)', r'\log \dot{f}$ $ ($Hz/s$)', r'\cos \iota', r'\log A']
m = 0
names = ['EclipticLongitude','EclipticLatitude','Frequency','FrequencyDerivative','Inclination','Amplitude']
samples = []
distributions = []
# for file_name in ['GW201457number of singal0LDC1-3overlap single signal', 'GW201457number of singal0LDC1-3overlap']:
#         df = pd.read_csv('/home/stefan/Repositories/ldc1_evaluation_data/submission/ETH_2/'+file_name+'.csv')
# for file_name in ['frequency1252567nHzLDC1-3', 'frequency1252567nHzLDC1-3fastGB']:
for file_name in ['frequency1666286nHzLDC1-3second', 'frequency1666286nHzLDC1-3fastGB']:
        df = pd.read_csv('/home/stefan/LDC/pictures/LDC1-3_v2/Chain/'+file_name+'.csv')
        df['Inclination'] = np.cos(df['Inclination'].values)
        df['EclipticLatitude'] = np.sin(df['EclipticLatitude'].values)
        df['FrequencyDerivative'] = np.log10(df['FrequencyDerivative'].values)
        df['Amplitude'] = np.log10(df['Amplitude'].values)
        df['Frequency'] *= 1000
        df = df.drop(labels= ['InitialPhase', 'Polarization'], axis=1)
        df = df.reindex(columns = names)
        samples.append(MCSamples(samples=df.to_numpy(), names = names, labels = lbls))
        distributions.append(df['Frequency'].to_numpy())
        samples[-1].updateSettings({'contours': [0.68, 0.95]})
        m += 1
print('Jansen Shannon divergence:', jensen_shannon_divergence(distributions)) 
g = plots.get_subplot_plotter(subplot_size=0.9)
g.settings.num_plot_contours = 2
g.triangle_plot(samples, shaded=True, legend_labels=[])

tr_s = np.zeros(len(parameters))
maxvalues = np.zeros(len(parameters))
i = 0
pGB = pGB_injected[2]
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
g.export('/home/stefan/LDC/LDC/pictures/corner overlap '+save_name+'fastGBsecond.png')
