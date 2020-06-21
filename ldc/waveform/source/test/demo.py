from ldc.waveform.source import SourceMaker
import numpy as np
import matplotlib.pyplot as plt
import sys

cfg_mbhb = {'approximant': 'IMRPhenomD',
            'catalogs': 'catalog_1_yrs_full_1.dat',
            'coalescence_time': [0.1, 0.95],
            'mass_ratio': [1, 10],
            'mass_total': [2, 10],
            'nsource': 3,
            'seed': 1234,
            'source_type': 'MBHB',
            'spin1': [0., 0.999],
            'spin2': [0., 0.999]}

## Verification Galactic Binairis
cfg_vgb = {'approximant': 'TD_fdot',
          'catalogs': 'VGB.npy',
          'nsource': 2,
          'seed': 1234,
          'source_type': 'GB'}

## Detached Galactic Binaries with gaussian randomization of mass, frequency and inclination
cfg_dgb = {'approximant': 'TD_fdot',
          'catalogs': '007_SeBa_r105_ag_wdwd_pop_highres_P025g70_N1e1.npy',
          'dfdt_GW' : True,
          'random_inclination': "gaussian_0.09",
          'random_mass' : "gaussian_0.01percent",
          'random_frequency' : "gaussian_0.01percent",
          'nsource': 2,
          'seed': 1234,
          'source_type': 'GB'}

cfg_dgb_rd = cfg_dgb.copy()
cfg_dgb_rd['nsource'] = 1
cfg_dgb_rd['catalogs'] = '007_SeBa_r105_ag_wdwd_pop_highres_P025g70_N1.npy'

## Detached Galactic Binaries without gaussian randomization for reference
cfg_dgb_ref = cfg_dgb_rd.copy()
cfg_dgb_ref['random_inclination'] = "0"
cfg_dgb_ref['random_mass'] = "0"
cfg_dgb_ref['random_frequency'] = "0"

## Interacting Galactic Binaries with gaussian randomization of mass and frequency
cfg_igb = {'approximant': 'TD_fdot',
          'catalogs': 'AMCVn_GWR_MLDC_bulgefix_opt_N1e1.npy',
          'dfdt_GW' : False,
          'random_inclination': "uniform",
          'random_mass' : "gaussian_0.01percent",
          'random_frequency' : "gaussian_0.01percent",
          'nsource': 2,
          'seed': 1234,
          'source_type': 'GB'}

cfg_igb_rd = cfg_igb.copy()
cfg_igb_rd['nsource'] = 1
cfg_igb_rd['catalogs'] = 'AMCVn_GWR_MLDC_bulgefix_opt_N1.npy'

## Interacting Galactic Binaries without gaussian randomization for reference
cfg_igb_ref = cfg_igb_rd.copy()
cfg_igb_ref['random_mass'] = "0"
cfg_igb_ref['random_frequency'] = "0"


### Choose sources from catalogs without randomization
for cfg in [cfg_mbhb, cfg_vgb]:
    source_maker = SourceMaker.type(cfg["source_type"],
                                    cfg["approximant"],
                                    catalogs=[cfg["catalogs"]]) 
    cat = source_maker.choose_from_cat(**cfg)

    source_maker.close_logger()
    print(cat)
    print(cat.dtype)
    

### Choose sources from catalogs with randomization
for cfg in [cfg_dgb,cfg_igb]:
    source_maker = SourceMaker.type(cfg["source_type"],
                                    cfg["approximant"],
                                    catalogs=[cfg["catalogs"]]) 
    cat = source_maker.draw_random_catalog(**cfg)

    source_maker.close_logger()
    print(cat)
    print(cat.dtype)

#sys.exit(0)

####### Check distribution

### Reference source for DWD GB
for cfg in [cfg_dgb_ref]:
    source_maker = SourceMaker.type(cfg["source_type"],
                                    cfg["approximant"],
                                    catalogs=[cfg["catalogs"]]) 
    catref = source_maker.draw_random_catalog(**cfg)

    source_maker.close_logger()
    print(catref)
    print(catref.dtype)

### Check distribution DWD GB 
Nrand = 10000
cats = np.zeros(Nrand,dtype=cat.dtype)
for i in range(Nrand):
    for cfg in [cfg_dgb_rd]:
        cfg['seed'] = cfg['seed'] + 1
        source_maker = SourceMaker.type(cfg["source_type"],
                                        cfg["approximant"],
                                        catalogs=[cfg["catalogs"]]) 
        cat = source_maker.draw_random_catalog(**cfg)

        source_maker.close_logger()
        cats[i] = cat
        #print(cat)
        #print(cat.dtype)

print(cats[:10])
print(catref)
print(catref.dtype)

for xR in catref.dtype.names[1:]: 
    plt.figure()
    plt.hist(cats[xR],bins=50)
    plt.title("DGB: "+xR+": "+str(catref[xR][0]))



### Reference source for DWD GB
for cfg in [cfg_igb_ref]:
    source_maker = SourceMaker.type(cfg["source_type"],
                                    cfg["approximant"],
                                    catalogs=[cfg["catalogs"]]) 
    catref = source_maker.draw_random_catalog(**cfg)

    source_maker.close_logger()
    print(catref)
    print(catref.dtype)

### Check distribution DWD GB 
Nrand = 10000
cats = np.zeros(Nrand,dtype=cat.dtype)
for i in range(Nrand):
    for cfg in [cfg_igb_rd]:
        cfg['seed'] = cfg['seed'] + 1
        source_maker = SourceMaker.type(cfg["source_type"],
                                        cfg["approximant"],
                                        catalogs=[cfg["catalogs"]]) 
        cat = source_maker.draw_random_catalog(**cfg)

        source_maker.close_logger()
        cats[i] = cat
        #print(cat)
        #print(cat.dtype)

print(cats[:10])
print(catref)
print(catref.dtype)

for xR in catref.dtype.names[1:]: 
    plt.figure()
    plt.hist(cats[xR],bins=50)
    plt.title("IGB: "+xR+": "+str(catref[xR][0]))

plt.show()

