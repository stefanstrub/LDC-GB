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

### Choose sources from catalogs without randomization
for cfg in [cfg_mbhb, cfg_vgb]:
    source_maker = SourceMaker.type(cfg["source_type"],
                                    cfg["approximant"],
                                    catalogs=[cfg["catalogs"]]) 
    cat = source_maker.choose_from_catalog(**cfg)

    source_maker.close_logger()
    print(cat)
    print(cat.dtype)
    
 
