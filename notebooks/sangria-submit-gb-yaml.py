import yaml
import numpy as np
import csv
import pandas as pd
import os
from os import listdir
from os.path import isfile, join

# get current directory
path = os.getcwd()
# parent directory
parent = os.path.dirname(path)
# grandparent directory
grandparent = os.path.dirname(parent)

entry = {'author': 'Stefan Strub',
         'e-mail': 'stefan.strub@erdw.ethz.ch',
         'date': '2024/3/21',           # standardize?
         'challenge': 'LDC2a',
         'dataset': 'LDC2_sangria_training_v2'}


folderpath_parent = grandparent+"/LDC/pictures/Sangria/"
folderpath_save = grandparent+"/LDC/pictures/Sangria/"
name = 'found_sources_not_anticorrelatedSangria_1_full_cut'
name = 'found_sources_not_anticorrelated30000to9261573Sangria_1_dynamic_noise'
name = 'found_sources_Sangria_12m_flat'
name = 'found_sources_not_anticorrelated_original_Sangria_12m_mbhb_SNR9_seed1'
name = 'found_sources_not_anticorrelated_original_Sangria_52w_mbhb_SNR9_seed1'
save_name = name
folderpath = folderpath_parent + name

# found_sources_flat = np.load(folderpath+'.npy', allow_pickle = True)
found_sources = np.load(folderpath+'.pkl', allow_pickle = True)
found_sources_flat = np.concatenate(found_sources)
# found_sources_flat = np.concatenate(found_sources)
# found_sources_flat = found_sources_flat[:100]
for i in range(len(found_sources_flat)):
    found_sources_flat[i]['Frequency'] *= 10**3
for parameter in found_sources_flat[0].keys():
    for i in range(len(found_sources_flat)):
        found_sources_flat[i][parameter] = float(found_sources_flat[i][parameter])
found_sources_flat = list(found_sources_flat)

units = {  
        'Amplitude' :1, 
        'Frequency' :'mHz',
        'FrequencyDerivative':'Hz',
        'EclipticLatitude': 'rad',
        'EclipticLongitude': 'rad',
        'Inclination': 'rad',
        'Polarization': 'rad', 
        'InitialPhase': 'rad' 
     }

entry['estimates'] = found_sources_flat
entry['units'] = units

open(folderpath_save+'ETH-LDC2-sangria-training-v2-training.yaml','w').write(yaml.dump(entry, default_flow_style=False))

print('end')