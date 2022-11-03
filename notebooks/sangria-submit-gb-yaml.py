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
         'date': '2022/11/2',           # standardize?
         'challenge': 'LDC2a',
         'dataset': 'LDC2_sangria_training_v2'}


folderpath_parent = grandparent+"/LDC/pictures/Sangria/"
folderpath_save = grandparent+"/LDC/pictures/Sangria/"
name = 'found_sources_not_anticorrelatedSangria_1_full_cut'
save_name = name
folderpath = folderpath_parent + name

found_sources = np.load(folderpath+'.npy', allow_pickle = True)
found_sources_flat = np.concatenate(found_sources)
# found_sources_flat = found_sources_flat[:100]
for i in range(len(found_sources_flat)):
    found_sources_flat[i]['Frequency'] *= 10**3
for parameter in found_sources_flat[0].keys():
    for i in range(len(found_sources_flat)):
        found_sources_flat[i][parameter] = float(found_sources_flat[i][parameter])
found_sources_flat = list(found_sources_flat)

gb = [{  'Amplitude' :4.065235140346656e-23, 
        'Frequency' :2.1869889108004115, 
        'FrequencyDerivative':4.837869734659327e-17, 
        'EclipticLatitude': -0.3982801352498764, 
        'EclipticLongitude': 4.4946831814774795, 
        'Inclination': 1.0208838997376564, 
        'Polarization': 3.737792079072995, 
        'InitialPhase': 5.6307387805247435, 
        'IntrinsicSNR': 9.42426286925173
     },{  'Amplitude' :4.065235140346656e-23, 
        'Frequency' :2.1869889108004115, 
        'FrequencyDerivative':4.837869734659327e-17, 
        'EclipticLatitude': -0.3982801352498764, 
        'EclipticLongitude': 4.4946831814774795, 
        'Inclination': 1.0208838997376564, 
        'Polarization': 3.737792079072995, 
        'InitialPhase': 5.6307387805247435, 
        'IntrinsicSNR': 9.42426286925173
     }]
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

open(folderpath_save+'ETH-LDC2-sangria-training-v2-training-gb-cut.yaml','w').write(yaml.dump(entry, default_flow_style=False))

print('end')