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
         'date': '2024/1/18',           # standardize?
         'challenge': 'LDC2a',
         'dataset': 'LDC2_sangria_training_v2'}


folderpath_parent = grandparent+"/LDC/MBHB/"
folderpath_save = grandparent+"/LDC/MBHB/"
name = 'Sangria_mbhbh_found_signals_12months_seed4'
save_name = name
folderpath = folderpath_parent + name

parameters = ['Mass1', 'Mass2', 'Spin1', 'Spin2', 'Distance', 'Phase', 'Inclination', 'EclipticLongitude', 'EclipticLatitude', 'Polarization', 'CoalescenceTime']
# found_sources_flat = np.load(folderpath+'.npy', allow_pickle = True)
found_sources = np.load(folderpath+'.pkl', allow_pickle = True)
found_sources_flat = []
# found_sources_flat = np.concatenate(found_sources)
# found_sources_flat = found_sources_flat[:100]
for i in range(len(found_sources)):
    found_sources_flat.append({})
    for j, parameter in enumerate(parameters):
        found_sources_flat[i][parameter] = found_sources[i][j]

PC_SI = 3.0856775814913674e+16
for i in range(len(found_sources_flat)):
    found_sources_flat[i]['Distance'] /= PC_SI*1e6

coalescence_time = []
for i in range(len(found_sources_flat)):
    coalescence_time.append(found_sources_flat[i]['CoalescenceTime'])
index_sorted = np.argsort(coalescence_time)
found_sources_flat = np.array(found_sources_flat)[index_sorted]

for parameter in found_sources_flat[0].keys():
    for i in range(len(found_sources_flat)):
        found_sources_flat[i][parameter] = float(found_sources_flat[i][parameter])
found_sources_flat = list(found_sources_flat)

units = {  
        'CoalescenceTime' :'s', 
        'Distance' :'Mpc',
        'EclipticLatitude': 'rad',
        'EclipticLongitude': 'rad',
        'Inclination': 'rad',
        'Mass1': 'M_sun',
        'Mass2': 'M_sun',
        'Phase': 'rad',
        'Polarization': 'rad', 
        'ProjectedSpin1': 'none',
        'ProjectedSpin2': 'none',
     }

entry['estimates'] = found_sources_flat
entry['units'] = units

open(folderpath_save+'ETH-LDC2-sangria-training-v2-training-mbhb.yaml','w').write(yaml.dump(entry, default_flow_style=False))

print('end')