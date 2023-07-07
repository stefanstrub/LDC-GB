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
         'date': '2022/10/12',           # standardize?
         'challenge': 'LDC2a',
         'dataset': 'LDC2_sangria_training_v2'}


folderpath_parent = grandparent+"/LDC/pictures/Sangria/Chain/"
folderpath_save = grandparent+"/LDC/pictures/Sangria/Chain_mHz/"
name = 'frequency799129nHzSangria_1_full_opt2'
save_name = 'Sangria' + name
folderpath = folderpath_parent + name

df = pd.read_csv(folderpath+'.csv', sep=',')
# df['Frequency'] = df['Frequency']*10**3

df.to_csv(folderpath_save+name+'.csv',index=False, float_format='%.8e')

print(df.values)