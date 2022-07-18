import numpy as np
import pandas as pd
import os

# import pyyaml module
import yaml
from yaml.loader import SafeLoader

# get current directory
path = os.getcwd()
 
# parent directory
parent = os.path.dirname(path)
# grandparent directory
grandparent = os.path.dirname(parent)

DATAPATH = grandparent+"/Repositories/ldc1_evaluation_data/submission"
SAVEPATH = grandparent+"/LDC/Radler/LDC1-4_evaluation/"


# Open the file and load the file
file_name = DATAPATH + '/MSFCMontanaU_Gal/msfc-montana-ldc1-4-gb.yaml'
with open(file_name) as f:
    data = yaml.load(f, Loader=SafeLoader)
    print(data)


##### Montana
i = 0
found_signals = []
for name in data.keys():
    i += 1
    if i < 7:
        continue
    if i == len(data.keys()):
        continue
    found_signals.append(data[name])
np.save(SAVEPATH+'Montana.npy', found_signals)



print('end')