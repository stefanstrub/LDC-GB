import numpy as np
import pandas as pd
import os
from os import listdir
from os.path import isfile, join
from copy import deepcopy

# get current directory
path = os.getcwd()
# parent directory
parent = os.path.dirname(path)
# grandparent directory
grandparent = os.path.dirname(parent)
folderpath_parent = grandparent+"/LDC/pictures/LDC1-4/lsf"

name = '_2_even3'
save_name = 'LDC1-4' + name
folderpath = folderpath_parent + name
onlyfiles = [f for f in listdir(folderpath) if isfile(join(folderpath, f))]
found_sources_mp_even_unsorted = []
frequencies = []
cpu_time = 0
run_time = 0
for i in range(len(onlyfiles)):
    print(i)
    with open(folderpath+'/'+onlyfiles[i]) as f:
        lines = f.readlines()
    if lines[18] == 'Successfully completed.\n':
        print('successfull')
        cpu_time += float(lines[22].split()[3])
        run_time += float(lines[30].split()[3])
print('CPU time: ',cpu_time, 's', cpu_time/3600, 'h')
print('run time: ',run_time, 's', run_time/3600, 'h')
print('efficiency: ',cpu_time/(run_time*16))
print('end')