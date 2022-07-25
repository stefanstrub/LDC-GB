import numpy as np
import os
from os import listdir
from os.path import isfile, join

# get current directory
path = os.getcwd()
# parent directory
parent = os.path.dirname(path)
# grandparent directory
grandparent = os.path.dirname(parent)
# folderpath_parent = grandparent+"/LDC/pictures/LDC1-4s3"
# folderpath_parent = grandparent+"/LDC/pictures/LDC1-4uneven"
# folderpath_parent = grandparent+"/LDC/pictures/LDC1-4_4mHz_Euler"
# folderpath_parent = grandparent+"/LDC/pictures/LDC1-4_half_even"
# folderpath_parent = grandparent+"/LDC/pictures/LDC1-4_4mHz/Found_signals"
# folderpath_parent = grandparent+"/LDC/pictures/LDC1-4/Found_signals"
# folderpath_parent = grandparent+"/LDC/pictures/LDC1-4/Found_signals_half_year_even3_T"
# folderpath_parent = grandparent+"/LDC/pictures/LDC1-4/Found_signals_half_year_odd_T"
# folderpath_parent = grandparent+"/LDC/pictures/LDC1-4/found_signals_half_year_even10"
# folderpath_parent = grandparent+"/LDC/pictures/LDC1-4/found_signals_half_year_full"
folderpath_parent = grandparent+"/LDC/pictures/LDC1-4/found_signals"

name = '_2_years_full'
save_name = 'LDC1-4' + name
folderpath = folderpath_parent + name
onlyfiles = [f for f in listdir(folderpath) if isfile(join(folderpath, f))]
found_sources_mp_even_unsorted = []
frequencies = []
for i in range(len(onlyfiles)):
    print(i)
    sources = np.load(folderpath+'/'+onlyfiles[i], allow_pickle = True).tolist()
    frequency = 0
    for j in range(len(sources)):
        try:
            frequency = sources[j][1][0][0][0]['Frequency']
            frequencies.append(frequency)
            found_sources_mp_even_unsorted.append(sources[j])
        except:
            pass
    # frequencies.append(frequency)
    # found_sources_mp_even_unsorted.append(sources)
sorted_indexes = np.argsort(frequencies)
found_sources_mp = []
for i in range(len(sorted_indexes)):
    found_sources_mp.append(found_sources_mp_even_unsorted[sorted_indexes[i]])

np.save(folderpath+'/found_sources' +save_name+'.npy', found_sources_mp)
