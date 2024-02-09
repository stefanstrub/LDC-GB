import numpy as np
import pandas as pd
from copy import deepcopy

recovered = []
matched = []
rates = []
with open('notebooks/results.txt') as f:
    lines = f.readlines()
    for line in lines:
        split = line.split()
        matched.append(int(split[1]))
        recovered.append(int(split[10]))
        rates.append(matched[-1]/recovered[-1])

mean_matched = np.mean(matched)
mean_recovered = np.mean(recovered)
mean_rates = np.mean(rates)
median_matched = np.median(matched)
median_recovered = np.median(recovered)
std_matched = np.std(matched)
std_recovered = np.std(recovered)
std_rates = np.std(rates)
print('Mean matched: ', mean_matched)
print('Mean recovered: ', mean_recovered)
print('Mean rates: ', mean_rates)
print('Median matched: ', median_matched)
print('Median recovered: ', median_recovered)
print('Std matched: ', std_matched)
print('Std recovered: ', std_recovered)
print('Std rates: ', std_rates)

groups = []
with open('notebooks/results_groups.txt') as f:
    lines = f.readlines()
    for line in lines:
        groups.append([])
        split = line.split()
        for item in split:
            groups[-1].append(int(item))
groups = np.array(groups)
groups = pd.DataFrame(groups)
print(groups)
print(groups.to_latex(index=False))

groups_extend = deepcopy(groups)
# groups = np.array(groups)
groups_extend.loc['mean'] = groups.mean()
groups_extend.loc['std'] = groups.std()

groups_extend = groups_extend.style.format(decimal='.', thousands='\,', precision=1)
print(groups_extend.to_latex())

recovered = 0