import numpy as np
import matplotlib.pyplot as plt
from ldc.lisa.noise import get_noise_model

fr = np.logspace(-5, 0, 10000)
model = 'SciRDv1'
wd = 1 # year

plt.figure()
plt.loglog(fr, get_noise_model(model, fr, wd=0).psd(option='A'),
           label=model)
plt.loglog(fr, get_noise_model(model, fr, wd=wd).psd(option='A'),
           label=model+" + %dyr wd"%wd)
plt.legend()
