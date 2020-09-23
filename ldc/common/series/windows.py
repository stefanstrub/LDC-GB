""" Provide time domain window functions. 
"""

import numpy as np

def window(tm, xl=1000.0, kap=0.005, show=False):
    """Return time domain window function to remove the first and last xl
    sec. 
    """
    ind_r = np.argwhere(tm[-1]-tm <= xl)[0][0]
    xr = tm[ind_r]
    winl = 0.5*(1.0 + np.tanh(kap*(tm-xl)))
    winr = 0.5*(1.0 - np.tanh(kap*(tm-xr)))
    if show:
        import matplotlib.pyplot as plt
        plt.plot(tm, winl)
        plt.plot(tm, winr)
        plt.grid(True)
        plt.show()
    return (winl*winr)
