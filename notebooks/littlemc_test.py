import numpy as np
import scipy
import littlemcmc as lmc

def logp_func(x, loc=0, scale=1):
    return np.log(scipy.stats.norm.pdf(x, loc=loc, scale=scale))


def dlogp_func(x, loc=0, scale=1):
    return -(x - loc) / scale


def logp_dlogp_func(x, loc=0, scale=1):
    return logp_func(x, loc=loc, scale=scale), dlogp_func(x, loc=loc, scale=scale)

# By default: 4 chains in 4 cores, 500 tuning steps and 1000 sampling steps.
trace, stats = lmc.sample(
    logp_dlogp_func=logp_dlogp_func,
    model_ndim=4,
    chains=4,
    cores=4,
    progressbar=None,  # HTML progress bars don't render well in RST.
)

print('end')