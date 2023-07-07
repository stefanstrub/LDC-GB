import aesara.tensor as at
import arviz as az
import numpy as np
import pymc as pm

print(f"Running on PyMC v{pm.__version__}")
az.style.use("arviz-darkgrid")

n = 8

mu1 = np.ones(n) * (1.0 / 2)
mu2 = -mu1

stdev = 0.1
sigma = np.power(stdev, 2) * np.eye(n)
isigma = np.linalg.inv(sigma)
dsigma = np.linalg.det(sigma)

w1 = 0.1  # one mode with 0.1 of the mass
w2 = 1 - w1  # the other mode with 0.9 of the mass


def two_gaussians(x):
    log_like1 = (
        -0.5 * n * at.log(2 * np.pi)
        - 0.5 * at.log(dsigma)
        - 0.5 * (x - mu1).T.dot(isigma).dot(x - mu1)
    )
    log_like2 = (
        -0.5 * n * at.log(2 * np.pi)
        - 0.5 * at.log(dsigma)
        - 0.5 * (x - mu2).T.dot(isigma).dot(x - mu2)
    )
    return pm.math.logsumexp([at.log(w1) + log_like1, at.log(w2) + log_like2])

with pm.Model() as model:
    X = pm.Uniform(
        "X",
        shape=n,
        lower=-2.0 * np.ones_like(mu1),
        upper=2.0 * np.ones_like(mu1),
        initval=-1.0 * np.ones_like(mu1),
    )
    llk = pm.Potential("llk", two_gaussians(X))
    step = pm.NUTS()
    idata_04 = pm.sample(2000, chains=10)

ax = az.plot_trace(idata_04, compact=True, kind="rank_vlines")
ax[0, 0].axvline(-0.5, 0, 0.9, color="k")
ax[0, 0].axvline(0.5, 0, 0.1, color="k")
f'Estimated w1 = {np.mean(idata_04.posterior["X"] < 0).item():.3f}'

print('end')