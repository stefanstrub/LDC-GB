import tensorflow.compat.v2 as tf
import tensorflow_probability as tfp
import numpy as np
import matplotlib.pyplot as plt

import corner

tf.enable_v2_behavior()

tfd = tfp.distributions
dtype = np.float32
true_mean = dtype([0, 0, 0])
true_cov = dtype([[1, 0.25, 0.25], [0.25, 1, 0.25], [0.25, 0.25, 1]])
num_results = 500
num_chains = 1

# Target distribution is defined through the Cholesky decomposition
chol = tf.linalg.cholesky(true_cov)
target = tfd.MultivariateNormalTriL(loc=true_mean, scale_tril=chol)
def target(x):
    print(x[:,0])
    return -x[:,0] - x[:,1]**2
# Here we define the volatility function to be non-constant
def volatility_fn(x):
  # Stack the input tensors together
  return 1. / (0.5 + 0.1 * tf.math.abs(x))

# Initial state of the chain
init_state = np.ones([num_chains, 3], dtype=dtype)

# Run MALA with normal proposal for `num_results` iterations for
# `num_chains` independent chains:
states = tfp.mcmc.sample_chain(
  num_results=1000,
  current_state=init_state,
  kernel=tfp.mcmc.RandomWalkMetropolis(target),
  num_burnin_steps=500,
  trace_fn=None,
  seed=42)



data = states.numpy()
data = data.reshape((np.shape(data)[0]*np.shape(data)[1],np.shape(data)[2]))
fig =  corner.corner(data,  bins=40, hist_kwargs={'density':True, 'lw':3}, plot_datapoints=False, fill_contours=False,  show_titles=True, \
                        color='#348ABD',  truth_color='k', use_math_test=True,\
                        levels=[0.9], title_kwargs={"fontsize": 12})
plt.show()
sample_mean = tf.reduce_mean(states, axis=[0, 1])
x = (states - sample_mean)[..., tf.newaxis]
sample_cov = tf.reduce_mean(
    tf.matmul(x, tf.transpose(x, [0, 1, 3, 2])), [0, 1])

print('sample mean', sample_mean.numpy())
print('sample covariance matrix', sample_cov.numpy())




tf.enable_v2_behavior()
tfd = tfp.distributions
dtype = np.float32
true_mean = dtype([0, 0, 0])
true_cov = dtype([[1, 0.25, 0.25], [0.25, 1, 0.25], [0.25, 0.25, 1]])
num_results = 5000
num_chains = 1
# Target distribution is defined through the Cholesky decomposition
chol = tf.linalg.cholesky(true_cov)
target = tfd.MultivariateNormalTriL(loc=true_mean, scale_tril=chol)
def target(x):
    return gpr.predict(x)
# Here we define the volatility function to be non-constant
def volatility_fn(x):
  # Stack the input tensors together
  return 1. / (0.5 + 0.1 * tf.math.abs(x))
init_state = np.ones([num_chains, 8], dtype=dtype)*0.5

states = tfp.mcmc.sample_chain(
    num_results=num_results,
    current_state=init_state,
    kernel=tfp.mcmc.MetropolisAdjustedLangevinAlgorithm(
        target_log_prob_fn=target,
        step_size=.1,
        volatility_fn=volatility_fn),
    num_burnin_steps=200,
    num_steps_between_results=1,
    trace_fn=None,
    seed=42)
data = states.numpy()
data = data.reshape((np.shape(data)[0]*np.shape(data)[1],np.shape(data)[2]))
fig =  corner.corner(data,  bins=40, hist_kwargs={'density':True, 'lw':3}, plot_datapoints=False, fill_contours=False,  show_titles=True, \
                        color='#348ABD',  truth_color='k', use_math_test=True,\
                        levels=[0.9], title_kwargs={"fontsize": 12})
plt.show()