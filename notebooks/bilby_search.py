#!/usr/bin/env python
"""
Tutorial to demonstrate running parameter estimation on a reduced parameter
space for an injected signal.

This example estimates the masses using a uniform prior in both component masses
and distance using a uniform in comoving volume prior on luminosity distance
between luminosity distances of 100Mpc and 5Gpc, the cosmology is Planck15.
"""

from math import log
import numpy as np
import bilby
import time
import matplotlib.pyplot as plt
from copy import deepcopy
import scipy
from pycbc.detector import Detector

d = Detector("V1")
tof = {}
tof['H1'] = d.light_travel_time_to_detector(Detector("H1"))
tof['L1'] = d.light_travel_time_to_detector(Detector("L1"))
# Set the duration and sampling frequency of the data segment that we're
# going to inject the signal into
duration = 4.
sampling_frequency = 2048.

# Specify the output directory and the name of the simulation.
outdir = 'outdir'
label = 'fast_tutorial'
bilby.core.utils.setup_logger(outdir=outdir, label=label)

# Set up a random seed for result reproducibility.  This is optional!
np.random.seed(88170235)

# We are going to inject a binary black hole waveform.  We first establish a
# dictionary of parameters that includes all of the different waveform
# parameters, including masses of the two black holes (mass_1, mass_2),
# spins of both black holes (a, tilt, phi), etc.
injection_parameters = dict(
    mass_1=36., mass_2=29., a_1=0.4, a_2=0.3, tilt_1=0.5, tilt_2=1.0,
    phi_12=1.7, phi_jl=0.3, luminosity_distance=2000., theta_jn=0.4, psi=2.659,
    phase=1.3, geocent_time=1126259642.413, ra=1.375, dec=-1.2108)
injection_parameters = dict(
    mass_1=35.6, mass_2=30.6, a_1=0.0, a_2=0.0, tilt_1=0.5, tilt_2=1.0,
    phi_12=1.7, phi_jl=0.3, luminosity_distance=440., theta_jn=0.4, psi=2.659,
    phase=1.3, geocent_time=1126259642.413, ra=1.375, dec=-1.2108) # GW150914
injection_parameters = dict(
    mass_1=43.5, mass_2=35.1, a_1=0.0, a_2=0.0, tilt_1=0.5, tilt_2=1.0,
    phi_12=1.7, phi_jl=0.3, luminosity_distance=570., theta_jn=0.4, psi=2.659,
    phase=1.3, geocent_time=1126259642.413, ra=1.375, dec=-1.2108) # GW190910

# Fixed arguments passed into the source model
waveform_arguments = dict(waveform_approximant='IMRPhenomPv2',
                          reference_frequency=50., minimum_frequency=20.)

# Create the waveform_generator using a LAL BinaryBlackHole source function
waveform_generator = bilby.gw.WaveformGenerator(
    duration=duration, sampling_frequency=sampling_frequency,
    frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
    waveform_arguments=waveform_arguments)

# Set up interferometers.  In this case we'll use two interferometers
# (LIGO-Hanford (H1), LIGO-Livingston (L1). These default to their design
# sensitivity
ifos = bilby.gw.detector.InterferometerList(['H1','L1'])
ifos.set_strain_data_from_power_spectral_densities(
    sampling_frequency=sampling_frequency, duration=duration,
    start_time=injection_parameters['geocent_time'] - 3)
ifos.inject_signal(waveform_generator=waveform_generator,
                   parameters=injection_parameters)

# Set up a PriorDict, which inherits from dict.
# By default we will sample all terms in the signal models.  However, this will
# take a long time for the calculation, so for this example we will set almost
# all of the priors to be equall to their injected values.  This implies the
# prior is a delta function at the true, injected value.  In reality, the
# sampler implementation is smart enough to not sample any parameter that has
# a delta-function prior.
# The above list does *not* include mass_1, mass_2, theta_jn and luminosity
# distance, which means those are the parameters that will be included in the
# sampler.  If we do nothing, then the default priors get used.
priors = bilby.gw.prior.BBHPriorDict()
priors.pop('mass_1')
priors.pop('mass_2')

priors['chirp_mass'] = bilby.prior.Uniform(
    name='chirp_mass', latex_label='$M$', minimum=10.0, maximum=100.0,
    unit='$M_{\\odot}$')

priors['mass_ratio'] = bilby.prior.Uniform(
    name='mass_ratio', latex_label='$q$', minimum=0.5, maximum=1.0)

priors['geocent_time'] = bilby.core.prior.Uniform(
    minimum=injection_parameters['geocent_time'] - 0.02,
    maximum=injection_parameters['geocent_time'] + 0.02,
    name='geocent_time', latex_label='$t_c$', unit='$s$')

# priors['luminosity_distance'].minimum = 100
# priors['luminosity_distance'].maximum = 4900
# priors['a_1'].minimum = 0
# priors['a_1'].maximum = 0.01
# priors['a_2'].minimum = 0
# priors['a_2'].maximum = 0.01

# Initialise the likelihood by passing in the interferometer data (ifos) and
# the waveform generator
likelihood = bilby.gw.GravitationalWaveTransient(
    interferometers=ifos, waveform_generator=waveform_generator)
likelihood.parameters = deepcopy(injection_parameters)

def funcMchirpofm1m2(m1, m2):
    return pow(m1*m2, 3./5) / pow(m1+m2, 1./5)
def funcm1m2ofMchirpq(Mc, q):
    m1 = (1+q)**(1/5)*Mc/q**(3/5)
    return m1, m1*q
print('target',likelihood.log_likelihood())

parameters_mass = ['mass_1','mass_2','a_1','a_2','tilt_1','tilt_2','phi_12','phi_jl','luminosity_distance','theta_jn','psi','phase','geocent_time','ra','dec']
parameters_wo_mass = ['a_1','a_2','tilt_1','tilt_2','phi_12','phi_jl','luminosity_distance','theta_jn','psi','phase','geocent_time','ra','dec']
parameters_chirp = ['chirp_mass','mass_ratio','a_1','a_2','tilt_1','tilt_2','phi_12','phi_jl','luminosity_distance','theta_jn','psi','phase','geocent_time','ra','dec']

def CoordinateMC(n):
    start = time.time()
    parameters_recorded = {}
    for parameter in parameters_chirp:
        parameters_recorded[parameter] = []
        parameters_recorded[parameter].append(priors[parameter].sample(1)[0])
    for parameter in parameters_wo_mass:
        likelihood.parameters[parameter] = parameters_recorded[parameter][0]
    likelihood.parameters['mass_1'], likelihood.parameters['mass_2'] = funcm1m2ofMchirpq(parameters_recorded['chirp_mass'][0], parameters_recorded['mass_ratio'][0])
    try:
        best_value
    except:
        best_value = likelihood.log_likelihood()
    likelihood2 = deepcopy(likelihood)
    previous_best = deepcopy(best_value)
    n_trials = 20
    for i in range(200):
        parameter1 = parameters_chirp[i % 15]
        parameter2 = parameters_chirp[np.random.randint(0, 14)]
        while parameter2 == parameter1:
            parameter2 = parameters_chirp[np.random.randint(0, 14)]
        changeableparameters = [parameter1, parameter2]
        for j in range(n_trials):
            change_chirp = False
            change_ratio = False
            if parameter1 == 'chirp_mass':
                change_chirp = True
                chirp_mass = priors[parameter1].sample(1)[0]
            elif parameter1 == 'mass_ratio':
                change_ratio = True
                mass_ratio = priors[parameter1].sample(1)[0]
            if parameter2 == 'chirp_mass':
                change_chirp = True
                chirp_mass = priors[parameter2].sample(1)[0]
            elif parameter2 == 'mass_ratio':
                change_ratio = True
                mass_ratio = priors[parameter2].sample(1)[0]
            if not(change_chirp):
                chirp_mass = parameters_recorded['chirp_mass'][-1]
            if not(change_ratio):
                mass_ratio = parameters_recorded['mass_ratio'][-1]
            if change_chirp or change_ratio:
                likelihood2.parameters['mass_1'], likelihood2.parameters['mass_2'] = funcm1m2ofMchirpq(chirp_mass, mass_ratio)
    
            for parameter in changeableparameters:
                if parameter in parameters_wo_mass:
                    likelihood2.parameters[parameter] = priors[parameter].sample(1)[0]

            suggestion = likelihood2.log_likelihood()
            if suggestion > previous_best:
                previous_best = suggestion
                for parameter in parameters_mass:
                    likelihood.parameters[parameter] = likelihood2.parameters[parameter]
        # if i in [30,40,50,60,70]:
        #     past_mean = 0
        #     sum_count = 0
        #     for l in range(n):
        #         try:
        #             past_mean += parameters_recorded[l][0]["Loglikelihood"][i]
        #             sum_count += 1
        #         except:
        #             pass
        #     try:
        #         past_mean = past_mean / sum_count
        #         if previous_best > past_mean:
        #             pass
        #         else:
        #             break
        #     except:
        #         pass

        # start = time.time()
        # print(n, i, previous_best, loglikelihood(maxpGB), maxpGB)
        # parameters_recorded["Loglikelihood"].append(loglikelihood(maxpGB))
        for parameter in parameters_wo_mass:
            parameters_recorded[parameter].append(likelihood.parameters[parameter])
        parameters_recorded['chirp_mass'].append(funcMchirpofm1m2(likelihood.parameters['mass_1'],likelihood.parameters['mass_2']))
        parameters_recorded['mass_ratio'].append(likelihood.parameters['mass_2']/likelihood.parameters['mass_1'])

        likelihood2 = deepcopy(likelihood)
    print('time',time.time()-start, n, i,previous_best )
    if previous_best > best_value:
        best_value = previous_best
    # parameters_recorded[n] = parameters_recorded1
    return parameters_recorded

def normalize_inv(value, min, max):
    return value*(max-min)+min
def normalize(value, min, max):
    return (value-min)/(max-min)

def optimize(pGBmodes):
    boundaries = []
    for parameter in parameters_chirp:
        boundaries.append((priors[parameter].minimum,priors[parameter].maximum))
    bounds = ((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1))
    for i in range(len(pGBmodes)):
        maxpGB = {}
        boundaries_reduced = []
        sample_parameters = pGBmodes[i]

        for j in range(1):
            x = []
            for parameter in parameters_chirp:
                if parameter in ['tilt_1','tilt_2','phi_12','phi_jl','theta_jn','psi','phase','ra','dec']:
                    length = priors[parameter].maximum - priors[parameter].minimum
                    x.append(normalize(sample_parameters[parameter],priors[parameter].minimum-length,priors[parameter].maximum+length))
                else:
                    x.append(normalize(sample_parameters[parameter],priors[parameter].minimum,priors[parameter].maximum))
            best_value_norm = -function(x,1)
            res = scipy.optimize.minimize(function, x,args=best_value_norm, method='SLSQP', bounds=bounds, tol=1e-10, options={'maxiter':400, 'disp': False})
            for k, parameter in enumerate(parameters_chirp):
                if parameter in ['tilt_1','tilt_2','phi_12','phi_jl','theta_jn','psi','phase','ra','dec']:
                    length = priors[parameter].maximum - priors[parameter].minimum
                    maxpGB[parameter] = normalize_inv(res.x[k],priors[parameter].minimum-length,priors[parameter].maximum+length)
                    if maxpGB[parameter] > priors[parameter].maximum:
                        maxpGB[parameter] -= length
                    elif maxpGB[parameter] < priors[parameter].minimum:
                        maxpGB[parameter] += length
                else:
                    maxpGB[parameter] = normalize_inv(res.x[k],priors[parameter].minimum,priors[parameter].maximum)      
            # print('optimized loglikelihood', loglikelihood(maxpGB),maxpGB)
            # print('boundaries reduced', boundaries_reduced)
        best_value = -function(res.x,1)
        print(best_value, res.fun,-function(res.x,best_value_norm),-function(res.x,best_value_norm)*best_value_norm,best_value_norm)
        if i == 0:
            current_best_value = best_value
            current_maxpGB = maxpGB
        try:
            if current_best_value < best_value:
                current_best_value = best_value
                current_maxpGB = maxpGB
        except:
            pass
        print(maxpGB)
    maxpGB = current_maxpGB
    print(current_best_value)
    return maxpGB

def function(sample_parameters, best_value):
    for k, parameter in enumerate(parameters_wo_mass):
        likelihood.parameters[parameter] = normalize_inv(sample_parameters[k+2],priors[parameter].minimum,priors[parameter].maximum)
        if parameter in ['tilt_1','tilt_2','phi_12','phi_jl','theta_jn','psi','phase','ra','dec']:
            length = priors[parameter].maximum - priors[parameter].minimum
            likelihood.parameters[parameter] = normalize_inv(sample_parameters[k+2],priors[parameter].minimum-length,priors[parameter].maximum+length)
            if likelihood.parameters[parameter] > priors[parameter].maximum:
                likelihood.parameters[parameter] -= length
            elif likelihood.parameters[parameter] < priors[parameter].minimum:
                likelihood.parameters[parameter] += length
    likelihood.parameters['mass_1'], likelihood.parameters['mass_2'] = funcm1m2ofMchirpq(normalize_inv(sample_parameters[0],priors['chirp_mass'].minimum,priors['chirp_mass'].maximum), normalize_inv(sample_parameters[1],priors['mass_ratio'].minimum,priors['mass_ratio'].maximum))
    return -likelihood.log_likelihood()/np.abs(best_value)

parameters_recorded = []
for n in range(5):
    parameters_recorded.append(CoordinateMC(n))
pGBmodes = []
for i in range(len(parameters_recorded)):
    pGBmodes.append({})
    for parameter in parameters_chirp:
        pGBmodes[i][parameter] = parameters_recorded[i][parameter][-1]

maxpGB = optimize(pGBmodes)
print(funcMchirpofm1m2(injection_parameters['mass_1'],injection_parameters['mass_2']), injection_parameters['mass_2']/injection_parameters['mass_1'])
print(maxpGB)
maxpGB2 = {'chirp_mass': 32.58854742193793, 'mass_ratio': 0.6829392403929653, 'a_1': 0.7706078009013267, 'a_2': 0.5581718691934835, 'tilt_1': 2.1083279826718924, 'tilt_2': 0.953148985323141, 'phi_12': 0.16229254241426716, 'phi_jl': 0.11727554557996811, 'luminosity_distance': 586.8833883991317, 'theta_jn': 0.06942314993813303, 'psi': 1.0856106634837044, 'phase': 0.8190299697919876, 'geocent_time': 1126259642.411288, 'ra': 1.6366446509765513, 'dec': -1.2514165244553155}

found_parameters = {}
for parameter in parameters_wo_mass:
    likelihood.parameters[parameter] = maxpGB[parameter]
    found_parameters[parameter] = maxpGB2[parameter]
likelihood.parameters['mass_1'], likelihood.parameters['mass_2'] = funcm1m2ofMchirpq(maxpGB['chirp_mass'],maxpGB['mass_ratio'])
found_parameters['mass_1'], found_parameters['mass_2'] = funcm1m2ofMchirpq(maxpGB2['chirp_mass'],maxpGB2['mass_ratio'])
print(likelihood.log_likelihood())
for parameter in parameters_mass:
    likelihood.parameters[parameter] = injection_parameters[parameter]
print(injection_parameters)
print(likelihood.log_likelihood())
# snr = likelihood.calculate_snrs(ifos[0])
for x_parameter, y_parameter in [['chirp_mass','mass_ratio'],['dec','ra'],['a_1','a_2'],['tilt_1','tilt_2'],['phi_12','phi_jl'],['luminosity_distance','theta_jn'],['geocent_time','psi'],['phase','psi']]:
    likelihood.parameters = deepcopy(found_parameters)
    print(likelihood.log_likelihood())
    fig = plt.figure()
    parameters = {}
    N = 1000
    loglikelihood = np.zeros(N)
    parameters[x_parameter] = priors[x_parameter].sample(N)
    parameters[y_parameter] = priors[y_parameter].sample(N)
    if x_parameter != 'chirp_mass':
        parameters[x_parameter][0] = found_parameters[x_parameter]
        parameters[y_parameter][0] = found_parameters[y_parameter]

    # likelihood.parameters['dec'] = -0.7
    # likelihood.parameters['ra'] = 4
    start = time.time()
    for i in range(N):
        if x_parameter == 'chirp_mass':
            likelihood.parameters['mass_1'], likelihood.parameters['mass_2'] = funcm1m2ofMchirpq(parameters[x_parameter][i], parameters[y_parameter][i])
        else:
            likelihood.parameters[x_parameter] = parameters[x_parameter][i]
            likelihood.parameters[y_parameter] = parameters[y_parameter][i]
        loglikelihood[i] = likelihood.log_likelihood()
    print('time:', time.time()-start)
    print('max',np.max(loglikelihood))
    plt.scatter(parameters[x_parameter],parameters[y_parameter], c= loglikelihood)
    plt.xlabel(x_parameter)
    plt.ylabel(y_parameter)
    cbar = plt.colorbar()
    if x_parameter != 'chirp_mass':
        plt.hlines(injection_parameters[y_parameter],xmin=np.min(parameters[x_parameter]),xmax=np.max(parameters[x_parameter]))
        plt.vlines(injection_parameters[x_parameter],ymin=np.min(parameters[y_parameter]),ymax=np.max(parameters[y_parameter]))
plt.show()


strain = waveform_generator.frequency_domain_strain(parameters=injection_parameters)

plt.plot(np.linspace(0,duration,len(ifos.frequency_array)*2-2),ifos[0].time_domain_strain)
plt.plot(np.linspace(0,duration,len(ifos.frequency_array)*2-2),strain_td['plus'])
plt.plot(np.linspace(0,duration,len(ifos.frequency_array)),np.fft.fft(strain['plus']))
plt.show()

plt.plot(ifos.frequency_array,strain['plus'])
plt.plot(ifos.frequency_array,np.fft.rfft(strain_td['plus']))
plt.plot(ifos.frequency_array,ifos[0].frequency_domain_strain)
plt.show()
plt.plot(np.linspace(0,duration,len(ifos.frequency_array)),np.fft.fft(strain['plus']))
plt.show()
# Run sampler.  In this case we're going to use the `dynesty` sampler
result = bilby.run_sampler(
    likelihood=likelihood, priors=priors, sampler='dynesty', npoints=10,
    injection_parameters=injection_parameters, outdir=outdir, label=label)

# Make a corner plot.
result.plot_corner()
