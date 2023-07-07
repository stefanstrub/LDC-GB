#%%
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Ellipse
from getdist import plots, MCSamples
import scipy
from scipy.optimize import differential_evolution
import numpy as np
import xarray as xr
# from getdist import plots, MCSamples
import time
from copy import deepcopy
import multiprocessing as mp
import pandas as pd
import os
import h5py
from KDEpy import FFTKDE
import sys
sys.path.append('/cluster/home/sstrub/Repositories/LDC/lib/lib64/python3.8/site-packages/ldc-0.1-py3.8-linux-x86_64.egg')

from ldc.lisa.noise import get_noise_model
from ldc.lisa.noise import AnalyticNoise
from ldc.common.series import TimeSeries
import ldc.waveform.fastGB as fastGB
from ldc.common.tools import window
# from ldc.waveform.fastGB import fastGB
# import ldc.waveform.fastGB as fastGB

# from ldc.common.tools import compute_tdi_snr

from fastkde import fastKDE
# from sklearn.metrics import mean_squared_error
# from sklearn.gaussian_process import GaussianProcessRegressor
# from sklearn.gaussian_process.kernels import RBF
# from chainconsumer import ChainConsumer


try:
    import cupy as xp

    gpu_available = True

except (ImportError, ModuleNotFoundError) as e:
    import numpy as xp

    gpu_available = False

from gbgpu.gbgpu import GBGPU

from gbgpu.utils.constants import *

from sources import *
import subprocess as sp
import os

gb_gpu = GBGPU(use_gpu=gpu_available)

class Posterior_computer():
    def __init__(self, tdi_fs, Tobs, frequencies, maxpGB, noise=None) -> None:
        self.tdi_fs = tdi_fs
        self.Tobs = Tobs
        self.frequencies = frequencies
        self.maxpGB = maxpGB
        self.search1 = Search(self.tdi_fs,self.Tobs, self.frequencies[0], self.frequencies[1], noise=noise, gb_gpu=gb_gpu)

    def reduce_boundaries(self, plot_confidance = False):
        fisher_information_matrix = self.search1.fisher_information(self.maxpGB)
        FIM = np.zeros((len(parameters),len(parameters)))
        for i,parameter1 in enumerate(parameters):
            for j,parameter2 in enumerate(parameters):
                FIM[i,j] = fisher_information_matrix[parameter1][parameter2]
        covariance_matrix = scipy.linalg.inv(FIM)
        maxpGB01 = scaleto01(self.maxpGB, self.search1.boundaries, self.search1.parameters, self.search1.parameters_log_uniform)

        # lambda_, v = scipy.linalg.eig(covariance_matrix)
        # transf = np.zeros((len(lambda_),len(lambda_)))
        # for i in range(len(lambda_)):
            # transf[i] = v[i] * np.sqrt(lambda_[i])
        # covariance2 = v @ np.diag(lambda_) @ scipy.linalg.inv(v)
        # transf = v @ np.diag(np.sqrt(lambda_))
        # scalematrix = np.max(np.abs(transf), axis=1)
        scalematrix = np.sqrt(np.diag(covariance_matrix))
        # print('scalematrix', scalematrix)
        maxpGB01_low = {}
        maxpGB01_high = {}
        self.boundaries_reduced = {}
        sigma_multiplier = 4
        for parameter in parameters:
            maxpGB01_low[parameter] = maxpGB01[parameter] - scalematrix[parameters.index(parameter)] * sigma_multiplier 
            maxpGB01_high[parameter] = maxpGB01[parameter] + scalematrix[parameters.index(parameter)] * sigma_multiplier 
            if parameter in [ 'Amplitude', 'Inclination']:
                maxpGB01_low[parameter] = maxpGB01[parameter] - scalematrix[parameters.index(parameter)] * sigma_multiplier
                maxpGB01_high[parameter] = maxpGB01[parameter] + scalematrix[parameters.index(parameter)] * sigma_multiplier
            if parameter in [ 'InitialPhase', 'Polarization']:
                maxpGB01_low[parameter] = maxpGB01[parameter] - 0.001
                maxpGB01_high[parameter] = maxpGB01[parameter] + 0.001
            if parameter in [ 'Frequency']:
                mulitplier = 1
                if np.abs(self.search1.boundaries['Frequency'][1] - self.search1.boundaries['Frequency'][0])  > 1e-4:
                    mulitplier = 0.2
                maxpGB01_low[parameter] = maxpGB01[parameter] - scalematrix[parameters.index(parameter)]  * mulitplier
                maxpGB01_high[parameter] = maxpGB01[parameter] + scalematrix[parameters.index(parameter)] * mulitplier
            if parameter == 'FrequencyDerivative':
                # print('scale',scalematrix[parameters.index(parameter)])
                # if scalematrix[parameters.index(parameter)] > 0.07:
                print('scale', scalematrix[parameters.index(parameter)], 'fd')
                range_fd = 1
                # if self.search1.boundaries['FrequencyDerivative'][1] > 1e-14:
                #     range_fd = 0.02
                # if self.search1.boundaries['FrequencyDerivative'][1] > 1e-17:
                #     range_fd = 0.4
                maxpGB01_low[parameter] = maxpGB01[parameter] - range_fd
                maxpGB01_high[parameter] = maxpGB01[parameter] + range_fd
            if maxpGB01_low[parameter] > maxpGB01_high[parameter]:
                placeholder = deepcopy(maxpGB01_low[parameter])
                maxpGB01_low[parameter] = deepcopy(maxpGB01_high[parameter])
                maxpGB01_high[parameter] = deepcopy(placeholder)
            if maxpGB01_low[parameter] < 0:
                maxpGB01_low[parameter] = 0
            if maxpGB01_high[parameter] > 1:
                maxpGB01_high[parameter] = 1
            maxpGB_fisher_low = (maxpGB01_low[parameter] * (self.search1.boundaries[parameter][1] - self.search1.boundaries[parameter][0])) + self.search1.boundaries[parameter][0]
            maxpGB_fisher_high = (maxpGB01_high[parameter] * (self.search1.boundaries[parameter][1] - self.search1.boundaries[parameter][0])) + self.search1.boundaries[parameter][0]
            self.boundaries_reduced[parameter] = [maxpGB_fisher_low, maxpGB_fisher_high]
        # self.boundaries_reduced = deepcopy(boundaries_reduced_fisher)

        # correct Frequency Derivative
        # split_fd = -17
        # if self.boundaries_reduced['FrequencyDerivative'][1] < split_fd+1:
        #     self.boundaries_reduced['FrequencyDerivative'][0] = -18.5
        #     self.boundaries_reduced['FrequencyDerivative'][1] = -16
        # elif self.boundaries_reduced['FrequencyDerivative'][0] < split_fd+0.5:
        #     self.boundaries_reduced['FrequencyDerivative'][0] = -18.5
        # correct Inclination and Amplitude
        # split_inclination = 0.95
        # if self.boundaries_reduced['Inclination'][1] > split_inclination or self.boundaries_reduced['Inclination'][0] < -split_inclination:
        #     if self.boundaries_reduced['Inclination'][1] > split_inclination and self.boundaries_reduced['Inclination'][0] < -split_inclination:
        #         if maxpGB01['Inclination'] > 0.5:
        #             self.boundaries_reduced['Inclination'][0] = 0.5
        #             self.boundaries_reduced['Inclination'][1] = 1
        #         else:
        #             self.boundaries_reduced['Inclination'][0] = -1
        #             self.boundaries_reduced['Inclination'][1] = -0.5
        #     if self.boundaries_reduced['Inclination'][1] > split_inclination:
        #         self.boundaries_reduced['Inclination'][0] = 0.5
        #         self.boundaries_reduced['Inclination'][1] = 1
        #     if self.boundaries_reduced['Inclination'][0] < -split_inclination:
        #         self.boundaries_reduced['Inclination'][0] = -1
        #         self.boundaries_reduced['Inclination'][1] = -0.5
        #     parameter = 'Amplitude'
        #     maxpGB01_low[parameter]  = maxpGB01[parameter] - 0.08
        #     maxpGB01_high[parameter]  = maxpGB01[parameter] + 0.08
        #     if maxpGB01_low[parameter] > maxpGB01_high[parameter]:
        #         placeholder = deepcopy(maxpGB01_low[parameter])
        #         maxpGB01_low[parameter] = deepcopy(maxpGB01_high[parameter])
        #         maxpGB01_high[parameter] = deepcopy(placeholder)
        #     if maxpGB01_low[parameter] < 0:
        #         maxpGB01_low[parameter] = 0
        #     if maxpGB01_high[parameter] > 1:
        #         maxpGB01_high[parameter] = 1
        #     maxpGB_fisher_low = (maxpGB01_low[parameter] * (self.search1.boundaries[parameter][1] - self.search1.boundaries[parameter][0])) + self.search1.boundaries[parameter][0]
        #     maxpGB_fisher_high = (maxpGB01_high[parameter] * (self.search1.boundaries[parameter][1] - self.search1.boundaries[parameter][0])) + self.search1.boundaries[parameter][0]
        #     self.boundaries_reduced[parameter] = [maxpGB_fisher_low, maxpGB_fisher_high]
        # print('boundaries reduced', self.boundaries_reduced)
    
    def reduce_boundaries_from_kde(self, resolution=10**6, proposal=None):
        test_x_m, probability = self.get_samples(resolution=10**6, proposal=proposal)
        min_parameters = np.min(test_x_m, axis=0)
        max_parameters = np.max(test_x_m, axis=0)
        self.boundaries_reduced_kde = deepcopy(self.boundaries_reduced)
        for i, parameter in enumerate(parameters):
            if parameter in ['InitialPhase', 'Polarization']:
                continue
            length = self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0] 
            self.boundaries_reduced_kde[parameter][0] = min_parameters[i]*length + self.boundaries_reduced[parameter][0]
            self.boundaries_reduced_kde[parameter][1] = max_parameters[i]*length + self.boundaries_reduced[parameter][0]
        # self.boundaries_reduced = self.boundaries_reduced_kde
        return self.boundaries_reduced_kde
    
    def reduce_boundaries_from_samples(self, samples):
        min_parameters = np.min(samples, axis=0)
        max_parameters = np.max(samples, axis=0)
        self.boundaries_reduced_samples = deepcopy(self.boundaries_reduced)
        for i, parameter in enumerate(parameters):
            if parameter in ['InitialPhase', 'Polarization']:
                continue
            length = self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0] 
            self.boundaries_reduced_samples[parameter][0] = min_parameters[i]*length + self.boundaries_reduced[parameter][0]
            self.boundaries_reduced_samples[parameter][1] = max_parameters[i]*length + self.boundaries_reduced[parameter][0]
        return self.boundaries_reduced_samples
    
    def train_model(self):
        rmse = 2
        train_size = 0
        test_size = 500
        added_trainig_size = 1000
        j = 0
        samples = np.random.rand(6000,8)
        samples[0] = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]
        samples[test_size] = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]
        while rmse > 0.6 and j < 5:
            j += 1
            train_size += added_trainig_size
            if j == 1:
                resolution = train_size + test_size
            else:
                resolution = added_trainig_size
            start = time.time()
            samples_likelihood = np.zeros(resolution)
            samples_likelihood2 = np.zeros(resolution)
            incl = np.zeros(resolution)
            for i in range(resolution):
                samples_p = scaletooriginal(samples[i+j*added_trainig_size], self.boundaries_reduced, self.search1.parameters, self.search1.parameters_log_uniform)
                incl[i] = samples_p['Inclination']
                samples_likelihood[i] = self.search1.loglikelihood([samples_p])
            print('sample time of', resolution, 'samples ',time.time() - start)

            samples_flat = np.zeros((resolution  , len(parameters)))
            i = 0
            boundary_ratio = 1# (boundaries_reduced1['FrequencyDerivative'][1]-boundaries_reduced1['FrequencyDerivative'][0])/(boundaries_reduced2['FrequencyDerivative'][1]-boundaries_reduced1['FrequencyDerivative'][0])
            for parameter in parameters:
                if parameter == 'FrequencyDerivative':
                    samples_flat[:, i] = samples[j*added_trainig_size:j*added_trainig_size+resolution,parameters.index(parameter)]*boundary_ratio
                else:
                    samples_flat[:, i] = samples[j*added_trainig_size:j*added_trainig_size+resolution,parameters.index(parameter)]
                i += 1
            # samples_flat = samples_flat*2-1
            if j == 1:
                train_y = samples_likelihood[test_size:]
                test_y = samples_likelihood[:test_size]
                train_x = samples_flat[test_size:]
                test_x = samples_flat[:test_size]
            else:
                train_y = np.append(train_y, samples_likelihood)
                train_x = np.append(train_x, samples_flat,axis=0)

            self.mu = np.mean(train_y)
            self.sigma = np.std(train_y)
            train_y_normalized = (train_y - self.mu) / self.sigma
            kernel = RBF(length_scale=[1,2,5,1,1,1,1,1],length_scale_bounds=[(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,10),(0.1,100),(0.1,100)])
            start = time.time()
            self.gpr = GaussianProcessRegressor(kernel=kernel, random_state=0).fit(train_x, train_y_normalized)
            print('train',time.time() - start)
            start = time.time()
            observed_pred_sk = self.gpr.predict(test_x)
            print('eval time of ', test_size, 'samples: ',time.time() - start)
            observed_pred_sk_scaled = observed_pred_sk*self.sigma + self.mu
            rmse = np.sqrt(mean_squared_error(test_y,observed_pred_sk_scaled))
            print("RMSE ",rmse,'with training size', len(train_y))
            if rmse > 30:
                print('high RMSE')

            # fig = plt.figure(figsize=(15,6))
            # plt.scatter(train_x[:,0], train_x[:,5], c=train_y, cmap='gray')
            # plt.show()
            if rmse > 5 and j == 5:
                j = 0
                samples = np.random.rand(6000,8)
                samples[0] = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]
                samples[test_size] = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]
                train_size = 0

    def evaluate(self, x):
        partial_length = 1*10**3
        # start = time.time()
        observed_pred_mean = np.zeros(len(x))
        observed_pred_sk = np.zeros(len(x))
        for n in range(int(len(x)/partial_length)):
            observed_pred_sk[n*partial_length:(n+1)*partial_length] = self.gpr.predict(x[(n)*partial_length:(n+1)*partial_length])
        try:
            observed_pred_sk[int(len(x)/partial_length)*partial_length:] = self.gpr.predict(x[int(len(x)/partial_length)*partial_length:])
        except:
            pass
        observed_pred_sk = np.asarray(observed_pred_sk)
        observed_pred_sk = observed_pred_sk.reshape(len(x))
        observed_pred_mean[:len(x)] = observed_pred_sk[:len(x)]*self.sigma + self.mu
        # print('eval time', time.time()-start)
        return observed_pred_mean
    
    def get_loglikelihood_gpu(self, x):
        partial_length = 1*10**3
        # start = time.time()
        observed_pred_mean = np.zeros(len(x))
        observed_pred_sk = np.zeros(len(x))
        for n in range(int(len(x)/partial_length)):
            observed_pred_sk[n*partial_length:(n+1)*partial_length] = self.gpr.predict(x[(n)*partial_length:(n+1)*partial_length])
        try:
            observed_pred_sk[int(len(x)/partial_length)*partial_length:] = self.gpr.predict(x[int(len(x)/partial_length)*partial_length:])
        except:
            pass
        observed_pred_sk = np.asarray(observed_pred_sk)
        observed_pred_sk = observed_pred_sk.reshape(len(x))
        observed_pred_mean[:len(x)] = observed_pred_sk[:len(x)]*self.sigma + self.mu
        # print('eval time', time.time()-start)
        return observed_pred_mean

    def sampler(self, resolution=1000, path_len=0.01, step_size=0.01):
        x0 = np.asarray([0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5])
        samples = hamiltonian_monte_carlo(n_samples=resolution, negative_log_prob=self.logp_func, grad_log_prob=self.dlogp_func, initial_position= x0, path_len=0.1, step_size=0.01)
        return samples

    def sample_dist(self, data0, data1, numPoints, resolution):
        data = np.array([data0, data1]).T
        grid, mypdf = FFTKDE(bw=0.01, kernel='gaussian').fit(data).evaluate(numPoints)
        axes = np.unique(grid[:, 0]), np.unique(grid[:, 1])
        mypdf = mypdf.reshape(numPoints, numPoints).T
        mypdf = np.clip(mypdf, a_min=0, a_max=None)

        dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
        data, pdfs = dist(resolution)
        return data, pdfs
    
    def sample_dist_fastkde(self, data0, data1, numPoints, resolution):
        ax = np.linspace(-0.15,1.15, numPoints)
        mypdf,axes = fastKDE.pdf(data0, data1, axes=[ax,ax])
        while np.min(mypdf) < 0:
            ax_expand = (np.random.rand(1)*0.2+0.1)[0]
            ax = np.linspace(-ax_expand,1+ax_expand,numPoints)
            mypdf,axes = fastKDE.pdf(data0, data1, axes=[ax,ax])

        dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
        data, pdfs = dist(resolution)
        return data, pdfs

    def get_samples(self, resolution=1000, proposal= None):
        numPoints = 2**5+1
        if proposal is None:
            test_x_m = np.random.uniform(size=(resolution,len(parameters)))
            probability = np.ones(resolution)
        else:
            probability = np.ones(resolution)
            test_x_m = np.zeros((resolution,len(parameters)))
            # [test_x_m[:,2],test_x_m[:,1]], pdfs = self.sample_dist(proposal[:,1],proposal[:,2], numPoints, resolution)
            # probability *= pdfs

            # [test_x_m[:,5],test_x_m[:,0]], pdfs = self.sample_dist(proposal[:,0],proposal[:,5], numPoints, resolution)
            # probability *= pdfs

            # [test_x_m[:,4],test_x_m[:,3]], pdfs = self.sample_dist(proposal[:,3],proposal[:,4], numPoints, resolution)
            # probability *= pdfs

            # plt.figure()
            # plt.hist2d(test_x_m[:,5],test_x_m[:,0])
            # plt.show()

            [test_x_m[:,2],test_x_m[:,1]], pdfs = self.sample_dist_fastkde(proposal[:,1],proposal[:,2], numPoints, resolution)
            probability *= pdfs

            [test_x_m[:,5],test_x_m[:,0]], pdfs = self.sample_dist_fastkde(proposal[:,0],proposal[:,5], numPoints, resolution)
            probability *= pdfs

            [test_x_m[:,4],test_x_m[:,3]], pdfs = self.sample_dist_fastkde(proposal[:,3],proposal[:,4], numPoints, resolution)
            probability *= pdfs

            # data = np.array([proposal[:,0],proposal[:,5]]).T
            # grid, mypdf = FFTKDE(bw=0.1, kernel='gaussian').fit(data).evaluate(numPoints)
            # axes = np.unique(grid[:, 0]), np.unique(grid[:, 1])
            # mypdf = mypdf.reshape(numPoints, numPoints).T
            # mypdf = np.clip(mypdf, a_min=0, a_max=None)
            # dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
            # data, pdfs = dist(resolution)
            # test_x_m[:,0] = data[1]
            # test_x_m[:,5] = data[0]
            # plt.figure()
            # plt.scatter(np.log10(test_x_m[:10000,0]),np.arccos(test_x_m[:10000,5]* (boundaries_reduced['Inclination'][1] - boundaries_reduced['Inclination'][0]) + boundaries_reduced['Inclination'][0]), c=pdfs[:10000])
            # probability *= pdfs

            ### fastkde
            # ax = np.linspace(-0.15,1.15,numPoints)
            # mypdf,axes = fastKDE.pdf(proposal[:,1],proposal[:,2], axes=[ax,ax])
            # #a pdf can not be negative
            # while np.min(mypdf) < 0:
            #     ax_expand = (np.random.rand(1)*0.2+0.1)[0]
            #     ax = np.linspace(-ax_expand,1+ax_expand,numPoints)
            #     mypdf,axes = fastKDE.pdf(data[0], data[1], axes=[ax,ax])

            # dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
            # data, pdfs = dist(resolution)
            # test_x_m[:,1] = data[1]
            # test_x_m[:,2] = data[0]
            # probability *= pdfs

            # ax = np.linspace(-0.15,1.15,numPoints)
            # mypdf,axes = fastKDE.pdf(proposal[:,0],proposal[:,5], axes=[ax,ax])
            # #a pdf can not be negative
            # while np.min(mypdf) < 0:
            #     ax_expand = (np.random.rand(1)*0.2+0.1)[0]
            #     ax = np.linspace(-ax_expand,1+ax_expand,numPoints)
            #     mypdf,axes = fastKDE.pdf(data[0], data[1], axes=[ax,ax])
            # dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
            # data, pdfs = dist(resolution)
            # test_x_m[:,0] = data[1]
            # test_x_m[:,5] = data[0]
            # probability *= pdfs

            # ax = np.linspace(-0.15,1.15,numPoints)
            # mypdf,axes = fastKDE.pdf(proposal[:,3],proposal[:,4], axes=[ax,ax])
            # #a pdf can not be negative
            # while np.min(mypdf) < 0:
            #     ax_expand = (np.random.rand(1)*0.2+0.1)[0]
            #     ax = np.linspace(-ax_expand,1+ax_expand,numPoints)
            #     mypdf,axes = fastKDE.pdf(data[0], data[1], axes=[ax,ax])
            # dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
            # data, pdfs = dist(resolution)
            # test_x_m[:,3] = data[1]
            # test_x_m[:,4] = data[0]
            # probability *= pdfs




            # mypdf,axes = fastKDE.pdf(proposal[:,3],proposal[:,4], axes=[ax,ax])
            # while not(np.all(mypdf>=0)):
            #     proposal = self.calculate_posterior(resolution = self.previous_resolution, proposal= self.previous_mcmc_samples, temperature= self.previous_T)
            #     mypdf,axes = fastKDE.pdf(proposal[:,3],proposal[:,4], axes=[ax,ax])
            # data = np.array([proposal[:,3],proposal[:,4]]).T
            # grid, mypdf = FFTKDE(bw=0.1, kernel='gaussian').fit(data).evaluate(numPoints)
            # axes = np.unique(grid[:, 0]), np.unique(grid[:, 1])
            # mypdf = mypdf.reshape(numPoints, numPoints).T
            # mypdf = np.clip(mypdf, a_min=0, a_max=None)
            # dist = Distribution(mypdf, transform=lambda i:i/len(mypdf[0,:])*(axes[0][-1]-axes[0][0])+axes[0][0])
            # data, pdfs = dist(resolution)
            # test_x_m[:,3] = data[1]
            # test_x_m[:,4] = data[0]
            # probability *= pdfs

            test_x_m[:,6] = np.random.uniform(size=resolution)
            test_x_m[:,7] = np.random.uniform(size=resolution)

            for n in range(8):
                index = np.where(test_x_m[:,n] > 0)
                test_x_m = test_x_m[index]
                probability = probability[index]
                index = np.where(test_x_m[:,n] < 1)
                test_x_m = test_x_m[index]
                probability = probability[index]
        return test_x_m, probability

    def calculate_posterior(self,resolution = 1*10**6, proposal= None, temperature = 1, pGB_true=None):
        test_x_m, probability = self.get_samples(resolution, proposal)
        start = time.time()

        # fig =  corner.corner(test_x_m,  bins=40, hist_kwargs={'density':True, 'lw':3}, plot_datapoints=False, fill_contours=False,  show_titles=True, \
        #                 color='#348ABD',  truth_color='k', use_math_test=True, \
        #                  levels=[0.9], title_kwargs={"fontsize": 12})
        # plt.show()

        use_gpu = True
        start = time.time()
        if use_gpu:
            # test_x_m[0] = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]
            test_x_m_scaled = scaletooriginal_array(test_x_m, self.boundaries_reduced, self.search1.parameters, self.search1.parameters_log_uniform)
            test_x_m_scaled = np.swapaxes(test_x_m_scaled, 0,1)

            params = np.array(
                [test_x_m_scaled[parameters.index('Amplitude')],
                test_x_m_scaled[parameters.index('Frequency')],
                test_x_m_scaled[parameters.index('FrequencyDerivative')],
                np.zeros_like(test_x_m_scaled[parameters.index('FrequencyDerivative')]),
                test_x_m_scaled[parameters.index('InitialPhase')],
                test_x_m_scaled[parameters.index('Inclination')],
                test_x_m_scaled[parameters.index('Polarization')],
                test_x_m_scaled[parameters.index('EclipticLongitude')],
                test_x_m_scaled[parameters.index('EclipticLatitude')]]
            )
            # for i, parameter in enumerate(['Amplitude', 'Frequency', 'FrequencyDerivative', 'InitialPhase', 'Inclination', 'Polarization', 'EclipticLongitude', 'EclipticLatitude']):
            #     j = 0
            #     if i > 2:
            #         j = 1
            #     params[i+j,0] = self.search1.pGBs[parameter]
            # params = np.array(
            #     [amp_in, f0_in, fdot_in, fddot_in, phi0_in, iota_in, psi_in, lam_in, beta_sky_in,]
            # )

            # params[1,0] = 0.00443235
            observed_pred_mean = self.search1.loglikelihood_gpu(params)
        else:
            observed_pred_mean = self.evaluate(test_x_m)

        # print('time loglikelihood for ', resolution, ' signals: ', time.time()-start)
        flatsamples = np.zeros(len(test_x_m))
        flatsamplesparameters = np.zeros((len(test_x_m),len(parameters)+1))
        i = 0
        flatsamples[:] = observed_pred_mean
        flatsamplesparameters[:,1:] = test_x_m
        flatsamplesparameters[:,0] = observed_pred_mean

        maxindx = np.unravel_index(flatsamplesparameters[:,0].argmax(), flatsamplesparameters[:,0].shape)
        max_parameters = flatsamplesparameters[maxindx[0],1:]
        # max_loglike = flatsamplesparameters[:,0].max()
        maxpGBpredicted = scaletooriginal(max_parameters, self.boundaries_reduced, self.search1.parameters, self.search1.parameters_log_uniform)
        # print(self.search1.loglikelihood_gpu_signal(pGB_dict_to_gpu_input([maxpGBpredicted])), self.search1.loglikelihood_gpu_signal(pGB_dict_to_gpu_input([self.maxpGB])))
        # print(self.search1.loglikelihood_gpu(pGB_dict_to_gpu_input([maxpGBpredicted])), self.search1.loglikelihood_gpu(pGB_dict_to_gpu_input([self.maxpGB])))
        # print(self.search1.loglikelihood([maxpGBpredicted]), self.search1.loglikelihood([self.maxpGB]))
        if self.search1.loglikelihood([maxpGBpredicted]) > self.search1.loglikelihood([self.maxpGB]):
            maxpGB = maxpGBpredicted
            # print('better maxpGB is found',maxpGB,self.search1.loglikelihood([maxpGBpredicted]), self.search1.loglikelihood([self.maxpGB]), self.maxpGB)
        
        best_value = np.max(observed_pred_mean)
        # print("pred", max_loglike, "true", self.search1.loglikelihood([scaletooriginal(max_parameters, self.boundaries_reduced)]), "max", self.search1.loglikelihood([self.maxpGB]), self.maxpGB)

        # np.random.shuffle(flatsamplesparameters)
        start = time.time()
        # normalizer = np.sum(np.exp(flatsamplesparameters[:,0]-best_value))
        flatsamples_normalized = np.exp(flatsamplesparameters[:,0]-best_value)
        # flatsamples_normalized = flatsamplesparameters[:,0]
        mcmc_samples = []
        mcmc_samples.append(flatsamplesparameters[0,1:])
        previous_p = flatsamples_normalized[0]
        if previous_p == 0:
            previous_p == 1e-300
        current_parameters = flatsamplesparameters[0,1:]
        # probability = scipy.stats.multivariate_normal.pdf(flatsamplesparameters[:,1:], mean=[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5], cov=[std_scale,std_scale,std_scale,std_scale,std_scale,std_scale,std_scale,std_scale])
        
        previous_probability = probability[0]
        accepted = 0
        for i in range(len(flatsamples_normalized)-1):
            if ((flatsamples_normalized[i+1] / previous_p) * (previous_probability/probability[i+1]))**(1/temperature) > np.random.uniform():
                previous_p = flatsamples_normalized[i+1]
                previous_probability = probability[i+1]
                current_parameters = flatsamplesparameters[i+1,1:]
                mcmc_samples.append(current_parameters)
                accepted += 1
            else:
                mcmc_samples.append(current_parameters)
        mcmc_samples = np.asarray(mcmc_samples)
        # print('time MHMC', time.time()-start)
        print('acceptance rate',np.round(accepted/len(probability)*100),'%')

        return mcmc_samples
    

    def predict(self,x,k=0):
        #x of shape (m)
        
        #returns the gp predictions where f is the true function and
        #df, ddf, If, IIf are its first and second derivate respectively antiderivates
        #the outputs are the predictions f_p,df_p,ddf_p,If_p,IIf_p where
        #f(x) = f_p(x), df(x) = df_p(x), ddf(x) = ddf_p(x), If(x) = If_p(x) + C1, 
        #IIf(x) = IIf_p(x) + C1*x + C2 with some constants C1,C2
        #set k = 0 for the normal prediction, K = 1,2 for the first or second derivates
        #and k = -1,-2 for the first or second antiderivates
    
        # x = x.reshape(-1,1)
    
        X = x - self.gpr.X_train_
        l = self.gpr.kernel_.length_scale
        A = self.gpr.alpha_

        K_trans = self.gpr.kernel_(x, self.gpr.X_train_)
        y_mean = K_trans @ self.gpr.alpha_

        f = np.prod(np.exp(-(X)**2 / (2*l**2)), axis=1)
        df = f * (-X / l ** 2).T
        
        if k == 0: 
            return f @ A
        elif k == 1: 
            return df @ A
        else:
            raise Exception('Unknown parameter k: {}'.format(k))

    def gradient(self,x):
        step_length = 0.01
        grad = np.ones_like(x)
        for i in range(len(x)):
            x_predict = np.copy(x)
            x_predict[i] = x[i]-step_length/2
            if x_predict[i] < 0:
                x_predict[i] = 0
            low = self.gpr.predict([x_predict])
            x_predict[i] = x[i]+step_length/2
            if x_predict[i] > 1:
                x_predict[i] = 1
            high = self.gpr.predict([x_predict])
            grad[i] = (high-low)/step_length
        return grad

    def predict2(self, x):
        # gets 'l' used in denominator of expected value of gradient for RBF kernel 
        k2_l = self.gpr.kernel_.length_scale

        # not necessary to do predict, but now y_pred has correct shape
        y_pred, sigma = self.gpr.predict(np.asanyarray(x).reshape(1,-1),  return_std=True)

        # allocate array to store gradient
        y_pred_grad = 0.0*y_pred

        # set of points where gradient is to be queried
        # x = np.atleast_2d(np.linspace(-5, 0.8, 1000)).T
        
        X = self.gpr.X_train_
        x_star = x

        # eval_gradient can't be true when eval site doesn't match X
        # this gives standard RBF kernel evaluations
        k_val= self.gpr.kernel_(X, np.atleast_2d(x_star), eval_gradient=False).ravel()

        # x_i - x_star / l^2
        x_diff_over_l_sq = ((X-x_star)/np.power(k2_l,2)).ravel()

        # pair-wise multiply
        intermediate_result = np.multiply(k_val, x_diff_over_l_sq)

        # dot product intermediate_result with the alphas
        final_result = np.dot(intermediate_result, self.gpr.alpha_)

        # store gradient at this point
        y_pred_grad[key] = final_result
        return y_pred_grad

    def logp_func(self,x):
        return -(self.evaluate([x])[0] * self.sigma + self.mu)

    def dlogp_func(self,x, loc=0, scale=1):
        return -(self.predict(x,k=1) * self.sigma)

    def plot_corner(self, mcmc_samples, pGB = {}, save_figure = False, save_chain = False, number_of_signal = 0, parameter_titles = False, rescaled = False):
        start = time.time()
        mcmc_samples_rescaled = []
        if not(rescaled):
            for i in range(len(mcmc_samples)):
                mcmc_samples_rescaled.append(scaletooriginal_array(mcmc_samples[i], self.boundaries_reduced, self.search1.parameters, self.search1.parameters_log_uniform))
            # mcmc_samples_rescaled = np.zeros(np.shape(mcmc_samples))
            # i = 0
            # for parameter in parameters:
            #     if parameter in ["EclipticLatitude"]:
            #         mcmc_samples_rescaled[:,i] = np.arcsin((mcmc_samples[:,parameters.index(parameter)] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
            #     elif parameter in ["Inclination"]:
            #         mcmc_samples_rescaled[:,i] = np.arccos((mcmc_samples[:,parameters.index(parameter)] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
            #     elif parameter in parameters_log_uniform:
            #         mcmc_samples_rescaled[:,i] = 10**((mcmc_samples[:,parameters.index(parameter)] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0])
            #     else:
            #         mcmc_samples_rescaled[:,i] = (mcmc_samples[:,parameters.index(parameter)] * (self.boundaries_reduced[parameter][1] - self.boundaries_reduced[parameter][0])) + self.boundaries_reduced[parameter][0]
            #     i += 1
            # print('time rescale', time.time()-start)
        else:
            mcmc_samples_rescaled = mcmc_samples
        
        save_frequency = self.maxpGB['Frequency']

        if save_chain:
            for i in range(len(mcmc_samples)):
                mcmc_samples_rescaled[i][:,parameters.index('Frequency')] *= 10**3 
            df = pd.DataFrame(data=mcmc_samples_rescaled[0], columns=parameters)
            for parameter in ['Amplitude','FrequencyDerivative','EclipticLatitude','EclipticLongitude','Inclination','InitialPhase','Polarization']:
                df[parameter] = df[parameter].astype('float32')
            df.to_csv(SAVEPATH+'Chains_gpu/frequency'+str(int(np.round(save_frequency*10**12)))+'pHz'+save_name+'.csv',index=False)
            for i in range(len(mcmc_samples)):
                mcmc_samples_rescaled[i][:,parameters.index('Frequency')] /= 10**3 
        # start = time.time()
        # df = pd.DataFrame(data=mcmc_samples_rescaled, columns=parameters)
        # df.to_csv('/home/stefan/Repositories/ldc1_evaluation_data/submission/Stefan_LDC14/GW'+str(int(np.round(maxpGB['Frequency']*10**8)))+'.csv',index=False)
        # print('saving time', time.time()-start)

        if save_figure:
            # print('full time', time.time()-first_start)

            lbls = [r'\lambda', r'\sin \beta', 'f$ $($mHz$)', r'\dot{f}$ $ ($Hz/s$)', r'\cos \iota', r'A', r'\phi', r'\Phi']

            if pGB:
                tr_s = np.zeros(len(parameters))
                i = 0
                for parameter in ['EclipticLongitude','EclipticLatitude','Frequency','FrequencyDerivative','Inclination','Amplitude','InitialPhase','Polarization']:
                    if parameter in parameters_log_uniform:
                        tr_s[i] = np.log10(pGB[parameter])
                    elif parameter in ['Frequency']:
                        tr_s[i] = pGB[parameter]*10**3
                    elif parameter in ['Inclination']:
                        tr_s[i] = np.cos(pGB[parameter])
                    elif parameter in ['EclipticLatitude']:
                        tr_s[i] = np.sin(pGB[parameter])
                    else:
                        tr_s[i] = pGB[parameter]
                    i += 1
            maxvalues = np.zeros(len(parameters))
            i = 0
            for parameter in ['EclipticLongitude','EclipticLatitude','Frequency','FrequencyDerivative','Inclination','Amplitude','InitialPhase','Polarization']:
                if parameter in parameters_log_uniform:
                    maxvalues[i] = np.log10(self.maxpGB[parameter])
                elif parameter in ['Frequency']:
                    maxvalues[i] = self.maxpGB[parameter]*10**3
                elif parameter in ['Inclination']:
                    maxvalues[i] = np.cos(self.maxpGB[parameter])
                elif parameter in ['EclipticLatitude']:
                    maxvalues[i] = np.sin(self.maxpGB[parameter])
                else:
                    maxvalues[i] = self.maxpGB[parameter]
                i += 1

            # rng = []
            # for i in range(len(lbls)):
            #     minrange = min(datS[:,i].min(), tr_s[i])
            #     maxrange = max(datS[:,i].max(), tr_s[i])
            #     range_width = np.abs(maxrange - minrange)
            #     oner = ( minrange- range_width/10, maxrange + range_width/10)
            #     rng.append(oner)
            # Get the getdist MCSamples objects for the samples, specifying same parameter
            # names and labels; if not specified weights are assumed to all be unity
            ndim = 2
            names = ['EclipticLongitude','EclipticLatitude','Frequency','FrequencyDerivative','Inclination','Amplitude']
            labels =  lbls[:ndim]
            samples = []
            for i in range(len(mcmc_samples_rescaled)):
                datS = np.zeros(np.shape(mcmc_samples_rescaled[i]))
                datS[:,0] = mcmc_samples_rescaled[i][:,2]
                datS[:,1] = np.sin(mcmc_samples_rescaled[i][:,1])
                datS[:,2] = mcmc_samples_rescaled[i][:,3]*10**3
                datS[:,3] = mcmc_samples_rescaled[i][:,4]
                datS[:,4] = np.cos(mcmc_samples_rescaled[i][:,5])
                datS[:,5] = np.log10(mcmc_samples_rescaled[i][:,0])
                datS[:,6] = mcmc_samples_rescaled[i][:,6]
                datS[:,7] = mcmc_samples_rescaled[i][:,7]
                samples.append(MCSamples(samples=datS[:,:ndim],names = names[:ndim], labels = labels[:ndim]))

                # g = plots.get_subplot_plotter(subplot_size=0.9)
                g = plots.get_subplot_plotter()
                samples[-1].updateSettings({'contours': [0.68, 0.95]})
                # g.settings.num_plot_contours = 3

            if parameter_titles:
                g.triangle_plot(samples, shaded=True, title_limit=2)
            else:
                g.triangle_plot(samples, shaded=False, legend_labels=['0.5 yr', '1 yr', '2 yr'])
            
            #markers vertical
            for i in range(ndim):
                for ax in g.subplots[i:,i]:
                    if pGB:
                        xlim = ax.get_xlim()
                        ax.set_xlim(np.min([xlim[0], tr_s[i]]),np.max([xlim[1], tr_s[i]]))
                        xlim = ax.get_xlim()
                        if xlim[0] + (xlim[1]-xlim[0])*0.1 > tr_s[i]:
                            ax.set_xlim(xlim[0] - (xlim[1]-xlim[0])*0.1, xlim[1])
                        if xlim[1] - (xlim[1]-xlim[0])*0.1 < tr_s[i]:
                            ax.set_xlim(xlim[0], xlim[1] + (xlim[1]-xlim[0])*0.1)
                        ax.axvline(tr_s[i], color='black', ls='--', lw = 1)
                    # ax.axvline(maxvalues[i], color='green', ls='--', lw = 1)
                # i += 1
            #markers horizontal
            for i in range(ndim):
                for ax in g.subplots[i,:i]:
                    if pGB:
                        ylim = ax.get_ylim()
                        ax.set_ylim(np.min([ylim[0], tr_s[i]]),np.max([ylim[1], tr_s[i]]))
                        ylim = ax.get_ylim()
                        if ylim[0] + (ylim[1]-ylim[0])*0.1 > tr_s[i]:
                            ax.set_ylim(ylim[0] - (ylim[1]-ylim[0])*0.1, ylim[1])
                        if ylim[1] - (ylim[1]-ylim[0])*0.1 < tr_s[i]:
                            ax.set_ylim(ylim[0], ylim[1] + (ylim[1]-ylim[0])*0.1)
                        # ax.axhline(tr_s[i], color='red', lw = 1)
                        ax.axhline(tr_s[i], color='black', ls='--', lw = 1)
                    # ax.axhline(maxvalues[i], color='green', ls='--', lw = 1)
                # i += 1
            g.export(SAVEPATH+'Posteriors_gpu/frequency'+ str(int(np.round(save_frequency*10**12)))+save_name+str(parameter_titles)+'.png')
            # plt.show()


# customized settings
plot_parameter = {  # 'backend': 'ps',
    "font.family": "DeJavu Serif",
    "font.serif": "Times",
    "font.size": 16,
    "mathtext.fontset": "cm",
    "axes.labelsize": "medium",
    "axes.titlesize": "medium",
    "legend.fontsize": "medium",
    "xtick.labelsize": "medium",
    "ytick.labelsize": "medium",
    "grid.color": "k",
    "grid.linestyle": ":",
    "grid.linewidth": 0.5,
    "savefig.dpi": 150,
}

# tell matplotlib about your param_plots
rcParams.update(plot_parameter)
# set nice figure sizes
fig_width_pt = 1.5*464.0  # Get this from LaTeX using \showthe\columnwidth
golden_mean = (np.sqrt(5.0) - 1.0) / 2.0  # Aesthetic ratio
ratio = golden_mean
inches_per_pt = 1.0 / 72.27  # Convert pt to inches
fig_width = fig_width_pt * inches_per_pt  # width in inches
fig_height = fig_width * ratio  # height in inches
fig_size = [fig_width, fig_height]
fig_size_squared = [fig_width, fig_width]
rcParams.update({"figure.figsize": fig_size})
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']


lbls = [r'\lambda', r'\sin \beta', 'f$ $($mHz$)', r'\log \dot{f}$ $ ($Hz/s$)', r'\cos \iota', r'A', r'\phi', r'\Phi']

# LDC1-4 ####################
def read_in_samples(chain_save_name):
    # if i != 5:
    #     continue
    # chain_save_name = SAVEPATH+'/Chain/397793number of singal0LDC1-4_4mHz.csv'
    # chain_save_name = SAVEPATH+'/Chain/frequency1666286nHzLDC1-3fastGB.csv'
    df = pd.read_csv(chain_save_name)
    # df['Inclination'] = df['Inclination'].values
    # df['EclipticLatitude'] = df['EclipticLatitude'].values
    # df['FrequencyDerivative'] = df['FrequencyDerivative'].values
    # df['Amplitude'] = df['Amplitude'].values
    mcmc_samples = df.to_numpy()
    names = ['EclipticLongitude','EclipticLatitude','Frequency','FrequencyDerivative','Inclination','Amplitude','InitialPhase','Polarizatoin']
    return mcmc_samples


parameters = [
    "Amplitude",
    "EclipticLatitude",
    "EclipticLongitude",
    "Frequency",
    "FrequencyDerivative",
    "Inclination",
    "InitialPhase",
    "Polarization",
]
# parameters_log_uniform = ['Amplitude','FrequencyDerivative']
parameters_log_uniform = ['Amplitude']
parameters_no_amplitude = parameters[1:]
intrinsic_parameters = ['EclipticLatitude','EclipticLongitude','Frequency', 'FrequencyDerivative']

# get current directory
path = os.getcwd()
 
# parent directory
parent = os.path.dirname(path)
# grandparent directory
grandparent = os.path.dirname(parent)
Radler = True
if Radler:
    DATAPATH = grandparent+"/LDC/Radler/data"
    SAVEPATH = grandparent+"/LDC/pictures/LDC1-4/"
    SAVEPATH = grandparent+"/LDC/Radler/LDC1-4_evaluation/"
else:
    DATAPATH = grandparent+"/LDC/Sangria/data"
    SAVEPATH = grandparent+"/LDC/pictures/Sangria/"

if Radler:
    sangria_fn = DATAPATH + "/LDC1-4_GB_v2.hdf5"
    # sangria_fn = DATAPATH + "/LDC1-3_VGB_v2_FD_noiseless.hdf5"
    # sangria_fn = DATAPATH + "/LDC1-3_VGB_v2.hdf5"
else:
    sangria_fn = DATAPATH + "/LDC2_sangria_training_v2.h5"
fid = h5py.File(sangria_fn)

reduction = 4
# get TDI 
if Radler:
    td = np.array(fid["H5LISA/PreProcess/TDIdata"])
    td = np.rec.fromarrays(list(td.T), names=["t", "X", "Y", "Z"])
    dt = float(np.array(fid['H5LISA/GWSources/GalBinaries']['Cadence']))
    Tobs = float(int(np.array(fid['H5LISA/GWSources/GalBinaries']['ObservationDuration']))/reduction)
else:
    td = fid["obs/tdi"][()]
    td = np.rec.fromarrays(list(td.T), names=["t", "X", "Y", "Z"])
    td = td['t']
    dt = td["t"][1]-td["t"][0]
    
    td_mbhb = fid["sky/mbhb/tdi"][()]
    # cat_mbhb = fid["sky/mbhb/cat"]
    td_mbhb  = np.rec.fromarrays(list(td_mbhb .T), names=["t", "X", "Y", "Z"])
    td_mbhb  = td_mbhb ['t']
    # tdi_ts_mbhb = dict([(k, TimeSeries(td_mbhb[k][:int(len(td_mbhb[k][:])/reduction)], dt=dt)) for k in ["X", "Y", "Z"]])
    # tdi_fs_mbhb = xr.Dataset(dict([(k, tdi_ts_mbhb[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))

    Tobs = float(int(np.array(fid['obs/config/t_max']))/reduction)
    for k in ["X", "Y", "Z"]:
        td[k] = td[k] - td_mbhb[k]


# Build timeseries and frequencyseries object for X,Y,Z
tdi_ts = dict([(k, TimeSeries(td[k][:int(len(td[k][:])/reduction)], dt=dt, t0=td.t[0])) for k in ["X", "Y", "Z"]])
tdi_fs = xr.Dataset(dict([(k, tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))
GB = fastGB.FastGB(delta_t=dt, T=Tobs)  # in seconds

f_Nyquist = 1/dt/2
search_range = [0.0003, f_Nyquist]
frequencies_search = create_frequency_windows(search_range, Tobs)

save_names = ['Radler_6m', 'Radler_12m', 'Radler_24m']
save_name = 'Radler_6m'
folderpath = SAVEPATH + '/Chain'

found_sources_in_flat_frequency_list = []
found_sources_in_flat_list = []
for i in range(len(save_names)):
    found_sources_in_flat = np.load(SAVEPATH+'found_sources_' +save_names[i]+'_flat.npy', allow_pickle = True)
    found_sources_in_flat_frequency = []
    for i in range(len(found_sources_in_flat)):
        found_sources_in_flat_frequency.append(found_sources_in_flat[i]['Frequency'])
    sorted_index = np.argsort(found_sources_in_flat_frequency)
    found_sources_in_flat_frequency = np.take(found_sources_in_flat_frequency,sorted_index)
    found_sources_in_flat = np.take(found_sources_in_flat,sorted_index)
    found_sources_in_flat_frequency_list.append(found_sources_in_flat_frequency)
    found_sources_in_flat_list.append(found_sources_in_flat)

end_string = '_SNR_scaled_03_injected_snr5'
pGB_injected_matched_flat_df = pd.read_pickle(SAVEPATH+'/injected_matched_windows_' +save_name+end_string+'_df')
pGB_injected_matched_flat_frequency = pGB_injected_matched_flat_df['Frequency']
pGB_injected_matched_flat = pGB_injected_matched_flat_df.to_records(index=False)

# onlyfiles = [f for f in os.listdir(chain_path) if os.path.isfile(os.path.join(chain_path, f))]
for k in range(17,18):
    samples_list = []
    for i in range(len(save_names)):
        print(i)
        # samples = read_in_samples(chain_path+'/'+onlyfiles[i])
        # end_index = onlyfiles[i].find('pHz')
        # frequency_max = int(onlyfiles[i][9:end_index])/10**12
        # frequency_max = found_sources_in_flat_frequency[j]

        chain_path = SAVEPATH + 'Chains_gpu_partial07_'+save_names[i]+'/'
        frequency_max = 0.004+0.00001*k
        # frequency_max = 0.004220352735
        # frequency_max = 0.009
        found_source_index = np.searchsorted(found_sources_in_flat_frequency_list[i], frequency_max)
        if found_source_index == len(found_sources_in_flat_frequency_list[i]):
            found_source_index -= 1
        if np.abs(found_sources_in_flat_frequency_list[i][found_source_index-1] - frequency_max) < np.abs(found_sources_in_flat_frequency_list[i][found_source_index] - frequency_max):
            print('s')
            found_source_index = found_source_index-1
        frequency_chain = found_sources_in_flat_list[i][found_source_index]['Frequency']
        # frequency_max = int(10049876211)/10**12
        samples = read_in_samples(chain_path+'/frequency'+str(int(np.round(frequency_chain*10**12)))+'pHz'+save_names[i]+'.csv')
        samples[:,3] /= 1000 
        for j in range(len(samples)):
            if samples[j,2] < 0:
                samples[j,2] += 2*np.pi 
        samples_list.append(samples)  

    injected_index = np.searchsorted(pGB_injected_matched_flat_frequency, frequency_max)
    if injected_index == len(pGB_injected_matched_flat_frequency):
        injected_index -= 1
    if np.abs(pGB_injected_matched_flat_frequency[injected_index-1] - frequency_max) < np.abs(pGB_injected_matched_flat_frequency[injected_index] - frequency_max):
        injected_index = injected_index-1
    frequencies_search = np.asarray(frequencies_search)
    frequencies_search_index = np.searchsorted(frequencies_search[:,0], frequency_max)
    if frequencies_search_index == len(frequencies_search[:,0]):
        frequencies_search_index -= 1
    if np.abs(frequencies_search[frequencies_search_index-1,0] - frequency_max) < np.abs(frequencies_search[frequencies_search_index,0] - frequency_max):
        frequencies_search_index = frequencies_search_index-1
    print(pGB_injected_matched_flat[injected_index]['Frequency'],pGB_injected_matched_flat[injected_index+1]['Frequency'],frequency_max)
    posterior1 = Posterior_computer(tdi_fs=tdi_fs, Tobs=Tobs, frequencies=frequencies_search[frequencies_search_index], maxpGB=pGB_injected_matched_flat[injected_index])
    # posterior1.reduce_boundaries()

    posterior1.plot_corner(samples_list, pGB=pGB_injected_matched_flat[injected_index], save_figure=True, rescaled=True)

print('end')