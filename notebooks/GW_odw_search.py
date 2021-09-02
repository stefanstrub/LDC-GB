from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import pylab
import bilby
from bilby.core.prior import Uniform
from bilby.gw.conversion import convert_to_lal_binary_black_hole_parameters, generate_all_bbh_parameters

from gwosc.datasets import event_gps
from gwpy.timeseries import TimeSeries
from pycbc.filter import highpass
from pycbc.waveform import get_td_waveform
import time

event = 'GW150914'
time_of_event = event_gps(event)

H1 = bilby.gw.detector.get_empty_interferometer("H1")
L1 = bilby.gw.detector.get_empty_interferometer("L1")
# Definite times in relation to the trigger time (time_of_event), duration and post_trigger_duration
post_trigger_duration = 1
duration = 2
analysis_start = time_of_event + post_trigger_duration - duration

# Use gwpy to fetch the open data
H1_analysis_data = TimeSeries.fetch_open_data(
    "H1", analysis_start, analysis_start + duration, sample_rate=4096, cache=True)

L1_analysis_data = TimeSeries.fetch_open_data(
    "L1", analysis_start, analysis_start + duration, sample_rate=4096, cache=True)
H1_analysis_data.plot()
plt.show()

H1.set_strain_data_from_gwpy_timeseries(H1_analysis_data)
L1.set_strain_data_from_gwpy_timeseries(L1_analysis_data)

psd_duration = duration * 32
psd_start_time = analysis_start - psd_duration

H1_psd_data = TimeSeries.fetch_open_data(
    "H1", psd_start_time, psd_start_time + psd_duration, sample_rate=4096, cache=True)

L1_psd_data = TimeSeries.fetch_open_data(
    "L1", psd_start_time, psd_start_time + psd_duration, sample_rate=4096, cache=True)

psd_alpha = 2 * H1.strain_data.roll_off / duration
H1_psd = H1_psd_data.psd(fftlength=duration, overlap=0, window=("tukey", psd_alpha), method="median")
L1_psd = L1_psd_data.psd(fftlength=duration, overlap=0, window=("tukey", psd_alpha), method="median")

H1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
    frequency_array=H1_psd.frequencies.value, psd_array=H1_psd.value)
L1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
    frequency_array=H1_psd.frequencies.value, psd_array=L1_psd.value)

H1.maximum_frequency = 1024
L1.maximum_frequency = 1024

fig, ax = plt.subplots()
idxs = H1.strain_data.frequency_mask
ax.loglog(H1.strain_data.frequency_array[idxs],
          np.abs(H1.strain_data.frequency_domain_strain[idxs]))
ax.loglog(H1.power_spectral_density.frequency_array[idxs],
          H1.power_spectral_density.asd_array[idxs])
ax.set_xlabel("Frequency [Hz]")
ax.set_ylabel("Strain [strain/$\sqrt{Hz}$]")
plt.show()

# specgram = H1_analysis_data.spectrogram2(fftlength=4, overlap=2, window='hann') ** (1/2.)
# plot = specgram.plot()
# ax = plot.gca()
# ax.set_yscale('log')
# ax.set_ylim(10, 1400)
# ax.colorbar(
#     clim=(1e-24, 1e-20),
#     norm="log",
#     label=r"Strain noise [$1/\sqrt{\mathrm{Hz}}$]",
# )


hq = H1_analysis_data.q_transform(frange=(20, 500), qrange=(20, 30))
plot = hq.plot()
ax = plot.gca()
ax.set_yscale('log')
ax.colorbar(label="Normalised energy")
plt.show()

start = time.time()
for i in range(1000):
    hp1, _ = get_td_waveform(approximant="SEOBNRv4_opt",
                         mass1=10,
                         mass2=10,
                         delta_t=1.0/H1.sampling_frequency,
                         f_lower=25)
print('time:',time.time()-start)


class Search():
    def __init__(self,signal_peak,tdi_fs,Tobs) -> None:
        self.tdi_fs = tdi_fs
        self.GB = fastGB.FastGB(delta_t=dt, T=Tobs)
        signal_peak_index = -signal_peak-1
        selected_frequency = f[peaks[indexes_peaks[signal_peak_index]]]
        amplitude = np.sqrt(psd[peaks[indexes_peaks[signal_peak_index]]])
        frequencyrange =  [selected_frequency - 2e-7, selected_frequency + 2e-7]
        indexes = np.argsort(p.get('Frequency'))
        index_low = np.searchsorted(p.get('Frequency')[indexes], frequencyrange[0])
        index_high = np.searchsorted(p.get('Frequency')[indexes], frequencyrange[1])
        frequencyrange =  [lower_frequency, upper_frequency]
        strongest_source_in_window = np.argmax(p.get('Amplitude')[indexes][index_low:index_high])
        index_closest = np.searchsorted(p.get('Frequency')[indexes], selected_frequency)
        self.pGB = {}
        for parameter in parameters:
            self.pGB[parameter] = p.get(parameter)[indexes][index_low:index_high][strongest_source_in_window]
        # pGB = deepcopy(pGBadded)
        # self.pGB = {'Amplitude': 3.676495e-22, 'EclipticLatitude': 0.018181, 'EclipticLongitude': 1.268061, 'Frequency': 0.01158392, 'FrequencyDerivative': 8.009579e-15, 'Inclination': 0.686485, 'InitialPhase': 4.201455, 'Polarization': 2.288223}
        print('pGB', self.pGB, signal_peak)

        self.boundaries = {
            "Amplitude": [np.log10(amplitude)-3,np.log10(amplitude)-1],
            # "Amplitude": [np.log10(self.pGB['Amplitude'])-2,np.log10(self.pGB['Amplitude'])+1],
            "EclipticLatitude": [-1.0, 1.0],
            "EclipticLongitude": [-np.pi, np.pi],
            # "Frequency": [self.pGB["Frequency"] * 0.99995, self.pGB["Frequency"] * 1.00015],
            # "Frequency": [self.pGB["Frequency"] - 3e-7, self.pGB["Frequency"] + 3e-7],
            "Frequency": frequencyrange,
            "FrequencyDerivative": [-20.0,-13.0],
            # "FrequencyDerivative": [np.log10(5e-6*self.pGB['Frequency']**(13/3)),np.log10(8e-8*self.pGB['Frequency']**(11/3))],
            "Inclination": [-1.0, 1.0],
            "InitialPhase": [0.0, 2.0 * np.pi],
            "Polarization": [0.0, 1.0 * np.pi],
        }
        if self.boundaries['FrequencyDerivative'][0] > self.boundaries['FrequencyDerivative'][1]:
            c = self.boundaries['FrequencyDerivative'][0]
            self.boundaries['FrequencyDerivative'][0] = self.boundaries['FrequencyDerivative'][1]
            self.boundaries['FrequencyDerivative'][1] = c

        print('amplitude boundaries',amplitude, 10**(self.boundaries['Amplitude'][0]), 10**(self.boundaries['Amplitude'][1]))
        previous_max = np.random.rand(8)
        previous_max[0] = np.random.rand(1)*0.1 +0.5
        previous_max[3] = np.random.rand(1)*0.1 +0.5
        i = 0
        self.pGBs = deepcopy(self.pGB)
        for parameter in parameters:
            if parameter in ["FrequencyDerivative"]:
                i -= 1
            elif parameter in ["EclipticLatitude"]:
                self.pGBs[parameter] = np.arcsin((previous_max[i] * (self.boundaries[parameter][1] - self.boundaries[parameter][0])) + self.boundaries[parameter][0])
            elif parameter in ["Inclination"]:
                self.pGBs[parameter] = np.arccos((previous_max[i] * (self.boundaries[parameter][1] - self.boundaries[parameter][0])) + self.boundaries[parameter][0])
            elif parameter in ['Amplitude',"FrequencyDerivative"]:
                self.pGBs[parameter] = 10**((previous_max[i] * (self.boundaries[parameter][1] - self.boundaries[parameter][0])) + self.boundaries[parameter][0])
            else:
                self.pGBs[parameter] = (previous_max[i] * (self.boundaries[parameter][1] - self.boundaries[parameter][0])) + self.boundaries[parameter][0]
            i += 1

        # cutoff_ratio = 1000
        # Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=self.pGBs, oversample=4, simulator="synthlisa")
        # psd_signal = np.abs(Xs.values) ** 2 + np.abs(Ys.values) ** 2 + np.abs(Zs.values) ** 2
        # highSNR = psd_signal > np.max(psd_signal) / cutoff_ratio
        # lowerindex = np.where(highSNR)[0][0] - 30
        # higherindex = np.where(highSNR)[0][-1] + 30
        # self.dataX = tdi_fs["X"].isel(f=slice(Xs.kmin, Xs.kmin + len(Xs)))[lowerindex:higherindex]
        # self.dataY = tdi_fs["Y"].isel(f=slice(Ys.kmin, Ys.kmin + len(Ys)))[lowerindex:higherindex]
        # self.dataZ = tdi_fs["Z"].isel(f=slice(Zs.kmin, Zs.kmin + len(Zs)))[lowerindex:higherindex]
        indexes = np.logical_and(tdi_fs['X'].f > frequencyrange[0]-padding, tdi_fs['X'].f < frequencyrange[1]+padding) 
        self.dataX = tdi_fs["X"][indexes]
        self.dataY = tdi_fs["Y"][indexes]
        self.dataZ = tdi_fs["Z"][indexes]

        self.DAf = (self.dataZ - self.dataX)/np.sqrt(2.0)
        self.DEf = (self.dataZ - 2.0*self.dataY + self.dataX)/np.sqrt(6.0)

        # Xs, Ys, Zs = (
        #     Xs[lowerindex:higherindex],
        #     Ys[lowerindex:higherindex],
        #     Zs[lowerindex:higherindex],
        # )
        spd_data = np.abs(self.dataX) ** 2 + np.abs(self.dataY) ** 2 + np.abs(self.dataZ) ** 2
        noise = (np.mean(spd_data[:2]) + np.mean(spd_data[-2:])).values / 2
        noise = 0  # (np.mean(spd_data).values)/2
        fmin, fmax = float(self.dataX.f[0]), float(self.dataX.f[-1] + self.dataX.attrs["df"])
        freq = np.array(self.dataX.sel(f=slice(fmin, fmax)).f)
        Nmodel = get_noise_model(noise_model, freq)
        self.Sn = Nmodel.psd(freq=freq, option="X")
        # diff = np.abs(self.dataX - Xs.values) ** 2 + np.abs(self.dataY - Ys.values) ** 2 + np.abs(self.dataZ - Zs.values) ** 2
        # p1 = float(np.sum(diff / (self.Sn + noise)) * Xs.df) / 2.0
        # p1 = -p1
        # diff = np.abs(self.dataX) ** 2 + np.abs(self.dataY) ** 2 + np.abs(self.dataZ) ** 2
        # null_pGBs = deepcopy(self.pGBs)
        # null_pGBs['Amplitude'] = 4*10**-25

    def plot(self, maxpGBs):
        plt.figure(figsize=fig_size)
        ax1 = plt.subplot(111)
        # plt.plot(dataX_training.f*1000,dataX_training.values, label='data')
        ax1.plot(self.dataX.f * 1000, self.dataX.values.real, label="data", marker="o", zorder=5)
        Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=self.pGB, oversample=4, simulator="synthlisa")
        index_low = np.searchsorted(Xs.f, self.dataX.f[0])
        Xs = Xs[index_low : index_low + len(self.dataX)]
        Ys = Ys[index_low : index_low + len(self.dataY)]
        Zs = Zs[index_low : index_low + len(self.dataZ)]
        ax1.plot(Xs.f * 1000, Xs.values.real, label="VGB", marker=".", zorder=5)
        Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=self.pGB, oversample=8, simulator="synthlisa")
        index_low = np.searchsorted(Xs.f, self.dataX.f[0])
        Xs = Xs[index_low : index_low + len(self.dataX)]
        Ys = Ys[index_low : index_low + len(self.dataY)]
        Zs = Zs[index_low : index_low + len(self.dataZ)]
        ax1.plot(Xs.f * 1000, Xs.values.real, label="VGB2", marker=".", zorder=5)
        Xs, Ys, Zs = self.GB.get_fd_tdixyz(template= maxpGBs, oversample=4, simulator="synthlisa")
        index_low = np.searchsorted(Xs.f, self.dataX.f[0])
        Xs = Xs[index_low : index_low + len(self.dataX)]
        Ys = Ys[index_low : index_low + len(self.dataY)]
        Zs = Zs[index_low : index_low + len(self.dataZ)]
        ax1.plot(Xs.f * 1000, Xs.values.real, label="start", marker=".", zorder=5)
        # ax1.plot(Xs_added2.f * 1000, Xs_added2.values.real, label="VGB2", marker=".", zorder=5)
        ax1.axvline(self.boundaries['Frequency'][0]* 1000, color= 'red')
        ax1.axvline(self.boundaries['Frequency'][1]* 1000, color= 'red')
        # ax1.plot(Xs.f * 1000, dataX.values.real - Xs.values.real, label="residual", alpha=0.8, color="red", marker=".")
        plt.legend()
        plt.show()
        # print("p true", self.loglikelihood([pGB]), "null hypothesis", self.loglikelihood([null_pGBs]))


    def loglikelihood(self, pGBs):
        for i in range(len(pGBs)):
            Xs, Ys, Zs = self.GB.get_fd_tdixyz(template=pGBs[i], oversample=4, simulator="synthlisa")
            index_low = np.searchsorted(Xs.f, self.dataX.f[0])
            if i == 0:
                Xs_total = Xs[index_low : index_low + len(self.dataX)]
                Ys_total = Ys[index_low : index_low + len(self.dataY)]
                Zs_total = Zs[index_low : index_low + len(self.dataZ)]
            else:
                Xs_total += Xs[index_low : index_low + len(self.dataX)]
                Ys_total += Ys[index_low : index_low + len(self.dataY)]
                Zs_total += Zs[index_low : index_low + len(self.dataZ)]
            if len(Xs_total) < len(self.dataX):
                a,Xs_total = xr.align(self.dataX, Xs, join='left',fill_value=0)
                a,Ys_total = xr.align(self.dataY, Ys, join='left',fill_value=0)
                a,Zs_total = xr.align(self.dataZ, Zs, join='left',fill_value=0)

        diff = np.abs(self.dataX - Xs_total.values) ** 2 + np.abs(self.dataY - Ys_total.values) ** 2 + np.abs(self.dataZ - Zs_total.values) ** 2

        p1 = float(np.sum(diff / self.Sn) * Xs_total.df) / 2.0

        return -p1


    def search(self):
        # np.random.seed(42)

        parameters_recorded = [None] * 50
        for n in range(len(parameters_recorded)):
            start = time.time()
            parameters_recorded[n] = CoordinateMC(n, self.pGBs, self.boundaries, parameters_recorded, self.loglikelihood)
            print('n',n+1,'time', int(time.time()-start), np.round(parameters_recorded[n][0]['Loglikelihood'][-1],2),len(parameters_recorded[n][0]['Loglikelihood']),np.round(self.loglikelihood([self.pGB]),3))
        # pbar = tqdm(total=len(parameters_recorded))
        # pool = mp.Pool(mp.cpu_count())
        # parameters_recorded = pool.map(CoordinateMC, [n for n in range(len(parameters_recorded))])
        # pool.close()
        # pool.join()
        # pbar.close()

        best_run = 0
        loglikelihoodofruns = np.zeros(len(parameters_recorded))
        for i in range(len(parameters_recorded)):
            loglikelihoodofruns[i] = parameters_recorded[i][0]['Loglikelihood'][-1]
        best_value = np.max(loglikelihoodofruns)
        best_run = np.argmax(loglikelihoodofruns)
        good_runs = loglikelihoodofruns > best_value*3
        indices = (-loglikelihoodofruns).argsort()[:len(good_runs)]
        pGBmodes = []
        for i in range(len(good_runs)):
            if good_runs[i]:
                pGBmodes.append({})
        indices = (-loglikelihoodofruns).argsort()[:len(pGBmodes)]
        pGBmodes = []
        for n in indices:
            pGBmodes.append([])
            for i in range(number_of_signals):
                pGBmodes[-1].append({})
                for parameter in parameters + ['Loglikelihood']:
                    if parameter == 'Loglikelihood':
                        pGBmodes[-1][i][parameter] = parameters_recorded[n][0][parameter][-1]
                    else:
                        pGBmodes[-1][i][parameter] = parameters_recorded[n][i][parameter][-1]
        if len(pGBmodes) > 10:
            for signal in range(number_of_signals):
                pGBmodes = pGBmodes[:10]
        return pGBmodes

    def optimize(self, pGBmodes):
        bounds = ()
        for signal in range(number_of_signals):
            bounds += ((0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1))
        for i in range(len(pGBmodes)):
            maxpGB = []
            boundaries_reduced = []
            pGBs01 = []

            for j in range(5):
                x = []
                for signal in range(number_of_signals):
                    if j == 0:
                        maxpGB.append({})
                        boundaries_reduced.append({})
                        for parameter in parameters:
                            maxpGB[signal][parameter] = pGBmodes[i][signal][parameter]
                    # print(maxpGB)
                    boundaries_reduced[signal] = Reduce_boundaries(maxpGB[signal], self.boundaries,ratio=0.1)
                    # if j == 0:
                    #     boundaries_reduced[signal] = deepcopy(self.boundaries)
                    pGBs01.append({})
                    for parameter in parameters:
                        if parameter in ["EclipticLatitude"]:
                            pGBs01[signal][parameter] = (np.sin(maxpGB[signal][parameter]) - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
                        elif parameter in ["Inclination"]:
                            pGBs01[signal][parameter] = (np.cos(maxpGB[signal][parameter]) - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
                        elif parameter in ['Amplitude',"FrequencyDerivative"]:
                            pGBs01[signal][parameter] = (np.log10(maxpGB[signal][parameter]) - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
                        else:
                            pGBs01[signal][parameter] = (maxpGB[signal][parameter] - boundaries_reduced[signal][parameter][0]) / (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])
                    for parameter in parameters:
                        x.append(pGBs01[signal][parameter])
                # print(loglikelihood(maxpGB))
                res = scipy.optimize.minimize(self.function, x, args=boundaries_reduced, method='SLSQP', bounds=bounds, tol=1e-10)
                for signal in range(number_of_signals):
                    maxpGB[signal] = scaletooriginal(res.x[signal*8:signal*8+8],boundaries_reduced[signal])
                # print('optimized loglikelihood', loglikelihood(maxpGB),maxpGB)
                # print('boundaries reduced', boundaries_reduced)

            best_value = self.loglikelihood(maxpGB)
            if i == 0:
                current_best_value = best_value
                current_maxpGB = maxpGB
            try:
                if current_best_value < best_value:
                    current_best_value = best_value
                    current_maxpGB = maxpGB
            except:
                pass
            print('optimized loglikelihood', self.loglikelihood(maxpGB),self.loglikelihood([self.pGB]))
        maxpGB = current_maxpGB
        print('final optimized loglikelihood', self.loglikelihood(maxpGB),self.loglikelihood([self.pGB]),maxpGB[0]['Frequency'])
        return maxpGB, self.pGB

    def function(self, pGBs01, boundaries_reduced):
        pGBs = []
        for signal in range(number_of_signals):
            pGBs.append({})
            i = 0
            for parameter in parameters:
                if parameter in ["EclipticLatitude"]:
                    pGBs[signal][parameter] = np.arcsin((pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
                elif parameter in ["Inclination"]:
                    pGBs[signal][parameter] = np.arccos((pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
                elif parameter in ['Amplitude',"FrequencyDerivative"]:
                    pGBs[signal][parameter] = 10**((pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0])
                else:
                    pGBs[signal][parameter] = (pGBs01[signal*8:signal*8+8][i] * (boundaries_reduced[signal][parameter][1] - boundaries_reduced[signal][parameter][0])) + boundaries_reduced[signal][parameter][0]
                i += 1
        p = -self.loglikelihood(pGBs)
        return p#/10**4


print('end')

def loglikelihood(likelihood,sample_parameters):
    for parameter in parameters_wo_mass:
        likelihood.parameters[parameter] = sample_parameters[parameter]
    likelihood.parameters['mass_1'], likelihood.parameters['mass_2'] = funcm1m2ofMchirpq(sample_parameters['chirp_mass'], sample_parameters['mass_ratio'])
    return(likelihood.log_likelihood())


def CoordinateMC(n):
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
    previous_best = deepcopy(best_value)
    likelihood2 = deepcopy(likelihood)
    n_trials = 50
    sample_parameters = {}
    for parameter in parameters_chirp:
        sample_parameters[parameter] = parameters_recorded[parameter][0]
    sample_parameters2 = deepcopy(sample_parameters)
    for i in range(100):
        parameter1 = parameters_chirp[i % 15]
        parameter2 = parameters_chirp[np.random.randint(0, 14)]
        while parameter2 == parameter1:
            parameter2 = parameters_chirp[np.random.randint(0, 14)]
        changeableparameters = [parameter1, parameter2]
        sample_parameters2 = deepcopy(sample_parameters)
        for j in range(n_trials):
            # change_chirp = False
            # change_ratio = False
            # if parameter1 == 'chirp_mass':
            #     change_chirp = True
            #     chirp_mass = priors[parameter1].sample(1)[0]
            # elif parameter1 == 'mass_ratio':
            #     change_ratio = True
            #     mass_ratio = priors[parameter1].sample(1)[0]
            # if parameter2 == 'chirp_mass':
            #     change_chirp = True
            #     chirp_mass = priors[parameter2].sample(1)[0]
            # elif parameter2 == 'mass_ratio':
            #     change_ratio = True
            #     mass_ratio = priors[parameter2].sample(1)[0]
            # if not(change_chirp):
            #     chirp_mass = parameters_recorded['chirp_mass'][-1]
            # if not(change_ratio):
            #     chirp_ratio = parameters_recorded['mass_ratio'][-1]
            # if change_chirp or change_ratio:
            #     likelihood2.parameters['mass_1'], likelihood2.parameters['mass_2'] = funcm1m2ofMchirpq(chirp_mass, chirp_ratio)
    
            # for parameter in changeableparameters:
            #     if parameter in parameters_wo_mass:
            #         likelihood2.parameters[parameter] = priors[parameter].sample(1)[0]
            for parameter in changeableparameters:
                sample_parameters2[parameter] = priors[parameter].sample(1)[0]

            suggestion = loglikelihood(likelihood,sample_parameters2)
            if suggestion > previous_best:
                previous_best = suggestion
                for parameter in changeableparameters:
                    sample_parameters[parameter] = sample_parameters2[parameter]
                for parameter in parameters_wo_mass:
                    likelihood.parameters[parameter] = sample_parameters[parameter]
                likelihood.parameters['mass_1'], likelihood.parameters['mass_2'] = funcm1m2ofMchirpq(sample_parameters['chirp_mass'], sample_parameters['mass_ratio'])
        print(previous_best)
        if i in [30,40,50,60,70]:
            past_mean = 0
            sum_count = 0
            for l in range(n):
                try:
                    past_mean += parameters_recorded[l][0]["Loglikelihood"][i]
                    sum_count += 1
                except:
                    pass
            try:
                past_mean = past_mean / sum_count
                if previous_best > past_mean:
                    pass
                else:
                    break
            except:
                pass

        # start = time.time()
        # print(n, i, previous_best, loglikelihood(maxpGB), maxpGB)
        # parameters_recorded1[0]["Loglikelihood"].append(loglikelihood(maxpGB))
        # for i in range(number_of_signals):
        #     for parameter in parameters:
        #         parameters_recorded1[i][parameter].append(maxpGB[i][parameter])

        # maxpGB2 = deepcopy(maxpGB)
    if previous_best > best_value:
        best_value = previous_best
    parameters_recorded[n] = parameters_recorded1
    return parameters_recorded1