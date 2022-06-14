
import matplotlib.font_manager

injected_frequencies = []
m = 0
for i in range(len(pGB_injected)):   
    for j in range(len( pGB_injected[i])):
        injected_frequencies.append(pGB_injected[i][j]['Frequency'])
        m += 1
pGB_injected_reshaped = np.reshape(pGB_injected, m)
def closest(list, Number):
    aux = []
    for valor in list:
        aux.append(abs(Number-valor))

    return aux.index(min(aux))

lbls = [ r'\log A', r'\sin \beta',r'\lambda', 'f - f_{True} $ $ ($nHz$\right)', '\log \dot{f} $ $ ($Hz/s$)', r'\cos \iota', r'\phi', r'\Phi']
g = plots.get_subplot_plotter(subplot_size_ratio=9/16*0.7, subplot_size=8)
g.settings.scaling_factor = 1
g.settings.line_styles = 'tab10'
g.settings.solid_colors='tab10'
boundaries = {
    "EclipticLatitude": [-1.0, 1.0],
    "EclipticLongitude": [-np.pi, np.pi],
    "Inclination": [-1.0, 1.0],
    "InitialPhase": [0.0, 2.0 * np.pi],
    "Polarization": [0.0, 1.0 * np.pi],
}
names = parameters
parameter_pairs = [['EclipticLongitude', 'EclipticLatitude'],['Inclination', 'Amplitude'],['Frequency', 'FrequencyDerivative']]
samples = []
pGB_injected_sorted_index = []
m = 0
for i in range(len(found_sources_mp)):
    for j in range(len(found_sources_in[i])):
        save_frequency = pGB_injected[i][j]['Frequency']
        df = pd.read_csv('/home/stefan/Repositories/ldc1_evaluation_data/submission/ETH_2/GW'+str(int(np.round(save_frequency*10**8)))+save_name+'.csv')
        df['Inclination'] = np.cos(df['Inclination'].values)
        df['EclipticLatitude'] = np.sin(df['EclipticLatitude'].values)
        df['FrequencyDerivative'] = np.log10(df['FrequencyDerivative'].values)
        df['Amplitude'] = np.log10(df['Amplitude'].values)
        pGB_injected_sorted_index.append(closest(injected_frequencies, df['Frequency'][0]))
        df['Frequency'] = (df['Frequency'] - pGB_injected_reshaped[pGB_injected_sorted_index[-1]]['Frequency'] + m*2e-9) * 1e9
        samples.append(MCSamples(samples=df.to_numpy(), names = names, labels = lbls))
        samples[-1].updateSettings({'contours': [0.68, 0.95]})
        m += 1
pGB_injected_sorted = []
for i in range(len(found_sources_mp)):
    pGB_injected_sorted.append(pGB_injected_reshaped[pGB_injected_sorted_index[i]])

g.settings.num_plot_contours = 2
# 3D (scatter) triangle plot
# you can adjust the scaling factor if font sizes are too small when
# making many subplots in a fixed size (default=2 would give smaller fonts)
g.settings.scaling_factor = 2
g.plots_2d(samples, param_pairs=parameter_pairs,legend_labels=[],lws=1.5)
for n, ax in enumerate(g.subplots[:,0]):
    parameter1, parameter2 = parameter_pairs[n]
    m = 0
    for i in range(len(pGB_injected_sorted)):   
        pGB_injected_scaled = deepcopy(pGB_injected_sorted[i])
        pGB_injected_scaled['Inclination'] = np.cos(pGB_injected_scaled['Inclination'])
        pGB_injected_scaled['EclipticLatitude'] = np.sin(pGB_injected_scaled['EclipticLatitude'])
        pGB_injected_scaled['FrequencyDerivative'] = np.log10(pGB_injected_scaled['FrequencyDerivative'])
        pGB_injected_scaled['Amplitude'] = np.log10(pGB_injected_scaled['Amplitude'])
        pGB_injected_scaled['Frequency'] = m * 2
        m += 1
        ax.plot(pGB_injected_scaled[parameter1],pGB_injected_scaled[parameter2],color='black', marker = '+',zorder=1, markersize = 10, label= 'true')
        ax.plot(pGB_injected_scaled[parameter1],pGB_injected_scaled[parameter2], marker = '+',zorder=1.1, markersize = 15,alpha = 0.5, label= 'true', linewidth = 4)
    try:
        ax.set_xlim(boundaries[parameter1])
    except:
        xlim = ax.get_xlim()
        x_length = xlim[1]-xlim[0]
        ax.set_xlim([xlim[0]-x_length*0.02, xlim[1]+x_length*0.02])
    try:
        ax.set_ylim(boundaries[parameter2])
    except:
        ylim = ax.get_ylim()
        y_length = ylim[1]-ylim[0]
        ax.set_ylim([ylim[0]-y_length*0.02, ylim[1]+y_length*0.02])
    if parameter2 in ['FrequencyDerivative']:
        ax.axhline(y=np.log10(1/Tobs**2/100), color='grey', linestyle = '--', zorder = 0.5)
    # if parameter2 in ['Amplitude', 'FrequencyDerivative']:
    #     ax.set_yscale('log')
g.export('/home/stefan/LDC/LDC/pictures/global fit all '+save_name+'log frequency.png')
