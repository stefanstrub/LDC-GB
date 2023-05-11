


### determine index of a specific frequency
index_of_interest_to_plot = np.searchsorted(np.asarray(frequencies_search)[:,0],  0.00399)+25
#plot strains
number_of_windows = 4
for i in range(len(frequencies_search)):
    if i != index_of_interest_to_plot:
        continue
    search1 = Search(tdi_fs,Tobs, frequencies_search[i][0], frequencies_search[i+number_of_windows-1][1], update_noise=False)
    vertical_lines = []
    for j in range(number_of_windows):
        vertical_lines.append(frequencies_search[i+j][0])
        vertical_lines.append(frequencies_search[i+j][1])
    found_extended = np.concatenate(found_sources_matched[i:i+number_of_windows])
    found_not_matched_extended = np.concatenate(found_sources_not_matched[i:i+number_of_windows])
    # matched_extended = np.concatenate(pGB_injected_matched[i:i+number_of_windows])
    # if len(pGB_injected_SNR[i]) > 0:
    pGB_injected_dict_list = []
    matched_extended = []
    not_matched_extended = []
    for j in range(number_of_windows):
        for k in range(len(pGB_injected_SNR_sorted_overlap[i+j])):
            pGB_injected_dict_list.append({})
            for parameter in parameters:
                pGB_injected_dict_list[-1][parameter] = pGB_injected_SNR_sorted_overlap[i+j][k][parameter]
            if k > 20:
                break
    for j in range(number_of_windows):
        for k in range(len(pGB_injected_matched[i+j])):
            matched_extended.append({})
            for parameter in parameters:
                matched_extended[-1][parameter] = pGB_injected_matched[i+j][k][parameter]
    
    # padding = (frequencies_search[i+number_of_windows-1][1] - frequencies_search[i][0])/2
    # pGB_injected_not_matched_flat_df_red = pGB_injected_not_matched_flat_df[pGB_injected_not_matched_flat_df['Frequency'] > frequencies_search[i][0]-padding]
    # pGB_injected_not_matched_flat_df_red = pGB_injected_not_matched_flat_df_red[pGB_injected_not_matched_flat_df_red['Frequency'] < frequencies_search[i+number_of_windows-1][1]+padding]
    # for k in range(len(pGB_injected_not_matched_flat_df_red)):
    #     not_matched_extended.append({})
    #     for parameter in parameters:
    #         not_matched_extended[-1][parameter] = pGB_injected_not_matched_flat_df_red.iloc[k][parameter]

    # found_not_matched_extended = []
    # found_not_matched_extended = found_not_matched_extended[:1]
    # not_matched_extended = []
    # found_extended = deepcopy(matched_extended)

    # found_extended = []
    # matched_extended = [matched_extended[3],matched_extended[4]]
    # matched_extended = []
    save_name_path = SAVEPATH+'/strain added Amplitude'+ str(int(np.round(frequencies_search[i][0]*10**8))) +save_name+str(int(len(matched_extended)))+'kristen.png'
    if len(pGB_injected_SNR_sorted_overlap[i]) > 20:
        search1.plot(found_sources_in=found_extended, found_sources_not_matched = found_not_matched_extended, pGB_injected= pGB_injected_dict_list[:20],  pGB_injected_matched= matched_extended, vertical_lines= vertical_lines, saving_label =save_name_path) 
    else:
        search1.plotA(found_sources_in=found_extended, found_sources_not_matched = found_not_matched_extended, pGB_injected= pGB_injected_dict_list,  pGB_injected_matched= matched_extended, vertical_lines= vertical_lines, saving_label =save_name_path) 
        # search1.plot(found_sources_in=found_sources_mp_best[i], pGB_injected=pGB_injected[i][:10], pGB_injected_matched= matched_extended, saving_label =SAVEPATH+'/strain added'+ str(int(np.round(frequencies_search[i][0]*10**8))) +save_name+'in.png') 
