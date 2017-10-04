#  function for plotting significance regions
def plot_significance_lines(data, time_points, offset, slope, p_value_cutoff = 0.05, pal = None):
    """plot_significance_lines takes , or regions, or something, and calculates cluster-based p-values against 0.
        data numpy.array, subjects by time by conditions
        offset float, offset in y position for lines in plot
        slope float, change in y position for consecutive lines in plot
        p_value_cutoff float, below which p_value to plot.
    """
    import matplotlib.pyplot as pl
    import numpy as np
    import seaborn as sn
    import mne

    if pal is None:
        pal = sn.dark_palette("green", data.shape[-1])

    for i in range(data.shape[-1]): # loop across regions
        clusters = mne.stats.permutation_cluster_1samp_test(data[...,i])
        for cluster_times, p_val in zip (clusters[1], clusters[2]):
            if p_val < p_value_cutoff:
                s = np.arange(time_points.shape[0])[cluster_times]
                pl.plot([time_points[s[0]], time_points[s[-1]]], [offset + slope * i, offset + slope * i], c = pal[i], linewidth = 3.0, alpha = 0.8)


def plot_deco_results(all_deco_files, subj_data, roi_name, interval = [-3,15], output_filename = '', 
    stop_FIR_amplitude_range = [0,0], stop_FIR_pe_range = [0,0]):
    import matplotlib.pyplot as pl
    import seaborn as sn
    import pandas as pd
    import numpy as np
    from IPython import embed as shell
    import statsmodels.api as sm


    all_data = np.array([np.loadtxt(df) for df in all_deco_files])
    timepoints = np.linspace(interval[0],interval[1],all_data.shape[-1])

    sets = [[0],[1],[2]]
    all_event_names = np.array([
                'correct', 'succesful_stop', 'failed_stop' ])
    colors = np.array(['g','r','orange'])

    color_dict = dict(zip(all_event_names, colors))

    sn.set_style('ticks')
    f = pl.figure(figsize = (5,11))
    s = f.add_subplot(3,1,1)
    s.set_title(roi_name + ' gain')
    s.axhline(0, c='k', lw = 0.25)
    s.axvline(0, c='k', lw = 0.25)
    s.set_xlabel('Time [s]')
    s.set_ylabel('BOLD % signal change')

    min_d = all_data.transpose((0,2,1))[:,:,[0,1,2]].mean(axis = 0).min()
    for x in range(len(sets)):
        sn.tsplot(all_data.transpose((0,2,1))[:,:,sets[x]], time = timepoints, condition = all_event_names[sets[x]], legend = True, ax = s, color = color_dict)
        plot_significance_lines(all_data.transpose((0,2,1))[:,:,sets[x]], time_points = timepoints, offset=0.0125+stop_FIR_amplitude_range[0], slope=0.025, p_value_cutoff = 0.05, pal = colors[sets[x]])
    s.set_ylim(stop_FIR_amplitude_range)
    s.set_xlim([timepoints[0], timepoints[-1]])
    sn.despine(offset = 10, ax = s)

    pl.legend()

    ##############################################################################################################
    #
    # Now, we compute the correlations with Beta
    #
    #
    ##############################################################################################################
    sig = -np.log10(0.0125)
    # if roi_name == 'Caudate':
    #     shell()

    # X = np.vstack([np.ones(len(all_deco_files)), np.array(subj_data['Beta'], dtype = float), np.array(subj_data['alphaL'], dtype = float), np.array(subj_data['alphaG'], dtype = float)]).T
    ssrt = np.array(subj_data['SSRT'], dtype = float)
    ssrt = (ssrt - ssrt.mean()) / ssrt.std()

    # X = np.vstack([np.ones(len(all_deco_files)), ssrt]).T

    # beta
    beta = np.array(subj_data['Beta'], dtype = float)
    beta = (beta - beta.mean()) / beta.std()

    X = np.vstack([np.ones(len(all_deco_files)), ssrt, beta]).T

    # shell()
    p_T_vals = np.zeros((len(all_event_names),all_data.shape[-1], X.shape[1]+1))
    tcs = np.zeros((len(all_event_names),all_data.shape[-1], X.shape[1]))
    for i, et in enumerate(all_event_names):
        for x in range(all_data.shape[-1]):
            model = sm.OLS(np.squeeze(all_data[:,i,x]),X)
            results = model.fit()
            p_T_vals[i,x,:X.shape[1]] = -np.log10(results.pvalues)
            p_T_vals[i,x,-1] = -np.log10(results.f_pvalue)
            tcs[i,x] = results.params

    for i, c in enumerate(['SSRT', 'Beta']): # , 'Beta'
        s = f.add_subplot(3,1,2+i)
        s.set_title(roi_name + ' corrs %s'%c)
        s.axhline(0, c='k', lw = 0.25)
        s.axvline(0, c='k', lw = 0.25)
        s.set_xlabel('Time [s]')
        s.set_ylabel('beta values')        
        for j, en in enumerate(all_event_names):
            data = tcs[j,:,i+1]
            pl.plot(timepoints, data, color = colors[j], label = en)

            sig_periods = p_T_vals[j,:,i+1] > sig

            # take care of start and end of deconvolution interval.
            if sig_periods[0] == True:
                sig_periods[0] = False
            if sig_periods[-1] == True:
                sig_periods[-1] = False

            nr_blocks = int(np.floor(np.abs(np.diff(sig_periods.astype(int))).sum() / 2.0))
            if nr_blocks > 0:
                print('# blocks found: %i'%nr_blocks)
                for b in range(nr_blocks):
                    time_sig = np.arange(timepoints.shape[0])[np.r_[False, np.abs(np.diff(sig_periods)) > 0]][b*2:(b*2)+2]
                    # time_sig = np.arange(timepoints.shape[0])[sig_periods]
                    pl.plot([timepoints[time_sig[0]]-0.5, timepoints[time_sig[-1]] + 0.5], [0.0125 + stop_FIR_pe_range[0]+0.0125*j, 0.0125 + stop_FIR_pe_range[0]+0.0125*j], color = colors[j], linewidth = 3.0, alpha = 0.8)
        s.set_ylim(stop_FIR_pe_range)
        s.set_xlim([timepoints[0], timepoints[-1]])
        pl.legend()
        sn.despine(offset = 10, ax = s)

    pl.tight_layout()
    # pl.show()
    pl.savefig(output_filename)

def plot_deco_results_covariates_per_event_type(all_deco_files, subj_data, 
            roi, interval = [-3,15], output_filename = '', 
            sj_covariates = {},
            stop_FIR_amplitude_range = [0,0], stop_FIR_pe_range = [0,0]):
    import matplotlib.pyplot as pl
    import seaborn as sn
    import pandas as pd
    import numpy as np
    from IPython import embed as shell
    import statsmodels.api as sm

    # shell()

    stats_threshold = 0.0125
    # all_data = np.array([np.loadtxt(df) for df in all_deco_files])
    all_data = [pd.read_csv(df, sep = '\t', index_col=0, header=0).T for df in all_deco_files]
    timepoints = np.array(all_data[0].index, dtype = float)

    all_data_np = np.array([np.array(ad[list(sj_covariates.keys())]) for ad in all_data])
    condition_names = list(sj_covariates.keys())

    ##############################################################################################################
    #
    # Now, we compute the correlations / betas
    #
    ##############################################################################################################

    # construct across subjects covariate design matrices
    X = []
    for i, nec in enumerate(condition_names):
        Xt = np.ones((len(sj_covariates[nec])+1, len(subj_data)))
        for j, cov_name in enumerate(sj_covariates[nec]):
            this_cov = np.array(subj_data[cov_name], dtype=float)
            # z-score
            Xt[j+1,:] = (this_cov - this_cov.mean()) / this_cov.std()
        X.append(Xt.T)

    # across subjects GLM
    p_T_dict = {}
    beta_dict = {}
    for i, nec in enumerate(condition_names):
        p_T_vals = np.zeros((all_data_np.shape[1], X[i].shape[1]+1))
        betas = np.zeros((all_data_np.shape[1], X[i].shape[1]))
        for x in range(all_data_np.shape[1]):
            model = sm.OLS(np.squeeze(all_data_np[:,x,i]),X[i])
            results = model.fit()
            p_T_vals[x,:X[i].shape[1]] = -np.log10(results.pvalues)
            p_T_vals[x,-1] = -np.log10(results.f_pvalue)
            betas[x] = results.params
        p_T_dict.update({nec:pd.DataFrame(p_T_vals, index = timepoints, columns = ['int'] + sj_covariates[nec] + ['all_F'])})
        beta_dict.update({nec:pd.DataFrame(betas, index = timepoints, columns = ['int'] + sj_covariates[nec])})

    ##############################################################################################################
    #
    # Plotting
    #
    ##############################################################################################################

    color_dict = dict(zip(['correct','succesful_stop','Failed_stop'], ['g','r','orange']))
    cv_color_dict = dict(zip(['int','SSRT','ac','med','Beta'],['k','r','g','orange','brown']))

    sig = -np.log10(stats_threshold)
    sn.set_style('ticks')
    f = pl.figure(figsize = (5,len(sj_covariates)*5+3))
    s = f.add_subplot(len(sj_covariates)+1,1,1)
    s.set_title(roi + ' gain')
    s.axhline(0, c='k', lw = 0.25)
    s.axvline(0, c='k', lw = 0.25)
    s.set_xlabel('Time [s]')
    s.set_ylabel('BOLD % signal change')

    min_d = all_data_np.min()
    sn.tsplot(all_data_np, time = timepoints, condition = condition_names, legend = True, ax = s, color = color_dict)
    # plot_significance_lines(all_data_np, time_points = timepoints, offset=0.0125+rl_test_FIR_amplitude_range[0], slope=0.025, p_value_cutoff = 0.05, pal = [color_dict[cn] for cn in condition_names])
    s.set_ylim(stop_FIR_amplitude_range)
    sn.despine(offset = 10, ax = s)
    pl.legend()

    # plotting betas and stats
    for i, nec in enumerate(condition_names):
        s = f.add_subplot(len(sj_covariates)+1,1,2+i)
        s.set_title(roi + ' corrs %s'%nec)
        s.axhline(0, c='k', lw = 0.25)
        s.axvline(0, c='k', lw = 0.25)
        s.set_xlabel('Time [s]')
        s.set_ylabel('beta values')        
        betas = beta_dict[nec]
        p_T_vals = p_T_dict[nec]
        for j, cn in enumerate(p_T_dict[nec].columns):
            if cn != 'all_F':
                # find out which color to use
                this_color = [cv_color_dict[ck] for ck in cv_color_dict.keys() if ck in cn][0]
                pl.plot(timepoints, betas[cn], c=this_color, label=cn)

                # take care of start and end of deconvolution interval.
                sig_periods = np.array(p_T_vals[cn]) > sig
                if sig_periods[0] == True:
                    sig_periods[0] = False
                if sig_periods[-1] == True:
                    sig_periods[-1] = False

                nr_blocks = int(np.floor(np.abs(np.diff(sig_periods.astype(int))).sum() / 2.0))
                if nr_blocks > 0:
                    print('# blocks found: %i'%nr_blocks)
                    for b in range(nr_blocks):
                        time_sig = np.arange(timepoints.shape[0])[np.r_[False, np.abs(np.diff(sig_periods)) > 0]][b*2:(b*2)+2]
                        pl.plot([timepoints[time_sig[0]]-0.5, timepoints[time_sig[-1]] + 0.5], [0.0125 + stop_FIR_pe_range[0]+0.0125*j, 0.0125 + stop_FIR_pe_range[0]+0.0125*j], color = this_color, linewidth = 3.0, alpha = 0.8)
        s.set_ylim(stop_FIR_pe_range)
        s.set_xlim([timepoints[0], timepoints[-1]])
        pl.legend()
        sn.despine(offset = 10, ax = s)

    pl.tight_layout()
    # shell()

    pl.savefig(output_filename)


def plot_deco_results_for_publication(all_deco_files, subj_data, 
            roi, interval = [-2,10], output_filename = '', 
            sj_covariates = {},
            stop_FIR_amplitude_range = [0,0], stop_FIR_pe_range = [0,0], 
            second_plot_covariate = 'SSRT', slow_fast_condition = 'succesful_stop'):
    import matplotlib.pyplot as pl
    import seaborn as sn
    import pandas as pd
    import numpy as np
    from IPython import embed as shell
    import statsmodels.api as sm
    import scipy.stats

    stats_threshold = 0.0125
    all_data = [pd.read_csv(df, sep = '\t', index_col=0, header=0).T for df in all_deco_files]
    timepoints = np.array(all_data[0].index, dtype = float)

    condition_names = ['correct', 'succesful_stop', 'Failed_stop'] # -> don't leave to chance but copied from roi.py where the deconv happens
    all_data_np = np.array([np.array(ad[condition_names]) for ad in all_data])

    ##############################################################################################################
    #
    # Now, we compute the correlations / betas
    #
    ##############################################################################################################

    # construct across subjects covariate design matrices
    X = []
    for i, nec in enumerate(condition_names):
        Xt = np.ones((len(sj_covariates[nec])+1, len(subj_data)))
        for j, cov_name in enumerate(sj_covariates[nec]):
            this_cov = np.array(subj_data[cov_name], dtype=float)
            # z-score
            Xt[j+1,:] = (this_cov - this_cov.mean()) / this_cov.std()
        X.append(Xt.T)

    # across subjects GLM
    p_T_dict = {}
    beta_dict = {}
    for i, nec in enumerate(condition_names):
        p_T_vals = np.zeros((all_data_np.shape[1], X[i].shape[1]+1))
        betas = np.zeros((all_data_np.shape[1], X[i].shape[1]))
        for x in range(all_data_np.shape[1]):
            model = sm.OLS(np.squeeze(all_data_np[:,x,i]),X[i])
            results = model.fit()
            p_T_vals[x,:X[i].shape[1]] = -np.log10(results.pvalues)
            p_T_vals[x,-1] = -np.log10(results.f_pvalue)
            betas[x] = results.params
        p_T_dict.update({nec:pd.DataFrame(p_T_vals, index = timepoints, columns = ['int'] + sj_covariates[nec] + ['all_F'])})
        beta_dict.update({nec:pd.DataFrame(betas, index = timepoints, columns = ['int'] + sj_covariates[nec])})

    ##############################################################################################################
    #
    # Plotting #1
    #
    ##############################################################################################################

    color_dict = dict(zip(condition_names, ['g','r','orange']))
    cv_color_dict = dict(zip(['int','SSRT','ac','med','Beta'],['k','r','g','orange','brown']))

    sig = -np.log10(stats_threshold)
    sn.set_style('ticks')
    f = pl.figure(figsize = (4,12))
    s = f.add_subplot(3,1,1)
    s.set_title(roi + ' gain')
    s.axhline(0, c='k', lw = 0.25)
    s.axvline(0, c='k', lw = 0.25)
    s.set_xlabel('Time [s]')
    s.set_ylabel('BOLD % signal change')

    min_d = all_data_np.min()
    sn.tsplot(all_data_np, time = timepoints, condition = condition_names, legend = True, ax = s, color = color_dict)
    # plot_significance_lines(all_data_np, time_points = timepoints, offset=0.0125+stop_FIR_amplitude_range[0], slope=0.025, p_value_cutoff = 0.05, pal = [color_dict[cn] for cn in condition_names])
    s.set_ylim(stop_FIR_amplitude_range)
    s.set_xlim([interval[0]+1, interval[1]-1])
    s.set_xticks([0,5,10])    
    sn.despine(offset = 10, ax = s)
    pl.legend()

    ##############################################################################################################
    #
    # Plotting #2, covariates
    #
    ##############################################################################################################
    s = f.add_subplot(3,1,2)
    s.set_title(roi + ' corrs ' + second_plot_covariate)
    s.axhline(0, c='k', lw = 0.25)
    s.axvline(0, c='k', lw = 0.25)
    s.set_xlabel('Time [s]')
    s.set_ylabel('beta values')        

    for j, bn in enumerate(condition_names):
        # find out which color to use
        this_color = color_dict[bn]
        betas = beta_dict[bn][second_plot_covariate]
        pl.plot(timepoints, betas, c=this_color, label=bn)

        # take care of start and end of deconvolution interval.
        sig_periods = np.array(p_T_dict[bn][second_plot_covariate]) > sig
        if sig_periods[0] == True:
            sig_periods[0] = False
        if sig_periods[-1] == True:
            sig_periods[-1] = False

        nr_blocks = int(np.floor(np.abs(np.diff(sig_periods.astype(int))).sum() / 2.0))
        if nr_blocks > 0:
            print('# blocks found: %i'%nr_blocks)
            for b in range(nr_blocks):
                time_sig = np.arange(timepoints.shape[0])[np.r_[False, np.abs(np.diff(sig_periods)) > 0]][b*2:(b*2)+2]
                pl.plot([timepoints[time_sig[0]]-0.5, timepoints[time_sig[-1]] + 0.5], [0.0125 + stop_FIR_pe_range[0]+0.0125*j, 0.0125 + stop_FIR_pe_range[0]+0.0125*j], color = this_color, linewidth = 3.0, alpha = 0.8)
    s.set_ylim(stop_FIR_pe_range)
    s.set_xlim([interval[0]+1, interval[1]-1])
    s.set_xticks([0,5,10])    
    pl.legend()
    sn.despine(offset = 10, ax = s)    
    pl.tight_layout()

    ##############################################################################################################
    #
    # Plotting #3, split on SSRT
    #
    ##############################################################################################################
    sf_names = ['SSRT_short', 'SSRT_long']
    sf_color_dict = dict(zip(sf_names,['k','gray']))
    sf_rename_dict = dict(zip(sf_names,['fast','slow']))

    s = f.add_subplot(3,1,3)
    s.set_title('Corr ' + slow_fast_condition)
    peak_timepoint = timepoints == np.argmax(p_T_dict[slow_fast_condition][second_plot_covariate])
    ssrt_pd = pd.DataFrame(np.array([all_data_np[:,peak_timepoint,1].squeeze(), np.array(subj_data[second_plot_covariate], dtype = float)]).T, index=np.arange(all_data_np.shape[0]), columns=['BOLD',second_plot_covariate])
    sn.regplot("BOLD", second_plot_covariate, data=ssrt_pd, color='r', ax=s)
    # s.set_ylim([0,600])
    s.set_xlim([-0.3,0.25])
    sn.despine(offset = 10, ax = s)    
    pl.tight_layout()

    ssrt_data = np.array(subj_data[second_plot_covariate], dtype = float)
    pl.savefig(output_filename)
    for i, cn in enumerate(condition_names):
        peak_time = timepoints == np.argmax(p_T_dict[cn][second_plot_covariate])
        print(cn, timepoints[peak_time], scipy.stats.pearsonr(all_data_np[:,peak_time,i].squeeze(), ssrt_data))

    shell()
