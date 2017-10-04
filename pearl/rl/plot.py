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
                print("Cluster p-value %1.5f between timepoints %2.2f - %2.2f, in color %s"%(p_val, time_points[s[0]], time_points[s[-1]], str(pal[i])))


def plot_deco_results_test(all_deco_files, subj_data, event_conditions, roi, event_conditions_for_covariates, sj_covariates = ['Beta', 'SSRT'], interval = [-3,15], output_filename = '', 
    rl_test_FIR_amplitude_range = [0,0], rl_test_FIR_pe_range = [0,0]):
    import matplotlib.pyplot as pl
    import seaborn as sn
    import pandas as pd
    import numpy as np
    from IPython import embed as shell
    import statsmodels.api as sm

    stats_threshold = 0.0125
    all_data = np.array([np.loadtxt(df) for df in all_deco_files])
    timepoints = np.linspace(interval[0],interval[1],all_data.shape[-1])
    new_event_conditions = [evc.replace('.','_') for evc in event_conditions]

    covariate_event_indices = [event_conditions.index(ecc) for ecc in event_conditions_for_covariates]
    cond_diffs = ['ll-ww', 'll-wl_u', 'ww-wl_u']
    all_event_names = new_event_conditions
    # shell()
    sets = [0,1,2,3]
    colors = np.array(['g','orange','orange','r'])

    color_dict = dict(zip(all_event_names, colors))


    roi_data_diffs = np.squeeze(np.array([  all_data[:,new_event_conditions=='ll'] - all_data[:,new_event_conditions=='ww'], 
                            all_data[:,new_event_conditions=='ll'] - all_data[:,new_event_conditions=='wl_u'], 
                            all_data[:,new_event_conditions=='ww'] - all_data[:,new_event_conditions=='wl_u'] ]
                            )).transpose(1,0,2)

    sig = -np.log10(stats_threshold)
    # if roi_name == 'Caudate':
    #     shell()
    sn.set_style('ticks')
    f = pl.figure(figsize = (5,len(sj_covariates)*4+3))
    s = f.add_subplot(len(sj_covariates)+1,1,1)
    s.set_title(roi + ' gain')
    s.axhline(0, c='k', lw = 0.25)
    s.axvline(0, c='k', lw = 0.25)
    s.set_xlabel('Time [s]')
    s.set_ylabel('BOLD % signal change')

    min_d = all_data.transpose((0,2,1))[:,:,[0,1,2]].mean(axis = 0).min()
    for x in range(len(sets)):
        sn.tsplot(all_data.transpose((0,2,1))[:,:,sets[x]], time = timepoints, condition = all_event_names[sets[x]], legend = True, ax = s, color = color_dict)
        plot_significance_lines(all_data.transpose((0,2,1))[:,:,[sets[x]]], time_points = timepoints, offset=0.0125+rl_test_FIR_amplitude_range[0], slope=0.025, p_value_cutoff = 0.05, pal = [colors[sets[x]]])
    s.set_ylim(rl_test_FIR_amplitude_range)
    s.set_xlim([timepoints[0], timepoints[-1]])
    sn.despine(offset = 10, ax = s)

    pl.legend()

    ##############################################################################################################
    #
    # Now, we compute the correlations with Beta
    #
    ##############################################################################################################

    # construct across subjects covariate design matrix
    X = [np.ones(len(all_deco_files))]
    for cov_name in sj_covariates:
        this_cov = np.array(subj_data[cov_name], dtype = float)
        # z-score
        this_cov = (this_cov - this_cov.mean()) / this_cov.std()
        X.append(this_cov)

    X = np.array(X).T

    p_T_vals = np.zeros((len(all_event_names),all_data.shape[-1], X.shape[1]+1))
    tcs = np.zeros((len(all_event_names),all_data.shape[-1], X.shape[1]))
    for i, et in enumerate(all_event_names):
        for x in range(all_data.shape[-1]):
            model = sm.OLS(np.squeeze(all_data[:,i,x]),X)
            results = model.fit()
            p_T_vals[i,x,:X.shape[1]] = -np.log10(results.pvalues)
            p_T_vals[i,x,-1] = -np.log10(results.f_pvalue)
            tcs[i,x] = results.params

    for i, c in enumerate(sj_covariates): # , 'Beta'
        s = f.add_subplot(len(sj_covariates)+1,1,2+i)
        s.set_title(roi + ' corrs %s'%c)
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
                    pl.plot([timepoints[time_sig[0]]-0.5, timepoints[time_sig[-1]] + 0.5], [0.0125 + rl_test_FIR_pe_range[0]+0.0125*j, 0.0125 + rl_test_FIR_pe_range[0]+0.0125*j], color = colors[j], linewidth = 3.0, alpha = 0.8)
        s.set_ylim(rl_test_FIR_pe_range)
        s.set_xlim([timepoints[0], timepoints[-1]])
        pl.legend()
        sn.despine(offset = 10, ax = s)

    pl.tight_layout()
    # pl.show()

    pl.savefig(output_filename)

def plot_deco_results_test_covariates_per_event_type(all_deco_files, 
                subj_data, event_conditions, roi, event_conditions_for_covariates, 
                sj_covariates = {'ll':['Beta', 'SSRT'],'ww':['Beta', 'SSRT'],'wl_u':['Beta', 'SSRT'],'wl_l':['Beta', 'SSRT']}, interval = [-3,15], output_filename = '', 
    rl_test_FIR_amplitude_range = [0,0], rl_test_FIR_pe_range = [0,0]):
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

    color_dict = dict(zip(['ww','wl_u','ll'], ['g','k','r']))
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
    s.set_ylim(rl_test_FIR_amplitude_range)
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
                        pl.plot([timepoints[time_sig[0]]-0.5, timepoints[time_sig[-1]] + 0.5], [0.0125 + rl_test_FIR_pe_range[0]+0.0125*j, 0.0125 + rl_test_FIR_pe_range[0]+0.0125*j], color = this_color, linewidth = 3.0, alpha = 0.8)
        s.set_ylim(rl_test_FIR_pe_range)
        s.set_xlim([timepoints[0], timepoints[-1]])
        pl.legend()
        sn.despine(offset = 10, ax = s)

    pl.tight_layout()

    pl.savefig(output_filename)

def extract_timepoints_results_test(all_deco_files, subj_data, event_conditions, roi, event_conditions_for_covariates, sj_covariates = ['Beta', 'SSRT'], interval = [-2,10]):
    import matplotlib.pyplot as pl
    import seaborn as sn
    import pandas as pd
    import numpy as np
    from IPython import embed as shell
    import statsmodels.api as sm
    import scipy.stats as st

    stats_threshold = 0.05
    sig = -np.log10(stats_threshold)

    all_data = np.array([np.loadtxt(df) for df in all_deco_files])
    timepoints = np.linspace(interval[0],interval[1],all_data.shape[-1])
    new_event_conditions = [evc.replace('.','_') for evc in event_conditions]
    all_event_names = new_event_conditions

    # construct across subjects covariate design matrix
    X = [np.ones(len(all_deco_files))]
    for cov_name in sj_covariates:
        this_cov = np.array(subj_data[cov_name], dtype = float)
        # z-score
        this_cov = (this_cov - this_cov.mean()) / this_cov.std()
        X.append(this_cov)

    X = np.array(X).T

    p_T_vals = np.zeros((len(all_event_names),all_data.shape[-1], X.shape[1]+1))
    tcs = np.zeros((len(all_event_names),all_data.shape[-1], X.shape[1]))
    for i, et in enumerate(all_event_names):
        for x in range(all_data.shape[-1]):
            model = sm.OLS(np.squeeze(all_data[:,i,x]),X)
            results = model.fit()
            p_T_vals[i,x,:X.shape[1]] = -np.log10(results.pvalues)
            p_T_vals[i,x,-1] = -np.log10(results.f_pvalue)
            tcs[i,x] = results.params

    append_df = {}
    for i, c in enumerate(sj_covariates): # , 'Beta'
        for j, en in enumerate(all_event_names):
            data = tcs[j,:,i+1]
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
                    print(c, str(b), en, timepoints[time_sig], np.argmax(p_T_vals[j,:,i+1]))
                    print(all_data[:,j,np.argmax(p_T_vals[j,:,i+1])])
                    print(10**-p_T_vals[j,:,i+1].max())
                    append_df.update({c + '_' + en: pd.Series(all_data[:,j,np.argmax(p_T_vals[j,:,i+1])], index = subj_data.index)})
    append_df = pd.DataFrame(append_df)
    all_stuff = subj_data.append(append_df)

    st.pearsonr(np.array(subj_data['ac_ww'], dtype=float), np.array(append_df['Beta_ww'], dtype=float))
    st.pearsonr(np.array(subj_data['medRT_ll'], dtype=float), np.array(append_df['SSRT_ll'], dtype=float))
    st.pearsonr(np.array(subj_data['medRT_ww'], dtype=float), np.array(append_df['SSRT_ww'], dtype=float))
    
    # shell()


def plot_deco_results_test_MR():

    rd_diff = np.squeeze(roi_data[:,:,conditions == 'll']) - np.squeeze(roi_data[:,:,conditions == 'ww'])
    betas, res, rank, sem = np.linalg.lstsq(X, rd_diff)
    bbruns = np.ones((nr_bs, X.shape[-1], rd_diff.shape[-1]))
    for bbr in range(nr_bs):
        inds = np.random.randint(0, X.shape[0], X.shape[0])
        bbs, _res, _rank, _sem = np.linalg.lstsq(X[inds], rd_diff[inds])
        bbruns[bbr] = bbs
    bs_sd = bbruns.std(axis = 0)     

    sig = -np.log10(0.05/8)
    p_T_vals = np.zeros((2, rd_diff.shape[-1], 4))
    for x in range(rd_diff.shape[-1]):
        model = sm.OLS(np.squeeze(rd_diff[:,x]),X)
        results = model.fit()
        p_T_vals[0,x,:3] = -np.log10(results.pvalues)
        p_T_vals[1,x,:3] = results.tvalues

        p_T_vals[0,x,3] = -np.log10(results.f_pvalue)
        p_T_vals[1,x,3] = results.fvalue



    f = pl.figure(figsize = (8,8))
    f.suptitle('STN')

    s = f.add_subplot(111)
    # rd = np.squeeze(roi_data[:,:,conditions == 'll'])
    # betas, res, rank, sem = np.linalg.lstsq(X, rd)

    bbruns = np.ones((nr_bs, X.shape[-1], rd_diff.shape[-1]))
    for bbr in range(nr_bs):
        inds = np.random.randint(0, X.shape[0], X.shape[0])
        bbs, _res, _rank, _sem = np.linalg.lstsq(X[inds], rd_diff[inds])
        bbruns[bbr] = bbs
    bs_sd = bbruns.std(axis = 0)        
    for x in range(betas.shape[0]):
        pl.plot(times, betas[x], colors[x], label = beta_names[x])
        pl.fill_between(times, betas[x] - bs_sd[x], betas[x] + bs_sd[x], color = colors[x], alpha = 0.2)

    # significance
    sig = -np.log10(0.05)
    p_T_vals = np.zeros((2, rd_diff.shape[-1], 4))
    for x in range(rd_diff.shape[-1]):
        model = sm.OLS(np.squeeze(rd_diff[:,x]),X)
        results = model.fit()
        p_T_vals[0,x,:3] = -np.log10(results.pvalues)
        p_T_vals[1,x,:3] = results.tvalues

        p_T_vals[0,x,3] = -np.log10(results.f_pvalue)
        p_T_vals[1,x,3] = results.fvalue

    # shell()
    sig_periods = p_T_vals[0,:,[0,1]] > sig
    for i in [0,1]:
        time_sig = np.arange(times.shape[0])[sig_periods[i]]
        pl.plot([times[time_sig[0]]-0.5, times[time_sig[-1]] + 0.5], [[-0.025, -0.025], [-0.035, -0.035]][i], color = ['red','green'][i], linewidth = 3.0, alpha = 0.8)


    s.set_title('Lose-Lose')
    s.axhline(0, color = 'k', alpha = 0.5, lw = 0.5)
    # s.axvline(0, color = 'k', alpha = 0.5, lw = 0.5)
    pl.legend()
    s.set_ylabel('Beta values')
    # s.set_xlabel('time [s]')
    sn.despine(offset=10)
    s.set_xlim(time_period)
    s.set_xticks([0,5,10])
    s.set_ylim([-0.04,0.06001])

    # shell()

    s = f.add_subplot(2,2,2)
    rd = np.squeeze(roi_data[:,:,conditions == 'ww'])
    betas, res, rank, sem = np.linalg.lstsq(X, rd)

    bbruns = np.ones((nr_bs, X.shape[-1], rd.shape[-1]))
    for bbr in range(nr_bs):
        inds = np.random.randint(0, X.shape[0], X.shape[0])
        bbs, _res, _rank, _sem = np.linalg.lstsq(X[inds], rd[inds])
        bbruns[bbr] = bbs
    bs_sd = bbruns.std(axis = 0)    
    for x in range(betas.shape[0]):
        pl.plot(times, betas[x], colors[x], label = beta_names[x])
        pl.fill_between(times, betas[x] - bs_sd[x], betas[x] + bs_sd[x], color = colors[x], alpha = 0.2)

    # # significance is never reached
    # sig = -np.log10(0.05/8)
    # p_T_vals = np.zeros((2, rd.shape[-1], 4))
    # for x in range(rd.shape[-1]):
    #     model = sm.OLS(np.squeeze(rd[:,x]),X)
    #     results = model.fit()
    #     p_T_vals[0,x,:3] = -np.log10(results.pvalues)
    #     p_T_vals[1,x,:3] = results.tvalues

    #     p_T_vals[0,x,3] = -np.log10(results.f_pvalue)
    #     p_T_vals[1,x,3] = results.fvalue
    # # shell()
    # sig_periods_ww = p_T_vals[0,:,[0,1]] > sig
    # for i in [0,1]:
    #     time_sig = np.arange(times.shape[0])[sig_periods_ww[i]]
    #     pl.plot([times[time_sig[0]]-0.5, times[time_sig[-1]] + 0.5], [[-0.025, -0.025], [-0.035, -0.035]][i], color = ['red','green'][i], linewidth = 3.0, alpha = 0.8)


    s.set_title('Win-Win')
    s.axhline(0, color = 'k', alpha = 0.5, lw = 0.5)
    # s.axvline(0, color = 'k', alpha = 0.5, lw = 0.5)
    pl.legend()
    s.set_ylabel('Beta values')
    s.set_xlabel('time [s]')
    sn.despine(offset=10)
    s.set_xlim(time_period)
    s.set_xticks([0,5,10])
    s.set_ylim([-0.04,0.06001])


    # group split
    evt = 'll'
    s = f.add_subplot(2,2,3)

    group_values = np.array([np.array(ssa.evts['SSRT'])[0] for ssa in self.ssas])
    group_median = np.median(group_values)

    group = group_values <= group_median


    idx = conditions == evt
    this_condition_data = (roi_data[group,:,idx],roi_data[-group,:,idx])

    s.set_title('Lose-Lose - SSRT')
    sn.tsplot(this_condition_data[0], time = times, condition = ['SSRT fast'], ci = 68, color = 'gray', ls = '--')
    sn.tsplot(this_condition_data[1], time = times, condition = ['SSRT slow'], ci = 68, color = 'gray')

    i = 0
    time_sig = np.arange(times.shape[0])[sig_periods[i]]
    pl.plot([times[time_sig[0]]-0.5, times[time_sig[-1]] + 0.5], [[-0.025, -0.025], [-0.035, -0.035]][i], color = ['red','green'][i], linewidth = 3.0, alpha = 0.8)

    s.axhline(0, color = 'k', alpha = 0.5, lw = 0.5)
    # s.axvline(0, color = 'k', alpha = 0.5, lw = 0.5)
    s.set_xlabel('time [s]')
    s.set_ylabel('Z-scored BOLD')
    sn.despine(offset=10)
    s.set_xlim(time_period)
    s.set_xticks([0,5,10])
    s.set_ylim([-0.05,0.09])


    s = f.add_subplot(2,2,4)

    group_values = np.array([np.array(ssa.evts['Beta'])[0] for ssa in self.ssas])
    group_median = np.median(group_values)

    group = group_values <= group_median


    idx = conditions == evt
    this_condition_data = (roi_data[group,:,idx],roi_data[-group,:,idx])

    s.set_title('Lose-Lose - Beta')
    sn.tsplot(this_condition_data[0], time = times, condition = ['Explore'], ci = 68, color = 'gray', ls = '--')
    sn.tsplot(this_condition_data[1], time = times, condition = ['Exploit'], ci = 68, color = 'gray')

    i = 1
    time_sig = np.arange(times.shape[0])[sig_periods[i]]
    pl.plot([times[time_sig[0]]-0.5, times[time_sig[-1]] + 0.5], [[-0.025, -0.025], [-0.025, -0.025]][i], color = ['red','green'][i], linewidth = 3.0, alpha = 0.8)

    s.axhline(0, color = 'k', alpha = 0.5, lw = 0.5)
    # s.axvline(0, color = 'k', alpha = 0.5, lw = 0.5)
    s.set_xlabel('time [s]')
    s.set_ylabel('Z-scored BOLD')
    sn.despine(offset=10)
    s.set_xlim(time_period)
    s.set_xticks([0,5,10])
    s.set_ylim([-0.05,0.09])

    pl.tight_layout()


def plot_deco_results_for_publication(all_deco_files, subj_data, 
            roi, interval = [-2,10], output_filename = '', 
            sj_covariates = {},
            rl_test_FIR_amplitude_range = [0,0], rl_test_FIR_pe_range = [0,0], 
            second_plot_covariate = 'SSRT', slow_fast_condition = 'll'):
    import matplotlib.pyplot as pl
    import seaborn as sn
    import pandas as pd
    import numpy as np
    from IPython import embed as shell
    import statsmodels.api as sm
    import scipy.stats

    stats_threshold = 0.0125
    # all_data = np.array([np.loadtxt(df) for df in all_deco_files])
    all_data = [pd.read_csv(df, sep = '\t', index_col=0, header=0).T for df in all_deco_files]
    timepoints = np.array(all_data[0].index, dtype = float)

    condition_names = ['ww','wl_u','ll']

    all_data_np = np.array([np.array(ad[condition_names]) for ad in all_data]) # -> don't leave to chance but copied from roi.py where the deconv happens

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

    color_dict = dict(zip(condition_names, ['g','k','r']))
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
    s.set_ylim(np.array(rl_test_FIR_amplitude_range)) # fix the scaling
    s.set_xlim([interval[0]+1, interval[1]-1])
    s.set_xticks([0,5,10])    
    pl.legend()
    sn.despine(offset = 10, ax = s)

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
                pl.plot([timepoints[time_sig[0]]-0.5, timepoints[time_sig[-1]] + 0.5], [0.0125 + rl_test_FIR_pe_range[0]+0.0125*j, 0.0125 + rl_test_FIR_pe_range[0]+0.0125*j], color = this_color, linewidth = 3.0, alpha = 0.8)
    s.set_ylim(rl_test_FIR_pe_range)
    s.set_xlim([interval[0]+1, interval[1]-1])
    s.set_xticks([0,5,10])    
    pl.legend()
    sn.despine(offset = 10, ax = s)    

    ##############################################################################################################
    #
    # Plotting #3, split on SSRT
    #
    ##############################################################################################################
    
    sf_color_dict = dict(zip(['SSRT_short', 'SSRT_long'],['k','gray']))
    sf_rename_dict = dict(zip(['SSRT_short', 'SSRT_long'],['fast','slow']))

    s = f.add_subplot(3,1,3)
    s.set_title('Corr ' + slow_fast_condition)

    bold = all_data_np[:,timepoints == np.argmax(p_T_dict[slow_fast_condition][second_plot_covariate]),condition_names.index(slow_fast_condition)].squeeze()
    ssrt_data = np.array(subj_data[second_plot_covariate], dtype = float)

    ssrt_pd = pd.DataFrame(np.array([bold,ssrt_data]).T, index=np.arange(all_data_np.shape[0]), columns=['BOLD',second_plot_covariate])
    sn.regplot("BOLD", second_plot_covariate, data=ssrt_pd, color='r', ax=s)
    # s.set_ylim([0,600])
    s.set_xlim([-0.3,0.25])

    sn.despine(offset = 10, ax = s)    
    pl.tight_layout()

    pl.savefig(output_filename)
    # shell()
    ssrt_data = np.array(subj_data[second_plot_covariate], dtype = float)

    for i, cn in enumerate(condition_names):
        peak_time = timepoints == np.argmax(p_T_dict[cn][second_plot_covariate])
        print(cn, timepoints[peak_time], scipy.stats.pearsonr(all_data_np[:,peak_time,i].squeeze(), ssrt_data))


    # shell()



