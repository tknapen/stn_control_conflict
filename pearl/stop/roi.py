from __future__ import division, print_function

# def prepare_mcf_data(in_mcf_data):

#     np.hstack(np.vstack([np.zeros(in_mcf_data.shape[0]), np.diff(in_mcf_data)])

#     return 


def fit_FIR_roi(experiment,
                h5_file,
                in_files,
                vol_regressor_list, 
                behavior_file_list, 
                mapper_file = 'zstat2_flirt',
                mask_threshold = 0.0,
                mask_direction = 'pos',
                fmri_data_type = 'psc',
                fir_frequency = 4,
                fir_interval = [-3.0,12.0],
                roi_list = ['maxSTN25exc','SThR25max','SST_GO_preSMA','SST_GO_rIFG', 'Caudate', 'PvmPFCNoventri'],
                TR = 2.0, 
                output_pdf_dir = '', 
                output_tsv_dir = ''):

    import nibabel as nib
    import numpy as np
    import numpy.linalg as LA
    import scipy as sp
    import os
    import os.path as op
    import pandas as pd
    import matplotlib.pyplot as pl
    import seaborn as sn
    from spynoza.nodes.utils import get_scaninfo
    from fir import FIRDeconvolution
    import tempfile
    from .behavior import process_tsv
    from ..utils.utils import roi_data_from_hdf
    from IPython import embed as shell

    run_durations = []
    for ifn in in_files:
        non_TR, dims, dyns, voxsize, affine = get_scaninfo(ifn)
        run_durations.append(TR*dyns)

    ################################################################################## 
    # behavior data generalizes across ROIs of course
    ##################################################################################
    all_event_df = process_tsv(behavior_file_list, run_durations)

    # first, just run a FIR on the image pairs
    stim_event_names = ['correct', 'succesful_stop', 'Failed_stop']
    stim_event_list = []
    for en in stim_event_names:
        stim_event_list.append(np.array(all_event_df['onset'])[np.array(all_event_df['Response'] == en)])

    all_event_names = stim_event_names

    ################################################################################## 
    # whole-brain nuisance data generalizes across ROIs of course
    ##################################################################################

    if vol_regressor_list != []:
        all_vol_regs = []
        for x in range(len(vol_regressor_list)):
            all_vol_regs.append(np.loadtxt(vol_regressor_list[x]))
        all_vol_regs = np.vstack(all_vol_regs)

    ################################################################################## 
    # per-roi data
    ##################################################################################
    for roi in roi_list:
        # shell()
        contrast_data = roi_data_from_hdf(data_types_wildcards = [roi], roi_name_wildcard = roi, hdf5_file = h5_file, folder_alias = 'rois')
        time_course_data = [roi_data_from_hdf(data_types_wildcards = [os.path.split(in_f)[-1][:-7]], roi_name_wildcard = roi, hdf5_file = h5_file, folder_alias = fmri_data_type) for in_f in in_files]

        time_course_data = np.hstack(time_course_data)

        # if mask_threshold < 0:
        #     mask_threshold = -mask_threshold
        #     contrast_data = -contrast_data

        over_mask_threshold = (contrast_data[:,0]>mask_threshold)
        iceberg_tip = contrast_data[over_mask_threshold, 0]

        projected_time_course = np.dot(time_course_data[over_mask_threshold].T, iceberg_tip) / np.sum(iceberg_tip)
        av_time_course = time_course_data.mean(axis = 0)

        # nuisance_regressors = np.nan_to_num(all_vol_reg)
        fd = FIRDeconvolution(
            signal = projected_time_course, 
            events = [stim_event_list[0], stim_event_list[1], stim_event_list[2]], # dictate order
            event_names = stim_event_names, 
            # durations = {'AB':stim_duration_list[0], 'CD':stim_duration_list[1], 'EF':stim_duration_list[2], 'fb':fb_durations},
            # events = [stim_events], # dictate order
            # event_names = ['stim'], 
            # durations = {'stim':stim_durations},
            # covariates = covariates,
            sample_frequency = 1.0/TR,
            deconvolution_frequency = fir_frequency,
            deconvolution_interval = fir_interval
            )

        fd.resampled_signal = np.nan_to_num(fd.resampled_signal)
        # we then tell it to create its design matrix
        fd.create_design_matrix()

        # resample mocos and so forth
        # all_nuisances = sp.signal.resample(nuisance_regressors, fd.resampled_signal_size, axis = -1)
        # fd.add_continuous_regressors_to_design_matrix(all_nuisances)

        # fit
        fd.regress(method = 'lstsq')
        # fd.ridge_regress(cv = 10)
        fd.calculate_rsq()

        # plot
        sn.set_style('ticks')
        f = pl.figure(figsize = (6,3))
        s = f.add_subplot(111)
        s.axhline(0, c='k', lw = 0.25)
        s.axvline(0, c='k', lw = 0.25)
        s.set_xlabel('Time [s]')
        s.set_ylabel('BOLD % signal change')
        for en in all_event_names:
            this_tc = np.squeeze(np.nan_to_num(fd.betas_for_cov(en).T))
            pl.plot(fd.deconvolution_interval_timepoints, this_tc, label = en)
        pl.legend()
        sn.despine(offset = 10, ax = s)
        pl.tight_layout()

        pl.savefig(op.join(output_pdf_dir, roi + '_deco.pdf'))

        f = pl.figure(figsize = (9,3))
        s = f.add_subplot(111)
        s.axhline(0, c='k', lw = 0.25)
        s.set_title('data and predictions, rsq %1.3f'%fd.rsq)
        s.set_xlabel('Time [s]')
        s.set_ylabel('BOLD % signal change')
        pl.plot(np.linspace(0,np.sum(run_durations), fd.resampled_signal.shape[1]), fd.resampled_signal.T, 'r', label = 'data')
        pl.plot(np.linspace(0,np.sum(run_durations), fd.resampled_signal.shape[1]), fd.predict_from_design_matrix(fd.design_matrix).T, 'k', label = 'model')
        pl.legend()
        sn.despine(offset = 10, ax = s)
        pl.tight_layout()
        pl.savefig(op.join(output_pdf_dir, roi + '_deco_tc.pdf'))

        op_df = pd.DataFrame(np.array([np.squeeze(np.nan_to_num(fd.betas_for_cov(en).T)) for en in all_event_names]), 
                            columns = fd.deconvolution_interval_timepoints, 
                            index = all_event_names)
        # np.savetxt(op.join(output_tsv_dir, roi + '_deco.tsv'), np.array(op_df), delimiter = '\t')
        op_df.to_csv(op.join(output_tsv_dir, roi + '_deco_stop.tsv'), sep = '\t')


