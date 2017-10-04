from __future__ import division, print_function


def fit_FIR_roi_test(experiment,
                h5_file,
                in_files,
                vol_regressor_list, 
                behavior_file_list, 
                mapper_file = 'zstat2_flirt',
                mask_threshold = 2.0,
                mask_direction = 'pos',
                fmri_data_type = 'psc',
                fir_frequency = 4,
                fir_interval = [-3.0,12.0],
                roi_list = ['maxSTN25exc','SThR25max','SST_GO_preSMA','PvmPFCNoventri','PstriatumNoVentri'],
                event_conditions = ['ww', 'wl.u', 'wl.l', 'll'],
                TR = 2.0, 
                output_pdf_dir = '', 
                output_tsv_dir = '',
                which_signal_selection = 'projection'):

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
    from .behavior import process_test_tsv
    from ..utils.utils import roi_data_from_hdf
    from IPython import embed as shell

    run_durations = []
    for ifn in in_files:
        non_TR, dims, dyns, voxsize, affine = get_scaninfo(ifn)
        run_durations.append(TR*dyns)

    ################################################################################## 
    # behavior data generalizes across ROIs of course
    ##################################################################################
    all_event_df = process_test_tsv(behavior_file_list, run_durations)

    event_types_times, event_types_durs = {}, {}
    for evc in event_conditions:
        event_types_times.update({evc.replace('.','_'): np.array(all_event_df[(all_event_df['Cond'] == evc)]['onset'])})
        event_types_durs.update({evc.replace('.','_'): np.array(all_event_df[(all_event_df['Cond'] == evc)]['RT'])})
    new_event_conditions = [evc.replace('.','_') for evc in event_conditions]

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
        contrast_data = roi_data_from_hdf(data_types_wildcards = [roi], roi_name_wildcard = roi, hdf5_file = h5_file, folder_alias = 'rois')
        # time_course_data = [roi_data_from_hdf(data_types_wildcards = [os.path.split(in_f)[-1][:-7]], roi_name_wildcard = roi, hdf5_file = h5_file, folder_alias = fmri_data_type) for in_f in in_files]
        time_course_data = []
        for in_f in in_files:
            time_course_data.append(roi_data_from_hdf(data_types_wildcards = [os.path.split(in_f)[-1][:-7]], roi_name_wildcard = roi, hdf5_file = h5_file, folder_alias = fmri_data_type))

        time_course_data = np.hstack(time_course_data)

        over_mask_threshold = (contrast_data[:,0]>mask_threshold)
        iceberg_tip = contrast_data[over_mask_threshold, 0]

        if which_signal_selection == 'projection':
            projected_time_course = np.dot(time_course_data[over_mask_threshold].T, iceberg_tip) / np.sum(iceberg_tip)
            this_timecourse = projected_time_course
        elif which_signal_selection == 'hard':
            av_time_course = time_course_data[over_mask_threshold,:].mean(axis = 0)
            this_timecourse = av_time_course

        nuisance_regressors = np.nan_to_num(all_vol_regs)

        # shell()

        fd = FIRDeconvolution(
            signal = this_timecourse, 
            events = [event_types_times[evt] for evt in new_event_conditions], # dictate order
            event_names = new_event_conditions, 
            durations = event_types_durs, #{evt: evd[evt] for evt, evd in zip(event_conditions, event_types_durs)},
            sample_frequency = 1.0/TR,
            deconvolution_frequency = fir_frequency,
            deconvolution_interval = fir_interval
            )

        fd.resampled_signal = np.nan_to_num(fd.resampled_signal)
        # we then tell it to create its design matrix
        fd.create_design_matrix()

        # resample mocos and so forth
        all_nuisances = sp.signal.resample(nuisance_regressors, fd.resampled_signal_size, axis = 0)
        fd.add_continuous_regressors_to_design_matrix(all_nuisances.T)

        # fit
        fd.regress(method = 'lstsq')
        # fd.ridge_regress(cv = 10)
        fd.calculate_rsq()

        # plot
        sn.set_style('ticks')
        f = pl.figure(figsize = (6,3))
        s = f.add_subplot(111)
        s.set_title(roi)
        s.axhline(0, c='k', lw = 0.25)
        s.axvline(0, c='k', lw = 0.25)
        s.set_xlabel('Time [s]')
        s.set_ylabel('BOLD % signal change')
        for en in new_event_conditions:
            this_tc = np.squeeze(np.nan_to_num(fd.betas_for_cov(en).T))
            pl.plot(fd.deconvolution_interval_timepoints, this_tc, label = en)
        pl.legend()
        sn.despine(offset = 10, ax = s)
        pl.tight_layout()

        pl.savefig(op.join(output_pdf_dir, roi + '_deco_%s.pdf'%which_signal_selection))

        f = pl.figure(figsize = (9,3))
        s = f.add_subplot(111)
        s.set_title(roi)
        s.axhline(0, c='k', lw = 0.25)
        s.set_title('data and predictions, rsq %1.3f, roi: %s'%(fd.rsq, roi))
        s.set_xlabel('Time [s]')
        s.set_ylabel('BOLD % signal change')
        pl.plot(np.linspace(0,np.sum(run_durations), fd.resampled_signal.shape[1]), fd.resampled_signal.T, 'r', label = 'data')
        pl.plot(np.linspace(0,np.sum(run_durations), fd.resampled_signal.shape[1]), fd.predict_from_design_matrix(fd.design_matrix).T, 'k', label = 'model')
        pl.legend()
        sn.despine(offset = 10, ax = s)
        pl.tight_layout()
        pl.savefig(op.join(output_pdf_dir, roi + '_deco_tc_%s.pdf'%which_signal_selection))

        op_df = pd.DataFrame(np.array([np.squeeze(np.nan_to_num(fd.betas_for_cov(en).T)) for en in new_event_conditions]), 
                            columns = fd.deconvolution_interval_timepoints, 
                            index = new_event_conditions)
        op_df.to_csv(op.join(output_tsv_dir, roi + '_deco_test_%s.tsv'%which_signal_selection), sep = '\t')






