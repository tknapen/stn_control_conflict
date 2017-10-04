def create_pearl_pp_workflow(analysis_info, name='pearl'):
    import os.path as op
    import nipype.pipeline as pe
    import tempfile
    import glob
    from nipype.interfaces import fsl
    from nipype.interfaces.utility import Function, Merge, IdentityInterface
    from spynoza.nodes.utils import get_scaninfo, dyns_min_1, topup_scan_params, apply_scan_params
    from nipype.interfaces.io import SelectFiles, DataSink
 
    # Importing of custom nodes from spynoza packages; assumes that spynoza is installed:
    # pip install git+https://github.com/spinoza-centre/spynoza.git@develop
    from spynoza.nodes.filtering import savgol_filter
    from spynoza.nodes.utils import get_scaninfo, pickfirst, percent_signal_change, average_over_runs, pickle_to_json, set_nifti_intercept_slope, non_uniformity_correct_4D_file
    from spynoza.workflows.topup_unwarping import create_topup_workflow
    from spynoza.workflows.B0_unwarping import create_B0_workflow
    from spynoza.workflows.motion_correction import create_motion_correction_workflow
    from spynoza.workflows.registration import create_registration_workflow
    from spynoza.workflows.retroicor import create_retroicor_workflow
    from spynoza.workflows.sub_workflows.masks import create_masks_from_surface_workflow
    from spynoza.nodes.fit_nuisances import fit_nuisances

    from .motion_correction import create_motion_correction_workflow
    from ..surf.masks import avg_label_to_subject_label
    from ..utils.utils import convert_mapper_data_to_session, mask_nii_2_hdf5

    ########################################################################################
    # nodes
    ########################################################################################

    input_node = pe.Node(IdentityInterface(
                fields=['raw_directory', 
                    'output_directory', 
                    'FS_ID', 
                    'FS_subject_dir',
                    'sub_id', 
                    'which_file_is_EPI_space',
                    'standard_file', 
                    'psc_func', 
                    'tr',
                    'exp_shorthand',
                    'masks']), name='inputspec')
    # get standard masks from MNI space
    input_node.inputs.masks = glob.glob(op.join(analysis_info['MNI_mask_folder'], '*.nii.gz'))

    # i/o node
    datasource_templates = dict(func='{sub_id}/func/*{exp_shorthand}*_bold.nii.gz', behavior='{sub_id}/func/*{exp_shorthand}*_events.tsv')
    datasource = pe.Node(SelectFiles(datasource_templates, sort_filelist = True, raise_on_empty = False), 
        name = 'datasource')

    output_node = pe.Node(IdentityInterface(fields=([
            'temporal_filtered_files', 
            'percent_signal_change_files'])), name='outputspec')

    bet_epi = pe.MapNode(interface=
        fsl.BET(frac=analysis_info['bet_f_value'], vertical_gradient = analysis_info['bet_g_value'], 
                functional=True, mask = True), name='bet_epi', iterfield=['in_file'])

    # node for converting pickle files to json
    sgfilter = pe.MapNode(Function(input_names=['in_file'],
                                    output_names=['out_file'],
                                    function=savgol_filter),
                      name='sgfilter', iterfield=['in_file'])

    # node for percent signal change
    psc = pe.MapNode(Function(input_names=['in_file', 'func'],
                                    output_names=['out_file'],
                                    function=percent_signal_change),
                      name='percent_signal_change', iterfield=['in_file'])

    # node for conversion of mapper data to session
    mapper_convert = pe.Node(Function(input_names=['workflow_output_directory', 'sub_id', 'hires_2_session_reg', 'example_func', 'str_repl'],
                                    output_names=['out_files'],
                                    function=convert_mapper_data_to_session),
                      name='mapper_convert')
    if analysis_info['experiment'] == 'rl':
        mapper_convert.inputs.str_repl = ['/rl/', '/map/']
    elif analysis_info['experiment'] == 'stop':
        mapper_convert.inputs.str_repl = ['/stop/', '/map/']
    

    hdf5_psc_masker = pe.Node(Function(input_names = ['in_files', 'mask_files', 'hdf5_file', 'folder_alias'], output_names = ['hdf5_file'],
                                    function = mask_nii_2_hdf5), 
                                    name = 'hdf5_psc_masker')
    hdf5_psc_masker.inputs.folder_alias = 'psc'
    hdf5_psc_masker.inputs.hdf5_file = op.join(tempfile.mkdtemp(), 'roi.h5')

    hdf5_stats_masker = pe.Node(Function(input_names = ['in_files', 'mask_files', 'hdf5_file', 'folder_alias'], output_names = ['hdf5_file'],
                                    function = mask_nii_2_hdf5), 
                                    name = 'hdf5_stats_masker')
    hdf5_stats_masker.inputs.folder_alias = 'stats'

    hdf5_roi_masker = pe.Node(Function(input_names = ['in_files', 'mask_files', 'hdf5_file', 'folder_alias'], output_names = ['hdf5_file'],
                                    function = mask_nii_2_hdf5), 
                                    name = 'hdf5_roi_masker')
    hdf5_roi_masker.inputs.folder_alias = 'rois'

    vol_trans_node = pe.MapNode(interface=fsl.ApplyXfm(apply_xfm = True, interp = 'sinc', padding_size = 0), name='vol_trans', iterfield = ['in_file'])
    thresh_node = pe.MapNode(fsl.Threshold(thresh = analysis_info['MNI_mask_threshold'], args = '-bin', output_datatype = 'int'), name='thresh', iterfield = ['in_file'])

    merge_masks = pe.Node(Merge(len(analysis_info['label_folders']) + 2), name='merge_masks')

    datasink = pe.Node(DataSink(), name='sinker')
    datasink.inputs.parameterization = False

    ########################################################################################
    # workflow
    ########################################################################################

    # the actual top-level workflow
    pearl_pp_workflow = pe.Workflow(name=name)

    # data source 
    pearl_pp_workflow.connect(input_node, 'raw_directory', datasource, 'base_directory')
    pearl_pp_workflow.connect(input_node, 'sub_id', datasource, 'sub_id')
    pearl_pp_workflow.connect(input_node, 'exp_shorthand', datasource, 'exp_shorthand')
    # and data sink
    pearl_pp_workflow.connect(input_node, 'output_directory', datasink, 'base_directory')

    pearl_pp_workflow.connect(datasource, 'behavior', datasink, 'behavior')

    if analysis_info['perform_mri'] == 1:
        # still have to decide on BET for correction. 
        # point for intern, to decide on topup and B0 correction
        # BET
        pearl_pp_workflow.connect(datasource, 'func', bet_epi, 'in_file')

        # motion correction
        motion_proc = create_motion_correction_workflow('moco', method = analysis_info['moco_method'])
        pearl_pp_workflow.connect(input_node, 'tr', motion_proc, 'inputspec.tr')
        pearl_pp_workflow.connect(input_node, 'output_directory', motion_proc, 'inputspec.output_directory')
        pearl_pp_workflow.connect(input_node, 'which_file_is_EPI_space', motion_proc, 'inputspec.which_file_is_EPI_space')

        pearl_pp_workflow.connect(bet_epi, 'out_file', motion_proc, 'inputspec.in_files')

        # registration
        reg = create_registration_workflow(analysis_info, name = 'reg')
        pearl_pp_workflow.connect(input_node, 'output_directory', reg, 'inputspec.output_directory')
        pearl_pp_workflow.connect(motion_proc, 'outputspec.EPI_space_file', reg, 'inputspec.EPI_space_file')
        pearl_pp_workflow.connect(input_node, 'FS_ID', reg, 'inputspec.freesurfer_subject_ID')
        pearl_pp_workflow.connect(input_node, 'FS_subject_dir', reg, 'inputspec.freesurfer_subject_dir')
        pearl_pp_workflow.connect(input_node, 'standard_file', reg, 'inputspec.standard_file')

        # temporal filtering
        pearl_pp_workflow.connect(motion_proc, 'outputspec.motion_corrected_files', sgfilter, 'in_file')

        # node for percent signal change
        pearl_pp_workflow.connect(input_node, 'psc_func', psc, 'func')
        pearl_pp_workflow.connect(sgfilter, 'out_file', psc, 'in_file')

        # connect filtering and psc results to output node 
        pearl_pp_workflow.connect(sgfilter, 'out_file', output_node, 'temporal_filtered_files')
        pearl_pp_workflow.connect(psc, 'out_file', output_node, 'percent_signal_change_files')

        # gather MNI mask files
        pearl_pp_workflow.connect(input_node, 'masks', vol_trans_node, 'in_file')
        pearl_pp_workflow.connect(reg, 'rename_standard2example_func.out_file', vol_trans_node, 'in_matrix_file')
        pearl_pp_workflow.connect(reg, 'rename_example_func.out_file', vol_trans_node, 'reference')

        pearl_pp_workflow.connect(vol_trans_node, 'out_file', thresh_node, 'in_file' )
        pearl_pp_workflow.connect(thresh_node, 'out_file', datasink, 'masks.MNI')
        pearl_pp_workflow.connect(vol_trans_node, 'out_file', datasink, 'masks.MNI_nt')

    ########################################################################################
    # masking stuff if doing mri analysis
    ########################################################################################
        # loop across different folders to mask
        # untested as yet.
        masking_list = []
        dilate_list = []
        # lllist = []
        for opd, label_directory in zip(['dc'] + analysis_info['label_folders'], [''] + analysis_info['label_folders']):
            # if label_directory != '':
            #     lllist.append(pe.Node(Function(input_names=['avg_sj', 'trg_sj', 'label_directory'],
            #                     output_names=['output_files'],
            #                     function=avg_label_to_subject_label),
            #       name='lbl2lbl_'+label_directory.replace('.', '_')))
            #     lllist[-1].inputs.avg_sj = analysis_info['avg_subjects_fsid']
            #     lllist[-1].inputs.label_directory = label_directory
            #     pearl_pp_workflow.connect(input_node, 'FS_ID', lllist[-1], 'trg_sj')
            #     pearl_pp_workflow.connect(lllist[-1], 'output_files', datasink, 'labels.'+label_directory)

            dilate_list.append(
                pe.MapNode(interface=fsl.maths.DilateImage(
                    operation = 'mean', kernel_shape = 'sphere', kernel_size = analysis_info['dilate_kernel_size']), 
                    name='dilate_'+label_directory.replace('.', '_'), iterfield=['in_file'])) 
            
            masking_list.append(create_masks_from_surface_workflow(name = 'masks_from_surface_'+label_directory.replace('.', '_')))

            masking_list[-1].inputs.inputspec.label_directory = label_directory
            masking_list[-1].inputs.inputspec.fill_thresh = 0.005
            masking_list[-1].inputs.inputspec.re = '*.label'
           
            pearl_pp_workflow.connect(motion_proc, 'outputspec.EPI_space_file', masking_list[-1], 'inputspec.EPI_space_file')
            pearl_pp_workflow.connect(input_node, 'output_directory', masking_list[-1], 'inputspec.output_directory')
            pearl_pp_workflow.connect(input_node, 'FS_subject_dir', masking_list[-1], 'inputspec.freesurfer_subject_dir')
            pearl_pp_workflow.connect(input_node, 'FS_ID', masking_list[-1], 'inputspec.freesurfer_subject_ID')
            pearl_pp_workflow.connect(reg, 'rename_register.out_file', masking_list[-1], 'inputspec.reg_file')

            pearl_pp_workflow.connect(masking_list[-1], 'outputspec.masks', dilate_list[-1], 'in_file')
            pearl_pp_workflow.connect(dilate_list[-1], 'out_file', datasink, 'masks.'+opd)

        # import stats for mapper GLM across sessions
        if analysis_info['experiment'] in ['rl', 'stop']:
            # we assume the mapper's already run
            pearl_pp_workflow.connect(input_node, 'output_directory', mapper_convert, 'workflow_output_directory')
            pearl_pp_workflow.connect(input_node, 'sub_id', mapper_convert, 'sub_id')
            pearl_pp_workflow.connect(reg, 'outputspec.T1_EPI_matrix_file', mapper_convert, 'hires_2_session_reg')
            pearl_pp_workflow.connect(reg, 'outputspec.EPI_space_file', mapper_convert, 'example_func')
            pearl_pp_workflow.connect(mapper_convert, 'out_files', datasink, 'mapper_stats')

            # to H5 file
            pearl_pp_workflow.connect(psc, 'out_file', hdf5_psc_masker, 'in_files')

            for i in range(len(analysis_info['label_folders'])+1):
                pearl_pp_workflow.connect(dilate_list[i], 'out_file', merge_masks, 'in'+str(i+1))
            # also add rois from MNI for hdf5 transplant
            pearl_pp_workflow.connect(thresh_node, 'out_file', merge_masks, 'in'+str(i+2))
            pearl_pp_workflow.connect(merge_masks, 'out', hdf5_psc_masker, 'mask_files')

            # the hdf5_file is created by the psc node, and then passed from masker to masker on into the datasink.
            pearl_pp_workflow.connect(hdf5_psc_masker, 'hdf5_file', hdf5_stats_masker, 'hdf5_file')

            # then mask the mapper stats
            pearl_pp_workflow.connect(mapper_convert, 'out_files', hdf5_stats_masker, 'in_files')
            pearl_pp_workflow.connect(merge_masks, 'out', hdf5_stats_masker, 'mask_files')
            pearl_pp_workflow.connect(hdf5_stats_masker, 'hdf5_file', hdf5_roi_masker, 'hdf5_file')

            # and the probability maps of the different rois
            pearl_pp_workflow.connect(vol_trans_node, 'out_file', hdf5_roi_masker, 'in_files')
            pearl_pp_workflow.connect(merge_masks, 'out', hdf5_roi_masker, 'mask_files')
            pearl_pp_workflow.connect(hdf5_roi_masker, 'hdf5_file', datasink, 'h5')

    ########################################################################################
    # wrapping up, sending data to datasink 
    ########################################################################################

        pearl_pp_workflow.connect(bet_epi, 'out_file', datasink, 'bet.epi')
        pearl_pp_workflow.connect(bet_epi, 'mask_file', datasink, 'bet.epimask')

        pearl_pp_workflow.connect(sgfilter, 'out_file', datasink, 'tf')
        pearl_pp_workflow.connect(psc, 'out_file', datasink, 'psc')

    return pearl_pp_workflow
