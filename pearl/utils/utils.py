from __future__ import division, print_function

def mask_nii_2_hdf5(in_files, mask_files, hdf5_file, folder_alias):
    """masks data in in_files with masks in mask_files,
    to be stored in an hdf5 file

    Takes a list of 3D or 4D fMRI nifti-files and masks the
    data with all masks in the list of nifti-files mask_files.
    These files are assumed to represent the same space, i.e.
    that of the functional acquisitions. 
    These are saved in hdf5_file, in the folder folder_alias.

    Parameters
    ----------
    in_files : list
        list of absolute path to functional nifti-files.
        all nifti files are assumed to have the same ndim
    mask_files : list
        list of absolute path to mask nifti-files.
        mask_files are assumed to be 3D
    hdf5_file : str
    	absolute path to hdf5 file.
   	folder_alias : str
   		name of the to-be-created folder in the hdf5 file.

    Returns
    -------
    hdf5_file : str
        absolute path to hdf5 file.
    """

    import nibabel as nib
    import os.path as op
    import numpy as np
    import tables

    success = True

    mask_data = [np.array(nib.load(mf).get_data(), dtype = bool) for mf in mask_files]
    nifti_data = [nib.load(nf).get_data() for nf in in_files]

    mask_names = [op.split(mf)[-1].split('_vol.nii.gz')[0] for mf in mask_files]
    nifti_names = [op.split(nf)[-1].split('.nii.gz')[0] for nf in in_files]

    h5file = tables.open_file(hdf5_file, mode = "a", title = hdf5_file)
    # get or make group for alias folder
    try:
        folder_alias_run_group = h5file.get_node("/", name = folder_alias, classname='Group')
    except tables.NoSuchNodeError:
        print('Adding group ' + folder_alias + ' to this file')
        folder_alias_run_group = h5file.create_group("/", folder_alias, folder_alias)

    for (roi, roi_name) in zip(mask_data, mask_names):
        # get or make group for alias/roi
        try:
            run_group = h5file.get_node(where = "/" + folder_alias, name = roi_name, classname='Group')
        except tables.NoSuchNodeError:
            print('Adding group ' + folder_alias + '_' + roi_name + ' to this file')
            run_group = h5file.create_group("/" + folder_alias, roi_name, folder_alias + '_' + roi_name)

        h5file.create_array(run_group, roi_name, roi, roi_name + ' mask file for reconstituting nii data from masked data')

        for (nii_d, nii_name) in zip(nifti_data, nifti_names):
            print('roi: %s, nifti: %s'%(roi_name, nii_name))
            n_dims = len(nii_d.shape)
            if n_dims == 3:
                these_roi_data = nii_d[roi]
            elif n_dims == 4:   # timeseries data, last dimension is time.
                these_roi_data = nii_d[roi,:]
            else:
                print("n_dims in data {nifti} do not fit with mask".format(nii_name))
                success = False

            h5file.create_array(run_group, nii_name, these_roi_data, roi_name + ' data from ' + nii_name)

    h5file.close()

    return hdf5_file

def roi_data_from_hdf(data_types_wildcards, roi_name_wildcard, hdf5_file, folder_alias):
    """takes data_type data from masks stored in hdf5_file

    Takes a list of 4D fMRI nifti-files and masks the
    data with all masks in the list of nifti-files mask_files.
    These files are assumed to represent the same space, i.e.
    that of the functional acquisitions. 
    These are saved in hdf5_file, in the folder folder_alias.

    Parameters
    ----------
    data_types_wildcards : list
        list of data types to be loaded.
        correspond to nifti_names in mask_2_hdf5
    roi_name_wildcard : str
        wildcard for masks. 
        corresponds to mask_name in mask_2_hdf5.
    hdf5_file : str
        absolute path to hdf5 file.
    folder_alias : str
        name of the folder in the hdf5 file from which data
        should be loaded.

    Returns
    -------
    output_data : list
        list of numpy arrays corresponding to data_types and roi_name_wildcards
    """
    import tables
    import itertools
    import fnmatch
    import numpy as np
    from IPython import embed as shell

    h5file = tables.open_file(hdf5_file, mode = "r")

    try:
        folder_alias_run_group = h5file.get_node(where = '/', name = folder_alias, classname='Group')
    except tables.NoSuchNodeError:
        # import actual data
        print('No group ' + folder_alias + ' in this file')
        # return None


    all_roi_names = h5file.list_nodes(where = '/' + folder_alias, classname = 'Group')
    roi_names = [rn._v_name for rn in all_roi_names if roi_name_wildcard in rn._v_name]
    if len(roi_names) == 0:
        print('No rois corresponding to ' + roi_name_wildcard + ' in group ' + folder_alias)
        # return None
    
    data_arrays = []
    for roi_name in roi_names:
        try:
            roi_node = h5file.get_node(where = '/' + folder_alias, name = roi_name, classname='Group')
        except tables.NoSuchNodeError:
            print('No data corresponding to ' + roi_name + ' in group ' + folder_alias)
            pass
        all_data_array_names = h5file.list_nodes(where = '/' + folder_alias + '/' + roi_name)
        data_array_names = [adan._v_name for adan in all_data_array_names]
        selected_data_array_names = list(itertools.chain(*[fnmatch.filter(data_array_names, dtwc) for dtwc in data_types_wildcards]))
        
        # if sort_data_types:
        selected_data_array_names = sorted(selected_data_array_names)
        if len(data_array_names) == 0:
            print('No data corresponding to ' + str(selected_data_array_names) + ' in group /' + folder_alias + '/' + roi_name)
            pass
        else:
            print('Taking data corresponding to ' + str(selected_data_array_names) + ' from group /' + folder_alias + '/' + roi_name)
            data_arrays.append([])
            for dan in selected_data_array_names:
                data_arrays[-1].append(eval('roi_node.__getattr__("' + dan + '").read()'))
            print('Taken data corresponding to ' + str(selected_data_array_names) + ' from group /' + folder_alias + '/' + roi_name)
            data_arrays[-1] = np.hstack(data_arrays[-1]) # stack across timepoints or other values per voxel
            if len(data_arrays[-1].shape) == 1:
                data_arrays[-1] = data_arrays[-1][:,np.newaxis]
    all_roi_data_np = np.vstack(data_arrays)    # stack across regions to create a single array of voxels by values (i.e. timepoints)

    h5file.close()

    return all_roi_data_np

def convert_mapper_data_to_session(workflow_output_directory, sub_id, hires_2_session_reg, example_func, str_repl = ['/rl/', '/map/'], stat_re = 'tf.feat/stats/*stat'):
    import os.path as op
    import glob
    import nipype.pipeline as pe
    from nipype.interfaces import fsl
    from nipype.interfaces import freesurfer
    from nipype.interfaces.utility import Function, IdentityInterface
    import nipype.interfaces.io as nio
    from IPython import embed as shell

    input_folder = workflow_output_directory.replace(str_repl[0], str_repl[1])
    input_files = glob.glob(op.join(input_folder, stat_re + '*.nii.gz'))
    input_files.append(op.join(input_folder, 'reg', 'example_func.nii.gz'))

    ### NODES
    input_node = pe.Node(IdentityInterface(
        fields=['input_files', 
        'output_folder', 
        'mapper_2_hires_reg', 
        'hires_2_session_reg',
        'template_file']), name='inputspec')

    output_node = pe.Node(IdentityInterface(
        fields=['output_files']), name='outputspec')

    input_node.inputs.input_files = input_files
    input_node.inputs.output_folder = workflow_output_directory # op.join(workflow_output_directory, 'mapper_stat')
    input_node.inputs.mapper_2_hires_reg = op.join(input_folder, 'reg', 'example_func2highres.mat')
    input_node.inputs.hires_2_session_reg = hires_2_session_reg # op.join(workflow_output_directory, 'reg', 'highres2example_func.mat')
    input_node.inputs.template_file = example_func # op.join(workflow_output_directory, 'reg', 'example_func.nii.gz')

    concat_N = pe.Node(fsl.ConvertXFM(concat_xfm = True), name = 'concat_Mapper')
    vol_trans_node = pe.MapNode(interface=fsl.ApplyXfm(apply_xfm = True, interp = 'sinc', padding_size = 0), name='vol_trans', iterfield = ['in_file'])
    
    datasink = pe.Node(nio.DataSink(), name='sinker')
    datasink.inputs.parameterization = False

    ### WORKFLOW
    convert_mapper_data_to_session_workflow = pe.Workflow(name='mapper2session')

    convert_mapper_data_to_session_workflow.connect(input_node, 'mapper_2_hires_reg', concat_N, 'in_file')
    convert_mapper_data_to_session_workflow.connect(input_node, 'hires_2_session_reg', concat_N, 'in_file2')

    convert_mapper_data_to_session_workflow.connect(concat_N, 'out_file', vol_trans_node, 'in_matrix_file')
    convert_mapper_data_to_session_workflow.connect(input_node, 'input_files', vol_trans_node, 'in_file')
    convert_mapper_data_to_session_workflow.connect(input_node, 'template_file', vol_trans_node, 'reference')

    convert_mapper_data_to_session_workflow.connect(input_node, 'output_folder', datasink, 'base_directory')
    convert_mapper_data_to_session_workflow.connect(vol_trans_node, 'out_file', datasink, 'mapper_stat')
    
    convert_mapper_data_to_session_workflow.connect(concat_N, 'out_file', datasink, 'mapper_stat.mat')
    convert_mapper_data_to_session_workflow.connect(vol_trans_node, 'out_file', output_node, 'output_files')

    convert_mapper_data_to_session_workflow.run('MultiProc', plugin_args={'n_procs': 24})

    out_files = glob.glob(op.join(workflow_output_directory, 'mapper_stat', '*.nii.gz'))

    return out_files

def natural_sort(l):
    import re
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

def import_MNI_masks(mask_folder, example_func, standard2example_func, output_folder):
    """Import MNI space masks to example_func space.
    """

    import os.path as op
    import glob
    import nipype.pipeline as pe
    from nipype.interfaces import fsl
    from nipype.interfaces.utility import Function, IdentityInterface
    import nipype.interfaces.io as nio
    from IPython import embed as shell

    ### NODES
    input_node = pe.Node(IdentityInterface(
        fields=['masks', 
        'output_folder', 
        'standard2example_func', 
        'template_file']), name='inputspec')

    output_node = pe.Node(IdentityInterface(
        fields=['output_files']), name='outputspec')

    input_node.inputs.masks = glob.glob(op.join(mask_folder, '*.nii.gz'))
    input_node.inputs.output_folder = output_folder # op.join(workflow_output_directory, 'mapper_stat')
    input_node.inputs.standard2example_func = standard2example_func
    input_node.inputs.template_file = example_func # op.join(workflow_output_directory, 'reg', 'example_func.nii.gz')

    datasink = pe.Node(nio.DataSink(), name='sinker')
    datasink.inputs.parameterization = False

    thresh_node = pe.MapNode(fsl.Threshold(thresh = 0.001, args = '-bin'), name='thresh', iterfield = ['in_file'])
    vol_trans_node = pe.MapNode(interface=fsl.ApplyXfm(apply_xfm = True, interp = 'sinc', padding_size = 0, datatype = 'int'), name='vol_trans', iterfield = ['in_file'])

    ### WORKFLOW
    import_MNI_masks_workflow = pe.Workflow(name='import_MNI_masks')

    import_MNI_masks_workflow.connect(input_node, 'masks', thresh_node, 'in_file')

    import_MNI_masks_workflow.connect(thresh_node, 'out_file', vol_trans_node, 'in_file')

    import_MNI_masks_workflow.connect(input_node, 'standard2example_func', vol_trans_node, 'in_matrix_file')
    import_MNI_masks_workflow.connect(input_node, 'template_file', vol_trans_node, 'reference')
    import_MNI_masks_workflow.connect(input_node, 'output_folder', datasink, 'base_directory')
    
    import_MNI_masks_workflow.connect(vol_trans_node, 'out_file', output_node, 'output_files')
    import_MNI_masks_workflow.connect(vol_trans_node, 'out_file', datasink, 'roi.MNI')

    import_MNI_masks_workflow.run('MultiProc', plugin_args={'n_procs': 24})

    out_files = glob.glob(op.join(output_folder, 'roi', 'MNI', '*.nii.gz'))

    return out_files

