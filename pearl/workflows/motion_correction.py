def curate_EPI_space(moco_target):
    """curate_EPI_space doubles the amount of timepoints of the 
    moco target if it has only one. This is mandatory for the AFNI moco workflow.
    Parameters
    ----------
    moco_target : string
        absolute path to moco file
    
    returns 
    out_file : string
        absolute path to moco file
    """
    from spynoza.nodes.utils import get_scaninfo
    import tempfile
    import nipype.pipeline as pe
    import nipype.interfaces.fsl as fsl
    import os.path as op

    TR, shape, dyns, voxsize, affine = get_scaninfo(moco_target)

    if len(shape) == 4: # 4D dataset
        out_file = moco_target

        print("CURATE MOCO: %s shape - %i"%(moco_target, len(shape)))
    elif len(shape) == 3: # 3D dataset
        tempdir = tempfile.mkdtemp()
        out_file = op.join(tempdir, op.split(moco_target)[-1][:-7]+'_dbl.nii.gz')

        print("CURATE MOCO: %s shape - %i"%(moco_target, len(shape)))

        mergenode = pe.Node(fsl.Merge(dimension='t'), name='merge')
        mergenode.inputs.in_files = [moco_target, moco_target]
        mergenode.inputs.merged_file = out_file
        mergenode.run()
        # out_file = mergenode.outputs.merged_file

    else:
        raise 'Moco target has %i dimensions, what am I to do?'%len(shape)

    return out_file

def create_motion_correction_workflow(name = 'moco', method = 'AFNI'):
    """uses sub-workflows to perform different registration steps.
    Requires fsl and freesurfer tools
    Parameters
    ----------
    name : string
        name of workflow
    
    Example
    -------
    >>> motion_correction_workflow = create_motion_correction_workflow('motion_correction_workflow')
    >>> motion_correction_workflow.inputs.inputspec.output_directory = '/data/project/raw/BIDS/sj_1/'
    >>> motion_correction_workflow.inputs.inputspec.in_files = ['sub-001.nii.gz','sub-002.nii.gz']
    >>> motion_correction_workflow.inputs.inputspec.which_file_is_EPI_space = 'middle'
 
    Inputs::
          inputspec.output_directory : directory in which to sink the result files
          inputspec.in_files : list of functional files
          inputspec.which_file_is_EPI_space : determines which file is the 'standard EPI space'
    Outputs::
           outputspec.EPI_space_file : standard EPI space file, one timepoint
           outputspec.motion_corrected_files : motion corrected files
           outputspec.motion_correction_plots : motion correction plots
           outputspec.motion_correction_parameters : motion correction parameters
    """
    import os.path as op
    import nipype.pipeline as pe
    import nipype.interfaces.fsl as fsl
    import nipype.interfaces.afni.preprocess as afni
    import nipype.interfaces.utility as util
    import nipype.interfaces.io as nio
    from nipype.interfaces.utility import Function, IdentityInterface
    import nipype.interfaces.utility as niu
    from spynoza.nodes import EPI_file_selector

    ### NODES
    input_node = pe.Node(IdentityInterface(fields=['in_files', 'output_directory', 'which_file_is_EPI_space',
                                                   'sub_id', 'tr']), name='inputspec')
    output_node = pe.Node(IdentityInterface(fields=([
                'motion_corrected_files', 
                'EPI_space_file', 
                'motion_correction_plots', 
                'motion_correction_parameters', 
                'extended_motion_correction_parameters', 
                'new_motion_correction_parameters'])), name='outputspec')

    ########################################################################################
    # Invariant nodes
    ########################################################################################

    EPI_file_selector_node = pe.Node(Function(input_names=['which_file', 'in_files'], output_names='raw_EPI_space_file',
                                       function=EPI_file_selector), name='EPI_file_selector_node')

    mean_bold = pe.Node(interface=fsl.maths.MeanImage(dimension='T'), name='mean_space')

    rename = pe.Node(niu.Rename(format_string='session_EPI_space',
                            keep_ext=True),
                    name='namer')  

    curate = pe.Node(Function(input_names=['moco_target'], output_names='out_file',
                                       function=curate_EPI_space), name='curate')

    ########################################################################################
    # Workflow
    ########################################################################################

    motion_correction_workflow = pe.Workflow(name=name)

    motion_correction_workflow.connect(input_node, 'which_file_is_EPI_space', EPI_file_selector_node, 'which_file')
    motion_correction_workflow.connect(input_node, 'in_files', EPI_file_selector_node, 'in_files')

    ########################################################################################
    # outputs via datasink
    ########################################################################################
    datasink = pe.Node(nio.DataSink(), name='sinker')
    datasink.inputs.parameterization = False

    # first link the workflow's output_directory into the datasink.
    motion_correction_workflow.connect(input_node, 'output_directory', datasink, 'base_directory')
    motion_correction_workflow.connect(input_node, 'sub_id', datasink, 'container')

    ########################################################################################
    # FSL MCFlirt
    ########################################################################################
    # new approach, which should aid in the joint motion correction of 
    # multiple sessions together, by pre-registering each run.
    # the strategy would be to, for each run, take the first TR
    # and FLIRT-align (6dof) it to the EPI_space file. 
    # then we can use this as an --infile argument to mcflirt.

    if method == 'FSL':
        motion_correct_EPI_space = pe.Node(interface=fsl.MCFLIRT(
                        cost = 'normcorr', 
                        interpolation = 'sinc',
                        mean_vol=True
                        ), name='realign_space')

        take_first_TR = pe.MapNode(fsl.ExtractROI(t_min=0, t_size=1), name='take_first_TR', iterfield = ['in_file'])

        # preregistration node is set up for rigid-body within-modality reg
        prereg_flirt_N = pe.MapNode(fsl.FLIRT(cost_func='normcorr', output_type = 'NIFTI_GZ', dof = 7, interp = 'sinc'), 
                            name = 'prereg_flirt_N', iterfield = ['in_file'])

        motion_correct_all = pe.MapNode(interface=fsl.MCFLIRT(
                        save_mats = True, 
                        save_plots = True, 
                        cost = 'normcorr', 
                        interpolation = 'sinc',
                        stats_imgs = True
                        ), name='realign_all',
                                    iterfield = ['in_file', 'init'])

        plot_motion = pe.MapNode(interface=fsl.PlotMotionParams(in_source='fsl'),
                                name='plot_motion',
                                iterfield=['in_file'])

        # extend_motion_pars = pe.MapNode(Function(input_names=['moco_par_file', 'tr'], output_names=['new_out_file', 'ext_out_file'],
        #                                function=_extend_motion_parameters), name='extend_motion_pars', iterfield = ['moco_par_file'])

        motion_correction_workflow.connect(EPI_file_selector_node, 'raw_EPI_space_file', curate, 'moco_target')
        motion_correction_workflow.connect(curate, 'out_file', motion_correct_EPI_space, 'in_file')

        motion_correction_workflow.connect(motion_correct_EPI_space, 'out_file', mean_bold, 'in_file')
        motion_correction_workflow.connect(mean_bold, 'out_file', motion_correct_all, 'ref_file')

        # the pre-registration
        motion_correction_workflow.connect(input_node, 'in_files', take_first_TR, 'in_file')
        motion_correction_workflow.connect(take_first_TR, 'roi_file', prereg_flirt_N, 'in_file')
        motion_correction_workflow.connect(mean_bold, 'out_file', prereg_flirt_N, 'reference')

        # motion correction across runs
        motion_correction_workflow.connect(prereg_flirt_N, 'out_matrix_file', motion_correct_all, 'init')
        motion_correction_workflow.connect(input_node, 'in_files', motion_correct_all, 'in_file')
            
        # output node, for later saving
        motion_correction_workflow.connect(mean_bold, 'out_file', output_node, 'EPI_space_file')
        motion_correction_workflow.connect(motion_correct_all, 'par_file', output_node, 'motion_correction_parameters')

        motion_correction_workflow.connect(motion_correct_all, 'out_file', output_node, 'motion_corrected_files')
        # motion_correction_workflow.connect(motion_correct_all, 'par_file', extend_motion_pars, 'moco_par_file')
        # motion_correction_workflow.connect(input_node, 'tr', extend_motion_pars, 'tr')
        # motion_correction_workflow.connect(extend_motion_pars, 'ext_out_file', output_node, 'extended_motion_correction_parameters')
        # motion_correction_workflow.connect(extend_motion_pars, 'new_out_file', output_node, 'new_motion_correction_parameters')


        ########################################################################################
        # Plot the estimated motion parameters
        ########################################################################################

        plot_motion.iterables = ('plot_type', ['rotations', 'translations'])
        motion_correction_workflow.connect(motion_correct_all, 'par_file', plot_motion, 'in_file')
        motion_correction_workflow.connect(plot_motion, 'out_file', output_node, 'motion_correction_plots')

        # and the output

        motion_correction_workflow.connect(mean_bold, 'out_file', rename, 'in_file')
        motion_correction_workflow.connect(rename, 'out_file', datasink, 'reg')

        motion_correction_workflow.connect(motion_correct_all, 'out_file', datasink, 'mcf')
        motion_correction_workflow.connect(motion_correct_all, 'par_file', datasink, 'mcf.motion_pars')
        motion_correction_workflow.connect(plot_motion, 'out_file', datasink, 'mcf.motion_plots')
        # motion_correction_workflow.connect(extend_motion_pars, 'ext_out_file', datasink, 'mcf.ext_motion_pars')
        # motion_correction_workflow.connect(extend_motion_pars, 'new_out_file', datasink, 'mcf.new_motion_pars')


    ########################################################################################
    # AFNI 3DVolReg
    ########################################################################################
    # for speed, we use AFNI's 3DVolReg brute-force.
    # this loses plotting of motion parameters but increases speed
    # we hold on to the same setup, first moco the selected run
    # and then moco everything to that image, but without the 
    # intermediate FLIRT step.

    if method == 'AFNI':
        motion_correct_EPI_space = pe.Node(interface=afni.Volreg(
                        outputtype = 'NIFTI_GZ', 
                        zpad = 5,
                        args = ' -cubic ' # -twopass -Fourier
                        ), name='realign_space')

        # take_first_TR = pe.MapNode(fsl.ExtractROI(t_min=0, t_size=1), name='take_first_TR', iterfield = ['in_file'])

        # prereg_all = pe.MapNode(interface=afni.Volreg(
        #                 outputtype = 'NIFTI_GZ', 
        #                 zpad = 5,
        #                 args = '-twopass -twodup  ' # -twopass 
        #                 ), name='realign_all',
        #                         iterfield = ['in_file'])

        motion_correct_all = pe.MapNode(interface=afni.Volreg(
                        outputtype = 'NIFTI_GZ', 
                        zpad = 5,
                        args = ' -cubic ' # -twopass 
                        ), name='realign_all',
                                iterfield = ['in_file'])
       

        # curate for moco between sessions
        motion_correction_workflow.connect(EPI_file_selector_node, 'raw_EPI_space_file', curate, 'moco_target')
        motion_correction_workflow.connect(curate, 'out_file', motion_correct_EPI_space, 'in_file')

        motion_correction_workflow.connect(motion_correct_EPI_space, 'out_file', mean_bold, 'in_file')
        # motion_correction_workflow.connect(input_node, 'in_files', take_first_TR, 'in_file')

        # # prereg files
        # motion_correction_workflow.connect(mean_bold, 'out_file', prereg_all, 'basefile')
        # motion_correction_workflow.connect(take_first_TR, 'out_file', prereg_all, 'in_file')


        # motion correction across runs
        motion_correction_workflow.connect(input_node, 'in_files', motion_correct_all, 'in_file')
        motion_correction_workflow.connect(mean_bold, 'out_file', motion_correct_all, 'basefile')
        # motion_correction_workflow.connect(mean_bold, 'out_file', motion_correct_all, 'rotparent')
        # motion_correction_workflow.connect(mean_bold, 'out_file', motion_correct_all, 'gridparent')
            
        # output node, for later saving
        motion_correction_workflow.connect(mean_bold, 'out_file', output_node, 'EPI_space_file')
        motion_correction_workflow.connect(motion_correct_all, 'md1d_file', output_node, 'max_displacement_info')
        motion_correction_workflow.connect(motion_correct_all, 'oned_file', output_node, 'motion_correction_parameter_info')
        motion_correction_workflow.connect(motion_correct_all, 'oned_matrix_save', output_node, 'motion_correction_parameter_matrix')

        motion_correction_workflow.connect(motion_correct_all, 'out_file', output_node, 'motion_corrected_files')

        # and the output
        motion_correction_workflow.connect(mean_bold, 'out_file', rename, 'in_file')
        motion_correction_workflow.connect(rename, 'out_file', datasink, 'reg')

        motion_correction_workflow.connect(motion_correct_all, 'out_file', datasink, 'mcf')
        motion_correction_workflow.connect(motion_correct_all, 'md1d_file', datasink, 'mcf.max_displacement_info')
        motion_correction_workflow.connect(motion_correct_all, 'oned_file', datasink, 'mcf.parameter_info')
        motion_correction_workflow.connect(motion_correct_all, 'oned_matrix_save', datasink, 'mcf.parameter_matrix')


    return motion_correction_workflow

