from __future__ import division, print_function

def avg_label_to_subject_label(avg_sj, trg_sj, label_directory, regmethod = 'surface' ):
    import os
    import glob

    cmd = 'mri_label2label'

    labels = glob.glob(os.path.join(os.environ['SUBJECTS_DIR'], avg_sj, 'label', label_directory, '*.label'))
    new_label_dir = os.path.join(os.environ['SUBJECTS_DIR'], trg_sj, 'label', label_directory)
    try:
        os.makedirs(new_label_dir)
    except:
        pass

    output_files = []

    for lbl in labels:
        label_file = os.path.split(lbl)[-1]
        hemisphere = label_file.split('.')[0]
        opf = os.path.join(new_label_dir, label_file)

        runcmd = cmd
        runcmd += ' --regmethod %s'%regmethod
        runcmd += ' --srclabel ' + lbl 
        runcmd += ' --srcsubject ' + avg_sj
        runcmd += ' --trgsubject ' + trg_sj
        runcmd += ' --hemi ' + hemisphere
        runcmd += ' --trglabel ' + opf

        os.system(runcmd)
        output_files.append(opf)

    return output_files



####
####    run this after preprocessing to 
####    incorporate novel definitions of surface based ROIs
####    into the preprocessed tree, after preprocessing.
####    folders are searched as defined in analysis_info['avg_subject_RS_label_folders']
####    Note that these actions on analysis_info['avg_subject_RS_label_folders'] are
####    automatically incorporated in the preprocessing pipeline
####    so is immediately run whenever the preprocessing is repeated.
####

# import os,glob
# import nipype.pipeline as pe
# from nipype.interfaces import fsl
# from spynoza.workflows.sub_workflows.masks import create_masks_from_surface_workflow

# experiment = 'rl'


# for s in range(1,50):#["sub-003",  "sub-005",  "sub-006",  "sub-012",  "sub-014",  "sub-016"]:
#     sub_id = 'sub-' + str(s).zfill(3)
#     execfile('pearl/parameters.py')
#     for label_directory in analysis_info['label_folders']:
#         mfs = create_masks_from_surface_workflow(name = 'masks_from_surface_nPRF')
#         mfs.inputs.inputspec.label_directory = str(label_directory)
#         mfs.inputs.inputspec.fill_thresh = 0.005
#         mfs.inputs.inputspec.re = '*.label'
#         mfs.inputs.inputspec.EPI_space_file = os.path.join(opd, 'reg', 'example_func.nii.gz')
#         mfs.inputs.inputspec.reg_file = os.path.join(opd, 'reg', 'register.dat')
#         mfs.inputs.inputspec.output_directory = os.path.join(opd, 'masks', str(label_directory))
#         mfs.inputs.inputspec.freesurfer_subject_dir = os.environ["SUBJECTS_DIR"]
#         mfs.inputs.inputspec.freesurfer_subject_ID = sub_id

#         mfs.run(plugin='MultiProc', plugin_args={'n_procs' : 6})

#         masks = glob.glob(os.path.join(opd, 'masks', str(label_directory), 'roi', '*.nii.gz'))

#         dilate_cortex = pe.MapNode(interface=
#             fsl.maths.DilateImage(operation = 'mean', kernel_shape = 'sphere', kernel_size = analysis_info['dilate_kernel_size']), 
#                         name='dilate_cortex', iterfield=['in_file']) 

#         dilate_cortex.inputs.in_file = masks
#         dilate_cortex.run()

#         os.system('cp -rf ' + os.path.join(dilate_cortex.output_dir(), 'mapflow', '_*', '*.nii.gz') + ' ' + os.path.join(opd, 'masks', str(label_directory)))

