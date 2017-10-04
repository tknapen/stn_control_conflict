# from nipype import config
# config.enable_debug_mode()
# Importing necessary packages
import os
import sys
import os.path as op
import glob
import json
import nipype
from nipype import config, logging
import matplotlib.pyplot as plt
import nipype.interfaces.fsl as fsl
import nipype.pipeline.engine as pe
#
#   run as in:
#
# source deactivate
# source activate cf
# export SUBJECTS_DIR=/home/shared/2017/reward/pearl_3T/FS_SJID
# for s in {001..049}
# do
#     echo sub-$s
#     # python preprocessing.py sub-$s stop &
#     python preprocessing.py sub-$s rl ;
# done

# export SUBJECTS_DIR=/home/shared/2017/reward/pearl_3T/FS_SJID
# for s in {013..026}
# do
#     echo sub-$s
#     # python preprocessing.py sub-$s stop &
#     python preprocessing.py sub-$s rl ;
# done

# export SUBJECTS_DIR=/home/shared/2017/reward/pearl_3T/FS_SJID
# for s in {027..036}
# do
#     echo sub-$s
#     # python preprocessing.py sub-$s stop &
#     python preprocessing.py sub-$s rl ;
# done

# export SUBJECTS_DIR=/home/shared/2017/reward/pearl_3T/FS_SJID
# for s in {037..049}
# do
#     echo sub-$s
#     # python preprocessing.py sub-$s stop &
#     python preprocessing.py sub-$s rl ;
# done



# older masks
#  ,"maxSTN25exc","SST_GO_preSMA","SST_GO_rIFG", "Caudate", "PvmPFCNoventri","PstriatumNoVentri"
# 
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio
import nibabel as nib
from nipype.interfaces.utility import Function, Merge, IdentityInterface
from nipype.interfaces.io import SelectFiles, DataSink

from pearl.workflows.pearl_pp_workflow import create_pearl_pp_workflow
from IPython import embed as shell

# the subject id and experiment vars are commandline arguments to this script.
sub_id = str(sys.argv[1])
experiment = str(sys.argv[2])

# from pearl.parameters import *
# execfile('pearl/parameters.py')
exec(open("pearl/parameters.py").read())

# we set up the folders and logging there.
if not op.isdir(preprocessed_data_dir):
    try:
        os.makedirs(preprocessed_data_dir)
    except OSError:
        pass

try:
    os.makedirs(op.join(opd, 'log'))
except OSError:
    pass

config.update_config({  'logging': {
                                    'log_directory': op.join(opd, 'log'),
                                    'log_to_file': True,
                                    'workflow_level': 'INFO',
                                    'interface_level': 'DEBUG'
                                  },
                        'execution': {
                                    'stop_on_first_crash': True
                                    }
                    })
logging.update_logging(config)

# the actual workflow
pearl_pp_workflow = create_pearl_pp_workflow(analysis_info, name = 'pearl_pp')

if experiment == 'rl':
    pearl_pp_workflow.inputs.inputspec.exp_shorthand = 'RL'
elif experiment == 'map':
    pearl_pp_workflow.inputs.inputspec.exp_shorthand = 'ppaffalo'
elif experiment == 'stop':
    pearl_pp_workflow.inputs.inputspec.exp_shorthand = 'Stopsignal'

# standard output variables
pearl_pp_workflow.inputs.inputspec.raw_directory = raw_data_dir
pearl_pp_workflow.inputs.inputspec.sub_id = sub_id
pearl_pp_workflow.inputs.inputspec.output_directory = opd

pearl_pp_workflow.inputs.inputspec.which_file_is_EPI_space = analysis_info['which_file_is_EPI_space']

# registration details
pearl_pp_workflow.inputs.inputspec.FS_ID = sub_id
pearl_pp_workflow.inputs.inputspec.FS_subject_dir = FS_subject_dir
pearl_pp_workflow.inputs.inputspec.standard_file = op.join(os.environ['FSL_DIR'], 'data/standard/MNI152_T1_1mm_brain.nii.gz')

# percent signal change and average-across-runs settings
pearl_pp_workflow.inputs.inputspec.psc_func = analysis_info['psc_func']
pearl_pp_workflow.inputs.inputspec.tr = acquisition_parameters['RepetitionTime']

# write out the graph and run
pearl_pp_workflow.write_graph(opd + '.png')
pearl_pp_workflow.run('MultiProc', plugin_args={'n_procs': 24})
