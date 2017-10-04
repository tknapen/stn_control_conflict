from __future__ import division, print_function
# freeview commands
# freeview -v /Applications/freesurfer_5.3/subjects/DE_110412/mri/T1.nii.gz -v example_func.nii.gz:reg=register.dat -f /Applications/freesurfer_5.3/subjects/DE_110412/surf/rh.white -f /Applications/freesurfer_5.3/subjects/DE_110412/surf/lh.white 
# -v ../deriv/all_rsq_s.nii.gz:reg=../reg/register.dat:colormap=heatscale:heatscale=-0.4,0,0.4

def make_complex_tiffs(surf_folder, fig_folder, sj_id, fs_id, condition = 'cv', pe_type = 'polar', mask = '_all', exit_when_ready = 1, delete_previous = False ):
    import os
    import glob
    import re
    import tempfile

    call_dir = os.getcwd()
    # find and create files and folders
    template_tcl_file = os.path.abspath('nPRF/surf/redraw_complex.tcl' )

    try:
        os.makedirs(fig_folder)
    except:
        pass
    if delete_previous:
        os.system('rm %s/*.tiff'%fig_folder)

    # make sure that polar is selected even if polar_sign is used as pe_type


    os.chdir(surf_folder)
    tempdir = tempfile.mkdtemp()

    for hemi in ['lh','rh']:
        REDict = {
        '---HEMI---': hemi,
        '---SJ---': sj_id, 
        '---FS_SJ---': fs_id, 
        '---FIGPATH---': fig_folder,
        '---CONDITION---': condition, 
        '---POLECC---': pe_type.split('_')[0], 
        '---POLECC_S---': pe_type, 
        '---MASK---': mask, 
        '---NAME---': fs_id,
        '---EXIT---': str(exit_when_ready),
        }

        sf = open(template_tcl_file,'r')
        working_string = sf.read()
        sf.close()
        for e in REDict:
            rS = re.compile(e)
            working_string = re.sub(rS, REDict[e], working_string)

        temporary_tcl_file = os.path.join(tempdir, 'temp_tcl.tcl')

        of = open(temporary_tcl_file, 'w')
        of.write(working_string)
        of.close()

        cmd = 'tksurfer {fs_id} {hemi} inflated -tcl {temporary_tcl_file}'.format(
            fs_id=fs_id, 
            hemi=hemi, 
            temporary_tcl_file=temporary_tcl_file)

        print(cmd)
        os.system(cmd)

    os.chdir(call_dir)


def make_scalar_tiffs(surf_folder, fig_folder, sj_id, fs_id, condition = 'cv', data_type = 'all_in_stim_ratio', mask = '_pos', exit_when_ready = 1, delete_previous = False ):
    import os
    import glob
    import re
    import tempfile

    call_dir = os.getcwd()
    # find and create files and folders
    template_tcl_file = os.path.abspath('nPRF/surf/redraw_scalar.tcl' )

    try:
        os.makedirs(fig_folder)
    except:
        pass
    if delete_previous:
        os.system('rm %s/*.tiff'%fig_folder)

    os.chdir(surf_folder)
    tempdir = tempfile.mkdtemp()

    for hemi in ['lh','rh']:
        REDict = {
        '---HEMI---': hemi,
        '---SJ---': sj_id, 
        '---FS_SJ---': fs_id, 
        '---FIGPATH---': fig_folder,
        '---CONDITION---': condition, 
        '---DATAYPE---': data_type,
        '---MASK---': mask, 
        '---NAME---': fs_id,
        '---EXIT---': str(exit_when_ready),
        }

        sf = open(template_tcl_file,'r')
        working_string = sf.read()
        sf.close()
        for e in REDict:
            rS = re.compile(e)
            working_string = re.sub(rS, REDict[e], working_string)

        temporary_tcl_file = os.path.join(tempdir, 'temp_tcl.tcl')

        of = open(temporary_tcl_file, 'w')
        of.write(working_string)
        of.close()

        cmd = 'tksurfer {fs_id} {hemi} inflated -tcl {temporary_tcl_file}'.format(
            fs_id=fs_id, 
            hemi=hemi, 
            temporary_tcl_file=temporary_tcl_file)

        print(cmd)
        os.system(cmd)

    os.chdir(call_dir)


def all2surf(all_deriv_nii_files, surf_folder, FS_subject_dir, fs_id, reg, pfr = None, target_subject = 'avg_6_selected', smooth = 0):
    import nipype.interfaces.freesurfer as fs
    import os
    import glob

    try:
        os.makedirs(surf_folder)
    except:
        pass

    if pfr == None:
        pfr = (-0.5,1.5,0.2) # 0.99

    # cmd = 'mri_vol2surf --mov {mov} --reg {reg} --trgsubject {fs_id} --hemi {hemi} --surf smoothwm --projfrac {pf} --o {o}'
    os.chdir(surf_folder)

    for dnf in all_deriv_nii_files:
        for hemi in ['lh','rh']:
            of = os.path.join(surf_folder, 
                                            os.path.split(dnf)[-1][:-7] + '-' + hemi + '.mgz')

            sampler = fs.SampleToSurface(hemi=hemi)
            sampler.inputs.subject_id = fs_id
            sampler.inputs.source_file = dnf
            sampler.inputs.reg_file = reg
            sampler.inputs.sampling_units = 'frac'
            sampler.inputs.sampling_method = 'average'
            sampler.inputs.sampling_range = pfr
            sampler.inputs.surface = 'smoothwm'
            sampler.inputs.out_type = 'mgz'
            sampler.inputs.out_file = of

            res = sampler.run() 

            # and to the average subject
            to_av = fs.SurfaceTransform(hemi=hemi)
            to_av.inputs.source_subject = fs_id
            to_av.inputs.target_subject = target_subject
            to_av.inputs.source_file = of
            to_av.inputs.source_type = 'mgz'
            to_av.inputs.subjects_dir = FS_subject_dir

            to_av.run()

            new_orig_of = os.path.join(surf_folder, 
                                            os.path.split(dnf)[-1][:-7] + '-' + hemi + '.' + fs_id + '.mgz')

            os.system('mv %s %s'%(of, new_orig_of))

def av_surf_across_sjs(all_surf_files, output_filename):
    import os
    import nibabel as nb
    import numpy as np

    all_data = []
    for sfn in all_surf_files:
        all_data.append(nb.load(sfn).get_data())
        mean_placeholder = nb.load(sfn)

    # now that we have all the files, let's average them together:
    mean_data = np.mean(all_data,axis=0)
    # and save
    mean_mgh = nb.MGHImage(mean_data,mean_placeholder.affine,mean_placeholder.header)
    nb.save(mean_mgh,output_filename)


def combine_inflated_tiffs(fig_folder, sub_id, fs_id, mask_signs, data_types, hemisphere, view_nr, base_rotation):
    import os
    import numpy as np
    import matplotlib.pyplot as pl
    import seaborn as sn
    import glob

    f = pl.figure(figsize = (6*len(mask_signs), 6*len(data_types)))
    for i, ms in enumerate(mask_signs):
        for j, dt in enumerate(data_types):
            ax = f.add_subplot(len(mask_signs), len(data_types), i*len(mask_signs)+j+1)

            tiff_file = glob.glob(os.path.join(fig_folder, dt, ms + '*' + dt + '*' + hemisphere + '-' + str(view_nr).zfill(3) + '-' + str(int(base_rotation)) + '-' + fs_id +'.tiff'))[0]

            ax.imshow(pl.imread(tiff_file))
            ax.set_title(ms + ' ' + dt)
            pl.axis('off')
    output_filename = os.path.join(fig_folder, sub_id + '_' + hemisphere + '.' + str(view_nr) + '.' + str(base_rotation) + '.' + '-'.join(mask_signs) + '_' + '-'.join(data_types) + '.pdf')
    pl.tight_layout()
    # pl.show()
    pl.savefig(output_filename)
