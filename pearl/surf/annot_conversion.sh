
for s in {001..049}
do
	echo sub-$s
	for h in 'lh' 'rh'
		do
		echo $h
		mri_surf2surf --srcsubject fsaverage --trgsubject sub-$s --hemi $h --sval-annot $FREESURFER_HOME/subjects/fsaverage/label/$h.Yeo2011_7Networks_N1000.annot --tval $SUBJECTS_DIR/sub-$s/label/$h.Yeo2011_7Networks_N1000.annot 
		mri_annotation2label --annotation $SUBJECTS_DIR/sub-$s/label/$h.Yeo2011_7Networks_N1000.annot --subject sub-$s --hemi $h --surface smoothwm --outdir $SUBJECTS_DIR/sub-$s/label/Yeo2011_7Networks
		mri_surf2surf --srcsubject fsaverage --trgsubject sub-$s --hemi $h --sval-annot $FREESURFER_HOME/subjects/fsaverage/label/$h.aparc.a2009s.annot --tval $SUBJECTS_DIR/sub-$s/label/$h.aparc.a2009s.annot 
		mri_annotation2label --annotation $SUBJECTS_DIR/sub-$s/label/$h.aparc.a2009s.annot --subject sub-$s --hemi $h --surface smoothwm --outdir $SUBJECTS_DIR/sub-$s/label/aparc.a2009s
		done
done
# for s in {001..049}
# do
# 	echo sub-$s
# 	for h in 'lh' 'rh'
# 		do
# 		echo $h
# 		mri_surf2surf --srcsubject fsaverage --trgsubject sub-$s --hemi $h --sval-annot $FREESURFER_HOME/subjects/fsaverage/label/$h.aparc.a2009s.annot --tval $SUBJECTS_DIR/sub-$s/label/$h.aparc.a2009s.annot 
# 		mri_annotation2label --annotation $SUBJECTS_DIR/sub-$s/label/$h.aparc.a2009s.annot --subject sub-$s --hemi $h --surface smoothwm --outdir $SUBJECTS_DIR/sub-$s/label/aparc.a2009s
# 		done
# done