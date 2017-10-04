# surface stuff
surface conversion done using freesurfer.

I have made an average subject, avg_6_selected, from the subjects that have performed this experiment to criterion. On the surface of this average subject, I have delineated regions of consistent negative PRF-like activations for further analysis. They are stored in a `label/nPRF` folder. 

For per-participant labels, we use mris_label2label to convert the average to single subject labels, which can then be converted to volume with the masking workflow already used in preprocessing. 
