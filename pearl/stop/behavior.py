from __future__ import division, print_function

def process_tsv(tsv_files, run_durations):
    import pandas as pd
    import numpy as np

    run_onset_offsets = np.r_[0, np.cumsum(run_durations)]
    trial_times = [pd.read_csv(tf, delimiter = '\t') for tf in tsv_files]
    for i, tt in enumerate(trial_times):
        trial_times[i]['onset'] += run_onset_offsets[i]

    tt_df = pd.concat(trial_times)
    
    return tt_df