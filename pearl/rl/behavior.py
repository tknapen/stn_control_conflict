from __future__ import division, print_function

def process_train_tsv(tsv_files, run_durations):
    import pandas as pd
    import numpy as np

    run_onset_offsets = np.r_[0, np.cumsum(run_durations)]
    trial_times = [pd.read_csv(tf, delimiter = '\t') for tf in tsv_files]
    for i, tt in enumerate(trial_times):
        trial_times[i]['onset'] += run_onset_offsets[i]
        trial_times[i]['Feedback_Onset'] += run_onset_offsets[i]

    tt_df = pd.concat(trial_times)
    
    return tt_df

def process_test_tsv(tsv_files, run_durations, event_conditions = ['ww', 'wl.u', 'wl.l', 'll']):
    import pandas as pd
    import numpy as np

    run_onset_offsets = np.r_[0, np.cumsum(run_durations)]
    trial_times = [pd.read_csv(tf, delimiter = '\t') for tf in tsv_files]

    for i, tt in enumerate(trial_times):
        trial_times[i]['onset'] += run_onset_offsets[i]

    tt_df = pd.concat(trial_times)
    tt_df = tt_df[tt_df['Response'] != 'omission']

    return tt_df