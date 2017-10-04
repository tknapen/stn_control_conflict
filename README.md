# pearl_3T

This is the analysis of the CIRAS project run by Sara Jahfari. The repo contains nipype preprocessing, as well as the FIR analysis on ROI level. 

To run, for example for the **stop** experiment:

```bash
export SUBJECTS_DIR=/home/shared/2017/reward/pearl_3T/FS_SJID
for s in {001..049}
do
    echo sub-$s
    python postprocessing.py sub-$s rl test &
done
python across.py stop Stop

```


