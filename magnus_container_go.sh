#!/bin/bash
module load singularity

# write this script and env to stdout
scontrol show job ${SLURM_JOBID} -ddd
env
cat timed_model.py

export TIME_LAUNCH_MPI=`date +%s%N | cut -b1-13`
srun -n ${NTASKS} singularity exec $IMAGE bash -c "TIME_LAUNCH_PYTHON=\`date +%s%N | cut -b1-13\` python3 timed_model.py"
