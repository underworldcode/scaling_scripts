#!/bin/bash
module load shifter

# write this script and env to stdout
scontrol show job ${SLURM_JOBID} -ddd
env
cat timed_model.py

srun -n ${NTASKS} shifter run --mpi $IMAGE  bash -c "TIME_LAUNCH_PYTHON=\`date +%s%N | cut -b1-13\` python3 timed_model.py"
