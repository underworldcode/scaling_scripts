#!/bin/bash
module load singularity

export UW_VIS_PORT=0

# write this script and env to stdout
scontrol show job ${SLURM_JOBID} -ddd
env
cat timed_model.py
echo ""
echo "---------- Running Job ----------"
echo ""
export TIME_LAUNCH_MPI=`date +%s%N | cut -b1-13`
srun -u -n ${NTASKS} singularity exec $IMAGE bash -c "TIME_LAUNCH_PYTHON=\`date +%s%N | cut -b1-13\` python3 timed_model.py"
