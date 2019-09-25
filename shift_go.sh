#!/bin/bash

##SBATCH --nodes=1
##SBATCH --ntasks-per-node=24
#SBATCH --time=00:10:00
#SBATCH --account=m18
##SBATCH --partition=debugq

# write this script and env to stdout
scontrol show job ${SLURM_JOBID} -ddd
env

export IMAGE=underworldcode/underworld2:v2.8_release

module load shifter
#export UW_RESOLUTION=256
#export NPROCS="$((${SLURM_NNODES} * ${SLURM_TASKS_PER_NODE}))"
export PREFIXSTRING=`head /dev/urandom | tr -dc A-Za-z0-9 | head -c 13 ; echo ''`
export UW_JOB_NAME=${SLURM_JOB_NAME}
export UW_JOB_ID=${SLURM_JOB_ID}
export TIME_LAUNCH_JOB=`date +%s%N | cut -b1-13`

srun -n ${SLURM_NTASKS} shifter run --mpi $IMAGE  bash -c "TIME_LAUNCH_PYTHON=\`date +%s%N | cut -b1-13\` python rt_timed.py"
