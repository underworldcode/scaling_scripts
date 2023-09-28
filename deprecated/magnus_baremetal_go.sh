#!/bin/bash
module swap PrgEnv-cray PrgEnv-gnu
module load /home/jmansour/software/cle60up05/modulefiles/petsc/3.9.4
module load python/3.6.3
module load numpy mpi4py
export PYTHONPATH=/home/jmansour/underworld2:/home/jmansour/underworld2/glucifer:$PYTHONPATH

# write this script and env to stdout
scontrol show job ${SLURM_JOBID} -ddd
env
cat timed_model.py

export TIME_LAUNCH_MPI=`date +%s%N | cut -b1-13`
srun -n ${NTASKS} bash -c "TIME_LAUNCH_PYTHON=\`date +%s%N | cut -b1-13\` python3 timed_model.py"
