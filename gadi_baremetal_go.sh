#!/bin/bash

# the follow load the full software stack and running environment on gadi
module load python3
module load openmpi
export PYTHONPATH=/home/565/jam565/petsc/build_opt_petsc-3.15.0_ompi41/lib/:/home/565/jam565/underworld3

cat timed_model.py

export TIME_LAUNCH_MPI=`date +%s%N | cut -b1-13`
mpiexec -n ${NTASKS} bash -c "TIME_LAUNCH_PYTHON=\`date +%s%N | cut -b1-13\` python3 timed_model.py"
