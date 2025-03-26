#!/bin/bash

module load petsc/3.21.3 python3/3.11.7 openmpi/4.1.7 python3-as-python hdf5/1.12.2p
 
UWENV=/g/data/m18/software/venv/uw216 # Project directory
export PATH=${UWENV}/bin:$PATH
export PYTHONPATH=${UWENV}/lib/python3.11/site-packages/:$PYTHONPATH

env
cat timed_model.py

export TIME_LAUNCH_MPI=`date +%s%N | cut -b1-13`
mpiexec -n ${NTASKS} bash -c "TIME_LAUNCH_PYTHON=\`date +%s%N | cut -b1-13\` python3 timed_model.py"
