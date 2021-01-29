#!/bin/bash

# the follow load the full software stack and running environment on gadi
source /g/data/m18/codes/gadi_setup.sh
#source /scratch/m18/codes/UWGeodynamics_2.10.sh
#source /g/data/m18/codes/UWGeodynamics_2.10.sh
#source /scratch/m18/codes/UWGeodynamics_2.10.2.sh

env
cat timed_model.py

export TIME_LAUNCH_MPI=`date +%s%N | cut -b1-13`
mpiexec -n ${NTASKS} bash -c "TIME_LAUNCH_PYTHON=\`date +%s%N | cut -b1-13\` python3 timed_model.py"
