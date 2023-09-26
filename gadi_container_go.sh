#!/bin/bash -l
 
# load Singularity
module load singularity
module load openmpi/4.1.4
 

export singularityDir=/g/data/m18/software/underworld/

# Define the container to use
export containerImage=$singularityDir/underworld-215.sif
 
env
cat timed_model.py
echo ""
echo "---------- Running Job ----------"
echo ""
export TIME_LAUNCH_MPI=`date +%s%N | cut -b1-13`
mpiexec -n $PBS_NCPUS singularity exec $containerImage bash -c "TIME_LAUNCH_PYTHON=\`date +%s%N | cut -b1-13\` python3 timed_model.py"
