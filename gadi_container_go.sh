#!/bin/bash -l
 
# load Singularity
module load singularity
module load openmpi/4.1.4
 

#export singularityDir=/home/565/jug565/scratch_m18/
export singularityDir=/scratch/m18/jug565/codes/

# Define the container to use
#export containerImage=$singularityDir/underworld-215.sif
export containerImage=$singularityDir/uw-rhel.sif
 
env
cat timed_model.py
echo ""
echo "---------- Running Job ----------"
echo ""
export TIME_LAUNCH_MPI=`date +%s%N | cut -b1-13`
#$PBS_NCPUS vs $NTASKS
mpiexec -n $NTASKS singularity exec -B /opt/pbs/default/lib/,/half-root/ --bind /lib64:/glib64 \
  $containerImage \
  bash -c "TIME_LAUNCH_PYTHON=\`date +%s%N | cut -b1-13\` LD_LIBRARY_PATH=/apps/openmpi/4.1.4/lib:/glib64:$LD_LIBRARY_PATH python3 timed_model.py"
#mpiexec -n $PBS_NCPUS singularity exec $containerImage bash -c "TIME_LAUNCH_PYTHON=\`date +%s%N | cut -b1-13\` LD_LIBRARY_PATH=/apps/openmpi/4.1.4/lib:$LD_LIBRARY_PATH python3 timed_model.py"
