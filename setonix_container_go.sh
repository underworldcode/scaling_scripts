#!/bin/bash -l
 
# load Singularity
module load singularity/3.8.6
 
export UW_VIS_PORT=0
 
# as per
# https://support.pawsey.org.au/documentation/pages/viewpage.action?pageId=116131367#UsewithSingularity-RunningPythonandR
# we unset all the host python-related ENV vars
unset $( env | grep ^PYTHON | cut -d = -f 1 | xargs )

scontrol show job ${SLURM_JOBID} -ddd
env
cat timed_model.py
echo ""
echo "---------- Running Job ----------"
echo ""
export TIME_LAUNCH_MPI=`date +%s%N | cut -b1-13`
srun --export=all -u -n ${SLURM_NTASKS} singularity exec -B ${PWD}:/work $IMAGE bash -c "cd work ; TIME_LAUNCH_PYTHON=\`date +%s%N | cut -b1-13\` python3 timed_model.py"
# execute
#srun --export=all -u -n $SLURM_NTASKS singularity exec -B ${PWD}:/work $containerImage bash -c "whoami; cd /work/; python3 $model"
