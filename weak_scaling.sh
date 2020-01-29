#!/bin/bash
source weak_params.sh
export NAME="${UW_NAME}_DIM_${UW_DIM}_BASE_${WEAK_SCALING_BASE}_ORDER_${UW_ORDER}_TOL_${UW_SOL_TOLERANCE}_PENALTY_${UW_PENALTY}_IO_${UW_ENABLE_IO}_MODEL_${UW_MODEL}"

## find the BATCH environment ##
#################################
if qstat --version &> /dev/null ; then
   BATCH_SYS="PBS"
   export NAME="${NAME}_Gadi"
elif squeue --version &> /dev/null ; then
   BATCH_SYS="SLURM"
   export NAME="${NAME}_Magnus"
else
   echo "Can't determine batch system"
   exit 1
fi

echo "Batch system is $BATCH_SYS"
#################################

for i in ${JOBS} 
do
   export UW_RESOLUTION="$((${WEAK_SCALING_BASE} * ${i}))"
   export NTASKS="$((${i}*${i}*${i}))"
   export TIME_LAUNCH_JOB=`date +%s%N | cut -b1-13`

   if [ $BATCH_SYS == "PBS" ] ; then
      export QUEUE="normal" # normal or express

      # memory requirement guess: 3GB * nprocs
      MEMORY="$((3*${NTASKS}))GB"
      PBSTASKS=`python2<<<"print((${NTASKS}/48 + (${NTASKS} % 48 > 0))*48)"`  # round up to nearest 48 as required by nci
      # -V to pass all env vars to PBS (raijin/nci) 
      CMD="qsub -V -N ${NAME} -l ncpus=${PBSTASKS},mem=${MEMORY},walltime=${WALLTIME},wd -P ${ACCOUNT} -q ${QUEUE} gadi_baremetal_go.sh"
      echo ${CMD}
      ${CMD}
   else
      export IMAGE=/group/m18/singularity/underworld/underworld2_v29.sif
      export QUEUE="workq" # workq or debugq

      CMD="sbatch --job-name=${NAME} --ntasks=${NTASKS} --time=${WALLTIME} --account=${ACCOUNT} --partition=${QUEUE} magnus_container_go.sh"
      echo ${CMD}
      ${CMD}
   fi

done

