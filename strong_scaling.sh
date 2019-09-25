#!/bin/bash
export STRONG_SCALING_BASE=64
export JNAME="Strong_Scaling_${STRONG_SCALING_BASE}_UW28"

## find the BATCH environment ##
#################################
if qstat --version > /dev/null ; then
   BATCH_SYS="PBS"
elif squeue --version > /dev/null ; then
   BATCH_SYS="SLURM"
else
   echo "Can't workout batch system"
   exit 1
fi

echo "Batch system is $BATCH_SYS"
#################################

for i in 1 2 4 8 12 
do
   export UW_RESOLUTION=${STRONG_SCALING_BASE}
   export NTASKS="$((${i}*${i}*${i}))"
   export UW_ENABLE_IO="0"
   export NAME=${JNAME}_IO_${UW_ENABLE_IO}
   if [ ${BATCH_SYS} = "PBS" ]; then
      # memory requirement guess: 2GB * nprocs
      MEMORY="$((4*${NTASKS}))GB"
      # -V passes all env vars to qsub job
      CMD="qsub -V -N ${NAME} -l ncpus=${NTASKS},mem=${MEMORY} pbs_go.sh"
      echo ${CMD}
      ${CMD}
   else
      CMD="sbatch --job-name=${NAME} --ntasks=${NTASKS} shift_go.sh"
      echo ${CMD}
      ${CMD}
   fi
done

