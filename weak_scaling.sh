#!/bin/bash
export WEAK_SCALING_BASE=32
export JNAME="Weak_Scaling_${WEAK_SCALING_BASE}_UW28"

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

#for i in 1 2 4
for i in 1 2 4 8 12 16 20 24
do
   export UW_RESOLUTION="$((${WEAK_SCALING_BASE} * ${i}))"
   export UW_ENABLE_IO="0"

   export NTASKS="$((${i}*${i}*${i}))"
   export NAME=${JNAME}_IO_${UW_ENABLE_IO}

   if [ $BATCH_SYS == "PBS" ] ; then
      # memory requirement guess: 2GB * nprocs
      MEMORY="$((2*${NTASKS}))GB"

      # -v do pass env vars to PBS (raijin/nci) 
      CMD="qsub -v UW_RESOLUTION=${UW_RESOLUTION},UW_ENABLE_IO=${UW_ENABLE_IO} -N ${NAME} -l ncpus=${NTASKS},mem=${MEMORY} pbs_go.sh"
      echo ${CMD}
      ${CMD}
   else
      CMD="sbatch --job-name=${NAME} --ntasks=${NTASKS} shift_go.sh"
      echo ${CMD}
      ${CMD}
   fi

done

