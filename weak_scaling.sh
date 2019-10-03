#!/bin/bash
export JOBS="4 8 10 12 14 16"
export UW_NAME="Weak_Dev"
export UW_ENABLE_IO="0"

export WALLTIME="00:60:00"
export ACCOUNT="m18"

export WEAK_SCALING_BASE=32
export UW_ORDER=1
export UW_DIM=3
export UW_SOL_TOLERANCE=1e-8
export UW_MODEL="SolH" # Use "SolH" for dim3,order1, use "SolDB3d" for dim3,order2.
export PICKLENAME="conv_test_results_high_res_tighter_tolerances_superconv_errors.pickle"  # set to "None" to disable conv testing

export NAME="${UW_NAME}_DIM_${UW_DIM}_BASE_${WEAK_SCALING_BASE}_ORDER_${UW_ORDER}_TOL_${UW_SOL_TOLERANCE}_IO_${UW_ENABLE_IO}_MODEL_${UW_MODEL}"

## find the BATCH environment ##
#################################
if qstat --version > /dev/null ; then
   BATCH_SYS="PBS"
   export NAME="${NAME}_Raijin"
elif squeue --version > /dev/null ; then
   BATCH_SYS="SLURM"
   export NAME="${NAME}_Magnus"
else
   echo "Can't workout batch system"
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
      PBSTASKS=`python2<<<"print((${NTASKS}/16 + (${NTASKS} % 16 > 0))*16)"`  # round up to nearest 16 as required by nci
      # -V to pass all env vars to PBS (raijin/nci) 
      CMD="qsub -V -N ${NAME} -l ncpus=${PBSTASKS},mem=${MEMORY},walltime=${WALLTIME},wd -P ${ACCOUNT} -q ${QUEUE} raijin_baremetal_go.sh"
      echo ${CMD}
      ${CMD}
   else
      export IMAGE=underworldcode/underworld2:dev
      export QUEUE="workq" # workq or debugq

      CMD="sbatch --job-name=${NAME} --ntasks=${NTASKS} --time=${WALLTIME} --account=${ACCOUNT} --partition=${QUEUE} magnus_shifter_go.sh"
      echo ${CMD}
      ${CMD}
   fi

done
