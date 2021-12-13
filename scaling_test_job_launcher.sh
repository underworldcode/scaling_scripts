#!/bin/bash
source params.sh
if (($SCALING_TYPE==1)); then
   TYPE="Weak"
elif (($SCALING_TYPE==2)); then
   TYPE="Strong"
else
   echo "Scaling type must be 1 (weak) or 2 (strong)."
   exit 1
fi

export NAME="Results_${TYPE}_${UW_NAME}_DIM_${UW_DIM}_BASE_${SCALING_BASE}_ORDER_${UW_ORDER}_MAXITS_${UW_MAX_ITS}_TOL_${UW_SOL_TOLERANCE}_PENALTY_${UW_PENALTY}_IO_${UW_ENABLE_IO}_MODEL_${UW_MODEL}"

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
mkdir -p ${NAME}
cp *.sh ${NAME}
cp *.py ${NAME}
cd ${NAME}


for i in ${JOBS} 
do
   if (($SCALING_TYPE==1)); then
      export UW_RESOLUTION="$((${SCALING_BASE} * ${i}))"
   else   
      export UW_RESOLUTION=${SCALING_BASE}
   fi
   export NTASKS="$((${i}*${i}*${i}))"
   export EXPORTVARS="UW_RESOLUTION,NTASKS,UW_ENABLE_IO,UW_ORDER,UW_DIM,UW_SOL_TOLERANCE,UW_MAX_ITS,UW_PENALTY,UW_MODEL,PICKLENAME"
   if [ $BATCH_SYS == "PBS" ] ; then
      # memory requirement guess: 3GB * nprocs
      MEMORY="$((4*${NTASKS}))GB"
      PBSTASKS=`python3 <<<"print((int(${NTASKS}/48) + (${NTASKS} % 48 > 0))*48)"`  # round up to nearest 48 as required by nci
      CMD="qsub -v ${EXPORTVARS} -N ${NAME} -l storage=gdata/m18,ncpus=${PBSTASKS},mem=${MEMORY},walltime=${WALLTIME},wd -P ${ACCOUNT} -q ${QUEUE} gadi_baremetal_go.sh"
      echo ${CMD}
      ${CMD}
   else
      export IMAGE=/group/m18/singularity/underworld/underworld2_2.11.0b.sif
      #export IMAGE=/group/m18/singularity/underworld/underworld2_2.10.0b_rc.sif
      #export IMAGE=/group/m18/singularity/underworld/underworld2_v29.sif
      if [[ "$QUEUE" == "express" ]] ; then
         export QUEUE="debugq"
      else
         export QUEUE="workq"
      fi
      export OUTNAME="Res_"${UW_RESOLUTION}"_Nproc_"${NTASKS}"_JobID_"%j".out"

      CMD="sbatch --export=IMAGE,${EXPORTVARS} --job-name=${NAME} --output=${OUTNAME} --ntasks=${NTASKS} --time=${WALLTIME} --account=${ACCOUNT} --partition=${QUEUE} magnus_container_go.sh"
      echo ${CMD}
      ${CMD}
   fi

done

