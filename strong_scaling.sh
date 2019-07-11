#!/bin/bash
export STRONG_SCALING_BASE=64
export JNAME="Strong_Scaling_${STRONG_SCALING_BASE}_UW28"

for i in 1 2 4 8 10 
do
   export UW_RESOLUTION=${STRONG_SCALING_BASE}
   export NTASKS="$((${i}*${i}*${i}))"
   export UW_ENABLE_IO="0"
   export NAME=${JNAME}_IO_${UW_ENABLE_IO}
   sbatch --job-name=${NAME} --ntasks=${NTASKS} shift_go.sh
done

