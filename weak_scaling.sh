#!/bin/bash
export WEAK_SCALING_BASE=32
export JNAME="Weak_Scaling_${WEAK_SCALING_BASE}_UW28"

for i in 1 2 4 8 10 
do
   export UW_RESOLUTION="$((${WEAK_SCALING_BASE} * ${i}))"
   export NTASKS="$((${i}*${i}*${i}))"
   export UW_ENABLE_IO="0"
   export NAME=${JNAME}_IO_${UW_ENABLE_IO}
   sbatch --job-name=${NAME} --ntasks=${NTASKS} shift_go.sh
done

