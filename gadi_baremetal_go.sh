#!/bin/bash
module purge
module load pbs openmpi/4.0.2 hdf5/1.10.5p python3/3.7.4
export PYTHONPATH=/scratch/m18/opt/:/scratch/m18/opt/underworld/lib:/scratch/m18/opt/lavavu/

export LD_PRELOAD=/apps/openmpi-mofed4.7-pbs19.2/4.0.2/lib/libmpi_usempif08_GNU.so.40:/apps/openmpi-mofed4.7-pbs19.2/4.0.2/lib/libmpi_usempi_ignore_tkr_GNU.so.40:/apps/openmpi-mofed4.7-pbs19.2/4.0.2/lib/libmpi_cxx.so.40

export OPENBLAS_NUM_THREADS=1

env
cat timed_model.py

export TIME_LAUNCH_MPI=`date +%s%N | cut -b1-13`
mpirun -n ${NTASKS} bash -c "TIME_LAUNCH_PYTHON=\`date +%s%N | cut -b1-13\` python3 timed_model.py"

