#!/bin/bash
#PBS -P m18
#PBS -q normal
#PBS -l walltime=00:10:00
##PBS -l ncpus=4
##PBS -l mem=6GB
#PBS -l wd
#PBS -M julian.giordani@unimelb.edu.au

module purge
module load pbs dot mpi4py/3.0.2-py36-ompi3

# uw-2.8 with petsc-3.11.3
export PYTHONPATH=/home/565/jug565/short/code/uw-2.8_petsc-3.11:/home/565/jug565/short/code/uw-2.8_petsc-3.11/libUnderworld/build/lib:/home/565/jug565/short/code/uw-2.8_petsc-3.11/glucifer:/apps/underworld/opt/h5py/2.9.0-py36-ompi3/lib/python3.6/site-packages/h5py-2.9.0-py3.6-linux-x86_64.egg/:/apps/underworld/opt/lavavu/1.4.2/:/apps/underworld/opt/pint/0.9_py36/lib/python3.6/site-packages/:/apps/mpi4py/3.0.2-py36-ompi3/lib/python3.6/site-packages/                                     
export LD_PRELOAD=/apps/openmpi/3.1.3/lib/libmpi.so

export PREFIXSTRING=`head /dev/urandom | tr -dc A-Za-z0-9 | head -c 13 ; echo ''`
export UW_JOB_NAME=${PBS_JOBNAME}
export UW_JOB_ID=${PBS_JOBID}
export TIME_LAUNCH_JOB=`date +%s%N | cut -b1-13`

#export OPENBLAS_NUM_THREADS=1 disabling for now

env
#mpiexec python3 -c "from mpi4py import MPI; print(MPI.COMM_WORLD.Get_size()) if MPI.COMM_WORLD.Get_rank() == 0 else None"
mpirun python3 rt_timed.py
