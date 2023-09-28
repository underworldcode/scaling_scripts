#!/bin/bash
#SBATCH --exclude=nid00[2024-2055],nid00[2792-2823] 

module load spack/0.19.0 petsc/3.19.4-3zidxfo py-numpy/1.23.4 py-h5py/3.7.0

# Use the python virtual environment indirectly, provent mishaps
export UWENV=/software/projects/pawsey0407/setonix/venv/py310/
export PATH=${UWENV}/bin/:$PATH
export PYTHONPATH=${UWENV}/lib/python3.10/site-packages/:$PYTHONPATH

export MPICH_OFI_STARTUP_CONNECT=1
export MPICH_OFI_VERBOSE=1
export OMP_NUM_THREADS=1
export MPICH_OFI_SKIP_NIC_SYMMETRY_TEST=1

export UW_IO_PATTERN=1 # for sequential I/O
#
# write this script and env to stdout
scontrol show job ${SLURM_JOBID} -ddd
env
cat timed_model.py

export TIME_LAUNCH_MPI=`date +%s%N | cut -b1-13`
srun -n ${NTASKS} bash -c "TIME_LAUNCH_PYTHON=\`date +%s%N | cut -b1-13\` python3 timed_model.py"
