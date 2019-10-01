#!/bin/bash
module purge
module load pbs dot mpi4py/3.0.2-py36-ompi3
export PYTHONPATH=/apps/mpi4py/3.0.2-py36-ompi3/lib/python3.6/site-packages/:/apps/underworld/opt/h5py/2.9.0-py36-ompi3/lib/python3.6/site-packages/h5py-2.9.0-py3.6-linux-x86_64.egg/:/apps/underworld/opt/lavavu/1.4.1_rc/:/apps/underworld/opt/pint/0.9_py36/lib/python3.6/site-packages/:/home/565/jam565/underworld:/home/565/jam565/underworld/underworld/libUnderworld/build/lib:/usr/local/lib/python2.7/site-packages
export LD_PRELOAD=/apps/openmpi/3.1.3/lib/libmpi.so

env
cat timed_model.py

mpirun -n ${NTASKS} bash -c "TIME_LAUNCH_PYTHON=\`date +%s%N | cut -b1-13\` python3 timed_model.py"

