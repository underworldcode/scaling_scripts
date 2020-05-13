Scaling Scripts
===============

This repo contains scripts to launch jobs for weak/strong scaling. 

To launch a set of weak scaling jobs, you will generally modify the parameters in `weak_params.sh` 
as necessary, and then execute `weak_scaling.sh` which will queue a set of jobs on the system. 

The executed jobs will generated JSON files which record timing information for the different execution
stages. A equivalent user readable TXT file will also be generated which records a summary. This data
is utilised by `scaling_graphs.ipynb` to generate graphs of the results. 


timed_model.py
--------------
This is an `underworld` script which executes a simple model and records timing information. 
The model covers most of the fundamental `underworld` constructs, including solvers, vis and data IO. 

weak_scaling.sh
---------------
This script will queue a set of jobs for testing underworld simulation weak scaling. It uses the settings 
in `weak_params.py` to determine run configuration. 

weak_params.sh
--------------
This script contains run parameters for weak scaling tests. The user should generally only need to 
modify values here. 

strong_scaling.sh
-----------------
This script will queue a set of jobs for testing underworld simulation strong scaling. Note that this 
script is old and requires updating. 

scaling_graphs.ipynb
--------------------
This notebook utilises simulation results to generate scaling graphs. 

magnus_container_go.sh
----------------------
This script configures an environment and launches the required script on Magnus using containers. 

magnus_baremetal_go.sh
----------------------
This script configures an environment and launches the required script on Magnus running directly on machine.

gadi_baremetal_go.sh
--------------------
This script configures an environment and launches the required script on Gadi running directly on machine.
