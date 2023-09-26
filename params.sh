#!/bin/bash
export JOBS="1 2 4 6 8"
#export JOBS="1 2 4 6 8 10 12 14 16 18 20 22"
#export JOBS="22 22"
export UW_NAME="2.15.0_WeakScaling"
export SCALING_TYPE=1  # 1=weak, 2=strong

export UW_MODEL="SolDB3d"  #  "SolDB3d" for dim3, though penalty method probably needed for q1.
export PICKLENAME=None   #"SolDB3d_Gadi_1e-11.pickle" #"conv_test_results_high_res_tighter_take2.pickle"  # set to "None" to disable conv testing
export UW_ENABLE_IO="0"

export WALLTIME="00:10:00"
export ACCOUNT="m18"
export QUEUE="normal" # normal or express  

export UW_ORDER=1
export UW_DIM=3
export SCALING_BASE=16

# Test style - UW_MAX_ITS (+ve, recommended >100): Fixed work, (-ve): Accuracy (UW_SOL_TERANCE is used)
export UW_MAX_ITS=50 # set to negative for accuracy test, positive for fixed iterative work irrepective of result fidelity
export UW_SOL_TOLERANCE=1e-9
export UW_PENALTY=-1.   # set to negative value to disable penalty 
