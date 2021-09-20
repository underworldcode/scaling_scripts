#!/bin/bash
export JOBS="4 5 6 7 8"
export UW_NAME="Strong_Scaling_Init_Test"
export SCALING_TYPE=2  # 1=weak, 2=strong

export UW_MODEL="SolDB3d"  #  "SolDB3d" for dim3, though penalty method probably needed for q1.
export PICKLENAME=None   #"SolDB3d_Gadi_1e-11.pickle" #"conv_test_results_high_res_tighter_take2.pickle"  # set to "None" to disable conv testing
export UW_ENABLE_IO="0"

export WALLTIME="00:15:00"
export ACCOUNT="director2186"
export QUEUE="normal" # normal or express  

export UW_ORDER=1
export UW_DIM=3
export SCALING_BASE=256

# Test style - UW_MAX_ITS (+ve, recommended >100): Fixed work, (-ve): Accuracy (UW_SOL_TERANCE is used)
export UW_MAX_ITS=50 # set to negative for accuracy test, positive for fixed iterative work irrepective of result fidelity
export UW_SOL_TOLERANCE=1e-9
export UW_PENALTY=-1.   # set to negative value to disable penalty 
