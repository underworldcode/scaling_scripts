#!/bin/bash
export JOBS="1 2 4 8 10 12 14 16" # 18 20 22 24"
export UW_NAME="V29_system_petsc_high_mem"
export UW_ENABLE_IO="0"

export WALLTIME="00:15:00"
export ACCOUNT="m18"

export WEAK_SCALING_BASE=16
export UW_ORDER=2
export UW_DIM=3
export UW_SOL_TOLERANCE=1e-6
export UW_PENALTY=-1.   # set to negative value to disable penalty 
export UW_MODEL="SolH" #  "SolDB3d" for dim3, though penalty method probably needed for q1.
export PICKLENAME=None  #"SolDB3d_Gadi_1e-11.pickle" #"conv_test_results_high_res_tighter_take2.pickle"  # set to "None" to disable conv testing
