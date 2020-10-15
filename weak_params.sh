#!/bin/bash
export JOBS="1 2 4 6 8 10"
export UW_NAME="bknight2d_test_take1"
export UW_ENABLE_IO="0"

export WALLTIME="00:10:00"
export ACCOUNT="m18"

export WEAK_SCALING_BASE=128
export UW_ORDER=1
export UW_SOL_TOLERANCE=1e-6
export UW_PENALTY=1e-3   # set to negative value to disable penalty 
export UW_SCRIPT="bknight_2d.py"
export UW_DIM=2