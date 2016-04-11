#!/bin/bash

# Uncoment the following to load modules (where available).
#module load intel/2016.0.109 gsl cuda/7.0 
module load intel/2013 gsl cuda/5.5

# The following environment variables must be specified.
export CUDA_HOME=$HPC_CUDA_DIR
export CUDAHOME=$HPC_CUDA_DIR
export MKL_HOME=$MKLROOT
export LIO_HOME=${HOME}/dev/lio
export LIOHOME=${HOME}/dev/lio
export GARCHAHOME=${HOME}/dev/lio
export QVA_HOME=${HOME}/dev/qumvia

echo $CUDA_HOME
echo $MKL_HOME

cd ${QVA_HOME}/src

# Uncomment to compile CPU version only
make cpu
# Uncomment to compile the LIO-compatible (GPU) version only
make lio
# Uncomment to compile the LIO-compatible (GPU) library (for qumvia-amber) only
#make lib
