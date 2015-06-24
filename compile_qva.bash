#!/bin/bash

# Uncoment the following to load modules (where available).
#module load intel gsl cuda mkl

# The following environment variables must be specified.
#export CUDA_HOME=$HPC_CUDA_DIR
#export MKL_HOME=/opt/intel/composerxe/mkl
#export LIO_HOME=${HOME}/progs/lio
export QVA_HOME=${HOME}/dev/qumvia

cd ${QVA_HOME}/src

# Uncomment to compile CPU version only
make cpu
# Uncomment to compile the LIO-compatible (GPU) version only
#make lio
