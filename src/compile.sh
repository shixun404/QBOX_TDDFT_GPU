#!/bin/bash


make clean

nvcc -c device_basis_mapping.cu -DAUTOTUNER -I/global/cfs/cdirs/m2956/adrianpd/DFTuning-master/dftuning-c/inc -L/global/cfs/cdirs/m2956/adrianpd/DFTuning-master/dftuning-c/lib -lutilities

make -j
