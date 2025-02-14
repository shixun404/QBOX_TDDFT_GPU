#-------------------------------------------------------------------------------
#
# module swap PrgEnv-gnu PrgEnv-nvhpc
# module load cray-fftw cray-libsci
# nvcc -c device_basis_mapping.cu
# -DOPTIMIZE_GPU
#-----------------------------------------------------------------------------
#
 PLT=Perlmutter
#-------------------------------------------------------------------------------

 PLTOBJECTS = readTSC.o

 CXX=CC
 LD=$(CXX)

 #XERCES=/global/cfs/cdirs/m2956/adrianpd/xerces-perl
 XERCES=/global/homes/s/sliu424/xerces-3.2.4
 #DFTUNING=/global/cfs/cdirs/m2956/adrianpd/DFTuning-master/dftuning-c


 OPT= -O3  -mp=multicore,gpu -gopt -traceback -gpu=cc80,cuda12.4 -gopt -traceback -cuda

 PLTFLAGS += $(OPT)  \
             -DOPTIMIZE_GPU -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES \
             -DMPICH_IGNORE_CXX_SEEK -DPARALLEL_FS -DOPTIMIZE_GPU#-DUSE_FFTW3

 INCLUDE = -I$(XERCES)/include -I$(DFTUNING)/inc

 CXXFLAGS= -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH = -L$(XERCES)/lib64 -L$(DFTUNING)/lib

 LIBS = -lfftw3_threads -lfftw3_omp -lfftw3 -lxerces-c -lcudart -lcufft -lcublas #-lutilities -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/cuda/lib64 

 LDFLAGS = $(LIBPATH) $(LIBS) -fopenmp -cudalib=cublas,cufft -cuda 

 NVCC = nvcc
 NVCCFLAGS = -arch=sm_80 # 对应 Perlmutter A100 GPU 
 %.o: %.cu
	$(CC) $(CXXFLAGS) -c $< -o $@                                                            
