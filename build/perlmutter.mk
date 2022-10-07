#-------------------------------------------------------------------------------
#
# module swap PrgEnv-gnu PrgEnv-nvhpc
# module load cray-fftw cray-libsci
#
#-------------------------------------------------------------------------------
#
 PLT=Perlmutter
#-------------------------------------------------------------------------------

 PLTOBJECTS = readTSC.o

 CXX=CC
 LD=$(CXX)

 XERCES=/global/cfs/cdirs/m2956/adrianpd/xerces-perl

 OPT= -O3  -gopt -traceback  

 PLTFLAGS += $(OPT)  \
             -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES \
             -DMPICH_IGNORE_CXX_SEEK -DPARALLEL_FS -DUSE_FFTW3

 INCLUDE = -I$(XERCES)/include

 CXXFLAGS= -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH = -L$(XERCES)/lib -L$(DFTUNING)/lib

 LIBS =  -lfftw3_threads -lfftw3_omp -lfftw3 -lxerces-c 

 LDFLAGS = $(LIBPATH) $(LIBS) -fopenmp
