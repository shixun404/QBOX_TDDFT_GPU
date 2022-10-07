#-------------------------------------------------------------------------------
#
#  cori.nersc.gov
#
# module load cray-fftw xerces
#
#-------------------------------------------------------------------------------
#
 PLT=Cori
#-------------------------------------------------------------------------------

 PLTOBJECTS = readTSC.o

 CXX=CC
 LD=$(CXX)

 OPT= -O3  -qopenmp -restrict

 PLTFLAGS += $(OPT)  \
             -DOPTIMIZE_TRANSPOSE -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES \
             -DMPICH_IGNORE_CXX_SEEK -DPARALLEL_FS -DUSE_FFTW3
 XERCES=/global/cfs/cdirs/m2956/adrianpd/xerces_cori
 INCLUDE = -I$(XERCES)/include

 CXXFLAGS= -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH = -L$(XERCES)/lib 

 LIBS =  -lfftw3 -lxerces-c

 LDFLAGS = $(LIBPATH) $(LIBS) -qopenmp 
