#-------------------------------------------------------------------------------
#
#  icc_mkl.mk
#
#-------------------------------------------------------------------------------
 MPIDIR=/cm/shared/apps/spack/cpu/opt/spack/linux-centos8-zen2/intel-19.1.1.217/mvapich2-2.3.4-bhl72eligau3gij652onqctsoau5x7tc
 XERCESCDIR=/home/p012chm/Utility/xerces-c-3.2.3
 COMPDIR=/cm/shared/apps/spack/cpu/opt/spack/linux-centos8-zen/gcc-8.3.1/intel-19.1.1.217-4d42ptjd6wsnh5bgbzcv6lp44vxpjwut
#COMPDIR=/cm/shared/apps/spack/cpu/opt/spack/linux-centos8-zen2/intel-19.1.1.217/intel-mkl-2019.1.144-inp2u4rjh6zrkdotn2fzyo35mksfyhgh
 FFTWDIR=/cm/shared/apps/spack/cpu/opt/spack/linux-centos8-zen2/intel-19.1.1.217/fftw-3.3.8-aeockdon5bukmh3bmztjsfvnpmtlluyx
 MKLDIR=/cm/shared/apps/spack/cpu/opt/spack/linux-centos8-zen2/intel-19.1.1.217/intel-mkl-2019.1.144-inp2u4rjh6zrkdotn2fzyo35mksfyhgh/mkl
 PLTOBJECTS = readTSC.o

#CXX=env OMPI_CXX=icc mpicxx
 CXX=icc
 LD=mpicxx
#LD=$(CXX)

 PLTFLAGS += -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES \
	     -DXERCES_3 -DMPICH_IGNORE_CXX_SEEK -DUSE_UUID \
	     -DUSE_FFTW3 -DFFTW3_2D

 INCLUDE = -I$(MKLDIR)/include -I$(XERCESCDIR)/include -I$(MPIDIR)/include -I$(FFTWDIR)/include -I$(COMPDIR)/include
#INCLUDE += -I$(FFTWDIR)/include -I$(MPIDIR)/include -I$(XERCESCDIR)/include -I$(COMPDIR)/include -I/usr/include/uuid

#CXXFLAGS=  -g -O3 -vec-report1 -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)
 CXXFLAGS=  -g -O3 $(INCLUDE) $(PLTFLAGS)
#CXXFLAGS=  -g -O3 -Dx86_64 $(INCLUDE) $(PLTFLAGS)

 LIBPATH += -L$(MPIDIR)/lib \
            -L$(MKLDIR)/lib/intel64  \
            -L$(XERCESCDIR)/lib \
	    -L$(COMPDIR)/lib/intel64 \
	    -L$(FFTWDIR)/lib \
	    -L/usr/lib64

 LIBS +=  -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 \
          -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
          -lmkl_lapack95_lp64 -lmkl_intel_thread \
          -lxerces-c -liomp5 -lpthread -luuid -lfftw3
#	  -lifcore -lirc -lsvml -lfftw3
#LIBS +=  -Wl,-Bstatic -Wl,--start-group \
#         -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 \
#         -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread \
#         -Wl,--end-group -Wl,-Bdynamic  \
#         -lxerces-c -liomp5 -lpthread -luuid -lfftw3

# Parallel libraries

 LDFLAGS = $(LIBPATH) $(LIBS)
#-------------------------------------------------------------------------------
