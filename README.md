# QBox-TDDFT
## Installation

1. Make sure you have installed xerces in your platform. Load fftw library too.
2. Go to ``build/`` and choose one of the existing templates. Customize it to your platform. Within the script, load the fftw/xerces modules or specify their binary paths.
3. If your build template was named ``template.mk``, then type ``make TARGET=../build/template -j 16`` from the src/ directory.

Note: If you want to execute the GPU version, please checkout to the corresponding branch and use a proper Programming Environment (e.g. ``PrgEnv-nvidia`` in Perlmutter). Also, add "DOPTIMIZE_GPU" to ``template.mk`` and comment the ``-DUSE_FFTW3``flag.

### Installing XERCES locally

0. Download the Xerces 3.2.4 library from their website ``wget https://dlcdn.apache.org//xerces/c/3/sources/xerces-c-3.2.4.zip`` on your local directory and unzip the file with ``unzip xerces-c-3.2.4.zip``
1. ``cd xerces-c-3.2.4`` and ``mkdir build`` 
2. ``cd build``
3. ``../configure --disable-static CC=cc CXX=CC CFLAGS=-O3 CXXFLAGS=-O3 --prefix=path_to_desired_installation_folder``
4. ``PREFIX=path_to_desired_installation_folder make install -j8``


## Run the example

1. Go to the ``example_RTTDDFT`` folder and choose the desired example. Modify ``Batch.sh`` with your platform requirements.


