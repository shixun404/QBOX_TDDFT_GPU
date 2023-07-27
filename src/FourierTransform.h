////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// FourierTransform.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FOURIERTRANSFORM_H
#define FOURIERTRANSFORM_H

#include <complex>
#include <vector>

#if !( defined(OPTIMIZE_GPU) || defined(USE_FFTW2) || defined(USE_FFTW3) || defined(USE_ESSL_FFT) || defined(FFT_NOLIB) )
#error "Must define USE_FFTW2, USE_FFTW3, USE_ESSL_FFT or FFT_NOLIB"
#endif

#if defined(USE_FFTW2) && defined(USE_FFTW3)
#error "Cannot define USE_FFTW2 and USE_FFTW3"
#endif



#if USE_FFTW2
#if USE_DFFTW
#include "dfftw.h"
#else
#include "fftw.h"
#endif
#endif

#if USE_FFTW3
#include "fftw3.h"
#if USE_FFTW3MKL
#include "fftw3_mkl.h"
#endif
#endif

#if AUTOTUNER
#include "utilities.hpp"

#endif

#if OPTIMIZE_GPU
#include <cufft.h>
#include "cublas_v2.h"
#endif 
//TODO: Review

#include "Timer.h"
#include "BasisMapping.h"

class Basis;

class FourierTransform
{
  private:

#if OPTIMIZE_GPU
  static int my_dev;
  static cudaStream_t* cuda_streams;
  static int nstreams; //THIS IS A PERFORMANCE PARAMETER
  //static cublasHandle_t handle;
  static int nbatches; //THISIS A PERFORMANCE PARAMETER
  double * c_device;
  double * zvec_device;
  double * f_device;
  BasisMapping bm_;

  
#else
  const BasisMapping bm_;
#endif	  

  const Basis& basis_;

  const int np0_,np1_,np2_;
  const int nvec_;

  int ntrans0_,ntrans1_,ntrans2_;

  std::vector<std::complex<double> > zvec_;
  
#if OPTIMIZE_TRANSPOSE
  bool optimization_;
#endif
  
  
  void init_lib(void);

#if OPTIMIZE_GPU
  cufftHandle planFWD;
  cufftHandle planBWD;
  cufftHandle planFWDbatch;
  cufftHandle planBWDbatch;
#elif USE_ESSL_FFT
#if USE_ESSL_2DFFT
  std::vector<double> aux1xyf,aux1zf;
  std::vector<double> aux1xyb,aux1zb;
  std::vector<double> aux2;
  int naux1xy,naux1z,naux2;
#else
  std::vector<double> aux1xf, aux1yf, aux1zf;
  std::vector<double> aux1xb, aux1yb, aux1zb;
  std::vector<double> aux2;
  int naux1x,naux1y,naux1z,naux2;
#endif
#elif USE_FFTW2
  fftw_plan fwplan0,fwplan1,fwplan2,bwplan0,bwplan1,bwplan2;
#elif USE_FFTW3
  //plans for np2_
  fftw_plan fwplan, bwplan;
#if defined(USE_FFTW3_2D) || defined(USE_FFTW3_THREADS)
  fftw_plan fwplan2d, bwplan2d;
#else
  fftw_plan fwplanx, fwplany, bwplanx, bwplany;
#endif
#elif defined(FFT_NOLIB)
  // no library
#else
#error "Must define USE_FFTW2, USE_FFTW3, USE_ESSL_FFT, OPTIMIZE_GPU or FFT_NOLIB"
#endif

  void fxy(std::complex<double>* val);
  void fxy_inv(std::complex<double>* val);
  void fz(void);
  void fz_inv(void);
  void fwd(std::complex<double>* val);
  void bwd(std::complex<double>* val);

#if OPTIMIZE_GPU 
  void cuda_do_fft3d( const int fsign, const int  *n, const double scale,  double *data, double* data2,cufftHandle &plan, cudaStream_t cuda_stream=0, const int batches = 1);
#endif

  public:
  FourierTransform (const Basis &basis, int np0, int np1, int np2,int nstloc=1);
  ~FourierTransform ();

  
#if OPTIMIZE_GPU 

    static void initializeConst() {
#if AUTOTUNER
	nstreams=atoi(DFTuning::getEnv("QBOX_NSTREAMS"));
	nbatches=atoi(DFTuning::getEnv("QBOX_NBATCHES"));
#else
   	nstreams=2;
	nbatches=4;	
#endif
  }


  static cudaStream_t& get_cuda_streams(int i){return cuda_streams[i];}
  static int get_nstreams(){return nstreams;}
  static int get_my_dev(){return my_dev;}
  static int get_nbatches(){return nbatches;}
  //CUBLAS IMPLEMENTATION
//  static cublasHandle_t & get_cublasHandle(){return handle;}
  
  double* get_f_device() {return f_device;}
  double* get_c_device(){return c_device;} 
#endif 
  
  
  int nstloc_;  
  // backward: Fourier synthesis, compute real-space function
  // forward:  Fourier analysis, compute Fourier coefficients
  // forward transform includes scaling by 1/np012
  // single transforms: c -> f, f -> c
  
  
  
  
  void backward (const std::complex<double>* c, std::complex<double>* f, cudaStream_t cuda_stream=0, bool enable_htod = true, bool enable_dtoh = true,const int band_id=0, const int nbaches=1);
  //void backward (const std::complex<double>* c, std::complex<double>* f);
 
#if OPTIMIZE_TRANSPOSE
  bool get_optimization() const {return optimization_;}
  void backward1 (const std::complex<double>* c, int band=-1);
  void backward2 ();
  void backward3 (std::complex<double>* f,int band=-1); 

  void forward1(std::complex<double>* f, int band=-1);
  void forward2();
  void forward3(std::complex<double>* c, int band=-1);

#endif


  // Note: forward transforms overwrite the array f
  void forward(std::complex<double>* f, std::complex<double>* c, cudaStream_t cuda_stream=0, bool enable_htod = true, bool enable_dtoh = true, const int band_id=0, const int batches=1);

  // double transforms: c1 + i*c2 -> f, f -> c1 + i*c2
  void backward (const std::complex<double>* c1,
                 const std::complex<double>* c2, std::complex<double>* f);
  // Note: forward transforms overwrite the array f
  void forward(std::complex<double>* f,
               std::complex<double>* c1, std::complex<double>* c2);

  int np0() const { return np0_; }
  int np1() const { return np1_; }
  int np2() const { return np2_; }
  int np2_loc(void) const { return bm_.np2_loc(); }
  int np2_loc(int iproc) const { return bm_.np2_loc(iproc); }
  int np2_first(void) const { return bm_.np2_first(); }
  int np2_first(int iproc) const { return bm_.np2_first(iproc); }
  long int np012() const { return ((long int)np0_) * np1_ * np2_; }
  int np012loc(int iproc) const { return np0_ * np1_ * np2_loc(iproc); }
  int np012loc(void) const { return np0_ * np1_ * np2_loc(); }
  int index(int i, int j, int k) const
  { return i + np0_ * ( j +  np1_ * k ); }
  int i(int ind) const { return ind % np0_; }
  int j(int ind) const { return (ind / np0_) % np1_; }
  int k(int ind) const { return (ind / np0_) / np1_ + np2_first(); }

  void reset_timers(void);
  Timer tm_fwd, tm_bwd, tm_map_fwd, tm_map_bwd, tm_trans_fwd, tm_trans_bwd,
        tm_fxy, tm_fxy_inv, tm_fz, tm_fz_inv, tm_init;
};
#endif
