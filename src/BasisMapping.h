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
// BasisMapping.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef BASISMAPPING_H
#define BASISMAPPING_H

#include <complex>
#include <vector>


class Basis;

class BasisMapping
{
  private:

  const Basis& basis_;
  int nprocs_, myproc_;
  int nstloc_;
  int np0_, np1_, np2_, np012loc_;
  int nvec_;

  std::vector<int> np2_loc_;   // np2_loc_[iproc], iproc=0, nprocs_-1
  std::vector<int> np2_first_; // np2_first_[iproc], iproc=0, nprocs_-1

  std::vector<int> scounts, sdispl, rcounts, rdispl;
  mutable std::vector<std::complex<double> > sbuf, rbuf;

  std::vector<int> ipack_, iunpack_, zvec_to_val_;


#if OPTIMIZE_GPU
  int * ip_device, * im_device;
  int * device_zvec_to_val;
  int * ip_, *im_;
#else
  std::vector<int>  ip_, im_;
#endif
  
  public:
  BasisMapping (const Basis &basis, int np0, int np1, int np2, int nstloc=1); 
 
#if OPTIMIZE_GPU
  ~BasisMapping();
#endif
   
  
  int np0(void) const { return np0_; }
  int np1(void) const { return np1_; }
  int np2(void) const { return np2_; }
  int np2_loc(void) const { return np2_loc_[myproc_]; }
  int np2_loc(int iproc) const { return np2_loc_[iproc]; }
  int np2_first(void) const { return np2_first_[myproc_]; }
  int np2_first(int iproc) const { return np2_first_[iproc]; }
  int np012loc(void) const { return np012loc_; }
  int nvec(void) const { return nvec_; }
  int zvec_size(void) const { return nvec_ * np2_; }


#if OPTIMIZE_GPU
  int allocate_device(cudaStream_t stream);
  void device_vector_to_zvec(const double *c, double *zvec, cudaStream_t stream, const int batch=1) const;
  void device_transpose_bwd(const double *zvec, double * ct, cudaStream_t stream, const int batch=1) const;
  void device_transpose_fwd(const double*ct, double * zvec, cudaStream_t stream, const int batch=1) const;
  void device_zvec_to_vector(const double * zvec, double * c, cudaStream_t stream, const int batch=1) const;
#endif



  // map a function c(G) to zvec_
  void vector_to_zvec(const std::complex<double> *c,
    std::complex<double> *zvec) const;
  // map two real functions c1(G) and c2(G) to zvec_
  void doublevector_to_zvec(const std::complex<double> *c1,
    const std::complex<double> *c2,std::complex<double> *zvec) const;
  // map zvec_ to a function c(G)
  void zvec_to_vector(const std::complex<double> *zvec,
    std::complex<double> *c) const;
  // map zvec_ to two real functions c1(G) and c2(G)
  void zvec_to_doublevector(const std::complex<double> *zvec,
    std::complex<double> *c1, std::complex<double> *c2) const;

  void transpose_bwd(const std::complex<double> *zvec,
                     std::complex<double> *ct) const;


#if OPTIMIZE_TRANSPOSE

  void transpose_fwd1(const std::complex<double> *ct,
                      int band=0) const;
  void transpose_fwd2() const;
  void transpose_fwd3(std::complex<double> *zvec, int band=0) const;
  void transpose_bwd1(const std::complex<double> *zvec,
  void transpose_bwd2() const;
  void transpose_bwd3(std::complex<double> *ct,int band=0) const;
#endif
  
  
void transpose_fwd(const std::complex<double> *ct,
                     std::complex<double> *zvec) const;
};
#endif
