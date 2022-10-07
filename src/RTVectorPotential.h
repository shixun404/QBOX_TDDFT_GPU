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
// RTVectorPotential.h
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RTVECTORPOTENTIAL_H
#define RTVECTORPOTENTIAL_H

#include <iostream>
#include <vector>
#include <valarray>
#include <complex>
#include "Basis.h"
#include "D3vector.h"
#include "Sample.h"
#include "UnitCell.h"

class RTVectorPotential 
{
  private:
  string calind_; // calculate induced vector potential T | F
  D3vector induced_;
  D3vector vecpot_; // vector potential
  D3vector velocity_;
  D3vector accel_;
  double nvecpot_; // norm of vector potential

  D3vector rtvpamp_;
  string rtvpeq_;
  double rtvpfreq_;
  double rtvpsigma_;
  double rtvpgauss_;
  double rtvplength_;
  double rtvpdelay_;

  public:

  RTVectorPotential(Sample& s);
  ~RTVectorPotential();

  double* get_kpgpa(const Basis& basis) const
  {
    const double* kpgpa2 = get_kpgpa2(basis);
    double* kpgpa = new double[basis.localsize()];
    for( int ig = 0; ig < basis.localsize(); ig++ )
    {
      kpgpa[ig] = sqrt(kpgpa2[ig]);
    }
    delete [] kpgpa2;
    return kpgpa;
  }

  double* get_kpgpa2(const Basis& basis) const
  {
    double* kpgpa2 = new double[basis.localsize()];
    for( int ig = 0; ig < basis.localsize(); ig++ )
    {
      kpgpa2[ig] = basis.kpg2_ptr()[ig] + (nvecpot());
      kpgpa2[ig] += 2 * vecpot_[0]*basis.kpgx_ptr(0)[ig];
      kpgpa2[ig] += 2 * vecpot_[1]*basis.kpgx_ptr(1)[ig];
      kpgpa2[ig] += 2 * vecpot_[2]*basis.kpgx_ptr(2)[ig];
    }
    return kpgpa2;
  }

  double* get_kpgpai(const Basis& basis) const
  {
    const double* kpgpa2 = get_kpgpa2(basis);
    double* kpgpai = new double[basis.localsize()];
    for( int ig = 0; ig < basis.localsize(); ig++ )
    {
      if( kpgpa2[ig] > 0.0 )
      {
        kpgpai[ig] = 1.0/sqrt(kpgpa2[ig]);
      }
      else
      {
        kpgpai[ig] = 0.0;
      }
    }
    delete [] kpgpa2;
    return kpgpai;
  }

  double* get_kpgpax(const Basis& basis, int j) const
  {
    if( j == 0 ) 
    {
      double* kpgpax = new double[3*basis.localsize()];
      for( int i = 0; i < 3; i++ )
      {
        for( int ig = 0; ig < basis.localsize(); ig++ )
	{
          kpgpax[i*basis.localsize()+ig] = basis.kpgx_ptr(i)[ig] + vecpot_[i];
	}
      }
      return kpgpax;
    }
    else
    {
      double* kpgpax = new double[basis.localsize()];
      for( int ig = 0; ig < basis.localsize(); ig++ )
      {
        kpgpax[ig] = basis.kpgx_ptr(j)[ig] + vecpot_[j];
      }
      return kpgpax;
    }
  }

  void calculate_acceleration(const double& dt, const D3vector& total_current, const UnitCell& cell);
  void vp_propagate(int step, const double& dt);
  const D3vector& vecpot() const { return vecpot_; }
  const double& nvecpot() const { return nvecpot_; }
  const string rtvpeq() const { return rtvpeq_; }
};
#endif
