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
// RTPosition.h
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RTPOSITION_H
#define RTPOSITION_H

#include <vector>
#include <valarray>
#include <complex>
#include "D3vector.h"
#include "FourierTransform.h"
#include "Basis.h"
#include "Wavefunction.h"
#include "AtomSet.h"

class RTPosition
{
  private:
  Wavefunction wf_;
  FourierTransform* ft_rho_;
  std::vector<FourierTransform*> ft_swfc_; // ft_wfc_[ikp_loc]
  std::vector<FourierTransform*> ft_dwfc_; // ft_wfc_[ikp_loc]
  std::vector<D3vector> rt_position_rho_; //rt_position_[np012loc][3]
  std::vector<std::vector<D3vector> > rt_position_swfc_; //rt_position_[ikp_loc][np012loc][3]
  std::vector<std::vector<D3vector> > rt_position_dwfc_; //rt_position_[ikp_loc][np012loc][3]
  int size_rho_;
  std::vector<int> size_swfc_; // size_wfc_[ikp_loc]
  std::vector<int> size_dwfc_; // size_wfc_[ikp_loc]

  public:
  RTPosition(Wavefunction& wf, AtomSet& atoms);
  ~RTPosition();

  D3vector rt_position_rho(int ir) { return rt_position_rho_[ir]; }
  D3vector rt_position_swfc(int ikp, int ir) { return rt_position_swfc_[ikp][ir]; }
  D3vector rt_position_dwfc(int ikp, int ir) { return rt_position_dwfc_[ikp][ir]; }
  int size_rho(void) const { return size_rho_; }
  int size_swfc(int ikp) const { return size_swfc_[ikp]; }
  int size_dwfc(int ikp) const { return size_dwfc_[ikp]; }
  FourierTransform* ft_rho(void) const { return ft_rho_; }
  FourierTransform* ft_swfc(int ikp) const { return ft_swfc_[ikp]; }
  FourierTransform* ft_dwfc(int ikp) const { return ft_dwfc_[ikp]; }
};
#endif
