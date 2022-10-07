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
// RTProjection.h
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RTPROJECTION_H
#define RTPROJECTION_H

#include "Wavefunction.h"
#include <valarray>
#include <iostream>

using namespace std;

class RTProjection
{
  private:
  Wavefunction pwf_;
  string rtprojtype_;
  int proj_bd_s_;
  int proj_bd_e_;
  int proj_kp_a_;
  int proj_kp_b_;
  int nbnd;

  std::vector<std::vector<int> > msize; // msize[ispin][ikp]
  std::vector<std::vector<std::vector<std::vector<std::complex<double> > > > > coll_init_; // coll_init_[ispin][ikp][ist][ig]
  std::vector<std::vector<std::vector<std::vector<std::complex<double> > > > > coll_band_; // coll_band_[ispin][ikp][ist][ig]
  std::vector<std::vector<std::vector<double> > > proj_band_; // proj_band_[ispin][ibnd][jbnd]
  std::vector<std::vector<std::vector<std::vector<double> > > > proj_; // proj_[ispin][ikp][ibnd][jbnd]

  public:
  double proj_band(int is, int ib, int jb) const { return proj_band_[is][ib][jb]; }
  double proj(int is, int ik, int ib, int jb) const { return proj_[is][ik][ib][jb]; }
  void init_collection();
  void init_collection2();
  void start_collection(Wavefunction& iwf);
  void band_collection(Wavefunction& cwf);
  void compute_projection(Wavefunction& ppwf);
  void set_collection();

  RTProjection(Wavefunction& wf, string rtprojtype, int proj_bd_s, int proj_bd_e, int proj_kp_a, int proj_kp_b);
  ~RTProjection();
};
#endif

