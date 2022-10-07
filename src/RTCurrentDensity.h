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
// RTCurrendDensity.h
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RTCURRENTDENSITY_H
#define RTCURRENTDENSITY_H

#include "Wavefunction.h"
#include "EnergyFunctional.h"
#include "RTPosition.h"
#include "ChargeDensity.h"
#include "Sample.h"
#include <valarray>

class RTCurrentDensity
{
  private:
  RTPosition rtp_;
  AtomSet atoms_;
  std::vector<std::vector<std::complex<double> > > obs_; // obs_[ispin][dir]
  std::vector<std::vector<std::complex<double> > > pobs_; // obs_[ispin][dir]
  std::vector<std::vector<std::complex<double> > > vobs_; // obs_[ispin][dir]
  std::vector<std::vector<std::complex<double> > > hobs_; // obs_[ispin][dir]

  protected:
  Sample& s_;

  public:

  complex<double> obs(int i, int j) const { return obs_[i][j]; }
  complex<double> pobs(int i, int j) const { return pobs_[i][j]; }
  complex<double> vobs(int i, int j) const { return vobs_[i][j]; }
  complex<double> hobs(int i, int j) const { return hobs_[i][j]; }
  void compute_current(Wavefunction& wf, EnergyFunctional& ef);
  void compute_momentum(Wavefunction& rwf, Wavefunction& pwf, EnergyFunctional& ef, int ipol);
  void compute_vnlr(Wavefunction& rwf, Wavefunction& vwf, EnergyFunctional& ef, int ipol);
  void compute_hr(Wavefunction& rwf, Wavefunction& vwf, EnergyFunctional& ef, int ipol);

  RTCurrentDensity(Sample& s, RTPosition& rtp, AtomSet& atoms);
  ~RTCurrentDensity();
};
#endif
