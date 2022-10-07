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
// BMDIonicStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
//
// IonicStepper is used in the following way
//
// input: r0,v0
//
// compute energy e0(r0) and forces f0(r0)
// for ( k=0; k<n; k++ )
// {
//   // r0,v0,e0,f0 known
//   stepper->compute_r(e0,f0)
//   {
//     computes rp using r0, v0 and f0
//     restores constraints on rp using rp and r0
//     updates rm<-r0, r0<-rp and update atomset positions
//   }
//
//   compute f0(r0)
//
//   stepper->compute_v(e0,f0)
//   {
//     computes v0 using r0,rm,f0
//     restores constraints on v0 using r0, v0
//     modifies velocities using the thermostat (rescaling)
//     updates atomset velocities
//   }
// }
// r0,v0,f0 consistent at this point
//

#ifndef BMDIONICSTEPPER_H
#define BMDIONICSTEPPER_H

#include "IonicStepper.h"

class BMDIonicStepper : public IonicStepper
{
  private:

  const double bmd_fac_;
  double e0_,em_;
  double ekin_;
  std::vector<std::vector< double> > fm_;
  void compute_ekin(void);

  public:

  BMDIonicStepper(Sample& s) : IonicStepper(s), bmd_fac_(0.025)
  {
    e0_ = 0.0;
    em_ = 0.0;
    ekin_ = 0.0;
    atoms_.get_positions(r0_);
    atoms_.get_velocities(v0_);
    fm_.resize(r0_.size());
    for ( int is = 0; is < fm_.size(); is++ )
      fm_[is].resize(r0_[is].size());
    compute_ekin();
  }

  void compute_r(double e0, const std::vector<std::vector< double> >& f0);
  void compute_v(double e0, const std::vector<std::vector< double> >& f0);
  double ekin(void) const { return ekin_; }
  double temp(void) const
  {
    const double boltz = 1.0 / ( 11605.0 * 2.0 * 13.6058 );
    if ( ndofs_ > 0 )
      return 2.0 * ( ekin_ / boltz ) / ndofs_;
    else
      return 0.0;
  }
};

#endif
