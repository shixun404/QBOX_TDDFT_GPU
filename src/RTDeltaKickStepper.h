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
// RTDeltaKickStepper.h
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RTDELTAKICKSTEPPER_H
#define RTDELTAKICKSTEPPER_H

#include "RTPosition.h"
#include "Sample.h"
#include "D3vector.h"
#include "ChargeDensity.h"
#include "Wavefunction.h"

class RTDeltaKickStepper
{
  private:
  RTPosition rtp_;
  D3vector rt_dipole_;

  protected:
  
  Sample& s_;

  public:

  RTPosition rtp(void) const { return rtp_; }
  D3vector rt_dipole(void) const { return rt_dipole_; }

  void update_dk(D3vector deltakick_amp);
  void compute_dip(ChargeDensity& cd);

  RTDeltaKickStepper(Sample& s, RTPosition& rtp);
  ~RTDeltaKickStepper();
};
#endif
