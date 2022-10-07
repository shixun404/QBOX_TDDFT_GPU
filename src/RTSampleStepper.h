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
// RTSampleStepper.h
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RTSAMPLESTEPPER_H
#define RTSAMPLESTEPPER_H

#include "SampleStepper.h"
#include "EnergyFunctional.h"
#include "Sample.h"
#include "ChargeDensity.h"
#include "Wavefunction.h"

class WavefunctionStepper;
class IonicStepper;

class RTSampleStepper : public SampleStepper
{
  private:

  Wavefunction wf_; // Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
  Wavefunction wf_prev_; // Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
  Wavefunction wf_next_; // Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
  Wavefunction wf_init_; // Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
  Wavefunction dwf;
  int rtitscf_;
  int rtite_;
  int rtas_;
  ChargeDensity cd_;
  EnergyFunctional ef_;

  WavefunctionStepper* wf_stepper;
  IonicStepper* ionic_stepper;

  bool update_density_first_;
  bool update_vh_, update_vxc_;

  // Do not allow construction of RTCNSampleStepper unrelated to a Sample
  RTSampleStepper(void);

  public:

  mutable TimerMap tmap;

  void step(int rtiter);
  void set_update_vh(bool update_vh) { update_vh_ = update_vh; }
  void set_update_vxc(bool update_vxc) { update_vxc_ = update_vxc; }
  void set_update_density_first(bool update_density_first)
    { update_density_first_ = update_density_first; }

  EnergyFunctional& ef(void) { return ef_; }
  ChargeDensity& cd(void) { return cd_; }

  RTSampleStepper(Sample& s, int rtitscf, int rtite, int rtas);
  ~RTSampleStepper();
};
#endif

