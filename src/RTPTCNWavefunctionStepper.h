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
// RTPTCNWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RTPTCNWAVEFUNCTIONSTEPPER_H
#define RTPTCNWAVEFUNCTIONSTEPPER_H

#include "EnergyFunctional.h"
#include "ChargeDensity.h"
#include "Sample.h"
#include "WavefunctionStepper.h"
#include "AndersonMixer.h"

class AndersonMixer;

class RTPTCNWavefunctionStepper : public WavefunctionStepper
{
  private:
  const double rtdt_;
  int rtitscf_;
  int iter_;
  std::vector<std::vector<std::vector<AndersonMixer*> > > mixer_;    // mixer[ispin][ikp][ist]

  protected:
  EnergyFunctional& ef_;
  ChargeDensity& cd_;
  Sample& s_;

 public:

  void update(Wavefunction& dwf);
  void get_iter(int& iter);
  bool conv_check(ChargeDensity& cd1, ChargeDensity& cd2, double& tsum);

  RTPTCNWavefunctionStepper(double rtdt, TimerMap& tmap, EnergyFunctional& ef, ChargeDensity& cd, Sample& s, int rtitscf);
  ~RTPTCNWavefunctionStepper();
};
#endif


