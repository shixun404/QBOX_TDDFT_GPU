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
// RTCNWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RTCNWAVEFUNCTIONSTEPPER_H
#define RTCNWAVEFUNCTIONSTEPPER_H

#include "EnergyFunctional.h"
#include "RTVectorPotential.h"
#include "ChargeDensity.h"
#include "Wavefunction.h"
#include "WavefunctionStepper.h"
#include "SampleStepper.h"
#include "UserInterface.h"
#include "FourierTransform.h"
#include "Basis.h"

typedef std::map<std::string,Timer> TimerMap;
class Wavefunction;

class RTCNWavefunctionStepper : public WavefunctionStepper
{
  private:
  double rtdt_;
  int rtitscf_;
  int rtite_;
  int iter_;

  protected:

  EnergyFunctional& ef_;
  ChargeDensity& cd_;
  Sample& s_;

  public:

  void update(Wavefunction& dwf);
  void get_iter(int iter);
  bool conv_check(ChargeDensity& cd1, ChargeDensity& cd2, double& tsum);

  RTCNWavefunctionStepper(double rtdt, TimerMap& tmap, EnergyFunctional& ef, ChargeDensity& cd, Sample& s, int rtitscf, int rtite);
  ~RTCNWavefunctionStepper() {};
};
#endif
