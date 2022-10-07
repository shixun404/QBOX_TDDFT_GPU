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
// RTVectorPotential.cpp
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <valarray>
#include <complex>
#include "Sample.h"
#include "D3vector.h"
#include "UnitCell.h"
#include "RTVectorPotential.h"

////////////////////////////////////////////////////////////////////////////////
RTVectorPotential::RTVectorPotential(Sample& s):
  calind_(s.rtctrl.rt_vp_ind), rtvpeq_(s.rtctrl.rt_vp_eq), rtvpamp_(s.rtctrl.rt_vp_amp), rtvpfreq_(s.rtctrl.rt_vp_freq), rtvpsigma_(s.rtctrl.rt_vp_sigma), rtvpgauss_(s.rtctrl.rt_vp_gauss), rtvplength_(s.rtctrl.rt_vp_length), rtvpdelay_(s.rtctrl.rt_vp_delay)
{
  induced_ = D3vector(0.0, 0.0, 0.0);
  if ( calind_ == "T" )
    vecpot_ += induced_;
  nvecpot_ = norm2(vecpot_);
  velocity_ = D3vector(0.0, 0.0, 0.0);
  accel_ = D3vector(0.0, 0.0, 0.0);
}

////////////////////////////////////////////////////////////////////////////////
RTVectorPotential::~RTVectorPotential(void)
{
}

////////////////////////////////////////////////////////////////////////////////
void RTVectorPotential::calculate_acceleration(const double& dt, const D3vector& total_current, const UnitCell& cell)
{
  velocity_ += 0.5*dt*accel_;
  if( calind_ == "T" )
  {
    accel_ = -4.0*M_PI*total_current / cell.volume();
  }
  else
  {
    accel_ = D3vector(0.0, 0.0, 0.0);
  }
  velocity_ += 0.5*dt*accel_;
}

////////////////////////////////////////////////////////////////////////////////
void RTVectorPotential::vp_propagate(int step, const double& dt)
{
  induced_ = dt*velocity_ + 0.5*dt*dt*accel_;
  double autoas = 24.18884254;
  double asdt = step * dt * autoas;
  double sol = 137.03599911;

  if( asdt <= rtvpdelay_ )
  {
    vecpot_ = D3vector(0.0, 0.0, 0.0);
    nvecpot_ = 0.0;
  }
  else
  {
    if( length(rtvpamp_) > 0.0 && rtvpeq_ == "AC" )
    {
      vecpot_ = sin(rtvpfreq_ * (asdt - rtvpdelay_) * 2.0 * M_PI) * rtvpamp_ / sol;
    }

    if( length(rtvpamp_) > 0.0 && rtvpeq_ == "GAUSSIAN" )
    {
      if ( asdt <= rtvpgauss_ + rtvpdelay_ )
      {
        vecpot_ = sin(rtvpfreq_ * (asdt - rtvpdelay_) * 2.0 * M_PI) * rtvpamp_ * exp(-(asdt - rtvpdelay_ - rtvpgauss_) * (asdt -rtvpdelay_ - rtvpgauss_) / (2.0 * rtvpsigma_ * rtvpsigma_)) / sol;
      }
      else if ( asdt >= rtvpgauss_ + rtvpdelay_ && asdt <= rtvpgauss_ + rtvplength_ + rtvpdelay_)
      {
        vecpot_ = sin(rtvpfreq_ * (asdt - rtvpdelay_) * 2.0 * M_PI) * rtvpamp_ / sol;
      }
      else
      {
        vecpot_ = sin(rtvpfreq_ * (asdt - rtvpdelay_) * 2.0 * M_PI) * rtvpamp_ * exp(-(asdt - rtvpdelay_ - rtvpgauss_ - rtvplength_) * (asdt -rtvpdelay_ - rtvpgauss_ - rtvplength_) / (2.0 * rtvpsigma_ * rtvpsigma_)) / sol;
      }
    }

    if( length(rtvpamp_) > 0.0 && rtvpeq_ == "CIRPOLAR" )
    {
      vecpot_[0] = cos(rtvpfreq_ * (asdt - rtvpdelay_) * 2.0 * M_PI) * rtvpamp_[2] / sol;
      vecpot_[1] = sin(rtvpfreq_ * (asdt - rtvpdelay_) * 2.0 * M_PI) * rtvpamp_[2] / sol;
      vecpot_[2] = 0.0; 
    }

    if( length(rtvpamp_) > 0.0 && rtvpeq_ == "PULSE" )
    {
      vecpot_ = sin(rtvpfreq_ * (asdt - rtvpdelay_) * 2.0 * M_PI) * rtvpamp_ * exp(-(asdt - rtvpdelay_ - rtvpgauss_) * (asdt -rtvpdelay_ - rtvpgauss_) / (2.0 * rtvpsigma_ * rtvpsigma_)) / sol;
    }

    if( length(rtvpamp_) > 0.0 && rtvpeq_ == "CIRPULSE" )
    {
      vecpot_[0] = rtvpamp_[2] * exp(-(asdt - rtvpdelay_ - rtvpgauss_) * (asdt - rtvpdelay_ - rtvpgauss_) / (2.0 * rtvpsigma_ * rtvpsigma_)) * cos(rtvpfreq_ * (asdt - rtvpdelay_) * 2.0 * M_PI) / sol;
      vecpot_[1] = rtvpamp_[2] * exp(-(asdt - rtvpdelay_ - rtvpgauss_) * (asdt - rtvpdelay_ - rtvpgauss_) / (2.0 * rtvpsigma_ * rtvpsigma_)) * sin(rtvpfreq_ * (asdt - rtvpdelay_) * 2.0 * M_PI) / sol;
      vecpot_[2] = 0.0;
    }
  }

  if( rtvpeq_ == "DELTA" )
  {
    vecpot_[0] = 0.0;
    vecpot_[1] = 0.0;
    vecpot_[2] = 0.0;
  }

  if ( calind_ == "T" )
    vecpot_ += induced_;
  nvecpot_ = norm2(vecpot_);
}
