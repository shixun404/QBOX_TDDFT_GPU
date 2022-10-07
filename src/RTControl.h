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
// RTControl.h:
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RTCONTROL_H
#define RTCONTROL_H

#include <string>
#include <vector>
#include <map>
#include "D3vector.h"

struct RTControl
{
  // rtcontrol variables
  std::string rt_propagator; // CN | ST2 | ST4
  double rt_anderson_coeff;
  int rt_anderson_dim;
  std::string rt_md; // ON | OFF
  double rt_scf_tol;
  double rt_rho_tol;
  int rt_as;
  // rt Vector Potential Variables
  std::string rt_vp; // ON | OFF
  std::string rt_vp_ind; // T | F
  std::string rt_vp_eq; // Formula for Vector Potential : AC | GAUSSIAN | PULSE | DELTA | CIRPULSE
  D3vector rt_vp_amp; // Amplitute of Vector Potential
  double rt_vp_freq; // Frequency of Vector Potential unit : PHz (10^15 Hz)
  double rt_vp_sigma; // Sigma for GAUSSIAN form Vector Potential
  double rt_vp_gauss; // t_exp for GAUSSIAN form Vector Potential
  double rt_vp_length; // Length for AC | GAUSSIAN  Vector Potential
  double rt_vp_delay; // Delay for Vector Potential
  // delta-kick
  std::string rt_delta_kick; // ON | OFF
  D3vector rt_delta_kick_amp;
  // projection
  std::string rt_proj; // ON | OFF
  std::string rt_proj_type; // BAND | SELKP | ALLKP
  int rt_proj_bd[2];
  int rt_proj_kp[2];

  double rt_dt; // Time step for RT-TDDFT Calculation
};
#endif
