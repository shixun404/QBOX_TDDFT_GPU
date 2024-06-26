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
// XCOperator.cpp
//
////////////////////////////////////////////////////////////////////////////////
#include "XCOperator.h"
#include "Sample.h"
#include "ChargeDensity.h"
#include "XCPotential.h"
#include "ExchangeOperator.h"
#include "HSEFunctional.h"
#include "RSHFunctional.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
XCOperator::XCOperator(Sample& s, const ChargeDensity& cd) :cd_(cd), s_(s)
{
  // set initial values
  xcp_ = 0;
  xop_ = 0;
  exc_ = 0.0 ;
  dxc_ = 0.0 ;
  hasHF_ = false;
  hasMeta_ = false;

  sigma_exc_.resize(6);

  string functional_name = s.ctrl.xc;

  // check the name of the functional
  if ( ( functional_name ==  "LDA" ) ||
       ( functional_name ==  "VWN" ) ||
       ( functional_name ==  "PBE" ) ||
       ( functional_name == "BLYP" ) )
  {
    // create only an xc potential
    xcp_ = new XCPotential(cd, functional_name, s_);
    hasPotential_ = true;
    hasGGA_ = xcp_->isGGA();
    hasHF_ = false;
  }
  else if ( functional_name == "HF" )
  {
    // create exchange operator with mixing coeff=1
    xop_ = new ExchangeOperator(s, 1.0, 1.0, 0.0);
    hasPotential_ = false;
    hasGGA_ = false;
    hasHF_ = true;
  }
  else if ( functional_name == "PBE0" )
  {
    // create an exchange potential
    xcp_ = new XCPotential(cd, functional_name, s_);

    // create the exchange operator with mixing coeff=0.25
    xop_ = new ExchangeOperator(s, s.ctrl.alpha_PBE0, s.ctrl.alpha_PBE0, 0.0);
    hasPotential_ = true;
    hasGGA_ = xcp_->isGGA();
    hasHF_ = true;
  }
  else if ( functional_name == "HSE" )
  {
    // create an exchange potential
    xcp_ = new XCPotential(cd, functional_name, s_);

    // create the exchange operator with mixing coeff=0.25
    xop_ = new ExchangeOperator(s, 0.0, 0.25, 0.11);
    hasPotential_ = true;
    hasGGA_ = xcp_->isGGA();
    hasHF_ = true;
  }
  else if ( functional_name == "RSH" )
  {
    // create an exchange potential
    xcp_ = new XCPotential(cd, functional_name, s_);

    // create the exchange operator with mixing coeff=beta_RSH
    xop_ = new ExchangeOperator(s, s.ctrl.alpha_RSH, s.ctrl.beta_RSH,
      s.ctrl.mu_RSH);
    hasPotential_ = true;
    hasGGA_ = xcp_->isGGA();
    hasHF_ = true;
  }
  else if ( functional_name == "B3LYP" )
  {
    // create an exchange potential
    xcp_ = new XCPotential(cd, functional_name, s_);

    // create the exchange operator with mixing coeff=0.20
    xop_ = new ExchangeOperator(s, 0.20, 0.20, 0.0);
    hasPotential_ = true;
    hasGGA_ = xcp_->isGGA();
    hasHF_ = true;
  }
  else if ( functional_name == "BHandHLYP" )
  {
    // create an exchange potential
    xcp_ = new XCPotential(cd, functional_name, s_);

    // create the exchange operator with mixing coeff=0.50
    xop_ = new ExchangeOperator(s, 0.50, 0.50, 0.0);
    hasPotential_ = true;
    hasGGA_ = xcp_->isGGA();
    hasHF_ = true;
  }
  else if ( functional_name == "SCAN" )
  {
    // create an exchange potential
    xcp_ = new XCPotential(cd, functional_name, s_);
    hasPotential_ = true;
    hasGGA_ = xcp_->isGGA();
    hasMeta_ = true;
  }
  else
  {
    throw XCOperatorException(
      "unknown functional name during exchange operator construction");
  }
}

XCOperator::~XCOperator()
{
  delete xcp_;
  delete xop_;
}

////////////////////////////////////////////////////////////////////////////////
void XCOperator::update(std::vector<std::vector<double> >& vr, bool compute_stress)
{
  // update xc potential and self-energy
  // used whenever the charge density and/or wave functions have changed
  // compute vxc potential and energy
  if ( hasPotential_ )
  {
    // update LDA/GGA/MetaGGA xc potential
    xcp_->update( vr );

    // LDA/GGA exchange energy
    exc_ = xcp_->exc();
    dxc_ = xcp_->dxc();

    if ( compute_stress )
      xcp_->compute_stress(sigma_exc_);
  }
  else
  {
    exc_ = 0.0;
    dxc_ = 0.0;
    sigma_exc_ = 0.0;
  }

  if ( hasHF() )
  {
    double ex_hf = xop_->update_operator(compute_stress);
    exc_ += ex_hf;
    dxc_ -= ex_hf;
    if ( compute_stress )
      xop_->add_stress(sigma_exc_);
  }
}

////////////////////////////////////////////////////////////////////////////////
void XCOperator::apply_self_energy(Wavefunction &dwf)
{
  if ( hasHF() )
    xop_->apply_operator(dwf);
  if ( hasMeta() )
    xcp_->apply_meta_operator(dwf);
}

////////////////////////////////////////////////////////////////////////////////
void XCOperator::rttd_apply_self_energy(Wavefunction& wf, Wavefunction &dwf)
{
  if ( hasHF() )
    xop_->rttd_apply_operator(wf, dwf);
//if ( hasMeta() )
//  xcp_->apply_meta_operator(dwf);
}

////////////////////////////////////////////////////////////////////////////////
void XCOperator::rt_apply_self_energy(int isp_loc, int ikp_loc, Wavefunction& wf, Wavefunction &dwf)
{
  if ( hasHF() )
    xop_->rt_apply_operator(isp_loc, ikp_loc, wf, dwf);
//if ( hasMeta() )
//  xcp_->apply_meta_operator(dwf);
}

////////////////////////////////////////////////////////////////////////////////
void XCOperator::compute_stress(std::valarray<double>& sigma)
{
  sigma = sigma_exc_;
}

////////////////////////////////////////////////////////////////////////////////
void XCOperator::cell_moved(void)
{
  if ( hasHF() )
    xop_->cell_moved();
}
