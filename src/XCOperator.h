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
// XCOperator.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef XCOPERATOR_H
#define XCOPERATOR_H

#include <valarray>
#include <vector>
#include <string>

class Sample;
class Wavefunction;
class ChargeDensity;
class XCPotential;
class ExchangeOperator;
class XCOperator
{
  private:

  XCPotential* xcp_;
  ExchangeOperator* xop_;

  const Sample& s_;
  const ChargeDensity& cd_;
  double exc_; // XC energy: includes local and HF terms
  double dxc_;

  std::valarray<double> sigma_exc_;

  bool hasPotential_;
  bool hasGGA_;
  bool hasHF_;
  bool hasMeta_;

  public:

  // constructor
  XCOperator( Sample& s, const ChargeDensity& cd);

  // destructor
  ~XCOperator();

  // return pointer to the exchange potential class
  const XCPotential* xcp() { return xcp_; }

  // return pointer to the ExchangeOperator class
  ExchangeOperator* xop() { return xop_; }

  bool hasGGA(void) { return hasGGA_; };
  bool hasHF(void) { return hasHF_; };
  bool hasMeta(void) { return hasMeta_; };

  void update(std::vector<std::vector<double> >& vr, bool compute_stress);
  void apply_self_energy(Wavefunction &dwf);
  void rttd_apply_self_energy(Wavefunction& wf, Wavefunction &dwf);
  void rt_apply_self_energy(int isp_loc, int ikp_loc, Wavefunction& wf, Wavefunction &dwf);
  void compute_stress(std::valarray<double>& sigma);
  void cell_moved(void);
  double exc(void) { return exc_ ; };
  double dxc(void) { return dxc_ ; };
};

class XCOperatorException
{
  public:
  std::string msg;
  XCOperatorException(std::string s) : msg(s) {}
};

#endif
