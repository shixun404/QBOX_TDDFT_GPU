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
// NonLocalPotential.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NONLOCALPOTENTIAL_H
#define NONLOCALPOTENTIAL_H

#include "AtomSet.h"
#include "Basis.h"
#include "SlaterDet.h"
#include "Context.h"
#include "Matrix.h"
#include "RTVectorPotential.h"
#include "RTPosition.h"
#include "D3vector.h"

class NonLocalPotential
{
  private:

  const Context& ctxt_;
  const AtomSet& atoms_;
  const SlaterDet& sd_;
  const Basis& basis_;
  const bool rtvp_;

  int nsp;   // number of species
  int nspnl; // number of non-local species

  std::vector<int>                  nop;      // nop[is]
  std::vector<int>                  lloc;     // lloc[is]
  std::vector<int>                  na;       // na[is]
  std::vector<int>                  npr;      // npr[is]
  std::vector<int>                  nprna;    // nprna[is]
  std::vector<std::vector<int> >    lproj;    // lproj[is][ipr]
  std::vector<std::vector<double> > wt;       // wt[is][ipr]
  std::vector<std::vector<double> > twnl;     // twnl[is][npr*ngwl]
  std::vector<std::vector<double> > dtwnl;    // dtwnl[is][6*npr*ngwl],ij=0,..,5

  std::vector<int>             nquad;    // nquad[is]
  // rquad[is][iquad], iquad = 0, nquad[is]-1
  std::vector<std::vector<double> > rquad;
  // wquad[is][iquad], iquad = 0, nquad[is]-1
  std::vector<std::vector<double> > wquad;

  mutable TimerMap tmap;
  void init(void);

  RTVectorPotential* vp_;

  public:

  NonLocalPotential(const AtomSet& as, const SlaterDet& sd, const bool rtvp, RTVectorPotential* vp);
  ~NonLocalPotential(void);

  void update_twnl(void);

  void update_twnl_dk(D3vector deltakick_amp, RTPosition rtp, int ikp_loc);

  void nlpsi(SlaterDet& psd, SlaterDet& hpsd);

  double energy(SlaterDet& rsd, bool compute_hpsi, SlaterDet& dsd,
    bool compute_forces, std::vector<std::vector<double> >& fion,
    bool compute_stress, std::valarray<double>& sigma_enl);

  double energy(bool compute_hpsi, SlaterDet& dsd,
    bool compute_forces, std::vector<std::vector<double> >& fion,
    bool compute_stress, std::valarray<double>& sigma_enl);
};
#endif
