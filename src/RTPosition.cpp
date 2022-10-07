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
// RTPosition.cpp
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#include "RTPosition.h"
#include <vector>
#include <valarray>
#include <complex>
#include "D3vector.h"
#include "MPIdata.h"
#include "FourierTransform.h"
#include "Basis.h"
#include "Wavefunction.h"
#include "AtomSet.h"
#include "Atom.h"
#include "Species.h"
#include "SlaterDet.h"

using namespace std;
////////////////////////////////////////////////////////////////////////////////
RTPosition::RTPosition(Wavefunction& wf, AtomSet& atoms):
  wf_(wf)
{
  // Calculate the Center of Charge
  double zval = 0;
  double coc[3];
  coc[0] = 0.0;
  coc[1] = 0.0;
  coc[2] = 0.0;
  for ( int is = 0; is < atoms.nsp(); is++ )
  {
    Species* ps = atoms.species_list[is];
    for ( int ia = 0; ia < atoms.na(is); ia++ )
    {
      Atom* pa = atoms.atom_list[is][ia];
      coc[0] = coc[0] + pa->position().x * ps->ztot();
      coc[1] = coc[1] + pa->position().y * ps->ztot();
      coc[2] = coc[2] + pa->position().z * ps->ztot();
      zval = zval + ps->ztot();
    }
  }
  coc[0] = coc[0]/(double)zval;
  coc[1] = coc[1]/(double)zval;
  coc[2] = coc[2]/(double)zval;
  //
  // Dense grid for rhor
  Basis* vbasis_rho_ = new Basis(MPIdata::g_comm(), D3vector(0,0,0));
  vbasis_rho_->set_real(true);
  vbasis_rho_->resize(wf.cell(), wf.refcell(), 4.0*wf.ecut());
  const Basis& vbr = *vbasis_rho_;
  int np0v = vbr.np(0)+2;
  int np1v = vbr.np(1)+2;
  int np2v = vbr.np(2)+2;
  while (!vbr.factorizable(np0v)) np0v += 2;
  while (!vbr.factorizable(np1v)) np1v += 2;
  while (!vbr.factorizable(np2v)) np2v += 2;
  ft_rho_ = new FourierTransform(vbr, np0v, np1v, np2v);
  if ( MPIdata::onpe0() )
    cout << "RTPosition DENST GRID : " << np0v << "  " << np1v << "  " << np2v << endl;

  int zidx = 0;
  size_rho_ = 0;
  double x_position[3];
  for ( int iz = 0; iz < ft_rho_->np2_loc(MPIdata::igb()); ++iz )
  {
    zidx = ft_rho_->np2_first(MPIdata::igb()) + iz;
    for ( int iy = 0; iy < ft_rho_->np1(); ++iy )
    {
      for ( int ix = 0; ix < ft_rho_->np0(); ++ix )
      {
        x_position[0] = ((double)ix/(double)ft_rho_->np0())*(wf.cell().a(0).x) + ((double)iy/(double)ft_rho_->np1())*(wf.cell().a(1).x) + ((double)zidx/(double)ft_rho_->np2())*(wf.cell().a(2).x);
        x_position[1] = ((double)ix/(double)ft_rho_->np0())*(wf.cell().a(0).y) + ((double)iy/(double)ft_rho_->np1())*(wf.cell().a(1).y) + ((double)zidx/(double)ft_rho_->np2())*(wf.cell().a(2).y);
        x_position[2] = ((double)ix/(double)ft_rho_->np0())*(wf.cell().a(0).z) + ((double)iy/(double)ft_rho_->np1())*(wf.cell().a(1).z) + ((double)zidx/(double)ft_rho_->np2())*(wf.cell().a(2).z);
        x_position[0] = x_position[0] - coc[0];
        x_position[1] = x_position[1] - coc[1];
        x_position[2] = x_position[2] - coc[2];
        rt_position_rho_.push_back(D3vector(x_position[0], x_position[1], x_position[2]));
        size_rho_++;
      }
    }
  }
  //
  // Dense grid for wavefunction
  if ( wf.nsp_loc() > 0 )
  {
    for ( int ikp_loc = 0; ikp_loc < wf.nkp_loc(); ++ikp_loc )
    {
      const Basis& vbw = wf.sd(0, ikp_loc)->basis();
      ft_dwfc_.push_back(new FourierTransform(vbw, np0v, np1v, np2v));
    }
  }

  zidx = 0;
  x_position[3];
  rt_position_dwfc_.resize(wf.nkp_loc());
  size_dwfc_.resize(wf.nkp_loc());
  for ( int ikp_loc = 0; ikp_loc < wf.nkp_loc(); ++ikp_loc )
  {
    size_dwfc_[ikp_loc] = 0;
    for ( int iz = 0; iz < ft_dwfc_[ikp_loc]->np2_loc(MPIdata::igb()); ++iz )
    {
      zidx = ft_dwfc_[ikp_loc]->np2_first(MPIdata::igb()) + iz;
      for ( int iy = 0; iy < ft_dwfc_[ikp_loc]->np1(); ++iy )
      {
        for ( int ix = 0; ix < ft_dwfc_[ikp_loc]->np0(); ++ix )
        {
          x_position[0] = ((double)ix/(double)ft_dwfc_[ikp_loc]->np0())*(wf.cell().a(0).x) + ((double)iy/(double)ft_dwfc_[ikp_loc]->np1())*(wf.cell().a(1).x) + ((double)zidx/(double)ft_dwfc_[ikp_loc]->np2())*(wf.cell().a(2).x);
          x_position[1] = ((double)ix/(double)ft_dwfc_[ikp_loc]->np0())*(wf.cell().a(0).y) + ((double)iy/(double)ft_dwfc_[ikp_loc]->np1())*(wf.cell().a(1).y) + ((double)zidx/(double)ft_dwfc_[ikp_loc]->np2())*(wf.cell().a(2).y);
          x_position[2] = ((double)ix/(double)ft_dwfc_[ikp_loc]->np0())*(wf.cell().a(0).z) + ((double)iy/(double)ft_dwfc_[ikp_loc]->np1())*(wf.cell().a(1).z) + ((double)zidx/(double)ft_dwfc_[ikp_loc]->np2())*(wf.cell().a(2).z);
          x_position[0] = x_position[0] - coc[0];
          x_position[1] = x_position[1] - coc[1];
          x_position[2] = x_position[2] - coc[2];
          rt_position_dwfc_[ikp_loc].push_back(D3vector(x_position[0], x_position[1], x_position[2]));
          size_dwfc_[ikp_loc]++;
        }
      }
    }
  }
  //
  // Smooth grid for wavefunction
  if ( wf.nsp_loc() > 0 )
  {
    for ( int ikp_loc = 0; ikp_loc < wf.nkp_loc(); ++ikp_loc )
    {
      const Basis& vbw = wf.sd(0, ikp_loc)->basis();
      np0v = vbw.np(0);
      np1v = vbw.np(1);
      np2v = vbw.np(2);
      ft_swfc_.push_back(new FourierTransform(vbw, np0v, np1v, np2v));
    }
  }

  if ( MPIdata::onpe0() )
    cout << "RTPosition SOFT GRID : " << np0v << "  " << np1v << "  " << np2v << endl;

  zidx = 0;
  x_position[3];
  rt_position_swfc_.resize(wf.nkp_loc());
  size_swfc_.resize(wf.nkp_loc());
  for ( int ikp_loc = 0; ikp_loc < wf.nkp_loc(); ++ikp_loc )
  {
    size_swfc_[ikp_loc] = 0;
    for ( int iz = 0; iz < ft_swfc_[ikp_loc]->np2_loc(MPIdata::igb()); ++iz )
    {
      zidx = ft_swfc_[ikp_loc]->np2_first(MPIdata::igb()) + iz;
      for ( int iy = 0; iy < ft_swfc_[ikp_loc]->np1(); ++iy )
      {
        for ( int ix = 0; ix < ft_swfc_[ikp_loc]->np0(); ++ix )
        {
          x_position[0] = ((double)ix/(double)ft_swfc_[ikp_loc]->np0())*(wf.cell().a(0).x) + ((double)iy/(double)ft_swfc_[ikp_loc]->np1())*(wf.cell().a(1).x) + ((double)zidx/(double)ft_swfc_[ikp_loc]->np2())*(wf.cell().a(2).x);
          x_position[1] = ((double)ix/(double)ft_swfc_[ikp_loc]->np0())*(wf.cell().a(0).y) + ((double)iy/(double)ft_swfc_[ikp_loc]->np1())*(wf.cell().a(1).y) + ((double)zidx/(double)ft_swfc_[ikp_loc]->np2())*(wf.cell().a(2).y);
          x_position[2] = ((double)ix/(double)ft_swfc_[ikp_loc]->np0())*(wf.cell().a(0).z) + ((double)iy/(double)ft_swfc_[ikp_loc]->np1())*(wf.cell().a(1).z) + ((double)zidx/(double)ft_swfc_[ikp_loc]->np2())*(wf.cell().a(2).z);
          x_position[0] = x_position[0] - coc[0];
          x_position[1] = x_position[1] - coc[1];
          x_position[2] = x_position[2] - coc[2];
          rt_position_swfc_[ikp_loc].push_back(D3vector(x_position[0], x_position[1], x_position[2]));
          size_swfc_[ikp_loc]++;
        }
      }
    }
  }
  //
  delete vbasis_rho_;
}

////////////////////////////////////////////////////////////////////////////////
RTPosition::~RTPosition(void)
{
  delete ft_rho_;
  for ( int ikp_loc = 0; ikp_loc < wf_.nkp_loc(); ++ikp_loc )
  {
    delete ft_swfc_[ikp_loc];
    delete ft_dwfc_[ikp_loc];
  }
}
