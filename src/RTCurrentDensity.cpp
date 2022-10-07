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
// RTCurrendDensity.cpp
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#include "RTCurrentDensity.h"
#include "Wavefunction.h"
#include "EnergyFunctional.h"
#include "RTPosition.h"
#include "NonLocalPotential.h"
#include "AtomSet.h"
#include "MPIdata.h"
#include "SlaterDet.h"
#include "FourierTransform.h"
#include "Basis.h"
#include "Sample.h"

#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
RTCurrentDensity::RTCurrentDensity(Sample& s, RTPosition& rtp, AtomSet& atoms):
rtp_(rtp), atoms_(atoms), s_(s)
{
  obs_.resize(s_.wf.nspin());
  pobs_.resize(s_.wf.nspin());
  vobs_.resize(s_.wf.nspin());
  hobs_.resize(s_.wf.nspin());
  for ( int ispin = 0; ispin < s_.wf.nspin(); ispin++ )
  {
    obs_[ispin].resize(3);
    pobs_[ispin].resize(3);
    vobs_[ispin].resize(3);
    hobs_[ispin].resize(3);
    for ( int ipol = 0; ipol < 3; ipol++ )
    {
      obs_[ispin][ipol] = complex<double>(0.0, 0.0);
      pobs_[ispin][ipol] = complex<double>(0.0, 0.0);
      vobs_[ispin][ipol] = complex<double>(0.0, 0.0);
      hobs_[ispin][ipol] = complex<double>(0.0, 0.0);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
RTCurrentDensity::~RTCurrentDensity(void)
{
}

////////////////////////////////////////////////////////////////////////////////
void RTCurrentDensity::compute_current(Wavefunction& wf, EnergyFunctional& ef)
{
  Wavefunction pwf(wf);
  Wavefunction vwf(wf);
  Wavefunction hwf(wf);

  ef.update_vhxc(false);
  for ( int ipol = 0; ipol < 3; ipol++ )
  {
    compute_momentum(wf, pwf, ef, ipol);
    compute_vnlr(wf, vwf, ef, ipol);
    compute_hr(wf, hwf, ef, ipol);
    std::vector<complex<double> > psum;
    std::vector<complex<double> > vsum;
    std::vector<complex<double> > hsum;
    psum.resize(wf.nspin());
    vsum.resize(wf.nspin());
    hsum.resize(wf.nspin());
    for ( int ispin = 0; ispin < wf.nspin(); ++ispin )
    {
      const int isp_loc = wf.isp_local(ispin);
      psum[ispin] = complex<double>(0.0, 0.0);
      vsum[ispin] = complex<double>(0.0, 0.0);
      hsum[ispin] = complex<double>(0.0, 0.0);
      if ( isp_loc >= 0 )
      {
        for ( int ikp_loc = 0; ikp_loc < wf.nkp_loc(); ++ikp_loc )
        {
          const int ikpg = wf.ikp_global(ikp_loc);
          for ( int ist_loc = 0; ist_loc < wf.sd(isp_loc, ikp_loc)->nstloc(); ++ist_loc )
          {
            const int istg = wf.sd(isp_loc, ikp_loc)->c().jglobal(ist_loc);
            const int mloc = wf.sd(isp_loc, ikp_loc)->c().mloc();
            int nmloc = ist_loc * mloc;
            double occ = wf.sd(isp_loc, ikp_loc)->occ(istg);
            int ngwloc = wf.sd(isp_loc, ikp_loc)->basis().localsize();

            for ( int ig = 0; ig < ngwloc; ig++ )
            {
              psum[ispin] += conj(wf.sd(isp_loc, ikp_loc)->c()[nmloc + ig]) * pwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] * occ * wf.weight(ikpg);
              //psum[ispin] += conj(wf.sd(isp_loc, ikp_loc)->c()[nmloc + ig]) * pwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] * wf.weight(ikpg);
              vsum[ispin] += conj(wf.sd(isp_loc, ikp_loc)->c()[nmloc + ig]) * vwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] * occ * wf.weight(ikpg);
              //vsum[ispin] += conj(wf.sd(isp_loc, ikp_loc)->c()[nmloc + ig]) * vwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] * wf.weight(ikpg);
              hsum[ispin] += conj(wf.sd(isp_loc, ikp_loc)->c()[nmloc + ig]) * hwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] * occ * wf.weight(ikpg);
              //hsum[ispin] += conj(wf.sd(isp_loc, ikp_loc)->c()[nmloc + ig]) * wf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] * occ * wf.weight(ikpg);
            }
          }
        }
      }
      MPI_Barrier(MPIdata::comm());
      pobs_[ispin][ipol] = complex<double>(0.0, 0.0);
      vobs_[ispin][ipol] = complex<double>(0.0, 0.0);
      hobs_[ispin][ipol] = complex<double>(0.0, 0.0);
      MPI_Allreduce(&psum[ispin], &pobs_[ispin][ipol], 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPIdata::comm());
      MPI_Allreduce(&vsum[ispin], &vobs_[ispin][ipol], 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPIdata::comm());
      MPI_Allreduce(&hsum[ispin], &hobs_[ispin][ipol], 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPIdata::comm());
      pobs_[ispin][ipol] += ef.vp->vecpot()[ipol]*wf.nel();
      obs_[ispin][ipol] = pobs_[ispin][ipol] + complex<double>(0.0, 1.0) * vobs_[ispin][ipol];
      hobs_[ispin][ipol] = complex<double>(0.0, 1.0) * hobs_[ispin][ipol];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void RTCurrentDensity::compute_momentum(Wavefunction& rwf, Wavefunction& pwf, EnergyFunctional& ef, int ipol)
{
  pwf = rwf;
  for ( int isp_loc = 0; isp_loc < rwf.nsp_loc(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < rwf.nkp_loc(); ++ikp_loc )
    {
      for ( int ist_loc = 0; ist_loc < rwf.sd(isp_loc, ikp_loc)->nstloc(); ++ist_loc )
      {
        int mloc = rwf.sd(isp_loc, ikp_loc)->c().mloc();
        int nmloc = ist_loc * mloc;
        int ngwloc = rwf.sd(isp_loc, ikp_loc)->basis().localsize();
        const double* kpg_ipol = rwf.sd(isp_loc, ikp_loc)->basis().kpgx_ptr(ipol);
        for ( int ig = 0; ig < ngwloc; ig++ )
        {
          complex<double> momentum = complex<double>(0.0, 0.0);
          momentum = kpg_ipol[ig];
          //momentum = kpg_ipol[ig] + ef.vp->vecpot()[ipol]*rwf.nel();
          pwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] = momentum * rwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig];
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void RTCurrentDensity::compute_vnlr(Wavefunction& rwf, Wavefunction& vwf, EnergyFunctional& ef, int ipol)
{
  Wavefunction awf(rwf);
  Wavefunction anwf(rwf);
  Wavefunction nwf(rwf);

  awf = rwf;
  anwf = rwf;
  nwf = rwf;
  vwf = rwf;

  vector<vector<double> > fion;
  valarray<double> sigma(6);
  fion.resize(atoms_.nsp());
  for ( int is = 0; is < atoms_.nsp(); is++ )
  {
    fion[is].resize(3*atoms_.na(is));
  }

  for ( int isp_loc = 0; isp_loc < rwf.nsp_loc(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < rwf.nkp_loc(); ++ikp_loc )
    {
      // nwf = Vnl|psi>
      nwf.sd(isp_loc, ikp_loc)->c().clear();
      //ef.nlpout(isp_loc, ikp_loc)->energy(*rwf.sd(isp_loc, ikp_loc), true, *nwf.sd(isp_loc, ikp_loc), false, fion, false, sigma);
      ef.nlpout(isp_loc, ikp_loc)->nlpsi(*rwf.sd(isp_loc, ikp_loc), *nwf.sd(isp_loc, ikp_loc));

      // awf = r|psi>
      valarray<complex<double> > wf_r(rtp_.ft_dwfc(ikp_loc)->np012loc());
      valarray<complex<double> > wf_n(rtp_.ft_dwfc(ikp_loc)->np012loc());
      int mloc = rwf.sd(isp_loc, ikp_loc)->c().mloc();

      for ( int ist_loc = 0; ist_loc < rwf.sd(isp_loc, ikp_loc)->nstloc(); ++ist_loc )
      {
        int nmloc = ist_loc * mloc;
        for ( int ig = 0; ig < wf_r.size(); ig++ )
        {
          wf_r[ig] = complex<double>(0.0, 0.0);
          wf_n[ig] = complex<double>(0.0, 0.0);
        }
        // awf = r|psi>
        rtp_.ft_dwfc(ikp_loc)->backward(awf.sd(isp_loc, ikp_loc)->c().cvalptr(nmloc), &wf_r[0]);
        // nwf = Vnl|psi>
        rtp_.ft_dwfc(ikp_loc)->backward(nwf.sd(isp_loc, ikp_loc)->c().cvalptr(nmloc), &wf_n[0]);

        for ( int ig = 0; ig < wf_r.size(); ig++ )
        {
          // wf_r = r|psi>
          wf_r[ig] = rtp_.rt_position_dwfc(ikp_loc, ig)[ipol] * wf_r[ig];
          // wf_rn = r*Vnl|psi>
          wf_n[ig] = rtp_.rt_position_dwfc(ikp_loc, ig)[ipol] * wf_n[ig];
        }
        rtp_.ft_dwfc(ikp_loc)->forward(&wf_r[0], awf.sd(isp_loc, ikp_loc)->c().valptr(nmloc));
        rtp_.ft_dwfc(ikp_loc)->forward(&wf_n[0], nwf.sd(isp_loc, ikp_loc)->c().valptr(nmloc));
      }
      // anwf = Vnl*r|psi>
      anwf.sd(isp_loc, ikp_loc)->c().clear();
      ef.nlpout(isp_loc, ikp_loc)->nlpsi(*awf.sd(isp_loc, ikp_loc), *anwf.sd(isp_loc, ikp_loc));
      for ( int ist_loc = 0; ist_loc < rwf.sd(isp_loc, ikp_loc)->nstloc(); ++ist_loc )
      {
        int nmloc = ist_loc * mloc;
        int ngwloc = rwf.sd(isp_loc, ikp_loc)->basis().localsize();
        // vwf = Vnl*r|psi> - r*Vnl|psi>
        for ( int ig = 0; ig < ngwloc; ig++ )
        {
          vwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] = anwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] - nwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig];
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void RTCurrentDensity::compute_hr(Wavefunction& rwf, Wavefunction& cwf, EnergyFunctional& ef, int ipol)
{
  Wavefunction awf(rwf);
  Wavefunction anwf(rwf);
  Wavefunction nwf(rwf);

  awf = rwf;

  vector<vector<double> > fion;
  valarray<double> sigma;


  for ( int isp_loc = 0; isp_loc < rwf.nsp_loc(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < rwf.nkp_loc(); ++ikp_loc )
    {
      // nwf = H|psi>
      ef.hpsi(isp_loc, ikp_loc, rwf, nwf);
      // awf = r|psi>
      valarray<complex<double> > wf_r(rtp_.ft_dwfc(ikp_loc)->np012loc());
      valarray<complex<double> > wf_n(rtp_.ft_dwfc(ikp_loc)->np012loc());
      int mloc = rwf.sd(isp_loc, ikp_loc)->c().mloc();

      for ( int ist_loc = 0; ist_loc < rwf.sd(isp_loc, ikp_loc)->nstloc(); ++ist_loc )
      {
        int nmloc = ist_loc * mloc;
        for ( int ig = 0; ig < wf_r.size(); ig++ )
        {
          wf_r[ig] = complex<double>(0.0, 0.0);
          wf_n[ig] = complex<double>(0.0, 0.0);
        }
        // awf = r|psi>
        rtp_.ft_dwfc(ikp_loc)->backward(awf.sd(isp_loc, ikp_loc)->c().cvalptr(nmloc), &wf_r[0]);
        // nwf = H|psi>
        rtp_.ft_dwfc(ikp_loc)->backward(nwf.sd(isp_loc, ikp_loc)->c().cvalptr(nmloc), &wf_n[0]);

        for ( int ig = 0; ig < wf_r.size(); ig++ )
        {
          // wf_r = r|psi>
          wf_r[ig] = rtp_.rt_position_dwfc(ikp_loc, ig)[ipol] * wf_r[ig];
          // wf_rn = r*Vnl|psi>
          wf_n[ig] = rtp_.rt_position_dwfc(ikp_loc, ig)[ipol] * wf_n[ig];
        }
        rtp_.ft_dwfc(ikp_loc)->forward(&wf_r[0], awf.sd(isp_loc, ikp_loc)->c().valptr(nmloc));
        rtp_.ft_dwfc(ikp_loc)->forward(&wf_n[0], nwf.sd(isp_loc, ikp_loc)->c().valptr(nmloc));
      }
      // anwf = H*r|psi>
      ef.hpsi(isp_loc, ikp_loc, awf, anwf);

      for ( int ist_loc = 0; ist_loc < rwf.sd(isp_loc, ikp_loc)->nstloc(); ++ist_loc )
      {
        int nmloc = ist_loc * mloc;
        int ngwloc = rwf.sd(isp_loc, ikp_loc)->basis().localsize();
        // hwf = H*r|psi> - r*H|psi>
        for ( int ig = 0; ig < ngwloc; ig++ )
        {
          cwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] = anwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] - nwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig];
        }
      }
    }
  }
}
