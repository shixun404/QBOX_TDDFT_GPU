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
// RTProjection.cpp
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#include "RTProjection.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "MPIdata.h"

#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
RTProjection::RTProjection(Wavefunction& wf, string rtprojtype, int proj_bd_s, int proj_bd_e, int proj_kp_a, int proj_kp_b):
pwf_(wf), rtprojtype_(rtprojtype), proj_bd_s_(proj_bd_s), proj_bd_e_(proj_bd_e), proj_kp_a_(proj_kp_a), proj_kp_b_(proj_kp_b), nbnd(0)
{
  nbnd = proj_bd_e_ - proj_bd_s_ + 1;
}

////////////////////////////////////////////////////////////////////////////////
RTProjection::~RTProjection(void)
{
}

////////////////////////////////////////////////////////////////////////////////
void RTProjection::set_collection(void)
{
  if ( rtprojtype_ == "BAND" )
  {
    proj_band_.resize(pwf_.nspin());
    for ( int ispin = 0; ispin < pwf_.nspin(); ++ispin )
    {
      proj_band_[ispin].resize(nbnd);
      for ( int ibnd = 0; ibnd < nbnd; ++ibnd )
      {
        proj_band_[ispin][ibnd].resize(nbnd);
      }
    }
  }

  if ( rtprojtype_ == "SELKP" )
  {
    proj_.resize(pwf_.nspin());
    for ( int ispin = 0; ispin < pwf_.nspin(); ++ispin )
    {
      if ( proj_kp_a_ == proj_kp_b_ )
      {
        proj_[ispin].resize(1);
	proj_[ispin][0].resize(nbnd);
	for ( int ibnd = 0; ibnd < nbnd; ++ibnd )
	{
	  proj_[ispin][0][ibnd].resize(nbnd);
        }
      }
      else
      {
        proj_[ispin].resize(2);
	proj_[ispin][0].resize(nbnd);
	proj_[ispin][1].resize(nbnd);
	for ( int ibnd = 0; ibnd < nbnd; ++ibnd )
	{
	  proj_[ispin][0][ibnd].resize(nbnd);
	  proj_[ispin][1][ibnd].resize(nbnd);
        }
      }
    }
  }

  //if ( rtprojtype_ == "ALLKP" )
  //{
  //  proj_.resize(pwf_.nspin());
  //  for ( int ispin = 0; ispin < pwf_.nspin(); ++ispin )
  //  {
  //    proj_[ispin].resize(pwf_.nkp());
  //    for ( int ikp = 0; ikp < pwf_.nkp(); ++ikp )
  //    {
  //      proj_[ispin][ikp].resize(nbnd);
  //      for ( int ibnd = 0; ibnd < nbnd; ++ibnd )
  //      {
  //        proj_[ispin][ikp][ibnd].resize(nbnd);
  //      }
  //    }
  //  }
  //}

  msize.resize(pwf_.nspin());
  for ( int ispin = 0; ispin < pwf_.nspin(); ++ispin )
  {
    msize[ispin].resize(pwf_.nkp());
    int isp_loc = pwf_.isp_local(ispin);
    for ( int ikp = 0; ikp < pwf_.nkp(); ++ikp )
    {
      int ikp_loc = pwf_.ikp_local(ikp);
      int m = 0;
      msize[ispin][ikp] = 0;
      MPI_Barrier(MPIdata::comm());
      if ( isp_loc >= 0 && ikp_loc >= 0 )
      {
        m = pwf_.sd(isp_loc, ikp_loc)->c().m();
      }
      MPI_Barrier(MPIdata::comm());
      MPI_Allreduce(&m, &msize[ispin][ikp], 1, MPI_INTEGER, MPI_MAX, MPIdata::kp_sp_comm());
    }
  }

  coll_init_.resize(pwf_.nspin());
  coll_band_.resize(pwf_.nspin());
  for ( int ispin = 0; ispin < pwf_.nspin(); ++ispin )
  {
    coll_init_[ispin].resize(pwf_.nkp());
    coll_band_[ispin].resize(pwf_.nkp());
    for ( int ikp = 0; ikp < pwf_.nkp(); ++ikp )
    {
      coll_init_[ispin][ikp].resize(pwf_.nst(ispin));
      coll_band_[ispin][ikp].resize(pwf_.nst(ispin));
      for ( int ist = 0; ist < pwf_.nst(ispin); ++ist )
      {
        coll_init_[ispin][ikp][ist].resize(msize[ispin][ikp]);
        coll_band_[ispin][ikp][ist].resize(msize[ispin][ikp]);
      }
    }
  }

  init_collection(); // coll_band_ initialize
  init_collection2(); // coll_init_ initialize
}

////////////////////////////////////////////////////////////////////////////////
void RTProjection::init_collection(void)
{
  for ( int ispin = 0; ispin < pwf_.nspin(); ++ispin )
  {
    for ( int ikp = 0; ikp < pwf_.nkp(); ++ikp )
    {
      for ( int ist = 0; ist < pwf_.nst(ispin); ++ist )
      {
	for ( int ig = 0; ig < msize[ispin][ikp]; ig++ )
          coll_band_[ispin][ikp][ist][ig] = complex<double>(0.0, 0.0);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void RTProjection::init_collection2(void)
{
  for ( int ispin = 0; ispin < pwf_.nspin(); ++ispin )
  {
    for ( int ikp = 0; ikp < pwf_.nkp(); ++ikp )
    {
      for ( int ist = 0; ist < pwf_.nst(ispin); ++ist )
      {
	for ( int ig = 0; ig < msize[ispin][ikp]; ig++ )
          coll_init_[ispin][ikp][ist][ig] = complex<double>(0.0, 0.0);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void RTProjection::start_collection(Wavefunction& iwf)
{
  for ( int ispin = 0; ispin < iwf.nspin(); ++ispin )
  {
    int isp_loc = iwf.isp_local(ispin);
    for ( int ikp = 0; ikp < iwf.nkp(); ++ikp )
    {
      int ikp_loc = iwf.ikp_local(ikp);
      std::vector<std::vector<complex<double> > > coll; // coll[nst][ig]
      coll.resize(iwf.nst(ispin));
      for ( int ist = 0; ist < iwf.nst(ispin); ++ist )
      {
	coll[ist].resize(msize[ispin][ikp], complex<double>(0.000, 0.000));
      }

      if ( isp_loc >= 0 && ikp_loc >= 0 )
      {
        for ( int ist_loc = 0; ist_loc < iwf.sd(isp_loc, ikp_loc)->nstloc(); ++ist_loc )
        {
          int mloc = iwf.sd(isp_loc, ikp_loc)->c().mloc();
	  int nmloc = mloc * ist_loc;
	  for ( int ig = 0; ig < mloc; ig++ )
	  {
	    int nidx = iwf.sd(isp_loc, ikp_loc)->c().jglobal(ist_loc);
	    int midx = iwf.sd(isp_loc, ikp_loc)->c().iglobal(ig);
	    coll[nidx][midx] = iwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig];
	  }
        }
      }

      MPI_Barrier(MPIdata::comm());
      for ( int ist = 0; ist < iwf.nst(ispin); ++ist )
      {
	MPI_Allreduce(&coll[ist][0], &coll_init_[ispin][ikp][ist][0], msize[ispin][ikp], MPI_DOUBLE_COMPLEX, MPI_SUM, MPIdata::comm());
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void RTProjection::band_collection(Wavefunction& cwf)
{
  for ( int ispin = 0; ispin < cwf.nspin(); ++ispin )
  {
    int isp_loc = cwf.isp_local(ispin);
    for ( int ikp = 0; ikp < cwf.nkp(); ++ikp )
    {
      int ikp_loc = cwf.ikp_local(ikp);
      std::vector<std::vector<complex<double> > > coll; // coll[nst][ig]
      coll.resize(cwf.nst(ispin));
      for ( int ist = 0; ist < cwf.nst(ispin); ++ist )
      {
	coll[ist].resize(msize[ispin][ikp], complex<double>(0.0, 0.0));
      }

      if ( isp_loc >= 0 && ikp_loc >= 0 )
      {
        for ( int ist_loc = 0; ist_loc < cwf.sd(isp_loc, ikp_loc)->nstloc(); ++ist_loc )
        {
          int mloc = cwf.sd(isp_loc, ikp_loc)->c().mloc();
	  int nmloc = mloc * ist_loc;
	  for ( int ig = 0; ig < mloc; ig++ )
	  {
	    int nidx = cwf.sd(isp_loc, ikp_loc)->c().jglobal(ist_loc);
	    int midx = cwf.sd(isp_loc, ikp_loc)->c().iglobal(ig);
	    coll[nidx][midx] = cwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig];
	  }
        }
      }

      MPI_Barrier(MPIdata::comm());
      for ( int ist = 0; ist < cwf.nst(ispin); ++ist )
      {
	MPI_Allreduce(&coll[ist][0], &coll_band_[ispin][ikp][ist][0], msize[ispin][ikp], MPI_DOUBLE_COMPLEX, MPI_SUM, MPIdata::comm());
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void RTProjection::compute_projection(Wavefunction& ppwf)
{
  init_collection();
  band_collection(ppwf);

  if ( rtprojtype_ == "BAND" )
  {
    for ( int ispin = 0; ispin < ppwf.nspin(); ++ispin )
    {
      int ibnd_idx = 0;
      for ( int ibnd = proj_bd_s_-1; ibnd <= proj_bd_e_-1; ++ibnd )
      {
        int jbnd_idx = 0;
        for ( int jbnd = proj_bd_s_-1; jbnd <= proj_bd_e_-1; ++jbnd )
        {
          proj_band_[ispin][ibnd_idx][jbnd_idx] = 0.0;
          for ( int ikp = 0; ikp < ppwf.nkp(); ++ikp )
          {
            complex<double> proj = complex<double>(0.0, 0.0);
            for ( int ig = 0; ig < msize[ispin][ikp]; ig++ )
            {
              proj += coll_band_[ispin][ikp][ibnd][ig] * conj(coll_init_[ispin][ikp][jbnd][ig]);
            }
            proj_band_[ispin][ibnd_idx][jbnd_idx] += ppwf.weight(ikp) * (proj * conj(proj)).real();
          }
          jbnd_idx = jbnd_idx + 1;
        }
        ibnd_idx = ibnd_idx + 1;
      }
    }
  }

  if ( rtprojtype_ == "SELKP" )
  {
    for ( int ispin = 0; ispin < ppwf.nspin(); ++ispin )
    {
      int ibnd_idx = 0;
      for ( int ibnd = proj_bd_s_-1; ibnd <= proj_bd_e_-1; ++ibnd )
      {
        int jbnd_idx = 0;
        for ( int jbnd = proj_bd_s_-1; jbnd <= proj_bd_e_-1; ++jbnd )
        {
          proj_[ispin][0][ibnd_idx][jbnd_idx] = 0.0;
          proj_[ispin][1][ibnd_idx][jbnd_idx] = 0.0;
          for ( int ig = 0; ig < min(msize[ispin][proj_kp_a_], msize[ispin][proj_kp_b_]); ig++ )
          {
            proj_[ispin][0][ibnd_idx][jbnd_idx] += ppwf.weight(proj_kp_b_) * (coll_band_[ispin][proj_kp_b_][ibnd][ig] * conj(coll_init_[ispin][proj_kp_a_][jbnd][ig])).real();
            proj_[ispin][1][ibnd_idx][jbnd_idx] += ppwf.weight(proj_kp_b_) * (coll_band_[ispin][proj_kp_b_][ibnd][ig] * conj(coll_init_[ispin][proj_kp_a_][jbnd][ig])).real();
          }
          jbnd_idx = jbnd_idx + 1;
        }
        ibnd_idx = ibnd_idx + 1;
      }
    }
  }

  //if ( rtprojtype_ == "ALLKP" )
  //{
  //  for ( int ispin = 0; ispin < wf.nspin(); ++ispin )
  //  {
  //    for ( int ikp = 0; ikp < wf.nkp(); ++ikp )
  //    {
  //      int ibnd_idx = 0;
  //      for ( int ibnd = proj_bd_s_; ibnd <= proj_bd_e; ++ibnd )
  //      {
  //        int jbnd_idx = 0;
  //        for ( int jbnd = proj_bd_s_; jbnd <= proj_bd_e_; ++jbnd )
  //        {
  //          proj_band_[ispin][ibnd_idx][jbnd_idx] = 0.0;
  //          for ( int ig = 0; ig < msize[ispin][ikp]; ig++ )
  //          {
  //            proj_band_[ispin][ibnd_idx][jbnd_idx] += wf.weight(ikp) * (coll_band_[ispin][ikp][ibnd][ig] * conj(coll_init_[ispin][ikp][jbnd][ig])).real();
  //          }
  //          jbnd_idx = jbnd_idx + 1;
  //        }
  //        ibnd_idx = ibnd_idx + 1;
  //      }
  //    }
  //  }
  //}
}

