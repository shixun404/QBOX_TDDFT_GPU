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
// RTCNWavefunctionStepper.cpp
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////
//
// Conjugate gradient square method (J. Comput. Appl. Math., 71, 125-146 (1996)
//
////////////////////////////////////////////////////////////////////////////////

#include "RTCNWavefunctionStepper.h"
#include "Wavefunction.h"
#include "EnergyFunctional.h"
#include "SlaterDet.h"
#include "Sample.h"
#include "MPIdata.h"
#include "Basis.h"
#include "FourierTransform.h"
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
RTCNWavefunctionStepper::RTCNWavefunctionStepper(double rtdt, TimerMap& tmap, 
  EnergyFunctional& ef, ChargeDensity& cd, Sample& s, int rtitscf, int rtite) 
  : rtdt_(rtdt), WavefunctionStepper(s.wf,tmap), ef_(ef), cd_(cd), s_(s), 
  rtitscf_(rtitscf), rtite_(rtite)
{
  tmap_["rt_cn"].reset();
  tmap_["rt_cg"].reset();
  tmap_["rt_cn_conv"].reset();
}

////////////////////////////////////////////////////////////////////////////////
void RTCNWavefunctionStepper::get_iter(int iter)
{
  iter_ = iter;
}

////////////////////////////////////////////////////////////////////////////////
void RTCNWavefunctionStepper::update(Wavefunction& rt_wf_next)
{
  const bool onpe0 = MPIdata::onpe0();
  cout.precision(18);

  bool scf_converged = false;
  bool ele_converged = false;
  vector<bool> ele_converged_tmp;
  int itscf;
  int itele;
  int itele_sum;
  complex<double> rho;
  complex<double> rhotmp;
  complex<double> rhos;
  complex<double> rhostmp;
  complex<double> rhor;
  complex<double> rhortmp;
  complex<double> alpha;
  complex<double> beta;
  vector<vector<double> > fion;
  valarray<double> sigma;
  double rhor_sum;
  int ngloc = cd_.vbasis()->localsize();

  Wavefunction owf(s_.wf); // original wavefunction |psi(t)>
  Wavefunction xwf2(s_.wf); // tmp for wavefunction |xwf(t+dt)>

  Wavefunction pwf(s_.wf);
  Wavefunction apwf(s_.wf);
  Wavefunction hpwf(s_.wf);

  Wavefunction bwf(s_.wf); 
  Wavefunction hbwf(s_.wf);

  Wavefunction xwf(s_.wf);
  Wavefunction hxwf(s_.wf);
  Wavefunction rwf(s_.wf);
  
  // initialize
  // owf = |psi(t)> --> FIX
  owf = s_.wf;
  // xwf = |psi(t+dt)> = x0
  xwf = rt_wf_next;
  // xwf2 = xwf previous step
  xwf2 = xwf;

  // Hamiltonian

  // cdori = rho(t)
#if OPTIMIZE_TRANSPOSE  
  ChargeDensity cdori(owf,true);
#else
  ChargeDensity cdori(owf);
#endif
  
  
  tmap_["charge"].start();
  cdori.update_density();
  tmap_["charge"].stop();

  // cdupd = rho(t+dt)
#if OPTIMIZE_TRANSPOSE
  ChargeDensity cdupd(xwf,true);
#else
  ChargeDensity cdupd(xwf);
#endif

  // cdpre = rho(t+dt) previous step
#if OPTIMIZE_TRANSPOSE 
  ChargeDensity cdpre(xwf2,true);
#else
  ChargeDensity cdpre(xwf2);
#endif  


  tmap_["charge"].start();
  cdpre.update_density();
  tmap_["charge"].stop();

  // efupd = H(t+dt)
  EnergyFunctional efupd(s_, xwf2, cdupd);
  if ( s_.rtctrl.rt_vp == "ON" )
    efupd.vp->vp_propagate(iter_, rtdt_);
  tmap_["update_vhxc"].start();
  efupd.update_vhxc(false);
  tmap_["update_vhxc"].stop();

  // Ax0 = b linear problem solving (CGS method)
  itscf = 0;
  if ( onpe0 )
    cout << "<RTCNWavefunctionStepper : CN iteration START : " << rtdt_ << " : " << rtitscf_ << " : " << rtite_ << " />" << endl;

  tmap_["rt_cn"].start();
  while ( !scf_converged && itscf < rtitscf_ )
  {
    /////////////////// Charge Density Mixing Start ///////////////////
    // rho(t+dt/2) = [rho(t)+rho(t+dt)]/2
    cdpre.update_density();
    for ( int ispin = 0; ispin < s_.wf.nspin(); ispin++ )
    {
      for ( int ig = 0; ig < ngloc; ig++ )
      {
        cdupd.rhog[ispin][ig] = 0.5*(cdori.rhog[ispin][ig] + cdpre.rhog[ispin][ig]);
      }
    }

    tmap_["charge"].start();
    cdupd.update_rhor();
    tmap_["charge"].stop();

    // H(t+dt/2) generate
    tmap_["update_vhxc"].start();
    efupd.update_vhxc(false);
    tmap_["update_vhxc"].stop();

    /////////////////// Charge Density Mixing End ///////////////////
    itele_sum = 0;
    rhor_sum = 0.0;


    /////////////////// Conjugate Gradient Initiation Start ///////////////////
    tmap_["rt_cg"].start();
    for ( int isp_loc = 0; isp_loc < s_.wf.nsp_loc(); ++isp_loc )
    {
      for ( int ikp_loc = 0; ikp_loc < s_.wf.nkp_loc(); ++ikp_loc )
      {

        int nstloc = s_.wf.sd(isp_loc, ikp_loc)->nstloc();
	ele_converged_tmp.resize(nstloc);
	hbwf.sd(isp_loc, ikp_loc)->c().clear();
	efupd.hpsi(isp_loc, ikp_loc, owf, hbwf);
	hxwf.sd(isp_loc, ikp_loc)->c().clear();
	efupd.hpsi(isp_loc, ikp_loc, xwf, hxwf);
	int mloc = s_.wf.sd(isp_loc, ikp_loc)->c().mloc();
	bwf.sd(isp_loc, ikp_loc)->c().clear();
	rwf.sd(isp_loc, ikp_loc)->c().clear();
	pwf.sd(isp_loc, ikp_loc)->c().clear();
	
	for ( int ist_loc = 0; ist_loc < nstloc; ++ist_loc )
	{
	  int nmloc = ist_loc * mloc;
	  ele_converged_tmp[ist_loc] = false;
	  for ( int ig = 0; ig < mloc; ig++ )
	  {
            bwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] = owf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] - complex<double>(0.0, 0.5 * rtdt_) * hbwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig];
	    rwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] = bwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] - (xwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] + complex<double>(0.0, 0.5 * rtdt_) * hxwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig]);
	    pwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] = rwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig];
	  }
	}

	hpwf.sd(isp_loc, ikp_loc)->c().clear();
	apwf.sd(isp_loc, ikp_loc)->c().clear();

	itele = 0;
	ele_converged = false;

        while ( !ele_converged && itele < rtite_ )
	{
          efupd.hpsi(isp_loc, ikp_loc, pwf, hpwf);
	  ele_converged = true;
	  for ( int ist_loc = 0; ist_loc < nstloc; ++ist_loc )
	  {
	    rho = complex<double>(0.0, 0.0);
	    rhos = complex<double>(0.0, 0.0);
	    rhotmp = complex<double>(0.0, 0.0);
	    rhostmp = complex<double>(0.0, 0.0);
	    rhor = complex<double>(0.0, 0.0);
	    rhortmp = complex<double>(0.0, 0.0);
	    alpha = complex<double>(0.0, 0.0);
	    beta = complex<double>(0.0, 0.0);
	    int nmloc = ist_loc * mloc;
            if ( !ele_converged_tmp[ist_loc] )
	    {

	    for ( int ig = 0; ig < mloc; ig++ )
	      {
	        apwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] = pwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] + complex<double>(0.0, 0.5 * rtdt_) * hpwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig];
	        rhotmp = rhotmp + conj( rwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] ) * ( rwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] );
	        rhostmp = rhostmp + conj( pwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] ) * ( apwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] );
	      }

	      MPI_Allreduce(&rhotmp, &rho, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPIdata::g_comm());
	      MPI_Allreduce(&rhostmp, &rhos, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPIdata::g_comm());
	      alpha = rho / rhos;


	      for ( int ig = 0; ig < mloc; ig++ )
	      {
	        xwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] = xwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] + alpha * pwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig];
	        rwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] = rwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] - alpha * apwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig];
	      }

	      for ( int ig = 0; ig < mloc; ig++ )
	      {
	        rhortmp = rhortmp + conj ( rwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] ) * ( rwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] );
	      }
	      MPI_Allreduce(&rhortmp, &rhor, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPIdata::g_comm());

	      beta = rhor / rho;
	      for ( int ig = 0; ig < mloc; ig++ )
	      {
	        pwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] = rwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] + beta * pwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig];
	      }
	      if ( abs(rhor) < s_.rtctrl.rt_scf_tol )
	      {
	        ele_converged_tmp[ist_loc] = true;
	      }
	      else
	      {
	        ele_converged_tmp[ist_loc] = false;
	      }
	      rhor_sum = rhor_sum + abs(rhor);
	    }


	    ele_converged = ( ele_converged && ele_converged_tmp[ist_loc] );
	  }
          if ( !ele_converged )
            itele++;
	}
	
      }
    }
    ///////////////////// CG LOOP END ///////////////////////////////////////
    itele_sum = itele_sum + itele;
    tmap_["rt_cg"].stop();
    /////////////////// Conjugate Gradient Algorithm End ///////////////////

    /////////////////// Hamiltonian Update Start ///////////////////
    double rhor_sum_tmp = 0.0;
    MPI_Allreduce(&rhor_sum, &rhor_sum_tmp, 1, MPI_DOUBLE, MPI_SUM, MPIdata::st_kp_sp_comm());
    rhor_sum = rhor_sum_tmp / (s_.wf.nspin() * s_.wf.nkp() * s_.wf.nst());
    int itele_sum_tmp = 0;

    MPI_Allreduce(&itele_sum, &itele_sum_tmp, 1, MPI_INT, MPI_SUM, MPIdata::st_kp_sp_comm());
    itele_sum = itele_sum_tmp / (s_.wf.nspin() * s_.wf.nkp() * s_.wf.nst());

    tmap_["charge"].start();
    cdupd.update_density();
    tmap_["charge"].stop();

    tmap_["charge"].start();
    cdpre.update_density();
    tmap_["charge"].stop();

    double tsum;
    tmap_["rt_cn_conv"].start();
    scf_converged = conv_check(cdpre, cdupd, tsum);
    tmap_["rt_cn_conv"].stop();

    if ( onpe0 )
      cout << " RTCN : " << itscf << " | DRHO : " << tsum << " | ITELE :" << itele_sum << " | RHOR : " << rhor_sum << endl;
    if ( !scf_converged )
    {
      itscf = itscf + 1;
      xwf2 = xwf;
    }
  }
  
  tmap_["rt_cn"].stop();
  s_.wf = xwf;
  if ( onpe0 )
    cout << "<RTCNWavefunctionStepper : CN iteration END : " << rtitscf_ << " : " << rtite_ << " />" << endl;
}

////////////////////////////////////////////////////////////////////////////////
bool RTCNWavefunctionStepper::conv_check(ChargeDensity& cd1, ChargeDensity& cd2, double& tsum)
{
  std::vector<std::vector<double> > drhor; 
  std::vector<std::vector<double> > srhor; 
  const double omega = cd_.vbasis()->cell().volume();
  double sumtmp;
  double ssumtmp;

  drhor.resize(s_.wf.nspin());
  srhor.resize(s_.wf.nspin());
  tsum = 0.0;
  double ssum = 0.0;
  sumtmp = 0.0;
  ssumtmp = 0.0;

  for ( int ispin = 0; ispin < s_.wf.nspin(); ispin++ )
  {
    drhor[ispin].resize(cd_.vft()->np012loc());
    srhor[ispin].resize(cd_.vft()->np012loc());
    double drhosum = 0.0;
    double srhosum = 0.0;
    for ( int ir = 0; ir < cd_.vft()->np012loc(); ir++ )
    {
      drhor[ispin][ir] = (cd2.rhor[ispin][ir] - cd1.rhor[ispin][ir]);
      srhor[ispin][ir] = cd1.rhor[ispin][ir];
      drhosum = drhosum + abs(drhor[ispin][ir]);
      srhosum = srhosum + abs(srhor[ispin][ir]);
    }
    double tsumtmp = 0.0;
    double tssumtmp = 0.0;
    MPI_Allreduce(&drhosum,&tsumtmp,1,MPI_DOUBLE,MPI_SUM,MPIdata::g_comm());
    MPI_Allreduce(&srhosum,&tssumtmp,1,MPI_DOUBLE,MPI_SUM,MPIdata::g_comm());
    tsumtmp = tsumtmp * omega / cd_.vft()->np012();
    tssumtmp = tssumtmp * omega / cd_.vft()->np012();
    tsumtmp = abs(tsumtmp) / s_.wf.nel();
    tssumtmp = abs(tssumtmp) / s_.wf.nel();
    sumtmp = sumtmp + tsumtmp;
    ssumtmp = ssumtmp + tssumtmp;
  }

  MPI_Allreduce(&sumtmp, &tsum, 1, MPI_DOUBLE, MPI_SUM, MPIdata::sp_comm());
  MPI_Allreduce(&ssumtmp, &ssum, 1, MPI_DOUBLE, MPI_SUM, MPIdata::sp_comm());
  tsum = tsum / (double) s_.wf.nspin();
  ssum = ssum / (double) s_.wf.nspin();

  if ( MPIdata::onpe0() )
	  cout << "RTCN RHO TEST : " << ssum << endl;
  if ( tsum < s_.rtctrl.rt_rho_tol )
  {
    return true;
  }
  else
  {
    return false;
  }
}
