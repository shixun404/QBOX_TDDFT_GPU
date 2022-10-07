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
// RTPTCNWavefunctionStepper.cpp
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

#include "RTPTCNWavefunctionStepper.h"
#include "Wavefunction.h"
#include "EnergyFunctional.h"
#include "SlaterDet.h"
#include "Sample.h"
#include "MPIdata.h"
#include "Basis.h"
#include "FourierTransform.h"
#include "RTVectorPotential.h"
#include "Context.h"
#include "AndersonMixer.h"
#include <iostream>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
RTPTCNWavefunctionStepper::RTPTCNWavefunctionStepper(double rtdt, TimerMap& tmap, 
  EnergyFunctional& ef, ChargeDensity& cd, Sample& s, int rtitscf) 
  : rtdt_(rtdt), WavefunctionStepper(s.wf,tmap), ef_(ef), cd_(cd), s_(s), 
  rtitscf_(rtitscf)
{
  MPI_Comm vcomm = MPIdata::g_comm();
  const int dim = s_.rtctrl.rt_anderson_dim;
  mixer_.resize(s_.wf.nsp_loc());
  for ( int isp_loc = 0; isp_loc < s_.wf.nsp_loc(); ++isp_loc )
  {
    mixer_[isp_loc].resize(s_.wf.nkp_loc());
    for ( int ikp_loc = 0; ikp_loc < s_.wf.nkp_loc(); ++ikp_loc )
    {
      mixer_[isp_loc][ikp_loc].resize(s_.wf.sd(isp_loc, ikp_loc)->nstloc());
      for ( int ist_loc = 0; ist_loc < s_.wf.sd(isp_loc, ikp_loc)->nstloc(); ++ist_loc)
      {
	int m = s_.wf.sd(isp_loc, ikp_loc)->c().m();
	int mloc = s_.wf.sd(isp_loc, ikp_loc)->c().mloc();
      //mixer_[isp_loc][ikp_loc][ist_loc] = new AndersonMixer(2*m, dim, 0);
	  ////////////////// PARALLELIZE /////////////////
        mixer_[isp_loc][ikp_loc][ist_loc] = new AndersonMixer(2*mloc, dim, &MPIdata::g_comm());
	  ////////////////////////////////////////////////
      }
    }
  }
  tmap_["rt_pt_cn"].reset();
  tmap_["rt_pt_cn_conv"].reset();
}

////////////////////////////////////////////////////////////////////////////////
RTPTCNWavefunctionStepper::~RTPTCNWavefunctionStepper(void)
{
  for ( int isp_loc = 0; isp_loc < s_.wf.nsp_loc(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < s_.wf.nkp_loc(); ++ikp_loc )
    {
      for ( int ist_loc = 0; ist_loc < s_.wf.sd(isp_loc, ikp_loc)->nstloc(); ++ist_loc )
      {
        delete mixer_[isp_loc][ikp_loc][ist_loc];
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void RTPTCNWavefunctionStepper::get_iter(int& iter)
{
  iter_ = iter;
}

////////////////////////////////////////////////////////////////////////////////
void RTPTCNWavefunctionStepper::update(Wavefunction& rt_wf)
{
  const bool onpe0 = MPIdata::onpe0();
  cout.precision(18);

  bool scf_converged = false;
  int itscf;
  vector<vector<double> > fion;
  valarray<double> sigma;

  const Context& sd_ctxt = rt_wf.sd_context();

  Wavefunction owf(rt_wf); // original wavefunction |psi(t)>
  Wavefunction howf(rt_wf); // original wavefunction H|psi(t)>
  Wavefunction rwf(rt_wf); // residual function rwf(t)>
  Wavefunction xwf(rt_wf); // initial step wavefunction |psi(t+dt)>
  Wavefunction fwf(rt_wf); // wavefunction |psi(t+dt)>
  Wavefunction hfwf(rt_wf); // wavefunction H|psi(t+dt)>
  Wavefunction fwf2(rt_wf); // previous step wavefunction |psi(t+dt)>

  // owf = |psi(t)> --> FIX
  owf = rt_wf;
  rwf = rt_wf;
  xwf = rt_wf;
  fwf = rt_wf;
  fwf2 = rt_wf;

  // Hamiltonian
  // cdori = rho(t)
  ChargeDensity cdori(owf);
  tmap_["charge"].start();
  cdori.update_density();
  tmap_["charge"].stop();
  // efori = H(t)
  EnergyFunctional efori(s_, cdori);
  if ( s_.rtctrl.rt_vp == "ON" )
    efori.vp->vp_propagate(iter_, rtdt_);
  tmap_["update_vhxc"].start();
  efori.update_vhxc(false);
  tmap_["update_vhxc"].stop();

  // howf = H(t)|psi(t)>
  efori.energy(true, howf, false, fion, false, sigma);

  for ( int isp_loc = 0; isp_loc < rt_wf.nsp_loc(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < rt_wf.nkp_loc(); ++ikp_loc )
    {
      int mloc = rt_wf.sd(isp_loc, ikp_loc)->c().mloc();
      rwf.sd(isp_loc, ikp_loc)->c().clear();
      xwf.sd(isp_loc, ikp_loc)->c().clear();
      fwf.sd(isp_loc, ikp_loc)->c().clear();
      fwf2.sd(isp_loc, ikp_loc)->c().clear();
      for ( int ist_loc = 0; ist_loc < rt_wf.sd(isp_loc, ikp_loc)->nstloc(); ++ist_loc )
      {
        int nmloc = ist_loc * mloc;
	for ( int ig = 0; ig < mloc; ig++ )
	{
          // Evaluate the initial residual R_n = H_n|psi_n> - |psi_n><psi_n|H_n|psi_n>
          rwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] = howf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] - owf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] * conj(owf.sd(isp_loc, ikp_loc)->c()[nmloc + ig]) * howf.sd(isp_loc, ikp_loc)->c()[nmloc + ig];
	  // Evaluate psi_n+1/2 = psi_n - (idt/2)R_n
	  xwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] = owf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] - complex<double>(0.0, 0.5 * rtdt_) * rwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig];
	  // , and let psi_f = psi_n+1/2
	  fwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] = xwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig];
	  fwf2.sd(isp_loc, ikp_loc)->c()[nmloc + ig] = fwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig];
	}
	mixer_[isp_loc][ikp_loc][ist_loc]->restart();
      }
    }
  }

  // cdupd = rho(t+dt)
  ChargeDensity cdupd(fwf);

  // efupd = H(t+dt)
  EnergyFunctional efupd(s_, cdupd);
  if ( s_.rtctrl.rt_vp == "ON" )
    efupd.vp->vp_propagate(iter_, rtdt_);

  // cdpre = rho(t+dt) previous step
  ChargeDensity cdpre(fwf2);
  tmap_["charge"].start();
  cdpre.update_density();
  tmap_["charge"].stop();

  itscf = 0;
  if ( onpe0 )
    cout << "<RTPTCNWavefunctionStepper : PT-CN iteration START : " << rtdt_ << " : " << rtitscf_ << " />" << endl;
 
  tmap_["rt_pt_cn"].start();
  while ( !scf_converged && itscf < rtitscf_ )
  {
    // Evaluete the electron density rho_f corresponding to psi_f
    tmap_["charge"].start();
    cdupd.update_density();
    tmap_["charge"].stop();

    // Update the potential and the Hamiltonian H_f
    tmap_["update_vhxc"].start();
    efupd.update_vhxc(false);
    tmap_["update_vhxc"].stop();
    efupd.energy(fwf, true, hfwf, false, fion, false, sigma);
 
    for ( int isp_loc = 0; isp_loc < rt_wf.nsp_loc(); ++isp_loc )
    {
      for ( int ikp_loc = 0; ikp_loc < rt_wf.nkp_loc(); ++ikp_loc )
      {
        rwf.sd(isp_loc, ikp_loc)->c().clear();
	int mloc = rt_wf.sd(isp_loc, ikp_loc)->c().mloc();
	int m = rt_wf.sd(isp_loc, ikp_loc)->c().m();
	for ( int ist_loc = 0; ist_loc < rt_wf.sd(isp_loc, ikp_loc)->nstloc(); ++ist_loc )
	{
	  int nmloc = ist_loc * mloc;
	  for ( int ig = 0; ig < mloc; ig++ )
	  {
	    // Evaluated the fixed point residual R_f = |psi_f> + (idt/2)*{H_f|psi_f> - |psi_f><psi_f|H_f|psi_f>} - |psi_n+1/2>
	    rwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] = fwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] + complex<double>(0.0, 0.5 * rtdt_) * (hfwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] - fwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] * conj(fwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig]) * hfwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig]) - xwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig];
	  }

	  //vector<complex<double> > fwf_tmp;
	  //vector<complex<double> > rwf_tmp;

	  vector<complex<double> > fwf_in;
	  vector<complex<double> > rwf_in;
	  vector<complex<double> > fwf_bar;
	  vector<complex<double> > rwf_bar;

	  //fwf_tmp.resize(m, complex<double>(0.0, 0.0));
	  //rwf_tmp.resize(m, complex<double>(0.0, 0.0));
	  //fwf_in.resize(m, complex<double>(0.0, 0.0));
	  //rwf_in.resize(m, complex<double>(0.0, 0.0));
	  //fwf_bar.resize(m);
	  //rwf_bar.resize(m);
	  ////////////////// PARALLELIZE /////////////////
	  fwf_in.resize(mloc, complex<double>(0.0, 0.0));
	  rwf_in.resize(mloc, complex<double>(0.0, 0.0));
	  fwf_bar.resize(mloc);
	  rwf_bar.resize(mloc);
	  ////////////////////////////////////////////////
	  for ( int ig = 0; ig < mloc; ig++ )
	  {
	  ////////////////// PARALLELIZE /////////////////
	    fwf_in[ig] = fwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig];
	    rwf_in[ig] = rwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig];
	  ////////////////////////////////////////////////
	  //int midx = rt_wf.sd(isp_loc, ikp_loc)->c().iglobal(ig);
	  //fwf_tmp[midx] = fwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig];
	  //rwf_tmp[midx] = rwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig];
	  }

	  MPI_Barrier(MPIdata::g_comm());
	  //MPI_Allreduce(&fwf_tmp[0], &fwf_in[0], m, MPI_DOUBLE_COMPLEX, MPI_SUM, MPIdata::g_comm());
	  //MPI_Allreduce(&rwf_tmp[0], &rwf_in[0], m, MPI_DOUBLE_COMPLEX, MPI_SUM, MPIdata::g_comm());

	  mixer_[isp_loc][ikp_loc][ist_loc]->update((double*)&fwf_in[0], (double*)&rwf_in[0], (double*)&fwf_bar[0], (double*)&rwf_bar[0]);

	  const double alpha = s_.rtctrl.rt_anderson_coeff;
	  ////////////////// PARALLELIZE /////////////////
	  for ( int ig = 0; ig < mloc; ig++ )
	  {
	    fwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] = fwf_bar[ig] + alpha * rwf_bar[ig];
	  }
	  ////////////////////////////////////////////////
	  //for ( int ig = 0; ig < m; ig++ )
          //  fwf_in[ig] = fwf_bar[ig] + alpha * rwf_bar[ig];

	  //for ( int ig = 0; ig < mloc; ig++ )
	  //{
	  //  int midx = rt_wf.sd(isp_loc, ikp_loc)->c().iglobal(ig);
	  //  fwf.sd(isp_loc, ikp_loc)->c()[nmloc + ig] = fwf_in[midx];
	  //}
	}
      }
    }

    // Evaluate the electron density rho_f corresponding to psi_f
    tmap_["charge"].start();
    cdupd.update_density();
    tmap_["charge"].stop();

    tmap_["charge"].start();
    cdpre.update_density();
    tmap_["charge"].stop();

    // If the change of the electron density is sufficiently small, exit the loop
    double tsum;
    tmap_["rt_pt_cn_conv"].start();
    scf_converged = conv_check(cdpre, cdupd, tsum);
    tmap_["rt_pt_cn_conv"].stop();

    if ( onpe0 )
      cout << " RTPTCN : " << itscf << " | DRHO : " << tsum << endl;
    if ( !scf_converged )
    {
      itscf = itscf + 1;
      fwf2 = fwf;
    }
  }
  tmap_["rt_pt_cn"].stop();
  // Orthogonalize |psi_f> to obtain |psi_n+1>
  fwf.gram();
  s_.wf = fwf;
  if ( onpe0 )
    cout << "<RTPTCNWavefunctionStepper : PT-CN iteration END : " << rtitscf_ << " />" << endl;
}

////////////////////////////////////////////////////////////////////////////////
bool RTPTCNWavefunctionStepper::conv_check(ChargeDensity& cd1, ChargeDensity& cd2, double& tsum)
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
	  cout << "RTPTCN RHO TEST : " << ssum << endl;
  if ( tsum < s_.rtctrl.rt_rho_tol )
  {
    return true;
  }
  else
  {
    return false;
  }
}


