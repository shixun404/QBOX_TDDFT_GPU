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
// RTSampleStepper.cpp
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#include "RTSampleStepper.h"
#include "EnergyFunctional.h"
#include "SlaterDet.h"
#include "Basis.h"
#include "RTCNWavefunctionStepper.h"
#include "RTPTCNWavefunctionStepper.h"
#include "WavefunctionStepper.h"
#include "UserInterface.h"
#include "Preconditioner.h"
#include "MDIonicStepper.h"
#include "D3tensor.h"
#include "cout0.h"
#include "NonLocalPotential.h"
#include "RTDeltaKickStepper.h"
#include "RTVectorPotential.h"
#include "RTCurrentDensity.h"
#include "RTProjection.h"

#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
RTSampleStepper::RTSampleStepper(Sample& s, int rtitscf, int rtite, int rtas) :
  SampleStepper(s), cd_(s.wf), ef_(s,cd_), wf_(s.wf), wf_prev_(s.wf),
  wf_next_(s.wf), wf_init_(s.wf),
  dwf(s.wf), rtitscf_(rtitscf), rtite_(rtite), rtas_(rtas),
  update_density_first_(true), update_vh_(true), update_vxc_(true) {}

////////////////////////////////////////////////////////////////////////////////
RTSampleStepper::~RTSampleStepper()
{
  for ( TimerMap::iterator i = tmap.begin(); i != tmap.end(); i++ )
  {
    double time = i->second.real();
    double tmin, tmax;
    MPI_Reduce(&time,&tmin,1,MPI_DOUBLE,MPI_MIN,0,MPIdata::comm());
    MPI_Reduce(&time,&tmax,1,MPI_DOUBLE,MPI_MAX,0,MPIdata::comm());
    if ( MPIdata::onpe0() && (tmax > 0.0) )
    {
      string s = "name=\"" + i->first + "\"";
      cout << "<timing " << left << setw(22) << s
           << " min=\"" << setprecision(3) << tmin << "\""
           << " max=\"" << setprecision(3) << tmax << "\"/>"
           << endl;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void RTSampleStepper::step(int rtiter)
{
  const bool compute_eigvec = false;
  const bool onpe0 = MPIdata::onpe0();
  const double gpa = 29421.5;
  double autoas = 24.18884254;
  if ( onpe0 )
    cout << "RTSampleStepper : step START" << endl;

  const double rtdt = s_.rtctrl.rt_dt;

  const string propagator = s_.rtctrl.rt_propagator;
  const string rtmd = s_.rtctrl.rt_md;
  const string rtdeltakick = s_.rtctrl.rt_delta_kick;
  const string rtvp = s_.rtctrl.rt_vp;
  const string rtvpeq = s_.rtctrl.rt_vp_eq;
  const string rtproj = s_.rtctrl.rt_proj;

  AtomSet& atoms = s_.atoms;
  wf_init_ = s_.wf;
  wf_ = s_.wf;
  const int nspin = s_.wf.nspin();
  const Context& sd_ctxt = s_.wf.sd_context();

  const bool atoms_move = ( rtiter > 0 && rtmd == "ON" );
  const bool compute_stress = ( s_.ctrl.stress == "ON" );
  const double force_tol = s_.ctrl.force_tol;
  const double stress_tol = s_.ctrl.stress_tol;

  Timer tm_iter;

  if ( onpe0 ) cout << "RTSampleStepper : RTPosition Preparing START" << endl;
  RTPosition rtp = RTPosition(s_.wf, atoms);
  if ( onpe0 ) cout << "RTSampleStepper : RTPosition Preparing END" << endl;

  RTDeltaKickStepper* wf_deltakick = 0;
  if ( rtdeltakick == "ON" || rtvpeq == "DELTA" )
  {
    if ( onpe0 ) cout << "RTSampleStepper : RTDeltaKickStepper Preparing START" << endl;
    wf_deltakick = new RTDeltaKickStepper(s_, rtp);
    if ( rtdeltakick == "ON" )
      wf_deltakick->update_dk(s_.rtctrl.rt_delta_kick_amp);
    if ( rtvpeq == "DELTA" )
    {
      wf_deltakick->update_dk(s_.rtctrl.rt_vp_amp);
    //ef_.update_nl_dk(s_.rtctrl.rt_vp_amp, rtp);
    }
    if ( onpe0 ) cout << "RTSampleStepper : RTDeltaKickStepper Preparing END" << endl;
  }

  RTCurrentDensity* rtc = 0;
  if ( rtvp == "ON" )
  {
    if ( onpe0 ) cout << "RTSampleStepper : RTCurrentDensity Preparing START" << endl;
    rtc = new RTCurrentDensity(s_, rtp, atoms);
    if ( onpe0 ) cout << "RTSampleStepper : RTCurrentDensity Preparing END" << endl;
  }

  RTProjection* rtpj = 0;
  if ( rtproj == "ON" )
  {
    if ( onpe0 ) cout << "RTSampleStepper : RTProjection Preparing START" << endl;
    rtpj = new RTProjection(s_.wf, s_.rtctrl.rt_proj_type, s_.rtctrl.rt_proj_bd[0], s_.rtctrl.rt_proj_bd[1], s_.rtctrl.rt_proj_kp[0], s_.rtctrl.rt_proj_kp[1]);
    if ( onpe0 ) cout << "RTSampleStepper : RTProjection Preparing 1" << endl;
    rtpj->set_collection();
    if ( onpe0 ) cout << "RTSampleStepper : RTProjection Preparing 2" << endl;
    rtpj->start_collection(wf_init_);
    if ( onpe0 ) cout << "RTSampleStepper : RTProjection Preparing END" << endl;
  }

  Preconditioner prec(s_.wf, ef_, s_.ctrl.ecutprec);

  WavefunctionStepper* wf_stepper = 0;
  if ( propagator == "CN" ) wf_stepper = new RTCNWavefunctionStepper(rtdt, tmap, ef_, cd_, s_, rtitscf_, rtite_);
  if ( propagator == "PTCN" ) wf_stepper = new RTPTCNWavefunctionStepper(rtdt, tmap, ef_, cd_, s_, rtitscf_);
//else if ( propagator == "ST2" )
//  wf_stepper = new RTSTSODWavefunctionStepper(wf,prec,tmap);
//else if ( propagator == "ST4" )
//  wf_stepper = new RTSTFOWavefunctionStepper(wf,prec,tmap);

  IonicStepper* ionic_stepper = 0;
  if ( rtmd == "ON" ) ionic_stepper = new MDIonicStepper(s_);
  if ( ionic_stepper ) ionic_stepper->setup_constraints();

///////////////////////////////////////////////ITERATION//////////////////////////////////////////
  for ( int iter = 0; iter < max(rtiter,1); ++iter )
  {
    tm_iter.start();

    double asdt = iter * rtdt * autoas;

    if ( onpe0 )
      cout << "<rt-iteration count=\"" << iter << "\">\n";

    // compute energy and ionic forces using existing wavefunction
    double maxforce = 0.0;
    double maxstress = 0.0;

    tmap["charge"].start();
    cd_.update_density();
    tmap["charge"].stop();
    tmap["update_vhxc"].start();
    ef_.update_vhxc(compute_stress);
    tmap["update_vhxc"].stop();
    const bool compute_forces = true;
    tmap["energy"].start();
    double energy = ef_.energy(false,dwf,compute_forces,fion,compute_stress,sigma_eks);
    tmap["energy"].stop();
    double enthalpy = ef_.enthalpy();

    if ( rtvp == "ON" )
    {
      ef_.vp->vp_propagate(iter, rtdt);
      if ( s_.rtctrl.rt_vp_ind == "T" )
      {
	D3vector curr;
	curr.x = rtc->obs(0, 0).real();
	curr.y = rtc->obs(0, 1).real();
	curr.z = rtc->obs(0, 2).real();
	ef_.vp->calculate_acceleration(rtdt, curr, atoms.cell());
      }
      rtc->compute_current(s_.wf, ef_);
      if ( onpe0 )
      {
        cout << "RTVectorPotential " << setprecision(5) << asdt << setprecision(16) << "  " << ef_.vp->vecpot() << endl;
        for ( int ispin = 0; ispin < s_.wf.nspin(); ispin++ )
        {
          cout << "RTCurrentDensity : " << ispin << "  " << setprecision(5) << asdt << setprecision(16) << "   " << rtc->obs(ispin, 0).real() << "  " << rtc->obs(ispin, 1).real() << "  " << rtc->obs(ispin, 2).real() << endl;
          cout << "RTCurrentPRDensity : " << ispin << "  " << setprecision(5) << asdt << setprecision(16) << "   " << rtc->pobs(ispin, 0).real() << "  " << rtc->pobs(ispin, 1).real() << "  " << rtc->pobs(ispin, 2).real() << endl;
        //cout << "RTCurrentPIDensity : " << ispin << "  " << setprecision(5) << asdt << setprecision(16) << "   " << rtc->pobs(ispin, 0).imag() << "  " << rtc->pobs(ispin, 1).imag() << "  " << rtc->pobs(ispin, 2).imag() << endl;
        //cout << "RTCurrentVRDensity : " << ispin << "  " << setprecision(5) << asdt << setprecision(16) << "   " << rtc->vobs(ispin, 0).real() << "  " << rtc->vobs(ispin, 1).real() << "  " << rtc->vobs(ispin, 2).real() << endl;
          cout << "RTCurrentVIDensity : " << ispin << "  " << setprecision(5) << asdt << setprecision(16) << "   " << rtc->vobs(ispin, 0).imag() << "  " << rtc->vobs(ispin, 1).imag() << "  " << rtc->vobs(ispin, 2).imag() << endl;
          cout << "RTCurrentHRDensity : " << ispin << "  " << setprecision(5) << asdt << setprecision(16) << "   " << rtc->hobs(ispin, 0).real() << "  " << rtc->hobs(ispin, 1).real() << "  " << rtc->hobs(ispin, 2).real() << endl;
          cout << "RTCurrentHIDensity : " << ispin << "  " << setprecision(5) << asdt << setprecision(16) << "   " << rtc->hobs(ispin, 0).imag() << "  " << rtc->hobs(ispin, 1).imag() << "  " << rtc->hobs(ispin, 2).imag() << endl;
        }
      }
    }

    if ( rtproj == "ON" )
    {
      rtpj->compute_projection(s_.wf);
      if ( onpe0 )
      {
        for ( int ispin = 0; ispin < s_.wf.nspin(); ispin++ )
        {
          cout << "RTProjection " << ispin << " / " << s_.wf.nspin()-1 << " time : " << setprecision(5) << asdt << endl;
          cout << "RTProjection " << "          ";
          for ( int ibnd = s_.rtctrl.rt_proj_bd[0]; ibnd <= s_.rtctrl.rt_proj_bd[1]; ibnd++ )
            cout << ibnd << "             ";
          cout << endl;
          int ibnd_idx = 0;
          for ( int ibnd = s_.rtctrl.rt_proj_bd[0]; ibnd <= s_.rtctrl.rt_proj_bd[1]; ibnd++ )
          {
            cout << "RTProjection  " << ibnd << "  "; 
            int jbnd_idx = 0;
            for ( int jbnd = s_.rtctrl.rt_proj_bd[0]; jbnd <= s_.rtctrl.rt_proj_bd[1]; jbnd++ )
            {
               cout << setprecision(8) << rtpj->proj_band(ispin, ibnd_idx, jbnd_idx) << "  ";
               jbnd_idx = jbnd_idx + 1;
            }
            ibnd_idx = ibnd_idx + 1;
            cout << endl;
          }
        }
      }
    }

    if ( rtdeltakick == "ON" )
    {
      wf_deltakick->compute_dip(cd_);
      if ( onpe0 )
        cout << " <rt_delta_dip> " << setprecision(18) << wf_deltakick->rt_dipole().x 
                           << "  " << setprecision(18) << wf_deltakick->rt_dipole().y 
      		           << "  " << setprecision(18) << wf_deltakick->rt_dipole().z 
      		           << " </rt_delta_dip> " << endl;
    }

    if ( force_tol > 0.0 )
    {
      maxforce = 0.0;
      for ( int is = 0; is < fion.size(); is++ )
        for ( int i = 0; i < fion[is].size(); i++ )
          maxforce = max(maxforce, fabs(fion[is][i]));
    }

    if ( stress_tol > 0.0 )
    {
      compute_sigma();
      for ( int i = 0; i < sigma.size(); i++ )
        maxstress = max(maxstress, gpa*fabs(sigma[i]));
    }
   
    if ( (iter%rtas_ == 0) && onpe0 )
    {
      cout << cd_;
      cout << ef_;
      if ( ef_.el_enth() )
        cout << *ef_.el_enth();
    }

    if ( iter > 0 && ionic_stepper )
    {
      ionic_stepper->compute_v(energy,fion);
    }
    double ekin_ion = 0.0, temp_ion = 0.0;
    if ( ionic_stepper )
    {
      ekin_ion = ionic_stepper->ekin();
      temp_ion = ionic_stepper->temp();
    }

    if ( (iter%rtas_ == 0) && onpe0 )
    {
      cout << "<atomset>" << endl;
      cout << atoms.cell();
      for ( int is = 0; is < atoms.atom_list.size(); is++ )
      {
        int i = 0;
        for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ )
        {
          Atom* pa = atoms.atom_list[is][ia];
          cout << "  <atom name=\"" << pa->name() << "\""
               << " species=\"" << pa->species()
               << "\">\n"
               << "    <position> " << pa->position() << " </position>\n"
               << "    <velocity> " << pa->velocity() << " </velocity>\n"
               << "    <force> "
               << fion[is][i] << " "
               << fion[is][i+1] << " "
               << fion[is][i+2]
               << " </force>\n";
          cout << "  </atom>" << endl;
          i += 3;
        }
      }
      cout << "</atomset>" << endl;
      cout << setprecision(6);
      cout << "<unit_cell_a_norm> " << atoms.cell().a_norm(0)
           << " </unit_cell_a_norm>" << endl;
      cout << "<unit_cell_b_norm> " << atoms.cell().a_norm(1)
           << " </unit_cell_b_norm>" << endl;
      cout << "<unit_cell_c_norm> " << atoms.cell().a_norm(2)
           << " </unit_cell_c_norm>" << endl;
      cout << setprecision(3) << "<unit_cell_alpha>  "
           << atoms.cell().alpha() << " </unit_cell_alpha>" << endl;
      cout << setprecision(3) << "<unit_cell_beta>   "
           << atoms.cell().beta() << " </unit_cell_beta>" << endl;
      cout << setprecision(3) << "<unit_cell_gamma>  "
           << atoms.cell().gamma() << " </unit_cell_gamma>" << endl;
      cout << setprecision(3) << "<unit_cell_volume> "
           << atoms.cell().volume() << " </unit_cell_volume>" << endl;

      // include the kinetic energy of the stepper
      // e.g. to include thermostat contributions
      double ekin_stepper;
      if ( ionic_stepper != 0 )
        ekin_stepper = ionic_stepper->ekin_stepper();
      cout << setprecision(8);
      cout << "  <econst> " << enthalpy+ekin_ion+ekin_stepper
           << " </econst>\n";
      cout << "  <ekin_ion> " << ekin_ion << " </ekin_ion>\n";
      cout << "  <temp_ion> " << temp_ion << " </temp_ion>\n";
    }

    if ( atoms_move )
    {
      if ( s_.constraints.size() > 0 )
      {
        s_.constraints.update_constraints(rtdt);
        s_.constraints.compute_forces(ionic_stepper->r0(), fion);
        if ( (iter%rtas_ == 0) && onpe0 )
        {
          s_.constraints.list_constraints(cout);
        }
      }
      ionic_stepper->compute_r(energy,fion);
      ef_.atoms_moved();
    }

    if ( compute_stress )
    {
      compute_sigma();
      if ( (iter%rtas_ == 0) )
        print_stress();
    }

    if ( wf_stepper != 0 )
    {
      wf_stepper->preprocess();
      if ( propagator == "CN" )
      {
        if ( onpe0 ) 
          cout << " RTSampleStpper : Start CN propagator " << iter << endl;
        if ( iter == 0 )
	{
	  // wf_prev = |psi(t=0)>
	  wf_prev_ = s_.wf;
	  wf_next_ = s_.wf;
	  //
	}
	else
	{
	  // wf_next = |psi(t)>
	  wf_next_ = s_.wf;

          // wf_next = 2|psi(t)>-|psi(t-dt)>
	  for ( int isp_loc = 0; isp_loc < s_.wf.nsp_loc(); ++isp_loc )
	  {
	    for ( int ikp_loc = 0; ikp_loc < s_.wf.nkp_loc(); ++ikp_loc )
	    {
              for ( int ist_loc = 0; ist_loc < s_.wf.sd(isp_loc, ikp_loc)->nstloc(); ++ist_loc )
	      {
	        int mloc = s_.wf.sd(isp_loc, ikp_loc)->c().mloc();
	        int nmloc = ist_loc * mloc;
	        for ( int ig = 0; ig < mloc; ig++ )
	        {
		  wf_next_.sd(isp_loc, ikp_loc)->c()[nmloc + ig] = complex<double>(2.0, 0.0) * wf_next_.sd(isp_loc, ikp_loc)->c()[nmloc + ig] - wf_prev_.sd(isp_loc, ikp_loc)->c()[nmloc + ig];
	        }
	      }
	    }
	  }
          // wf_prev = |psi(t)>
	  wf_prev_ = s_.wf;
	}
	wf_stepper->get_iter(iter);
	wf_stepper->update(wf_next_);
        if ( onpe0 ) 
          cout << " RTSampleStpper : End CN propagator " << endl;
	wf_ = s_.wf;
      }
      else if ( propagator == "PTCN" )
      {
        if ( onpe0 ) 
          cout << " RTSampleStpper : Start PTCN propagator " << iter << endl;
	wf_next_ = s_.wf;
	wf_stepper->get_iter(iter);
	wf_stepper->update(wf_next_);
        if ( onpe0 ) 
          cout << " RTSampleStpper : End PTCN propagator " << endl;
	wf_ = s_.wf;
      }
    }
    double ehart, ehart_m;
    double etotal = 0.0, etotal_m = 0.0, etotal_mm = 0.0;

    tmap["charge"].start();
    cd_.update_density();
    tmap["charge"].stop();
    tmap["update_vhxc"].start();
    ef_.update_vhxc(false);
    tmap["update_vhxc"].stop();
    tmap["energy"].start();
    ef_.energy(wf_, true, dwf, false, fion, false, sigma_eks);
    tmap["energy"].stop();

    double eigenvalue_sum = 0.0;
    double etotal_int;
    eigenvalue_sum = real(wf_.dot(dwf));
    double fac = ( wf_.nspin() == 1 ) ? 2.0 : 1.0 ;
    double w_eigenvalue_sum = fac * eigenvalue_sum;
    if ( (iter%rtas_ == 0) && onpe0 )
      cout << "  <rt_eigenvalue_sum>  " << setprecision(8)
           << eigenvalue_sum << " </rt_eigenvalue_sum>" << endl;
    wf_.diag(dwf, compute_eigvec);

    // print eigenvalues
    if ( (iter%rtas_ == 0) && onpe0 )
      cout << "<rt_eigenset>" << endl;
    if ( iter%rtas_ == 0 )
    {
      for ( int ispin = 0; ispin < wf_.nspin(); ++ispin )
      {
        const int isp_loc = wf_.isp_local(ispin);
        for ( int ikp = 0; ikp < wf_.nkp(); ++ikp )
        {
          const int ikp_loc = wf_.ikp_local(ikp);
          ostringstream ostr;
          ostr.setf(ios::fixed,ios::floatfield);
          ostr.setf(ios::right,ios::adjustfield);
          int isrc = -1;
          if ( ( isp_loc >= 0 ) && ( ikp_loc >= 0 ) )
          {
            if ( MPIdata::sd_rank() == 0 )
            {
              ostr.str("");
              isrc = MPIdata::rank();
              const int nst = wf_.sd(isp_loc,ikp_loc)->nst();
              const double eVolt = 2.0 * 13.6058;
              ostr <<    "  <rt_eigenvalues spin=\"" << ispin
                   << "\" kpoint=\""
                   << setprecision(8)
                   << wf_.sd(isp_loc,ikp_loc)->kpoint()
                   << "\" weight=\""
                   << setprecision(8)
                   << wf_.weight(ikp)
                   << "\" n=\"" << nst << "\">" << endl;
              for ( int i = 0; i < nst; i++ )
              {
                ostr << setw(12) << setprecision(5)
                     << wf_.sd(isp_loc,ikp_loc)->eig(i)*eVolt;
                if ( i%5 == 4 ) ostr << endl;
              }
              if ( nst%5 != 0 ) ostr << endl;
              ostr << "  </rt_eigenvalues>" << endl;
            }
          }
          cout0(ostr.str(),isrc);
          MPI_Barrier(MPIdata::comm());
        }
      }
    }
    if ( (iter%rtas_ == 0) && onpe0 )
      cout << "</rt_eigenset>" << endl;
    MPI_Barrier(MPIdata::comm());
    //

    w_eigenvalue_sum = 0.0;
    for ( int isp_loc = 0; isp_loc < wf_.nsp_loc(); ++isp_loc )
    {
      for ( int ikp_loc = 0; ikp_loc < wf_.nkp_loc(); ++ikp_loc )
      {
        const int nst = wf_.sd(isp_loc,ikp_loc)->nst();
        const int ikpg = wf_.ikp_global(ikp_loc);
        const double wkp = wf_.weight(ikpg);
        for ( int n = 0; n < nst; n++ )
        {
          const double occ = wf_.sd(isp_loc,ikp_loc)->occ(n);
          w_eigenvalue_sum += wkp * occ * wf_.sd(isp_loc,ikp_loc)->eig(n);
        }
      }
    }
    double tsum = 0.0;
    MPI_Allreduce(&w_eigenvalue_sum, &tsum, 1, MPI_DOUBLE, MPI_SUM, MPIdata::kp_sp_comm());
    w_eigenvalue_sum = tsum;
    
    // Harris-Foulkes estimate of the total energy
    etotal_int = w_eigenvalue_sum - ef_.ehart_e() + ef_.ehart_p() +
                 ef_.esr() - ef_.eself() + ef_.dxc() + ef_.ets();
    if ( (iter%rtas_ == 0) && onpe0 )
    {
      cout << setprecision(8);
      cout << "rt_w_eigenvalue_sum = " << setw(15) << w_eigenvalue_sum << endl;
      cout << "rt_ef.dxc()         = " << setw(15) << ef_.dxc() << endl;
      cout << "rt_ef.ehart()       = " << setw(15) << ef_.ehart() << endl;
      cout << "rt_ef.ehart_e()     = " << setw(15) << ef_.ehart_e() << endl;
      cout << "rt_ef.ehart_ep()    = " << setw(15) << ef_.ehart_ep() << endl;
      cout << "rt_ef.esr()         = " << setw(15) << ef_.esr() << endl;
      cout.setf(ios::fixed,ios::floatfield);
      cout.setf(ios::right,ios::adjustfield);
      cout << "  <rt_etotal_int>  " << asdt << " as " <<  setprecision(15) << setw(18)
           << etotal_int << " </rt_etotal_int>\n";
    }

    if ( force_tol > 0.0 )
    {
      if ( (iter%rtas_ == 0) && onpe0 )
        cout << "  maxforce: " << scientific
             << setprecision(4) << maxforce << fixed << endl;
    }
    if ( stress_tol > 0.0 )
    {
      if ( (iter%rtas_ == 0) && onpe0 )
        cout << "  maxstress: " << scientific
             << setprecision(4) << maxstress << fixed << endl;
    }

    double time = tm_iter.real();
    double tmin, tmax;
    MPI_Reduce(&time,&tmin,1,MPI_DOUBLE,MPI_MIN,0,MPIdata::comm());
    MPI_Reduce(&time,&tmax,1,MPI_DOUBLE,MPI_MAX,0,MPIdata::comm());
    if ( (iter%rtas_ == 0) && onpe0 )
    {
      string s = "name=\"rt-iteration\"";
      cout << "<timing " << left << setw(22) << s
           << " min=\"" << setprecision(3) << tmin << "\""
           << " max=\"" << setprecision(3) << tmax << "\"/>"
           << endl;
    }

    if ( (iter%rtas_ == 0) && onpe0 )
      cout << "</rt-iteration>" << endl;
  } // for iter

  if ( atoms_move )
  {
    // compute ionic forces at last position to update velocities
    // consistently with last position
    //tmap["charge"].start();
    cd_.update_density();
    //tmap["charge"].stop();

    //tmap["update_vhxc"].start();
    ef_.update_vhxc(compute_stress);
    //tmap["update_vhxc"].stop();
    const bool compute_forces = true;
    //tmap["energy"].start();
    double energy = ef_.energy(false,dwf,compute_forces,fion,compute_stress,sigma_eks);
    //tmap["energy"].stop();

    ionic_stepper->compute_v(energy,fion);
    // positions r0 and velocities v0 are consistent
  }

  // delete steppers
  if ( onpe0 ) cout << "RTSAMPLESTEPPER_TEST1" << endl;
  delete wf_stepper;
  if ( onpe0 ) cout << "RTSAMPLESTEPPER_TEST2" << endl;
  delete ionic_stepper;
  if ( onpe0 ) cout << "RTSAMPLESTEPPER_TEST3" << endl;
  delete rtc;
  if ( onpe0 ) cout << "RTSAMPLESTEPPER_TEST4" << endl;
  delete wf_deltakick;
  if ( onpe0 ) cout << "RTSAMPLESTEPPER_TEST5" << endl;
  delete rtpj;
  if ( onpe0 ) cout << "RTSAMPLESTEPPER_END" << endl;

}
