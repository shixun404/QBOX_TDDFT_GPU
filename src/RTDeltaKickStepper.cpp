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
// RTDeltaKickStepper.cpp
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#include "RTDeltaKickStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Sample.h"
#include "Basis.h"
#include "FourierTransform.h"
#include "MPIdata.h"
#include "ChargeDensity.h"

#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
RTDeltaKickStepper::RTDeltaKickStepper(Sample& s, RTPosition& rtp) :
  s_(s), rtp_(rtp)
{
//if ( MPIdata::onpe0() )
//  cout << "RTDKS GRID : " << rtp_.ft_wfc()->np0() << " " << rtp_.ft_wfc()->np1() << " " << rtp_.ft_wfc()->np2() << endl;
}

////////////////////////////////////////////////////////////////////////////////
RTDeltaKickStepper::~RTDeltaKickStepper(void)
{
}

////////////////////////////////////////////////////////////////////////////////
void RTDeltaKickStepper::update_dk(D3vector deltakick_amp)
{
  double sol = 137.03599911;
  double alpha = 1.0/sol;
  for ( int isp_loc = 0; isp_loc < s_.wf.nsp_loc(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < s_.wf.nkp_loc(); ++ikp_loc )
    {
      vector<complex<double> > wfr(rtp_.ft_dwfc(ikp_loc)->np012loc());
      int nstloc = s_.wf.sd(isp_loc, ikp_loc)->nstloc();
      int mloc = s_.wf.sd(isp_loc, ikp_loc)->c().mloc();
      for ( int ist_loc = 0; ist_loc < nstloc; ++ist_loc )
      {
        // psi = psi * exp(iK0.dot(r))
        rtp_.ft_dwfc(ikp_loc)->backward(s_.wf.sd(isp_loc, ikp_loc)->c().cvalptr(ist_loc*mloc), &wfr[0]);
        for ( int ir = 0; ir < rtp_.ft_dwfc(ikp_loc)->np012loc(); ir++ )
	{
          double phase = - deltakick_amp.x * rtp_.rt_position_dwfc(ikp_loc, ir).x - deltakick_amp.y * rtp_.rt_position_dwfc(ikp_loc, ir).y - deltakick_amp.z * rtp_.rt_position_dwfc(ikp_loc, ir).z;
	  phase = phase * alpha;
          wfr[ir] *= complex<double>(cos(phase), sin(phase));
	}
	rtp_.ft_dwfc(ikp_loc)->forward(&wfr[0], s_.wf.sd(isp_loc, ikp_loc)->c().valptr(ist_loc*mloc));
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void RTDeltaKickStepper::compute_dip(ChargeDensity& cd)
{
  D3vector dipole = D3vector(0.0, 0.0, 0.0);
  rt_dipole_ = D3vector(0.0, 0.0, 0.0);
  const double omega = cd.vbasis()->cell().volume();

  for ( int isp_loc = 0; isp_loc < s_.wf.nsp_loc(); ++isp_loc )
  {
    for ( int ir = 0; ir < rtp_.ft_rho()->np012loc(); ir++ )
    {
      dipole.x = dipole.x + rtp_.rt_position_rho(ir).x * cd.rhor[isp_loc][ir];
      dipole.y = dipole.y + rtp_.rt_position_rho(ir).y * cd.rhor[isp_loc][ir];
      dipole.z = dipole.z + rtp_.rt_position_rho(ir).z * cd.rhor[isp_loc][ir];
    }
  }
  MPI_Allreduce(&dipole.x, &rt_dipole_.x, 1, MPI_DOUBLE, MPI_SUM, MPIdata::g_comm());
  MPI_Allreduce(&dipole.y, &rt_dipole_.y, 1, MPI_DOUBLE, MPI_SUM, MPIdata::g_comm());
  MPI_Allreduce(&dipole.z, &rt_dipole_.z, 1, MPI_DOUBLE, MPI_SUM, MPIdata::g_comm());
  rt_dipole_.x = rt_dipole_.x * omega / rtp_.ft_rho()->np012();
  rt_dipole_.y = rt_dipole_.y * omega / rtp_.ft_rho()->np012();
  rt_dipole_.z = rt_dipole_.z * omega / rtp_.ft_rho()->np012();
}  
