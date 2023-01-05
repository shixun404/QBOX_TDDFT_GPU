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
// RTTDCmd.cpp
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#include "RTTDCmd.h"

#include "Wavefunction.h"

#include "SlaterDet.h"
#include "FourierTransform.h"
#include "Context.h"
#include "blas.h" // daxpy
#include "Base64Transcoder.h"
#include "Timer.h"
#include "MPIdata.h"

#include <iostream>
using namespace std;
#include "BOSampleStepper.h"
#include "RTSampleStepper.h"

#include <ctime>
#include <cassert>
#include <vector>
#include <valarray>
#include <complex>

int RTTDCmd::action(int argc, char **argv)
{
  if ( argc < 2 || argc > 7 )
  {
    if ( ui->onpe0() )
    {
      cout << " use: rttd [-atomic_density] rtiter [-auto-save rtas]" << endl;
      cout << "      rttd [-atomic_density] rtiter [-auto-save rtas] rtitscf" << endl;
      cout << "      rttd [-atomic_density] rtiter [-auto-save rtas] rtitscf rtite" << endl;
    }
    return 1;
  }

  if ( s->wf.nst() == 0 )
  {
    if ( ui->onpe0() )
      cout << " RTTDCmd: no states, cannot run" << endl;
    return 1;
  }
  if ( s->wf.ecut() == 0.0 )
  {
    if ( ui->onpe0() )
      cout << " RTTDCmd: ecut = 0.0, cannot run" << endl;
    return 1;
  }
  if ( s->wf.cell().volume() == 0.0 )
  {
    if ( ui->onpe0() )
      cout << " RTTDCmd: volume = 0.0, cannot run" << endl;
    return 1;
  }

  int iarg = 1;
  bool atomic_density = false;
  if ( !strcmp(argv[iarg],"-atomic_density") )
  {
    atomic_density = true;
    iarg++;
    argc--;
  }

  int rtiter = atoi(argv[iarg]);
  int rtitscf = 1;
  int rtite = 1;

  int rtas = 1;
  bool rt_autosave = false;
  if ( !strcmp(argv[iarg+1],"-auto-save") )
  {
    rt_autosave = true;
    rtas = atoi(argv[iarg+2]);
    iarg = iarg + 2;
    argc = argc - 2;
  }

  if ( argc == 3 )
  {
    // rttd niter nitscf
    rtitscf = atoi(argv[iarg+1]);
  }
  else if ( argc == 4 )
  {
    // rttd niter nitscf nite
    rtitscf = atoi(argv[iarg+1]);
    rtite = atoi(argv[iarg+2]);
  }

  s->ctrl.dt = s->rtctrl.rt_dt;
  s->rtctrl.rt_as = rtas;
  if ( ui->onpe0() )  cout << " RTTDCmd : " << rtiter << " | " << rtitscf << " | " << rtite << " | " << rtas << endl;

  s->extforces.setup(s->atoms);

  //if ( ui->onpe0() )
  //  cout << " RTTDCmd nempty before : " << s->wf.nempty() << endl;

  //if ( s->wf.nempty() == 0 )
  //{
  //  s->wf.set_nempty(max((int)ceil(s->wf.nst() * 0.5),4));
  //  if ( s->wfv != 0 )
  //  {
  //    s->wfv->set_nempty(max((int)ceil(s->wf.nst() * 0.5),4));
  //  }
  //}

  //if ( ui->onpe0() )
  //  cout << " RTTDCmd nempty after : " << s->wf.nempty() << endl;

  SampleStepper* stepper;

  stepper = new BOSampleStepper(*s,rtitscf,rtite);

  assert(stepper!=0);
  stepper->set_iter_cmd(s->ctrl.iter_cmd);
  stepper->set_iter_cmd_period(s->ctrl.iter_cmd_period);

  //if ( atomic_density )
  stepper->initialize_density();

//ADRIAN ADDED JUST FOR TESTING  

  s->wf.info(cout,"wavefunction");
  stepper->step(0);

//ADRIAN ADDED JUST FOR TESTING

  if ( s->wfv != 0) delete s->wfv;
  s->wfv = 0;

  delete stepper;

  SampleStepper* rtstepper;

  rtstepper = new RTSampleStepper(*s, rtitscf, rtite, rtas);

  rtstepper->set_iter_cmd(s->ctrl.iter_cmd);
  rtstepper->set_iter_cmd_period(s->ctrl.iter_cmd_period);

  if ( ui->onpe0() )  cout << "RTTDCmd TEST1" << endl;
  s->wf.info(cout,"wavefunction");
  rtstepper->step(rtiter);
  if ( ui->onpe0() )  cout << "RTTDCmd TEST2" << endl;

  delete rtstepper;
  if ( ui->onpe0() )  cout << "RTTDCmd TEST3" << endl;

  return 0;
}
