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
// RTTDCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RTTDCMD_H
#define RTTDCMD_H

#include <iostream>

#include "UserInterface.h"

class Sample;
class RTTDCmd : public Cmd
{
  private:

  public:

  Sample *s;

  RTTDCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "rttd"; }
  const char *help_msg(void) const
  {
    return
    "\n rttd\n\n";
    " syntax: rttd [-atomic_density] rtiter [-auto-save rtas] [rtitscf [rtite]]\n\n"
    "   The rttd command runs rtiter steps of RT-TDDFT simulation. Each step\n"
    "   consists of one or more (rtitscf) scf steps, each consisting\n"
    "   of one or more (rtite) electronic steps.\n"
    "   If the -atomic_density option is used, the initial charge\n"
    "   density is a sum of atomic charge densities.\n"
    "   If the -auto-save option is used, the wavefunction is saved for\n"
    "   each rtas step.\n\n";
  }

  int action(int argc, char **argv);

};
#endif
