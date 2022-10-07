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
// HelpCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef HELPCMD_H
#define HELPCMD_H

#include <iostream>
#include <iomanip>
#include <string.h>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class HelpCmd : public Cmd
{
  public:

  Sample *s;

  HelpCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "help"; }
  const char *help_msg(void) const
  {
    return
    "\n help\n\n"
    " syntax: help [command_name]\n\n"
    "   The help command gives a short description of a command. If used\n"
    "   without arguments, help prints a list of valid commands.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( ui->onpe0() )
    {
      if ( argc == 1 )  // no arguments
      {
        cout << endl << " valid commands are: " << endl << endl;
        list<Cmd*>::iterator cmd = ui->cmdlist.begin();
        int n = 0;
        while ( cmd != ui->cmdlist.end() )
        {
          n++;
          cout.setf(ios::left,ios::adjustfield);
          cout << " " << setw(16) << (*cmd)->name();
          cout.setf(ios::right,ios::adjustfield);
          if ( n%4 == 0 ) cout << endl;
          cmd++;
        }
        if ( n%4 != 0 ) cout << endl;
        cout << endl;
      }
      else if ( argc == 2 ) // one argument
      {
        // search command list
        Cmd *cmdptr = ui->findCmd(argv[1]);

        if ( cmdptr )
        {
          cout << cmdptr->help_msg();
        }
        else
        {
          cout << " help: " << argv[1] << " is not a valid command" << endl;
        }
      }
      else
      {
        cout << " use: help [command_name]" << endl;
      }
    }
    return 0;
  }
};
#endif
