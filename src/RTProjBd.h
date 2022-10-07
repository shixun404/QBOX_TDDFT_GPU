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
// RTProjBd.h
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RTPROJBD_H
#define RTPROJBD_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>

#include "Sample.h"

class RTProjBd : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "rt_proj_bd"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 3 )
    {
      if ( ui->onpe0() )
      cout << " rt_proj_bd must be specified with 2 integer values (start, end)" << endl;
      return 1;
    }

    int v0 = atoi(argv[1]);
    int v1 = atoi(argv[2]);

    if ( v0 > v1 )
    {
      if ( ui->onpe0() )
        cout << " rt_proj_bd start value must be less or equal than end value. arg[1] <= arg[2]" << endl;
      return 1;
    }

    s->rtctrl.rt_proj_bd[0] = v0;
    s->rtctrl.rt_proj_bd[1] = v1;
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) 
	<< s->rtctrl.rt_proj_bd[0] << " " 
	<< s->rtctrl.rt_proj_bd[1] << " ";
     return st.str();
  }

  RTProjBd(Sample *sample) : s(sample) 
  { 
    s->rtctrl.rt_proj_bd[0] = 0; 
    s->rtctrl.rt_proj_bd[1] = 0; 
  }
};
#endif



