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
// RTProjKp.h
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RTPROJKP_H
#define RTPROJKP_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>

#include "Sample.h"

class RTProjKp : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "rt_proj_kp"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 3 )
    {
      if ( ui->onpe0() )
      cout << " rt_proj_kp must be specified with 2 integer values" << endl;
      return 1;
    }

    int v0 = atoi(argv[1]);
    int v1 = atoi(argv[2]);

    s->rtctrl.rt_proj_kp[0] = v0;
    s->rtctrl.rt_proj_kp[1] = v1;
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) 
	<< s->rtctrl.rt_proj_kp[0] << " " 
	<< s->rtctrl.rt_proj_kp[1] << " ";
     return st.str();
  }

  RTProjKp(Sample *sample) : s(sample) 
  { 
    s->rtctrl.rt_proj_kp[0] = 0; 
    s->rtctrl.rt_proj_kp[1] = 0; 
  }
};
#endif




