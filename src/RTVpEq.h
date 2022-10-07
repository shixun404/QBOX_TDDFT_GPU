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
// RTVpEq.h
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RTVPEQ_H
#define RTVPEQ_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>

#include "Sample.h"

class RTVpEq : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "rt_vp_eq"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " rt_vp_eq takes only one value" << endl;
      return 1;
    }

    string v = argv[1];
    if ( !( v == "AC" || v == "GAUSSIAN" || v == "PULSE" || v == "CIRPOLAR" || v == "CIRPULSE" || v == "DELTA" ) )
    {
      if ( ui->onpe0() )
        cout << " rt_vp_eq must be AC or GAUSSIAN or PULSE or CIRPOLAR or CIRPULSE or DELTA" << endl;
      return 1;
    }

    s->rtctrl.rt_vp_eq = v;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->rtctrl.rt_vp_eq;
     return st.str();
  }

  RTVpEq(Sample *sample) : s(sample) { s->rtctrl.rt_vp_eq = "NONE"; }
};
#endif

