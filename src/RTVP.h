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
// RTVP.h
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RTVP_H
#define RTVP_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>

#include "Sample.h"

class RTVP : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "rt_vp"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " rt_vp takes only one value" << endl;
      return 1;
    }

    string v = argv[1];
    if ( !( v == "ON" || v == "OFF" ) )
    {
      if ( ui->onpe0() )
        cout << " rt_vp must be ON or OFF" << endl;
      return 1;
    }

    s->rtctrl.rt_vp = v;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->rtctrl.rt_vp;
     return st.str();
  }

  RTVP(Sample *sample) : s(sample) { s->rtctrl.rt_vp = "OFF"; }
};
#endif
