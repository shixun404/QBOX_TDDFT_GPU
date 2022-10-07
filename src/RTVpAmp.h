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
// RTVpAmp.h
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RTVPAMP_H
#define RTVPAMP_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>

#include "Sample.h"

class RTVpAmp : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "rt_vp_amp"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 4 )
    {
      if ( ui->onpe0() )
      cout << " rt_vp_amp must be specified with 3 values (x, y, z)" << endl;
      return 1;
    }

    D3vector v(atof(argv[1]), atof(argv[2]), atof(argv[3]));

    s->rtctrl.rt_vp_amp = v;
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->rtctrl.rt_vp_amp;
     return st.str();
  }

  RTVpAmp(Sample *sample) : s(sample) { s->rtctrl.rt_vp_amp = D3vector(0.0,0.0,0.0); }
};
#endif


