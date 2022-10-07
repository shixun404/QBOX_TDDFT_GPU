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
// RTVpSigma.h
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RTVPSIGMA_H
#define RTVPSIGMA_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>

#include "Sample.h"

class RTVpSigma : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "rt_vp_sigma"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " rt_vp_sigma takes only one value" << endl;
      return 1;
    }

    double v = atof(argv[1]);
    if ( v == 0.0 )
    {
      if ( ui->onpe0() )
        cout << " rt_vp_sigma must be non-zero" << endl;
      return 1;
    }

    s->rtctrl.rt_vp_sigma = v;
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->rtctrl.rt_vp_sigma;
     return st.str();
  }

  RTVpSigma(Sample *sample) : s(sample) { s->rtctrl.rt_vp_sigma = 0.0; }
};
#endif


