////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2014 The Regents of the University of California
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
// RTRhoTol.h
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////
#ifndef RTRHOTOL_H
#define RTRHOTOL_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class RTRhoTol : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "rt_rho_tol"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " rt_rho_tol takes only one value" << endl;
      return 1;
    }

    double v = atof(argv[1]);
    if ( v < 0.0 )
    {
      if ( ui->onpe0() )
        cout << " rt_rho_tol must be non-negative" << endl;
      return 1;
    }

    s->rtctrl.rt_rho_tol = v;
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->rtctrl.rt_rho_tol;
     return st.str();
  }

  RTRhoTol(Sample *sample) : s(sample) { s->rtctrl.rt_rho_tol = 1.0e-15; }
};
#endif


