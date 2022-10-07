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
// RTPropagator.h
//
////////////////////////////////////////////////////////////////////////////////
//
// Developed by Dr. Min Choi and Prof. Bryan Wong in UCR
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RTPROPAGATOR_H
#define RTPROPAGATOR_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>

#include "Sample.h"

class RTPropagator : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "rt_propagator"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " rt_propagator takes only one value" << endl;
      return 1;
    }

    string v = argv[1];
    if ( !( v == "CN" || v == "PTCN" || v == "ST2" || v == "ST4" ) )
    {
      if ( ui->onpe0() )
        cout << " rt_propagator must be CN or PTCN or ST2 or ST4"
             << endl;
      return 1;
    }

    s->rtctrl.rt_propagator = v;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->rtctrl.rt_propagator;
     return st.str();
  }

  RTPropagator(Sample *sample) : s(sample) { s->rtctrl.rt_propagator = "NONE"; };
};
#endif
