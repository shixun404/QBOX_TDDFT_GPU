////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2009 The Regents of the University of California
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
// PeriodicTable.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef PERIODICTABLE_H
#define PERIODICTABLE_H

#include <map>
#include <vector>
#include <string>
using namespace std;

struct Element
{
  int z;
  string symbol;
  string config;
  double mass;
  Element(int zz, string s, string c, double m) : z(zz), symbol(s), config(c),
    mass(m) {}
};

class PeriodicTable
{
  private:

  vector<Element> ptable;
  map<string,int> zmap;

  public:

  PeriodicTable(void);
  int z(string symbol) const;
  string symbol(int zval) const;
  string configuration(int zval) const;
  string configuration(string symbol) const;
  double mass(int zval) const;
  double mass(string symbol) const;
  int size(void) const;

};
#endif
