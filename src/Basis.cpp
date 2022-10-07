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
// Basis.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "Basis.h"
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <vector>
#include <set>
#include <algorithm>
#include <functional>
#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
double Basis::localmemsize(void) const
{
  return
  5.0 * (npes_*nrods_*sizeof(int)) // x[ipe][irod]
  + localsize_[mype_] * (3.0*sizeof(int) + 10 * sizeof(double));
}
double Basis::memsize(void) const { return npes_*localmemsize(); }

MPI_Comm Basis::comm(void) const { return comm_; }

const UnitCell& Basis::cell() const { return cell_; }
const UnitCell& Basis::refcell() const { return refcell_; }
int Basis::idxmin(int i) const { return idxmin_[i]; }
int Basis::idxmax(int i) const { return idxmax_[i]; }
double Basis::ecut() const { return ecut_; }

int Basis::size() const { return size_; }
int Basis::localsize() const { return localsize_[mype_]; }
int Basis::localsize(int ipe) const { return localsize_[ipe]; }
int Basis::maxlocalsize() const { return maxlocalsize_; }
int Basis::minlocalsize() const { return minlocalsize_; }

int Basis::nrods() const { return nrods_; }
int Basis::nrod_loc() const { return nrod_loc_[mype_]; }
int Basis::nrod_loc(int ipe) const { return nrod_loc_[ipe]; }

int Basis::rod_h(int irod) const
{ return rod_h_[mype_][irod]; }
int Basis::rod_h(int ipe, int irod) const
{ return rod_h_[ipe][irod]; }

int Basis::rod_k(int irod) const { return rod_k_[mype_][irod]; }
int Basis::rod_k(int ipe, int irod) const { return rod_k_[ipe][irod]; }

int Basis::rod_lmin(int irod) const { return rod_lmin_[mype_][irod]; }
int Basis::rod_lmin(int ipe, int irod) const { return rod_lmin_[ipe][irod]; }

// size of rod irod on current process
int Basis::rod_size(int irod) const { return rod_size_[mype_][irod]; }
int Basis::rod_size(int ipe, int irod) const { return rod_size_[ipe][irod]; }

// position of first elem. of rod irod in the local list of g vectors
int Basis::rod_first(int irod) const { return rod_first_[mype_][irod]; }
int Basis::rod_first(int ipe, int irod) const { return rod_first_[ipe][irod]; }

int    Basis::idx(int i) const   { return idx_[i]; }
double Basis::g(int i) const     { return g_[i]; }
double Basis::kpg(int i) const   { return kpg_[i]; }
double Basis::gi(int i) const    { return gi_[i]; }
double Basis::kpgi(int i) const  { return kpgi_[i]; }
double Basis::g2(int i) const    { return g2_[i]; }
double Basis::kpg2(int i) const  { return kpg2_[i]; }
double Basis::g2i(int i) const   { return g2i_[i]; }
double Basis::kpg2i(int i) const { return kpg2i_[i]; }
double Basis::gx(int i) const    { return gx_[i]; }
double Basis::kpgx(int i) const  { return kpgx_[i]; }
int    Basis::isort(int i) const { return isort_loc[i]; }

const int*    Basis::idx_ptr(void) const   { return &(idx_[0]); }
const double* Basis::g_ptr(void)  const    { return &(g_[0]); }
const double* Basis::kpg_ptr(void)  const  { return &(kpg_[0]); }
const double* Basis::gi_ptr(void) const    { return &(gi_[0]); }
const double* Basis::kpgi_ptr(void) const  { return &(kpgi_[0]); }
const double* Basis::g2_ptr(void) const    { return &(g2_[0]); }
const double* Basis::kpg2_ptr(void) const  { return &(kpg2_[0]); }
const double* Basis::g2i_ptr(void) const   { return &(g2i_[0]); }
const double* Basis::kpg2i_ptr(void) const { return &(kpg2i_[0]); }
const double* Basis::gx_ptr(int j) const
{ return &(gx_[j*localsize_[mype_]]); }
const double* Basis::kpgx_ptr(int j) const
{ return &(kpgx_[j*localsize_[mype_]]); }

////////////////////////////////////////////////////////////////////////////////
bool Basis::factorizable(int n) const
{
  // next lines: use AIX criterion for all platforms (AIX and fftw)

//#if AIX

  // Acceptable lengths for FFTs in the ESSL library:
  // n = (2^h) (3^i) (5^j) (7^k) (11^m) for n <= 37748736
  // where:
  //  h = 1, 2, ..., 25
  //  i = 0, 1, 2
  //  j, k, m = 0, 1
  if ( n % 11 == 0 ) n /= 11;
  if ( n % 7 == 0 ) n /= 7;
  if ( n % 5 == 0 ) n /= 5;
  if ( n % 3 == 0 ) n /= 3;
  if ( n % 3 == 0 ) n /= 3;
  // the bound h <= 25 is not tested since 2^25 would cause
  // memory allocation problems
  while ( ( n % 2 == 0 ) ) n /= 2;
  return ( n == 1 );

// #else
//   while ( n % 5 == 0 ) n /= 5;
//   while ( n % 3 == 0 ) n /= 3;
//   while ( n % 2 == 0 ) n /= 2;
//   return ( n == 1 );
// #endif
}

////////////////////////////////////////////////////////////////////////////////
int Basis::np(int i) const { return np_[i]; }

////////////////////////////////////////////////////////////////////////////////
bool Basis::fits_in_grid(int np0, int np1, int np2) const
{ return
  ( idxmax_[0] < np0/2 ) && ( idxmin_[0] >= -np0/2 ) &&
  ( idxmax_[1] < np1/2 ) && ( idxmin_[1] >= -np1/2 ) &&
  ( idxmax_[2] < np2/2 ) && ( idxmin_[2] >= -np2/2 );
}

////////////////////////////////////////////////////////////////////////////////
const D3vector Basis::kpoint(void) const { return kpoint_; }

////////////////////////////////////////////////////////////////////////////////
bool Basis::real(void) const { return real_; }
void Basis::set_real(bool isreal) { real_ = isreal; }
inline double sqr( double x ) { return x*x; }
inline void swap(int &x, int &y) { int tmp = x; x = y; y = tmp; }

////////////////////////////////////////////////////////////////////////////////
class Rod
{
  // z-column of non-zero reciprocal lattice vectors
  int h_, k_, lmin_, size_;

  public:
  Rod(int h, int k, int lmin, int size) : h_(h), k_(k),
  lmin_(lmin), size_(size) {}

  int h(void) const { return h_; }
  int k(void) const { return k_; }
  int lmin(void) const { return lmin_; }
  int size(void) const { return size_; }
  bool operator< (const Rod& x) const
  {
    return size_ < x.size();
  }
};

////////////////////////////////////////////////////////////////////////////////
class Node
{
  int id_, nrods_, size_;

  public:
  Node() : id_(0), nrods_(0), size_(0) {}
  Node(int id) : id_(id), nrods_(0), size_(0) {}

  int id(void) const { return id_; }
  int nrods(void) const { return nrods_; }
  int size(void) const { return size_; }

  void addrod(const Rod& r)
  {
    nrods_++;
    size_ += r.size();
  }
  bool operator< (const Node& x) const
  {
    return size_ < x.size();
  }
};

////////////////////////////////////////////////////////////////////////////////
template <class T>
struct ptr_less
{
  public:
  bool operator() ( T* x, T* y ) { return *x < *y; }
};
template <class T>
struct ptr_greater
{
  public:
  bool operator() ( T* x, T* y ) { return *y < *x; }
};

////////////////////////////////////////////////////////////////////////////////
template <class T>
struct VectorLess
{
  // function object for indirect comparison of vector elements
  public:
  vector<T>& a_;
  VectorLess<T>(vector<T>& a) : a_(a) {};
  bool operator() (int i, int j) const
  {
    return a_[i] < a_[j];
  }
};

////////////////////////////////////////////////////////////////////////////////
Basis::Basis(MPI_Comm comm, D3vector kpoint) : comm_(comm)
{
  // Construct the default empty basis
  // cell and refcell are (0,0,0)
  MPI_Comm_rank(comm_,&mype_);
  MPI_Comm_size(comm_,&npes_);

  ecut_ = 0.0;
  kpoint_ = kpoint;
  //real_ = ( kpoint == D3vector(0.0,0.0,0.0) );
  real_ = false;

  localsize_.resize(npes_);
  nrod_loc_.resize(npes_);
  rod_h_.resize(npes_);
  rod_k_.resize(npes_);
  rod_lmin_.resize(npes_);
  rod_size_.resize(npes_);
  rod_first_.resize(npes_);

  // resize with zero cutoff to initialize empty Basis
  resize(cell_,refcell_,0.0);
}

////////////////////////////////////////////////////////////////////////////////
Basis::~Basis(void) {}

////////////////////////////////////////////////////////////////////////////////
bool Basis::resize(const UnitCell& cell, const UnitCell& refcell,
  double ecut)
{
  assert(ecut>=0.0);
  assert(cell.volume() >= 0.0);
  assert(refcell.volume() >= 0.0);

  if ( ecut == ecut_ && refcell == refcell_ && refcell_.volume() != 0.0 )
  {
    cell_ = cell;
    // only the cell changes, ecut and the refcell remain unchanged
    update_g();
    return true;
  }

  ecut_ = ecut;
  cell_ = cell;
  refcell_ = refcell;

  if ( ecut == 0.0 || cell.volume() == 0.0)
  {
    idxmax_[0] = 0;
    idxmax_[1] = 0;
    idxmax_[2] = 0;
    idxmin_[0] = 0;
    idxmin_[1] = 0;
    idxmin_[2] = 0;

    size_ = 0;
    nrods_ = 0;
    for ( int ipe = 0; ipe < npes_; ipe++ )
    {
      localsize_[ipe] = 0;
      nrod_loc_[ipe] = 0;
    }
    maxlocalsize_ = minlocalsize_ = 0;
    np_[0] = np_[1] = np_[2] = 0;
    idx_.resize(3*localsize_[mype_]);
    g_.resize(localsize_[mype_]);
    kpg_.resize(localsize_[mype_]);
    gi_.resize(localsize_[mype_]);
    kpgi_.resize(localsize_[mype_]);
    g2_.resize(localsize_[mype_]);
    kpg2_.resize(localsize_[mype_]);
    g2i_.resize(localsize_[mype_]);
    kpg2i_.resize(localsize_[mype_]);
    gx_.resize(3*localsize_[mype_]);
    kpgx_.resize(3*localsize_[mype_]);
    isort_loc.resize(localsize_[mype_]);
    return true;
  }

  const double two_ecut = 2.0 * ecut;
  const double twopi = 2.0 * M_PI;

  const double kpx = kpoint_.x;
  const double kpy = kpoint_.y;
  const double kpz = kpoint_.z;

  UnitCell defcell;
  // defcell: cell used to define which vectors are contained in the Basis
  // if refcell is defined, defcell = refcell
  // otherwise, defcell = cell
  if ( norm2(refcell.b(0)) + norm2(refcell.b(1)) + norm2(refcell.b(2)) == 0.0 )
  {
    defcell = cell;
  }
  else
  {
    defcell = refcell;
  }

  const D3vector b0 = defcell.b(0);
  const D3vector b1 = defcell.b(1);
  const D3vector b2 = defcell.b(2);

  const double normb2 = norm2(b2);
  const double b2inv2 = 1.0 / normb2;

  const D3vector kp = kpx*b0 + kpy*b1 + kpz*b2;

  if ( !cell.in_bz(kp) )
  {
    if ( mype_ == 0 )
      cout << " Basis::resize: warning: " << kpoint_
           << " out of the BZ: " << kp << endl;
  }

  const double fac = sqrt(two_ecut) / twopi;

  // define safe enclosing domain for any k-point value in the BZ
  const int hmax = (int) ( 1.5 + fac * ( length(defcell.a(0) ) ) );
  const int hmin = - hmax;

  const int kmax = (int) ( 1.5 + fac * ( length(defcell.a(1) ) ) );
  const int kmin = - kmax;

  const int lmax = (int) ( 1.5 + fac * ( length(defcell.a(2) ) ) );
  const int lmin = - lmax;

  multiset<Rod> rodset;

  // build rod set

  int hmax_used = hmin;
  int hmin_used = hmax;
  int kmax_used = kmin;
  int kmin_used = kmax;
  int lmax_used = lmin;
  int lmin_used = lmax;

  if ( real_ )
  {
    // build basis at kpoint (0,0,0)

    // rod(0,0,0)
    // length of rod(0,0,0) is lend+1
    int lend = (int) ( sqrt(two_ecut * b2inv2) );
    size_ = lend + 1;
    rodset.insert(Rod(0,0,0,lend+1));
    nrods_ = 1;

    hmax_used = 0;
    hmin_used = 0;
    kmin_used = 0;
    kmax_used = 0;
    lmin_used = 0;
    lmax_used = lend;

    // rods (0,k,l)
    for ( int k = 1; k <= kmax; k++ )
    {
      int lstart=lmax,lend=lmin;
      bool found = false;
      for ( int l = lmin; l <= lmax; l++ )
      {
        const double two_e = norm2(k*b1+l*b2);
        if ( two_e < two_ecut )
        {
          lstart = min(l,lstart);
          lend = max(l,lend);
          found = true;
        }
      }
      if ( found )
      {
        // non-zero intersection was found
        const int rodsize = lend - lstart + 1;
        size_ += rodsize;
        rodset.insert(Rod(0,k,lstart,rodsize));
        nrods_++;
        kmax_used = max(k,kmax_used);
        kmin_used = min(k,kmin_used);
        lmax_used = max(lend,lmax_used);
        lmin_used = min(lstart,lmin_used);
      }
    }
    // rods (h,k,l) h>0
    for ( int h = 1; h <= hmax; h++ )
    {
      for ( int k = kmin; k <= kmax; k++ )
      {
        int lstart=lmax,lend=lmin;
        bool found = false;
        for ( int l = lmin; l <= lmax; l++ )
        {
          const double two_e = norm2(h*b0+k*b1+l*b2);
          if ( two_e < two_ecut )
          {
            lstart = min(l,lstart);
            lend = max(l,lend);
            found = true;
          }
        }
        if ( found )
        {
          // non-zero intersection was found
          const int rodsize = lend - lstart + 1;
          size_ += rodsize;
          rodset.insert(Rod(h,k,lstart,rodsize));
          nrods_++;
          hmax_used = max(h,hmax_used);
          hmin_used = min(h,hmin_used);
          kmax_used = max(k,kmax_used);
          kmin_used = min(k,kmin_used);
          lmax_used = max(lend,lmax_used);
          lmin_used = min(lstart,lmin_used);
        }
      }
    }
  }
  else
  {
    // build Basis for k != (0,0,0)

    size_ = 0;
    nrods_ = 0;
    // rods (h,k,l)
    for ( int h = hmin; h <= hmax; h++ )
    {
      for ( int k = kmin; k <= kmax; k++ )
      {
        int lstart=lmax,lend=lmin;
        bool found = false;
        for ( int l = lmin; l <= lmax; l++ )
        {
          const double two_e = norm2((kpx+h)*b0 + (kpy+k)*b1 + (kpz+l)*b2);
          if ( two_e < two_ecut )
          {
            lstart = min(l,lstart);
            lend = max(l,lend);
            found = true;
          }
        }
        if ( found )
        {
          // non-zero intersection was found
          const int rodsize = lend - lstart + 1;
          size_ += rodsize;
          rodset.insert(Rod(h,k,lstart,rodsize));
          nrods_++;
          hmax_used = max(h,hmax_used);
          hmin_used = min(h,hmin_used);
          kmax_used = max(k,kmax_used);
          kmin_used = min(k,kmin_used);
          lmax_used = max(lend,lmax_used);
          lmin_used = min(lstart,lmin_used);
        }
      }
    }
  }

#if DEBUG
  cout << " hmin/hmax: " << hmin << " / " << hmax << endl;
  cout << " kmin/kmax: " << kmin << " / " << kmax << endl;
  cout << " lmin/lmax: " << lmin << " / " << lmax << endl;
  cout << " hmin/hmax used: " << hmin_used << " / " << hmax_used << endl;
  cout << " kmin/kmax used: " << kmin_used << " / " << kmax_used << endl;
  cout << " lmin/lmax used: " << lmin_used << " / " << lmax_used << endl;
#endif

  idxmax_[0] = hmax_used;
  idxmin_[0] = hmin_used;

  idxmax_[1] = kmax_used;
  idxmin_[1] = kmin_used;

  idxmax_[2] = lmax_used;
  idxmin_[2] = lmin_used;

  assert(hmax_used - hmin_used + 1 <= 2 * hmax);
  assert(kmax_used - kmin_used + 1 <= 2 * kmax);
  assert(lmax_used - lmin_used + 1 <= 2 * lmax);

  // compute good FFT sizes
  // use values independent of the kpoint
  int n;
  n = 2 * hmax;
  while ( !factorizable(n) ) n+=2;
  np_[0] = n;
  n = 2 * kmax;
  while ( !factorizable(n) ) n+=2;
  np_[1] = n;
  n = 2 * lmax;
  while ( !factorizable(n) ) n+=2;
  np_[2] = n;

  // Distribute the basis on npes_ processors

  // build a min-heap of Nodes

  vector<Node*> nodes(npes_);

  for ( int ipe = 0; ipe < npes_; ipe++ )
  {
    nodes[ipe] = new Node(ipe);
    localsize_[ipe] = 0;
    nrod_loc_[ipe] = 0;
    rod_h_[ipe].resize(0);
    rod_k_[ipe].resize(0);
    rod_lmin_[ipe].resize(0);
    rod_size_[ipe].resize(0);
  }

  // nodes contains a valid min-heap of zero-size Nodes

  // insert rods into the min-heap
  // keep track of where rod(0,0,0) goes
  int pe_rod0 = -1, rank_rod0 = -1;
  multiset<Rod>::iterator p = rodset.begin();
  while ( p != rodset.end() )
  {
    // pop smallest element
    pop_heap(nodes.begin(), nodes.end(), ptr_greater<Node>());

    // add rod size to smaller element
    nodes[npes_-1]->addrod(*p);
    int ipe = nodes[npes_-1]->id();

    // update info about rod on process ipe
    nrod_loc_[ipe]++;
    rod_h_[ipe].push_back(p->h());
    rod_k_[ipe].push_back(p->k());
    rod_lmin_[ipe].push_back(p->lmin());
    rod_size_[ipe].push_back(p->size());
    localsize_[ipe] += p->size();
    if ( p->h() == 0 && p->k() == 0 )
    {
      pe_rod0 = ipe;
      rank_rod0 = nodes[npes_-1]->nrods()-1;
    }

    // push modified element back in the heap
    push_heap(nodes.begin(), nodes.end(), ptr_greater<Node>());

    p++;
  }

  maxlocalsize_ = (*max_element(nodes.begin(), nodes.end(),
    ptr_less<Node>()))->size();
  minlocalsize_ = (*min_element(nodes.begin(), nodes.end(),
    ptr_less<Node>()))->size();

  for ( int ipe = 0; ipe < npes_; ipe++ )
  {
    delete nodes[ipe];
  }

  // swap node pe_rod0 with node 0 in order to have rod(0,0,0) on node 0
  swap(nrod_loc_[0], nrod_loc_[pe_rod0]);
  rod_h_[pe_rod0].swap(rod_h_[0]);
  rod_k_[pe_rod0].swap(rod_k_[0]);
  rod_lmin_[pe_rod0].swap(rod_lmin_[0]);
  rod_size_[pe_rod0].swap(rod_size_[0]);
  swap(localsize_[0], localsize_[pe_rod0]);
  //Node *tmpnodeptr = nodes[0]; nodes[0] = nodes[pe_rod0];
  //  nodes[pe_rod0]=tmpnodeptr;

  // reorder rods on node 0 so that rod(0,0,0) comes first
  swap(rod_h_[0][rank_rod0], rod_h_[0][0]);
  swap(rod_k_[0][rank_rod0], rod_k_[0][0]);
  swap(rod_lmin_[0][rank_rod0], rod_lmin_[0][0]);
  swap(rod_size_[0][rank_rod0], rod_size_[0][0]);

  // compute position of first element of rod (ipe,irod)
  for ( int ipe = 0; ipe < npes_; ipe++ )
  {
    rod_first_[ipe].resize(nrod_loc_[ipe]);
    if ( nrod_loc_[ipe] > 0 )
      rod_first_[ipe][0] = 0;
    for ( int irod = 1; irod < nrod_loc_[ipe]; irod++ )
    {
      rod_first_[ipe][irod] = rod_first_[ipe][irod-1] + rod_size_[ipe][irod-1];
    }
  }

  // local arrays idx, g, gi, g2i, g2, gx
  idx_.resize(3*localsize_[mype_]);
  int i = 0;
  for ( int irod = 0; irod < nrod_loc_[mype_]; irod++ )
  {
    for ( int l = 0; l < rod_size_[mype_][irod]; l++ )
    {
      idx_[3*i]   = rod_h_[mype_][irod];
      idx_[3*i+1] = rod_k_[mype_][irod];
      idx_[3*i+2] = rod_lmin_[mype_][irod] + l;

      i++;
    }
  }

  g_.resize(localsize_[mype_]);
  kpg_.resize(localsize_[mype_]);
  gi_.resize(localsize_[mype_]);
  kpgi_.resize(localsize_[mype_]);
  g2_.resize(localsize_[mype_]);
  kpg2_.resize(localsize_[mype_]);
  g2i_.resize(localsize_[mype_]);
  kpg2i_.resize(localsize_[mype_]);
  gx_.resize(3*localsize_[mype_]);
  kpgx_.resize(3*localsize_[mype_]);
  isort_loc.resize(localsize_[mype_]);

  update_g();

  // basis set construction is complete

  return true;
}

////////////////////////////////////////////////////////////////////////////////
void Basis::update_g(void)
{
  // compute the values of g, kpg, gi, g2i, g2, kpg2, gx
  // N.B. use the values of cell (not defcell)

  const int locsize = localsize_[mype_];
  for ( int i = 0; i < locsize; i++ )
  {
    D3vector gt = idx_[3*i+0] * cell_.b(0) +
                  idx_[3*i+1] * cell_.b(1) +
                  idx_[3*i+2] * cell_.b(2);

    D3vector kpgt = (kpoint_.x + idx_[3*i+0]) * cell_.b(0) +
                    (kpoint_.y + idx_[3*i+1]) * cell_.b(1) +
                    (kpoint_.z + idx_[3*i+2]) * cell_.b(2);

    gx_[i] = gt.x;
    gx_[locsize+i] = gt.y;
    gx_[locsize+locsize+i] = gt.z;
    kpgx_[i] = kpgt.x;
    kpgx_[locsize+i] = kpgt.y;
    kpgx_[locsize+locsize+i] = kpgt.z;

    g2_[i] = norm2(gt);
    g_[i] = sqrt( g2_[i] );

    kpg2_[i] = norm2(kpgt);
    kpg_[i] = sqrt( kpg2_[i] );

    gi_[i] = g_[i] > 0.0 ? 1.0 / g_[i] : 0.0;
    kpgi_[i] = kpg_[i] > 0.0 ? 1.0 / kpg_[i] : 0.0;
    g2i_[i] = gi_[i] * gi_[i];
    kpg2i_[i] = kpgi_[i] * kpgi_[i];
    isort_loc[i] = i;
  }

  VectorLess<double> g2_less(g2_);
  sort(isort_loc.begin(), isort_loc.end(), g2_less);
#if DEBUG
  for ( int i = 0; i < locsize; i++ )
  {
    cout << mype_ << " sorted " << i << " " << g2_[isort_loc[i]] << endl;
  }
#endif
}

////////////////////////////////////////////////////////////////////////////////
void Basis::print(ostream& os)
{
  os << mype_ << ": ";
  os << " Basis.kpoint():   " << kpoint() << endl;
  os << mype_ << ": ";
  os << " Basis.kpoint():   " << kpoint().x << " * b0 + "
                              << kpoint().y << " * b1 + "
                              << kpoint().z << " * b2" << endl;
  os << mype_ << ": ";
  os << " Basis.kpoint():   " << kpoint().x * cell().b(0) +
                                 kpoint().y * cell().b(1) +
                                 kpoint().z * cell().b(2) << endl;
  os << mype_ << ": ";
  os << " Basis.cell():     " << endl << cell() << endl;
  os << mype_ << ": ";
  os << " Basis.ref cell(): " << endl << refcell() << endl;
  os << mype_ << ": ";
  os << " Basis.ecut():     " << ecut() << endl;
  os << mype_ << ": ";
  os << " Basis.np(0,1,2):  " << np(0) << " "
       << np(1) << " " << np(2) << endl;
  os << mype_ << ": ";
  os << " Basis.idxmin:       " << idxmin(0) << " "
       << idxmin(1) << " " << idxmin(2) << endl;
  os << mype_ << ": ";
  os << " Basis.idxmax:       " << idxmax(0) << " "
       << idxmax(1) << " " << idxmax(2) << endl;
  os << mype_ << ": ";
  os << " Basis.size():     " << size() << endl;
  os << mype_ << ": ";
  os << " Basis.localsize(): " << localsize() << endl;
  os << mype_ << ": ";
  os << " Basis.nrods():    " << nrods() << endl;
  os << mype_ << ": ";
  os << " Basis.real():     " << real() << endl;
  os << mype_ << ": ";
  os << " Basis total mem size (MB): " << memsize() / 1048576 << endl;
  os << mype_ << ": ";
  os << " Basis local mem size (MB): " << localmemsize() / 1048576 << endl;

  os << mype_ << ": ";
  os << "   ig      i   j   k        gx      gy      gz       |k+g|^2"
     << endl;
  os << mype_ << ": ";
  os << "   --      -   -   -        --      --      --       -------"
     << endl;
  for ( int i = 0; i < localsize(); i++ )
  {
    os << mype_ << ": ";
    os << setw(5) << i << "   "
       << setw(4) << idx(3*i)
       << setw(4) << idx(3*i+1)
       << setw(4) << idx(3*i+2)
       << "    "
       << setw(8) << setprecision(4) << gx(i)
       << setw(8) << setprecision(4) << gx(i+localsize())
       << setw(8) << setprecision(4) << gx(i+2*localsize())
       << setw(12) << setprecision(4) << 0.5 * kpg2(i)
       << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream& os, Basis& b)
{
  b.print(os);
  return os;
}
