/**
*
* Copyright (C) 2004, 2005, 2008 by the Gascoigne 3D authors
*
* This file is part of Gascoigne 3D
*
* Gascoigne 3D is free software: you can redistribute it
* and/or modify it under the terms of the GNU General Public
* License as published by the Free Software Foundation, either
* version 3 of the License, or (at your option) any later
* version.
*
* Gascoigne 3D is distributed in the hope that it will be
* useful, but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* Please refer to the file LICENSE.TXT for further information
* on this license.
*
**/


#ifndef __nvector_h
#define __nvector_h

#include  <vector>
#include  <iterator>
#include  <iostream>
#include  <iterator>
#include  <climits>
#include  <numeric>
#include  <cassert>
#include  <cstdlib>

#include  "gascoignemath.h"

/*----------------------------------------------*/

namespace Gascoigne
{
template<class T>
class nvector : public std::vector<T>
{
private:

public:

  typedef typename std::vector<T>::iterator       iterator;
  typedef typename std::vector<T>::const_iterator const_iterator;

  ~nvector()  {}
  nvector()                    : std::vector<T>()    {}
  nvector(size_t n)               : std::vector<T>(n)   {}
  nvector(size_t n, const T& d)   : std::vector<T>(n,d) {}
  nvector(const std::vector<T>& v) : std::vector<T>(v)   {}

  friend std::ostream& operator<<(std::ostream &s, const nvector<T>& A)
  {
    std::ostream_iterator<T>  os(s, " ");
    copy(A.begin(),A.end(),os);
    return s;
  }
  friend std::istream& operator>>(std::istream &s, nvector<T>& A)
    {
      iterator p = A.begin();
      while(p!=A.end())
	s >> *p++;
      return s;
    }

  void write_data(std::ostream& s) const
  {
    s << std::vector<T>::size() << std::endl << *this;
  }

  void read_data(std::istream& s)
  {
    size_t n;
    s >> n;
    reservesize(n);
    s >> *this;
  }

  const T& secure_access(int i) const {
    assert(i<std::vector<T>::size());
    assert(i>=0);
    return std::vector<T>::operator[](i);
  }
  T& secure_access(int i) {
    assert(i<std::vector<T>::size());
    assert(i>=0);
    return std::vector<T>::operator[](i);
  }


  T operator*(const nvector& v) const;
  nvector<T>&   operator=(const T&);
  nvector<T>&   operator=(const std::vector<T>&);
  nvector<T>&   operator*=(const T& d);
  nvector<T>&   operator+=(const T& d) { add(d); return *this; }
  nvector<T>&   operator+=(const nvector& v) { add(1,v); return *this; }
  nvector<T>&   operator-=(const nvector& v) { add(-1,v); return *this; }

  void   zero ();
  void   equ  (const T&);
  void   equ  (const T&, const nvector&);
  void   equ  (const T&, const nvector&, const T&, const nvector&);
  void   equ  (const T&, const nvector&, const T&, const nvector&, 
	       const T&, const nvector&);
  void   sequ (const T&, const T&, const nvector&);
  void   add  (const T&);
  void   add  (const T&, const nvector&);
  void   add  (const T&, const nvector&, const T&, const nvector&);
  void   add  (const T&, const nvector&, const T&, const nvector&, 
	       const T&, const nvector&);
  void   sadd (const T&, const T&, const nvector&);
  double max  ()     const;
  double min  ()     const;
  T      sum()     const;
  double      norm()     const;
  double      norm_l1()     const;
  double      norm_l8()     const;

  void ReInit(size_t n)
    {
      std::vector<T>::reserve(n); std::vector<T>::resize(n);
    }
  void memory(size_t n)
    {
      ReInit(n);
    }
  void reservesize(size_t n)
    {
      ReInit(n);
    }
  void reservesize(size_t n, const T& s)
    {
      std::vector<T>::reserve(n); std::vector<T>::resize(n,s);
    }
  void reservesize(const nvector<T>& v) {reservesize(v.size());}

  void BinWrite(std::ostream& out) const;
  void BinRead (std::istream& in);
  int  find(const T& x) const;
};


/**************************************************/

template<class T>
inline double nvector<T>::norm() const
{
  const_iterator first  = std::vector<T>::begin();
  const_iterator last   = std::vector<T>::end();

  T n(0);
  while( first != last)
    {
      n += ((*first)) * ((*first));
      first++;
    }
  return sqrt(static_cast<double>(n));
}

/**************************************************/

template<class T>
inline T nvector<T>::sum() const
{
  const_iterator first  = std::vector<T>::begin();
  const_iterator last   = std::vector<T>::end();

  T n(0);
  while( first != last)
    {
      n += (*first++);
    }
  return n;
}

/**************************************************/

template<class T>
inline double nvector<T>::norm_l1() const
{
  const_iterator first  = std::vector<T>::begin();
  const_iterator last   = std::vector<T>::end();

  double n(0);
  while( first != last)
    {
      n += fabs((*first++));
    }
  return n;
}

/**************************************************/

template<class T>
inline double nvector<T>::norm_l8() const
{
  const_iterator first  = std::vector<T>::begin();
  const_iterator last   = std::vector<T>::end();

  double n(0);
  while( first != last)
    {
      n = Gascoigne::max(n,fabs(*first));
      first++;
    }
  return n;
}

/**************************************************/

template<class T>
inline nvector<T>& nvector<T>::operator=(const T& d)
{
  iterator  first  = std::vector<T>::begin();
  const_iterator last   = std::vector<T>::end();

  while( first != last)
    {
      *first++ = d;
    }
  return *this;
}

/**************************************************/

template<class T>
inline nvector<T>& nvector<T>::operator=(const std::vector<T>& v)
{
  assert(std::vector<T>::size()==v.std::template vector<T>::size());
  iterator  first  = std::vector<T>::begin();
  const_iterator last   = std::vector<T>::end();
  const_iterator vfirst = v.std::template vector<T>::begin();

  while( first != last)
    {
      *first++ = *vfirst++;
    }
  return *this;
}

/**************************************************/

template<class T>
inline nvector<T>& nvector<T>::operator*=(const T& d)
{
  iterator  first  = std::vector<T>::begin();
  const_iterator last   = std::vector<T>::end();

  while(first != last)
    {
      (*first++) *= d;
   }
  return *this;
}

/**************************************************/

template<class T>
inline void nvector<T>::zero()
{
  iterator  first  = std::vector<T>::begin();
  const_iterator last   = std::vector<T>::end();

  while( first != last)
    {
      *first++ = 0;
    }
}

/**************************************************/

template<class T>
inline double nvector<T>::max() const
{
  double d = 0;//std::numeric_limits<double>::min();
/*   double d = std::numeric_limits<double>::min(); */
  const_iterator first  = std::vector<T>::begin();
  const_iterator last   = std::vector<T>::end();

  while( first != last)
    {
      d = Gascoigne::max( d, fabs((*first)));
      first++;
    }
  return d;
}

/**************************************************/

template<class T>
inline double nvector<T>::min() const
{
  double d = 100000.;//std::numeric_limits<double>::max();
/*   double d = std::numeric_limits<double>::max(); */
  const_iterator first  = std::vector<T>::begin();
  const_iterator last   = std::vector<T>::end();

  while( first != last)
    {
      d = Gascoigne::min( d, fabs((*first)));
      first++;
    }
  return d;
}

/**************************************************/

template<class T>
inline T nvector<T>::operator* (const nvector<T>& v) const
{
  const_iterator first  = std::vector<T>::begin();
  const_iterator last   = std::vector<T>::end();
  const_iterator first2 = v.begin();

  T d(0);
  while(first != last)
    {
      d += (*first++) * (*first2++);
    }
  return d;

  //return inner_product(first,last,first2,0.);
}

/**************************************************/

template<class T>
inline void nvector<T>::equ (const T& d)
{
  iterator  first  = std::vector<T>::begin();
  const_iterator last   = std::vector<T>::end();

  while(first != last)
    {
      (*first++) = d;
    }
}

/**************************************************/

template<class T>
inline void nvector<T>::equ (const T& d, const nvector<T>& v)
{
  iterator  first  = std::vector<T>::begin();
  const_iterator last   = std::vector<T>::end();
  const_iterator first2 = v.begin();

  while(first != last)
    {
      (*first++) = d*(*first2++);
    }
}

/**************************************************/

template<class T>
inline void nvector<T>::equ (const T& d, const nvector<T>& v, 
			     const T& e, const nvector<T>& w)
{
  iterator  first  = std::vector<T>::begin();
  const_iterator last   = std::vector<T>::end();
  const_iterator first2 = v.begin();
  const_iterator first3 = w.begin();

  while(first != last)
    {
      (*first++) = d*(*first2++) + e*(*first3++);
    }
}

/**************************************************/

template<class T>
inline void nvector<T>::equ (const T& d, const nvector<T>& v, 
			     const T& e, const nvector<T>& w,
			     const T& f, const nvector<T>& x)
{
  iterator  first  = std::vector<T>::begin();
  const_iterator last   = std::vector<T>::end();
  const_iterator first2 = v.begin();
  const_iterator first3 = w.begin();
  const_iterator first4 = x.begin();

  while(first != last)
    {
      (*first++) = d*(*first2++) + e*(*first3++) + f*(*first4++);
    }
}

/**************************************************/

template<class T>
inline void nvector<T>::sequ (const T& s, const T& d, const nvector<T>& v)
{
  iterator  first  = std::vector<T>::begin();
  const_iterator last   = std::vector<T>::end();
  const_iterator first2 = v.begin();

  while(first != last)
    {
      (*first) = s*(*first) + d*(*first2++);
      first++;
    }
}

/**************************************************/

template<class T>
inline void nvector<T>::add (const T& d)
{
  iterator  first  = std::vector<T>::begin();
  const_iterator last   = std::vector<T>::end();

  while(first != last)
    {
      (*first++) += d;
   }
}

/**************************************************/

template<class T>
inline void nvector<T>::add (const T& d, const nvector<T>& v)
{
  iterator  first  = std::vector<T>::begin();
  const_iterator last   = std::vector<T>::end();
  const_iterator first2 = v.begin();

  while(first != last)
    {
      (*first++) += d*(*first2++);
   }
}

/**************************************************/

template<class T>
inline void nvector<T>::add (const T& d, const nvector<T>& v,
			     const T& e, const nvector<T>& w)
{
  iterator  first  = std::vector<T>::begin();
  const_iterator last   = std::vector<T>::end();
  const_iterator first2 = v.begin();
  const_iterator first3 = w.begin();

  while(first != last)
    {
      (*first++) += d*(*first2++) + e*(*first3++);
   }
}

/**************************************************/

template<class T>
inline void nvector<T>::add (const T& d, const nvector<T>& v,
			     const T& e, const nvector<T>& w,
			     const T& f, const nvector<T>& x)
{
  iterator  first  = std::vector<T>::begin();
  const_iterator last   = std::vector<T>::end();
  const_iterator first2 = v.begin();
  const_iterator first3 = w.begin();
  const_iterator first4 = x.begin();

  while(first != last)
    {
      (*first++) += d*(*first2++) + e*(*first3++) + f*(*first4++);
   }
}

/**************************************************/

template<class T>
inline void nvector<T>::sadd (const T& a, const T& d, const nvector<T>& v)
{
  iterator  first  = std::vector<T>::begin();
  const_iterator last   = std::vector<T>::end();
  const_iterator first2 = v.begin();

  while(first != last)
    {
      (*first) = a*(*first) + d*(*first2++);
      first++;
   }
}

/**************************************************/

template<class T>
inline void nvector<T>::BinWrite(std::ostream& out) const
{
  out << std::vector<T>::size() << std::endl << "[";
  
  int sizeT = sizeof(T);
  for(int i=0; i<std::vector<T>::size(); i++)
    {
      out.write (reinterpret_cast<const char*>(&(std::vector<T>::operator[](i))),sizeT);
    }
  out << "]"; 
}

/**********************************************************/

template<class T>
inline void nvector<T>::BinRead(std::istream& in)
{
  char c;
  int  n;
  in >> n >> c;
  std::vector<T>::resize(n);
  
  int sizeT = sizeof(T);
  for(int i=0; i<std::vector<T>::size(); i++)
    {
      in.read(reinterpret_cast<char*>(&(nvector<T>::operator[](i))),sizeT);
    }
  in >> c;
}

/**********************************************************/

template<class T>
inline int nvector<T>::find(const T& x) const
{
  for (int i=0; i<std::vector<T>::size(); i++)
    {
      if ((*this)[i]==x) return i;
    }
  return -1;
}

/**********************************************************/
}

#endif

