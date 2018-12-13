/**
*
* Copyright (C) 2004 by the Gascoigne 3D authors
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


#ifndef __triple_h
#define __triple_h

namespace Gascoigne
{
template <class T1, class T2, class T3>
struct triple {
  typedef T1 first_type;
  typedef T2 second_type;
  typedef T3 third_type;

  T1 first;
  T2 second;
  T3 third;
  triple() : first(T1()), second(T2()), third(T3()) {}
  triple(const T1& a, const T2& b, const T3& c) : first(a), second(b), third(c) {}

#ifdef __STL_MEMBER_TEMPLATES
  template <class U1, class U2, class U3>
  triple(const triple<U1, U2, U3>& p) : first(p.first), second(p.second), third(p.third) {}
#endif
};

template <class T1, class T2, class T3>
inline bool operator==(const triple<T1, T2, T3>& x, const triple<T1, T2, T3>& y) { 
  return x.first == y.first && x.second == y.second && x.third == y.third; 
}

template <class T1, class T2, class T3>
inline bool operator<(const triple<T1, T2, T3>& x, const triple<T1, T2, T3>& y) { 
  return x.first < y.first || (!(y.first < x.first) && x.second < y.second)
    || ( !(y.first < x.first) && !(y.second < x.second) && x.third < y.third)  ; 
}

template <class T1, class T2, class T3>
inline triple<T1, T2, T3> make_triple(const T1& x, const T2& y, const T3& z) {
  return triple<T1, T2, T3>(x, y,z);
}
}
#endif
