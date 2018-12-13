/**
*
* Copyright (C) 2004, 2005, 2008, 2011 by the Gascoigne 3D authors
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


#ifndef  __hanglist_h
#define  __hanglist_h

#include  <string>
#include  "edgearray.h"
#include  "hang.h"

#ifdef __OLDCOMPILER__
#include  <hash_map>
#define HANGMAP  hash_map<EdgeArray<N>,Hang,EdgeHash>
#else
#ifdef __NEWER_THAN_GCC_4_2__
#include <tr1/unordered_map>
#define HANGMAP   std::tr1::unordered_map<EdgeArray<N>,Hang,EdgeHash> 
#else
#include  <ext/hash_map>
#define HANGMAP  __gnu_cxx::hash_map<EdgeArray<N>,Hang,EdgeHash> 
#endif
#endif

/*------------------------------------------------------*/

namespace Gascoigne
{

//
/// This hash function has to be consistent with the operator "=="
/// for EdgeArrays, i.e. permutated fixarrays
//

class EdgeHash
{
 public:
  template<int N>
    int operator()(const EdgeArray<N>& h) const { return h.sum();}
};


/*------------------------------------------------------*/

template<int N>
class HangList : public HANGMAP
{
 protected:

 public:

  typedef typename HANGMAP::iterator        iterator;
  typedef typename HANGMAP::const_iterator  const_iterator;

  void update(const std::vector<int>&);
  void update(const std::vector<int>&, const std::vector<int>&);
  void make_consistent(HangList<N>&);
  void move(HangList<N>& src, iterator& p);
  HangList<N>& operator=(const HangList<N>& A);
  void BinWrite(std::ostream &s) const;
  void BinRead(std::istream &s);
};

/*------------------------------------------------------*/

template<int N>
std::ostream& operator<<(std::ostream &s, const HangList<N>& A);

template<int N>
std::istream& operator>>(std::istream &s, HangList<N>& A);
}

/*------------------------------------------------------*/

#undef HANGMAP

#endif
