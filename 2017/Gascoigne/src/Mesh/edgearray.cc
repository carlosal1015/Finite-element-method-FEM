/**
*
* Copyright (C) 2004, 2005, 2011 by the Gascoigne 3D authors
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


#include "edgearray.h"

/*--------------------------------------------------------------*/

namespace Gascoigne
{

/*--------------------------------------------------------------*/

template<int N>
EdgeArray<N>::EdgeArray(const fixarray<N,int>& e) : fixarray<N,int>(e) { }

/*--------------------------------------------------------------*/

template<>
bool EdgeArray<2>::operator==(const fixarray<2,int> &A) const
{
  // the correctness of this method relies on all vertices (in *this) being different

  if ( ((*this)[0]!=A[0]) && ((*this)[0]!=A[1]) ) return 0;
  if ( ((*this)[1]!=A[0]) && ((*this)[1]!=A[1]) ) return 0;
  return 1;
}

/*--------------------------------------------------------------*/

template<>
int EdgeArray<2>::sum() const 
{ 
  return (*this)[0]+(*this)[1];
}

/*--------------------------------------------------------------*/

template<>
bool EdgeArray<4>::operator==(const fixarray<4,int> &A) const
{
  // the correctness of this method relies on all vertices (in *this) being different
  // and that we can identify a face by looking at only three vertices

  if ( ((*this)[0]!=A[0]) && ((*this)[0]!=A[1]) && 
       ((*this)[0]!=A[2]) && ((*this)[0]!=A[3]) ) return 0;
  if ( ((*this)[1]!=A[0]) && ((*this)[1]!=A[1]) && 
       ((*this)[1]!=A[2]) && ((*this)[1]!=A[3]) ) return 0;
  if ( ((*this)[2]!=A[0]) && ((*this)[2]!=A[1]) && 
       ((*this)[2]!=A[2]) && ((*this)[2]!=A[3]) ) return 0;
  return 1;
}

/*--------------------------------------------------------------*/

template<>
int EdgeArray<4>::sum() const 
{ 
  return (*this)[0]+(*this)[1]+(*this)[2]+(*this)[3];
}

/*--------------------------------------------------------------*/

template class EdgeArray<2>;
template class EdgeArray<4>;
}

/*--------------------------------------------------------------*/

