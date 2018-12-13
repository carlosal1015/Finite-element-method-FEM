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


#ifndef __find_in_linehang_h
#define __find_in_linehang_h

/*---------------------------------------------------*/

namespace Gascoigne
{
template <int N>
std::pair<typename HangList<N>::iterator,bool> 
find_in_linehang(HangList<N>& LineHang, const fixarray<N,int>& lineglob)
{
  // sort(lineglob.begin(),lineglob.end());
  typename HangList<N>::iterator p = LineHang.find(lineglob);
  bool b=0;
  if(p!=LineHang.end()) b = 1;
  return std::make_pair(p,b);
}
}

#endif
