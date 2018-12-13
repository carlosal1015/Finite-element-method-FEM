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


#ifndef  __hang_h
#define  __hang_h

#include  <vector>
#include  "fixarray.h"
#include  <string>

/*------------------------------------------------------*/

namespace Gascoigne
{
class Hang : public fixarray<3,int>
{
 public:

  Hang(); 
  Hang(const Hang& h);
  Hang(int nh, int nr, int nc) ;
  
  int  hanging   () const { return (*this)[0]; }
  int& hanging   ()       { return (*this)[0]; }
  int  rneighbour() const { return (*this)[1]; }
  int& rneighbour()       { return (*this)[1]; }
  int  cneighbour() const { return (*this)[2]; }
  int& cneighbour()       { return (*this)[2]; }

  const fixarray<3,int>& operator()() const { return (*this);}

  friend std::ostream& operator<<(std::ostream &s, const Hang& A);
  friend std::istream& operator>>(std::istream &s, Hang& A);
};
}

#endif
