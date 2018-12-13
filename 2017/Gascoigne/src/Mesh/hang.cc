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


#include "hang.h"

using namespace std;

/*********************************************************************/

namespace Gascoigne
{
Hang::Hang() : fixarray<3,int>(-1)
{}

/*********************************************************************/

Hang::Hang(const Hang& h) 
  : fixarray<3,int>(h) {}

/*********************************************************************/

Hang::Hang(int nh, int nr, int nc)
{
  hanging() = nh;
  rneighbour() = nr;
  cneighbour() = nc;
}
  
/*********************************************************************/

ostream& operator<<(ostream &s, const Hang& A)
{
  s << A.hanging()    << " : ";
  s << A.rneighbour() << " ";
  s << A.cneighbour() << " ";
  return s;
}

/*********************************************************************/

istream& operator>>(istream &s, Hang& A)
{
  char symbol;
  s >> A.hanging() >> symbol;
  s >> A.rneighbour() >> A.cneighbour();
  return s;
}
}
