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


#ifndef __hangsort_h
#define __hangsort_h

/*---------------------------------------------------*/

namespace Gascoigne
{
class HangEdgeSort{
protected:
  const EdgeManager& HR;
public:
  HangEdgeSort(const EdgeManager& H) : HR(H) {}
  bool operator() (int i, int j) const
    {
      return !HR.EdgeIsHanging(i) && HR.EdgeIsHanging(j);
    }
};

/*---------------------------------------------------*/

class HangEdgeSort2{
protected:
  const EdgeManager& HR;
public:
  HangEdgeSort2(const EdgeManager& H) : HR(H) {}
  bool operator() (const Edge& i, const Edge& j) const
    {
      return !HR.EdgeIsHanging(i) && HR.EdgeIsHanging(j);
    }
};
}

#endif
