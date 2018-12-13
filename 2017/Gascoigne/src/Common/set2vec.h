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


#ifndef  __set2vec_h
#define  __set2vec_h

#include  <vector>
#include  <set>

/*---------------------------------------------------*/

namespace Gascoigne
{
void Set2Vec(std::vector<int>& v, const std::set<int>& h);
void Vec2Set(std::set<int>& h, const std::vector<int>& v);
}

#endif
