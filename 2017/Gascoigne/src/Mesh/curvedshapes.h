/**
*
* Copyright (C) 2004, 2005 by the Gascoigne 3D authors
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


#ifndef __curvedshapes_h
#define __curvedshapes_h

#include <map>
#include <set>
#include <vector>
#include  <string>
#include "boundaryfunction.h"
#include "boundaryline.h"
#include "boundaryquad.h"
#include "vertex.h"

/******************************************************/

namespace Gascoigne
{
template<int DIM>
class CurvedShapes : public std::map<int,BoundaryFunction<DIM>* >
{
 public:

  typedef typename std::map<int,BoundaryFunction<DIM>* >::iterator iterator;

  ~CurvedShapes() {
    for (iterator p=std::map<int,BoundaryFunction<DIM>* >::begin();p!=std::map<int,BoundaryFunction<DIM>* >::end();++p)
    if (p->second) {
      std::cerr<< "not deleting shape: "<< p->second->GetName() << std::endl;
    }
  }

  const BoundaryFunction<DIM>& GetShape(int col) const { return *this->find(col)->second;}

  void AddShape(int col, BoundaryFunction<DIM>* f) {
    (*this)[col] = f;
  }
  
  void newton(int col,Vertex<DIM>& V) { GetShape(col).newton(V); }

    int Curved(int col) const { return (this->find(col)!=std::map<int,BoundaryFunction<DIM>* >::end());}

    bool empty() const {return (std::map<int,BoundaryFunction<DIM>* >::size()==0);}
};
}

/******************************************************/

#endif
