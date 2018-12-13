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


#ifndef __vertex_h
#define __vertex_h

#include "numfixarray.h"

/*------------------------------------------------------------*/
namespace Gascoigne
{

template<int N>
class Vertex : public numfixarray<N,double>
{
public:

  Vertex<N>() : numfixarray<N,double>() {}
  Vertex<N>(const Vertex& c) : numfixarray<N,double>(c) {}
  Vertex<N>(const double& x0) : numfixarray<N,double>(x0) {}
  Vertex<N>(const double& x0, const double& y0) : numfixarray<N,double>() 
    {x()=x0; y()=y0;}
  Vertex<N>(const double& x0, const double& y0, const double& z0) : numfixarray<N,double>() 
    {x()=x0; y()=y0; z()=z0;}
 
  Vertex<N>& operator=(const Vertex<N>& c) 
    {
      numfixarray<N,double>::operator=(c);
      return *this;
    }
  Vertex<N>& operator=(double d) 
    {
      numfixarray<N,double>::operator=(d);
      return *this;
    }

  const double& x() const  { return (*this)[0]; }
  const double& y() const  { return (*this)[1]; }
  const double& z() const  { return (*this)[2]; }
  double&       x()        { return (*this)[0]; }
  double&       y()        { return (*this)[1]; }
  double&       z()        { return (*this)[2]; }
};

/*------------------------------------------------------------*/

typedef  Vertex<1>   Vertex1d; 
typedef  Vertex<2>   Vertex2d; 
typedef  Vertex<3>   Vertex3d; 
}

#endif
