/**
*
* Copyright (C) 2007 by the Gascoigne 3D authors
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


#include "nsface.h"
#include "filescanner.h"


using namespace std;

namespace Gascoigne
{
  
  
  NavierStokesFace2d::NavierStokesFace2d(const ParamFile* pf) : FaceEquation()
  {
    DataFormatHandler DFH;
    DFH.insert("visc"   , &__visc , 1.);
    DFH.insert("alpha"  , &__alpha0,0.5);
    DFH.insert("delta"  , &__delta0, 0.);
    
    FileScanner FS(DFH, pf, "Equation");
  }

  // --------------------------------------------------

  void NavierStokesFace2d::point_face(double h, const FemFunction& U1,const FemFunction& U2, const Vertex2d& v,const Vertex2d& n) const 
  {
    __n=n;
    __alpha = __alpha0 / (1+__visc) * h*h;
  }
  
  // --------------------------------------------------

  void NavierStokesFace2d::FaceForm(VectorIterator b, const FemFunction& U1,const FemFunction& U2, const TestFunction& N1,const TestFunction& N2) const
  {
    double un=0,pn=0;
    for (int i=0;i<2;++i)
      {
	un += __n[i]*(U1[0][i+1]-U2[0][i+1]);
	pn += __n[i]*(N1[i+1]   -N2[i+1]);
      }
    
    b[0] += __alpha*un*pn;
  }
  
  
  void NavierStokesFace2d::FaceMatrix(EntryMatrix& A, const FemFunction& U1,const FemFunction& U2, const TestFunction& M1,const TestFunction& M2, const TestFunction& N1, const TestFunction& N2) const
  {
    double un=0,pn=0;
    for (int i=0;i<2;++i)
      {
	un += __n[i]*(M1[i+1]-M2[i+1]);
	pn += __n[i]*(N1[i+1]-N2[i+1]);
      }
    
    A(0,0) += __alpha*un*pn;
  }
  

}
