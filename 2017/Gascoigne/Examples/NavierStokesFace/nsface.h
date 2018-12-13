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


/*----------------------------   nsface.h     ---------------------------*/
/*      $Id$                 */
#ifndef __nsface_H
#define __nsface_H
/*----------------------------   nsface.h     ---------------------------*/

#include "faceequation.h"
#include "paramfile.h"


using namespace std;

namespace Gascoigne
{
  
  class NavierStokesFace2d : public FaceEquation
    {
    protected:
      double __alpha0,__delta0,__visc;
      mutable double __alpha,__delta;
      
      mutable Vertex2d __n;
      
      
    public:
      
      NavierStokesFace2d(const ParamFile* pf);

      std::string GetName() const { return "NavierStokesFace2d";}
      int  GetNcomp() const { return 3; }
      
      void point_face(double h, const FemFunction& U1,const FemFunction& U2, const Vertex2d& v,const Vertex2d& n) const;
            
      void FaceForm(VectorIterator b, const FemFunction& U1,const FemFunction& U2, const TestFunction& N1,const TestFunction& N2) const;
      void FaceMatrix(EntryMatrix& A, const FemFunction& U1,const FemFunction& U2, const TestFunction& M1,const TestFunction& M2, const TestFunction& N1, const TestFunction& N2) const;
      
      
    };
  
}




/*----------------------------   nsface.h     ---------------------------*/
/* end of #ifndef __nsface_H */
#endif
/*----------------------------   nsface.h     ---------------------------*/
