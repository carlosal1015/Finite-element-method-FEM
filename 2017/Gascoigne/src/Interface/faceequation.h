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


/*----------------------------   faceequation.h     ---------------------------*/
/*      $Id$                 */
#ifndef __faceequation_H
#define __faceequation_H
/*----------------------------   faceequation.h     ---------------------------*/

#include  "entrymatrix.h"
#include  "vertex.h"
#include  "application.h"


/*-------------------------------------------------------------------------*/

namespace Gascoigne
{
  
  //////////////////////////////////////////////
  ///
  ///@brief
  /// Interface class for Equation

  ///
  ///
  //////////////////////////////////////////////

  class FaceEquation : public virtual Application
    {
    private:
      
    protected:
      
    public:
      //
      // Constructors
      //
      FaceEquation() : Application() {}
      virtual ~FaceEquation() {}


      virtual void point_face(double h, const FemFunction& U1,const FemFunction& U2, const Vertex2d& v,const Vertex2d& n) const {}
      virtual void point_face(double h, const FemFunction& U1,const FemFunction& U2, const Vertex3d& v,const Vertex3d& n) const {}
     
      virtual void pointmatrix_face(double h, const FemFunction& U1,const FemFunction& U2, const Vertex2d& v,const Vertex2d& n) const
	{
	  point_face(h,U1,U2,v,n);
	}
      virtual void pointmatrix_face(double h, const FemFunction& U1,const FemFunction& U2, const Vertex3d& v,const Vertex3d& n) const
	{
	  point_face(h,U1,U2,v,n);
	}
      
      virtual int GetNcomp() const=0;
      virtual void FaceForm(VectorIterator b, const FemFunction& U1, const FemFunction& U2, const TestFunction& N1, const TestFunction& N2) const=0;
      virtual void FaceMatrix(EntryMatrix& A, const FemFunction& U1, const FemFunction& U2, const TestFunction& M1, const TestFunction& M2, const TestFunction& N1, const TestFunction& N2) const=0;
  };
}

/*----------------------------   faceequation.h     ---------------------------*/
/* end of #ifndef __faceequation_H */
#endif
/*----------------------------   faceequation.h     ---------------------------*/
