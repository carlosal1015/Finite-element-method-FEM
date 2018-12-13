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


#ifndef  __GlsEquation_h
#define  __GlsEquation_h

#include "equation.h"

namespace Gascoigne
{
  
  //////////////////////////////////////////////
  //
  ///@brief
  /// Interface class for  Elements

  ///
  ///
  //////////////////////////////////////////////
  
  class GlsEquation : public virtual Equation
  {
    public:
    GlsEquation() {};
      ~GlsEquation() {};

      //
      /// computation of stabilization parameters
      //
      virtual void glspoint(double h, const FemFunction& U, const Vertex2d& v) const {
        std::cerr << "\"GlsEquation::glspoint\" not written!" << std::endl;
        abort();
      } 
      
      //
      /// computation of stabilization parameters
      //
      virtual void glspoint(double h, const FemFunction& U, const Vertex3d& v) const {
        std::cerr << "\"GlsEquation::glspoint\" not written!" << std::endl;
        abort();
      }
     
      //
      /// computation of stabilization parameters for the matrix
      //
      virtual void glspointmatrix(double h, const FemFunction& U, const Vertex2d& v) const {
        glspoint(h,U,v);
      }
      
      //
      /// computation of stabilization parameters for the matrix
      //
      virtual void glspointmatrix(double h, const FemFunction& U, const Vertex3d& v) const {
        glspoint(h,U,v);
      }
      
      //
      /// describes the strong form of the PDE
      //
      virtual void L(DoubleVector& b, const FemFunction& U) const=0;
      
      //
      /// describes the stabilization term of the PDE;
      /// can be chosen as -L^
      //
      virtual void S(DoubleMatrix& A, const FemFunction& U, const TestFunction& N) const=0;

      //
      /// describes the strong derivative of the PDE
      //
      virtual void LMatrix(DoubleMatrix& A, const FemFunction& U, const TestFunction& M) const=0;

      //
      /// describes the derivative of the stabilization term S;
      //
      virtual void SMatrix(DoubleVector& b, const FemFunction& U, const FemFunction& M, const FemFunction& N) const {};
  };
}

/*-----------------------------------------*/

#endif
