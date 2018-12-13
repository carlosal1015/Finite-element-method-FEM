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


#ifndef  __LpsEquation_h
#define  __LpsEquation_h

#include "equation.h"

namespace Gascoigne
{

/*-----------------------------------------*/

  //////////////////////////////////////////////
  //
  ///@brief
  /// Interface class for Lps Elements
  ///
  ///
  //////////////////////////////////////////////

  class LpsEquation : public virtual Equation
  {
    private:

    protected:
      
    public:
      LpsEquation() {}
      ~LpsEquation() {}

      virtual void lpspoint(double h, const FemFunction& U, const Vertex2d& v) const {
        std::cerr << "\"LpsEquation::lpspoint\" not written!" << std::endl;
        abort();
      } 
      virtual void lpspoint(double h, const FemFunction& U, const Vertex3d& v) const {
        std::cerr << "\"LpsEquation::lpspoint\" not written!" << std::endl;
        abort();
      } 
   
      virtual void lpspointmatrix(double h, const FemFunction& U, const Vertex2d& v) const {
        lpspoint(h,U,v);
      }
      virtual void lpspointmatrix(double h, const FemFunction& U, const Vertex3d& v) const {
        lpspoint(h,U,v);
      }

      virtual void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const=0;
      virtual void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const=0;
  };
}

#endif
