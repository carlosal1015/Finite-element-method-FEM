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


/*----------------------------   faceintegrator.h     ---------------------------*/
/*      $Id$                 */
#ifndef __faceintegrator_H
#define __faceintegrator_H
/*----------------------------   faceintegrator.h     ---------------------------*/

#include "faceequation.h"
#include "feminterface.h"
#include "integrationformula.h"

namespace Gascoigne
{
  class FaceIntegratorInterface
    {

    protected:
      IntegrationFormulaInterface*  __IFF;
      
      mutable TestFunction          __NN;
      
      IntegrationFormulaInterface*& FaceFormulaPointer() { return __IFF;}
      const IntegrationFormulaInterface* FaceFormula() const { assert(__IFF); return __IFF; }

      // --------------------------------------------------

      virtual void universal_point_face(const FemInterface& FEM1,const FemInterface& FEM2, FemFunction& U1, FemFunction& U2, const LocalVector& U) const=0;


    public:

      FaceIntegratorInterface() : __IFF(0) {}
      virtual ~FaceIntegratorInterface() { }
      
      
      virtual void FaceForm(const FaceEquation& EQ, LocalVector& F, const FemInterface& FEM1, const FemInterface& FEM2, const LocalVector& U) const=0;
      virtual void FaceMatrix(const FaceEquation& EQ, EntryMatrix& E, const FemInterface& FEM1, const FemInterface& FEM2, const LocalVector& U) const=0;

    };

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////

  template<int DIM, int Q>
    class FaceIntegrator : public FaceIntegratorInterface
    {
    protected:
      void universal_point_face(const FemInterface& FEM1,const FemInterface& FEM2, FemFunction& U1, FemFunction& U2, const LocalVector& U) const;
      
    public:
      FaceIntegrator();
      ~FaceIntegrator()
	{
	  if (__IFF) { delete __IFF; __IFF=0; }
	}
      
      void FaceForm(const FaceEquation& EQ, LocalVector& F, const FemInterface& FEM1, const FemInterface& FEM2, const LocalVector& U) const;
      void FaceMatrix(const FaceEquation& EQ, EntryMatrix& E, const FemInterface& FEM1, const FemInterface& FEM2, const LocalVector& U) const;
    };
  
  
}




/*----------------------------   faceintegrator.h     ---------------------------*/
/* end of #ifndef __faceintegrator_H */
#endif
/*----------------------------   faceintegrator.h     ---------------------------*/
