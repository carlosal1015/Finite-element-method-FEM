/**
*
* Copyright (C) 2007, 2011 by the Gascoigne 3D authors
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


/*----------------------------   facediscretization.h     ---------------------------*/
/*      $Id$                 */
#ifndef __facediscretization_H
#define __facediscretization_H
/*----------------------------   facediscretization.h     ---------------------------*/


#include  <string>
#include  "meshinterface.h"
#include  "gascoigne.h"
#include  "equation.h"
#include  "matrixinterface.h"
#include  "paramfile.h"
#include  "feminterface.h"
#include  "faceintegrator.h"

namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments DiscretizationInterface

  ///
  ///
  /////////////////////////////////////////////

  class FaceDiscretization
    {
    protected:
      
      const MeshInterface*  __MP;

      nvector<nvector<int> > __faces;

      FemInterface*            __FEM1;
      FemInterface*            __FEM2;
      FaceIntegratorInterface* __INT;

      mutable EntryMatrix __E;
      mutable LocalVector __U,__F;
      
      // --------------------------------------------------
      
      virtual void TransformationFace(FemInterface::Matrix& T1,FemInterface::Matrix& T2,int f) const
	{std::cerr << "FaceDiscretization::Transformation not written!" << std::endl; abort();}
      virtual void LocalToGlobalFace(GlobalVector& f, const LocalVector& F, int iq, double s) const;
      virtual void GlobalToLocalFace(LocalVector& U, const GlobalVector& u, int iq) const;
      virtual void LocalToGlobalFace(MatrixInterface& A, EntryMatrix& E, int iq, double s) const;
      
   public:
      FaceDiscretization() : __MP(0),__FEM1(0),__FEM2(0),__INT(0) {}
      virtual ~FaceDiscretization() {}


      // init
      virtual std::string GetName() const { return "FaceDiscretizationInterface"; }
      
      virtual void BasicInit(const ParamFile* pf)
	{ std::cerr << "\"FaceDiscretization::BasicInit not written\"" << std::endl;  abort(); } 
      virtual void ReInit   (const MeshInterface* M)
	{ std::cerr << "\"FaceDiscretization::ReInit not written\"" << std::endl;  abort(); } 
      virtual void build_faces()
	{ std::cerr << "\"FaceDiscretization::build_faces not written\"" << std::endl;  abort(); } 
      virtual void Structure(SparseStructureInterface* SI) const;
	

      // acc
      int nfaces() const { return __faces.size(); }
      const nvector<int>& GetFace(int i) const { return __faces[i]; }
      const MeshInterface* GetMesh() const { return __MP; }

      virtual const FemInterface* GetFem1() const {assert(__FEM1); return __FEM1;}
      virtual const FemInterface* GetFem2() const {assert(__FEM2); return __FEM2;}
      virtual FemInterface*& GetFem1Pointer() {return __FEM1;}
      virtual FemInterface*& GetFem2Pointer() {return __FEM2;}

      virtual const FaceIntegratorInterface* GetFaceIntegrator() const {assert(__INT); return __INT;}
      virtual FaceIntegratorInterface*&      GetFaceIntegratorPointer() {return __INT;}

      //
      //// Functions called from the Solver
      //
      
      virtual void FaceForm(GlobalVector& f, const GlobalVector& u, const FaceEquation& EQ, double d) const;
      virtual void FaceMatrix(MatrixInterface& A, const GlobalVector& u, const FaceEquation& EQ, double) const;
    };
}


/*----------------------------   facediscretization.h     ---------------------------*/
/* end of #ifndef __facediscretization_H */
#endif
/*----------------------------   facediscretization.h     ---------------------------*/
