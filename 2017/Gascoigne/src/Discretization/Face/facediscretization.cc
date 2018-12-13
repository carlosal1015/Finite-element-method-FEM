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



#include "facediscretization.h"
#include "faceequation.h"
#include "sparsestructure.h"

using namespace std;

namespace Gascoigne
{


  void FaceDiscretization::Structure(SparseStructureInterface* SI) const
  {
    SparseStructure* S = dynamic_cast<SparseStructure*>(SI);
    assert(S);
    for(int iq=0;iq<nfaces();++iq)
      {
	IntVector indices = GetFace(iq);
	S->build_add(indices.begin(), indices.end());
      }
    S->build_end();
  }

  /* ----------------------------------------- */  

  void FaceDiscretization::GlobalToLocalFace(LocalVector& U, const GlobalVector& u, int iq) const
  {
    const nvector<int> ei = GetFace(iq);
    U.ReInit(u.ncomp(),ei.size());
    for(int ii=0; ii<ei.size(); ii++) 
      {
	int i = ei[ii];
	U.equ_node(ii,i,u);
      }
  } 

  /* ----------------------------------------- */
  
  void FaceDiscretization::LocalToGlobalFace(GlobalVector& f, const LocalVector& F, int iq, double s) const
  {
    const nvector<int> ei = GetFace(iq);
    for(int ii=0; ii<ei.size(); ii++) 
      {
	int i = ei[ii];
	f.add_node(i,s,ii,F);
      }
  }

  /* ----------------------------------------- */
  
  void FaceDiscretization::LocalToGlobalFace(MatrixInterface& A, EntryMatrix& E, int iq, double s) const
  {
    const nvector<int> indices = GetFace(iq);
    IntVector::const_iterator  start = indices.begin();
    IntVector::const_iterator  stop  = indices.end();
    A.entry(start,stop,E,s);
  }

  
  //////////////////////////////////////////////////

  void FaceDiscretization::FaceForm(GlobalVector& f, const GlobalVector& u, const FaceEquation& FEQ, double d) const
  {
    nmatrix<double> T1,T2;
  
    for(int iq=0;iq<nfaces();++iq)
      {
	TransformationFace(T1,T2,iq);
	GetFem1()->ReInit(T1);
	GetFem2()->ReInit(T2);
	GlobalToLocalFace(__U,u,iq);
	
	GetFaceIntegrator()->FaceForm(FEQ,__F,*GetFem1(),*GetFem2(),__U);
	
	LocalToGlobalFace(f,__F,iq,d);
      }
  }

  /* ----------------------------------------- */

  void FaceDiscretization::FaceMatrix(MatrixInterface& A, const GlobalVector& u, const FaceEquation& FEQ, double d) const
  {

    nmatrix<double> T1,T2;  
    
    for(int iq=0;iq<nfaces();++iq)
      {
	TransformationFace(T1,T2,iq);
	GetFem1()->ReInit(T1);
	GetFem2()->ReInit(T2);
	GlobalToLocalFace(__U,u,iq);
	
	GetFaceIntegrator()->FaceMatrix(FEQ,__E,*GetFem1(),*GetFem2(),__U);
	
	LocalToGlobalFace(A,__E,iq,d);
      }
  }



}
