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


#include "faceq2.h"
#include "gascoignemesh2d.h"
#include "baseq22d.h"
#include "baseq23d.h"
#include "finiteelement.h"
#include "transformation2d.h"
#include "transformation3d.h"

using namespace std;


namespace Gascoigne
{
  
  //////////////////////////////////////////////////

  template<int DIM>
  void FaceQ2<DIM>::BasicInit(const ParamFile* pf)
  {
    if (!GetFem1Pointer())
      {
	if (DIM==2)
	  GetFem1Pointer() = new FiniteElement<2,1,Transformation2d<BaseQ22d>,BaseQ22d >;
	if (DIM==3)
	  GetFem1Pointer() = new FiniteElement<3,2,Transformation3d<BaseQ23d>,BaseQ23d >;
      }
    if (!GetFem2Pointer())
      {
	if (DIM==2)
	  GetFem2Pointer() = new FiniteElement<2,1,Transformation2d<BaseQ22d>,BaseQ22d >;
	if (DIM==3)
	  GetFem2Pointer() = new FiniteElement<3,2,Transformation3d<BaseQ23d>,BaseQ23d >;
      }

    if (!GetFaceIntegratorPointer())
      GetFaceIntegratorPointer() = new FaceIntegrator<DIM,2>;
  }

  //////////////////////////////////////////////////
  
  template<int DIM>
  void FaceQ2<DIM>::ReInit   (const MeshInterface* M)
  {
    // build face/edge-structur
    __MP = M;
    assert(__MP);
    build_faces();
  }
  
  //////////////////////////////////////////////////
  
  template<>
  void FaceQ2<2>::build_faces() 
  {
    this->__faces.clear();
    const GascoigneMesh2d* GM = dynamic_cast<const GascoigneMesh2d*> (__MP);
    assert(GM->HasPatch());
    
    // list of neighbors (middle-node-of-line-of-patch)
    std::map<int ,fixarray<2,int> > fm;
    int ind[4] = {1,5,7,3};
    for (int c=0;c<GM->npatches();++c)
      {
	const nvector<int>& iop = *(GM->IndicesOfPatch(c));
	for (int i=0;i<4;++i)
	  {
	    int f = iop[ind[i]];
	    if (fm.find(f)==fm.end()) { fm[f][0]=c; fm[f][1]=-1; }
	    else
	      {
		assert(fm[f][1]==-1);
		assert(fm[f][0]!=-1);
		fm[f][1]=c;
	      }
	  }
      }

    // create list of internal faces
    for (std::map<int ,fixarray<2,int> >::const_iterator it = fm.begin();it!=fm.end();++it)
      {
	if (it->second[1]==-1) continue;
	assert(it->second[0]!=-1);
	int c0=it->second[0];
	int c1=it->second[1];
	const nvector<int>& I0 = *(GM->IndicesOfPatch(c0));
	const nvector<int>& I1 = *(GM->IndicesOfPatch(c1));

	int f = it->first;
	
	nvector<int> ei(15,-1);
	// left
	int n0=-1;
	for (n0=0;n0<4;++n0) if (I0[ind[n0]]==f) break;
	assert(n0<4);
	for (int iy=0;iy<3;++iy)
	  for (int ix=0;ix<3;++ix)
	    if (n0==0)      ei[5*iy+ix]=I0[6-3*ix+iy];
	    else if (n0==1) ei[5*iy+ix]=I0[3*iy+ix];
	    else if (n0==2) ei[5*iy+ix]=I0[2+3*ix-iy];
	    else if (n0==3) ei[5*iy+ix]=I0[8-ix-3*iy];
	
	// right
	n0=-1;
	for (n0=0;n0<4;++n0) if (I1[ind[n0]]==f) break;
	assert(n0<4);

	if (n0==0)      assert((ei[2]==I1[2]));
	else if (n0==1) assert((ei[2]==I1[8]));
	else if (n0==2) assert((ei[2]==I1[6]));
	else if (n0==3) assert((ei[2]==I1[0]));
	for (int iy=0;iy<3;++iy)
	  for (int ix=1;ix<3;++ix)
	    if (n0==0)      ei[5*iy+ix+2]=I1[2+3*ix-iy];
	    else if (n0==1) ei[5*iy+ix+2]=I1[8-3*iy-ix];
	    else if (n0==2) ei[5*iy+ix+2]=I1[6-3*ix+iy];
	    else if (n0==3) ei[5*iy+ix+2]=I1[3*iy+ix];
	
	this->__faces.push_back(ei);
      }
    cout << "FaceQ2: build " << this->__faces.size() << " edges" << endl;
  }

  
  template<>
  void FaceQ2<3>::build_faces() 
  {
  }

  //////////////////////////////////////////////////

  template<int DIM>
  void FaceQ2<DIM>::TransformationFace(FemInterface::Matrix& T1,FemInterface::Matrix& T2,int f) const
  {
    const PatchMesh* PM = dynamic_cast<const PatchMesh*> (GetMesh());
    assert(PM);
    
    int dim = PM->dimension();
    int ne  = PM->nodes_per_patch();
  
    const nvector<int>& ei = this->GetFace(f);
    assert(ei.size()==15+30*(DIM-2));
  
    T1.memory(dim,ne);
    T2.memory(dim,ne);
    assert(ne==18*(DIM-2)+9);
  
    if(dim==2)
      {
	Vertex2d v;
	for (int iy=0;iy<3;++iy)
	  for (int ix=0;ix<3;++ix)
	    {
	      v = this->GetMesh()->vertex2d(ei[5*iy+ix]);
	      T1(0,3*iy+ix) = v.x();               
	      T1(1,3*iy+ix) = v.y();
	      v = this->GetMesh()->vertex2d(ei[5*iy+ix+2]);
	      T2(0,3*iy+ix) = v.x();               
	      T2(1,3*iy+ix) = v.y();	      
	    }
      }
    else if(dim==3)
      {
	Vertex3d v;
	for (int iz=0;iz<3;++iz)
	  for (int iy=0;iy<3;++iy)
	    for (int ix=0;ix<3;++ix)
	    {
	      v = this->GetMesh()->vertex3d(ei[15*iz+5*iy+ix]);
	      T1(0,9*iz+3*iy+ix) = v.x();
	      T1(1,9*iz+3*iy+ix) = v.y();
	      T1(2,9*iz+3*iy+ix) = v.z();
	      v = this->GetMesh()->vertex3d(ei[15*iz+5*iy+ix+1]);
	      T2(0,9*iz+3*iy+ix) = v.x();
	      T2(1,9*iz+3*iy+ix) = v.y();
	      T2(2,9*iz+3*iy+ix) = v.z();
	    }
      }
  }
  
  //////////////////////////////////////////////////

  template class FaceQ2<2>;
  template class FaceQ2<3>;
  

}
