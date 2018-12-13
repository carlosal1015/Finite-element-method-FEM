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


#include "faceq1.h"
#include "gascoignemesh2d.h"
#include "baseq12d.h"
#include "baseq13d.h"
#include "finiteelement.h"
#include "transformation2d.h"
#include "transformation3d.h"

using namespace std;


namespace Gascoigne
{


  

  
  //////////////////////////////////////////////////

  template<int DIM>
  void FaceQ1<DIM>::BasicInit(const ParamFile* pf)
  {
    if (!GetFem1Pointer())
      {
	if (DIM==2)
	  GetFem1Pointer() = new FiniteElement<2,1,Transformation2d<BaseQ12d>,BaseQ12d >;
	if (DIM==3)
	  GetFem1Pointer() = new FiniteElement<3,2,Transformation3d<BaseQ13d>,BaseQ13d >;
      }
    if (!GetFem2Pointer())
      {
	if (DIM==2)
	  GetFem2Pointer() = new FiniteElement<2,1,Transformation2d<BaseQ12d>,BaseQ12d >;
	if (DIM==3)
	  GetFem2Pointer() = new FiniteElement<3,2,Transformation3d<BaseQ13d>,BaseQ13d >;
      }    

    if (!GetFaceIntegratorPointer())
      GetFaceIntegratorPointer() = new FaceIntegrator<DIM,1>;
  }

  //////////////////////////////////////////////////
  
  template<int DIM>
  void FaceQ1<DIM>::ReInit   (const MeshInterface* M)
  {
    // build face/edge-structur
    __MP = M;
    assert(__MP);

    build_faces();
  }
  
  //////////////////////////////////////////////////
  
  template<>
  void FaceQ1<2>::build_faces() 
  {
    this->__faces.clear();
    const GascoigneMesh2d* GM = dynamic_cast<const GascoigneMesh2d*> (__MP);

    // list of neighbors
    std::map<pair<int,int> ,fixarray<2,int> > fm;
    for (int c=0;c<GM->ncells();++c)
      {
	for (int i=0;i<4;++i)
	  {
	    pair<int,int> f=make_pair<int,int> (GM->vertex_of_cell(c,i),GM->vertex_of_cell(c,(i+1)%4));
	    if (f.first>f.second) std::swap(f.first,f.second);
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
    for (std::map<pair<int,int> ,fixarray<2,int> >::const_iterator it = fm.begin();it!=fm.end();++it)
      {
	if (it->second[1]==-1) continue;
	assert(it->second[0]!=-1);
	int c0=it->second[0];
	int c1=it->second[1];
	
	nvector<int> ei(6,-1);
	// left
	int n0=-1,n1=-1,dir = 1;
	for (n0=0;n0<4;++n0) if (GM->vertex_of_cell(c0,n0)==it->first.first) break;
	assert(n0<4);
	for (n1=0;n1<4;++n1) if (GM->vertex_of_cell(c0,n1)==it->first.second) break;
	assert(n1<4);
	if ((n0+dir)%4!=n1) dir=-1;
	assert((4+n0+dir)%4==n1);
	ei[0]=GM->vertex_of_cell(c0,(4+n0-dir)%4);
	ei[1]=GM->vertex_of_cell(c0,n0);
	ei[3]=GM->vertex_of_cell(c0,(4+n1+dir)%4);
	ei[4]=GM->vertex_of_cell(c0,n1);

	n0=-1;n1=-1;dir=1;
	for (n0=0;n0<4;++n0) if (GM->vertex_of_cell(c1,n0)==it->first.first) break;
	assert(n0<4);
	for (n1=0;n1<4;++n1) if (GM->vertex_of_cell(c1,n1)==it->first.second) break;
	assert(n1<4);
	if ((n0+dir)%4!=n1) dir=-1;
	assert((4+n0+dir)%4==n1);
	assert(ei[1]==GM->vertex_of_cell(c1,n0));
	assert(ei[4]==GM->vertex_of_cell(c1,n1));
	ei[2]=GM->vertex_of_cell(c1,(4+n0-dir)%4);
	ei[5]=GM->vertex_of_cell(c1,(4+n1+dir)%4);

	this->__faces.push_back(ei);
      }

    cout << "FaceQ1: build " << this->__faces.size() << " edges" << endl;
  }

  
  template<>
  void FaceQ1<3>::build_faces() 
  {
  }

  //////////////////////////////////////////////////

  template<int DIM>
  void FaceQ1<DIM>::TransformationFace(FemInterface::Matrix& T1,FemInterface::Matrix& T2,int f) const
  {
    int dim = this->GetMesh()->dimension();
    int ne  = this->GetMesh()->nodes_per_cell(0);
  
    const nvector<int>& ei = GetFace(f);
    assert(ei.size()==6*(DIM-1));
  
    T1.memory(dim,ne);
    T2.memory(dim,ne);
    assert(ne==(1<<DIM));
  
    if(dim==2)
      {
	Vertex2d v;
	for (int iy=0;iy<2;++iy)
	  for (int ix=0;ix<2;++ix)
	    {
	      v = this->GetMesh()->vertex2d(ei[3*iy+ix]);
	      T1(0,2*iy+ix) = v.x();               
	      T1(1,2*iy+ix) = v.y();
	      v = this->GetMesh()->vertex2d(ei[3*iy+ix+1]);
	      T2(0,2*iy+ix) = v.x();               
	      T2(1,2*iy+ix) = v.y();	      
	    }
      }
    else if(dim==3)
      {
	Vertex3d v;
	for (int iz=0;iz<2;++iz)
	  for (int iy=0;iy<2;++iy)
	    for (int ix=0;ix<2;++ix)
	    {
	      v = this->GetMesh()->vertex3d(ei[6*iz+3*iy+ix]);
	      T1(0,4*iz+2*iy+ix) = v.x();
	      T1(1,4*iz+2*iy+ix) = v.y();
	      T1(2,4*iz+2*iy+ix) = v.z();
	      v = this->GetMesh()->vertex3d(ei[6*iz+3*iy+ix+1]);
	      T2(0,4*iz+2*iy+ix) = v.x();
	      T2(1,4*iz+2*iy+ix) = v.y();
	      T2(2,4*iz+2*iy+ix) = v.z();
	    }
      }
  }
  


  template class FaceQ1<2>;
  template class FaceQ1<3>;
  

}
