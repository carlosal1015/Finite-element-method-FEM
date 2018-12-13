/**
*
* Copyright (C) 2004, 2005, 2006 by the Gascoigne 3D authors
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


#include "edgeinfocontainer.h"

using namespace std;

/**********************************************************/

namespace Gascoigne
{
template<int DIM>
EdgeInfoContainer<DIM>::~EdgeInfoContainer<DIM>()
{
  for (int i=0; i<nvector<EdgeInfo<DIM>*>::size(); i++)
    {
      if ((*this)[i])
	{
	  delete (*this)[i];
	  (*this)[i] = NULL;
	}
    }
  nvector<EdgeInfo<DIM>*>::resize(0);
}

/**********************************************************/

template<int DIM>
void EdgeInfoContainer<DIM>::BasicInit(const HierarchicalMesh* HM, int ncomp)
{
  nvector<EdgeInfo<DIM>*>::resize(0);
  _HMP   = HM;
  _ncomp = ncomp;

  nvector<EdgeInfo<DIM>*>::resize(_HMP->nedges());
  for (int i=0; i<nvector<EdgeInfo<DIM>*>::size(); i++)
    {
      (*this)[i]=NULL;
    }
}

/**********************************************************/

template<>
void EdgeInfoContainer<2>::ModifyHanging()
{
  const QuadLawAndOrder& QLAO = dynamic_cast<const HierarchicalMesh2d*>(_HMP)->QuadLawOrder();
  fixarray<2,int>        vertexes;
  LocalVector            lu,lul,lur;
  
  lu.ReInit(_ncomp,2);
  lu.zero();
  lul.ReInit(_ncomp,2);
  lul.zero();
  lur.ReInit(_ncomp,2);
  lur.zero();
  
  for (int i=0; i<size(); i++)
    {
      if ((*this)[i]!=NULL && (*this)[i]->GetCount()==1 && (*this)[i]->GetEdge().slave()!=-1)
	{
	  const Edge& edge = (*this)[i]->GetEdge();
	  vertexes = (*this)[i]->GetVertex();
	  
	  int left  = QLAO.GlobalChildEdge(vertexes,edge.slave(),0);
	  int right = QLAO.GlobalChildEdge(vertexes,edge.slave(),1);
	  
	  for (int c=0; c<_ncomp; c++)
	    {
	      lu(0,c) = ((*this)[right]->GetValue())(1,c);
	      lu(1,c) = ((*this)[left]->GetValue())(0,c);
	      
	      lul(0,c) = ((*this)[i]->GetValue())(1,c);
	      lul(1,c) = 0.5 * (((*this)[i]->GetValue())(0,c) + ((*this)[i]->GetValue())(1,c));
	      
	      lur(0,c) = 0.5 * (((*this)[i]->GetValue())(0,c) + ((*this)[i]->GetValue())(1,c));
	      lur(1,c) = ((*this)[i]->GetValue())(0,c);
	    }
	  (*this)[i]->AddNodes(lu);
	  (*this)[left]->AddNodes(lul);
	  (*this)[right]->AddNodes(lur);
	}
    }
}

/**********************************************************/

template<>
void EdgeInfoContainer<3>::ModifyHanging()
{
  const HexLawAndOrder& HLAO = dynamic_cast<const HierarchicalMesh3d*>(_HMP)->HexLawOrder();
  fixarray<4,int>       vertexes;
  LocalVector           lugr,lulu,luru,lulo,luro;
  
  lugr.ReInit(_ncomp,4);
  lugr.zero();
  lulu.ReInit(_ncomp,4);
  lulu.zero();
  luru.ReInit(_ncomp,4);
  luru.zero();
  lulo.ReInit(_ncomp,4);
  lulo.zero();
  luro.ReInit(_ncomp,4);
  luro.zero();
  
  for (int i=0; i<size(); i++)
    {
      if ((*this)[i]!=NULL && (*this)[i]->GetCount()==1 && (*this)[i]->GetEdge().slave()!=-1)
	{
	  const Edge& quad = (*this)[i]->GetEdge();
	  vertexes = (*this)[i]->GetVertex();
	  
	  int lu = HLAO.GlobalChildFace(vertexes,quad.master(),0);
	  int ru = HLAO.GlobalChildFace(vertexes,quad.master(),1);
	  int lo = HLAO.GlobalChildFace(vertexes,quad.master(),3);
	  int ro = HLAO.GlobalChildFace(vertexes,quad.master(),2);
	  
	  LocalVector help = (*this)[i]->GetValue();
	  
	  for (int c=0; c<_ncomp; c++)
	    {
	      lugr(0,c) = ((*this)[lu]->GetValue())(0,c);
	      lugr(1,c) = ((*this)[ru]->GetValue())(1,c);
	      lugr(2,c) = ((*this)[ro]->GetValue())(2,c);
	      lugr(3,c) = ((*this)[lo]->GetValue())(3,c);
	      
	      lulu(0,c) = help(0,c);
	      lulu(1,c) = 0.5 * (help(0,c) + help(1,c));
	      lulu(2,c) = 0.25 * (help(0,c) + help(1,c) + help(2,c) + help(3,c));
	      lulu(3,c) = 0.5 * (help(0,c) + help(3,c));
	      
	      luru(0,c) = 0.5 * (help(0,c) + help(1,c));
	      luru(1,c) = help(1,c);
	      luru(2,c) = 0.5 * (help(1,c) + help(2,c));
	      luru(3,c) = 0.25 * (help(0,c) + help(1,c) + help(2,c) + help(3,c));
	      
	      lulo(0,c) = 0.5 * (help(0,c) + help(3,c));
	      lulo(1,c) = 0.25 * (help(0,c) + help(1,c) + help(2,c) + help(3,c));
	      lulo(2,c) = 0.5 * (help(2,c) + help(3,c));
	      lulo(3,c) = help(3,c);
	      
	      luro(0,c) = 0.25 * (help(0,c) + help(1,c) + help(2,c) + help(3,c));
	      luro(1,c) = 0.5 * (help(1,c) + help(2,c));
	      luro(2,c) = help(2,c);
	      luro(3,c) = 0.5 * (help(2,c) + help(3,c));
	    }
	  (*this)[i]->AddNodes(lugr);
	  (*this)[lu]->AddNodes(lulu);
	  (*this)[ru]->AddNodes(luru);
	  (*this)[lo]->AddNodes(lulo);
	  (*this)[ro]->AddNodes(luro);
	}
    }
}

/**********************************************************/

template<int DIM>
void EdgeInfoContainer<DIM>::ShowStatistics() const
{
  for (int i=0; i<nvector<EdgeInfo<DIM>*>::size(); i++)
    {
      if ((*this)[i]!=NULL)
	{
	  cout << "Edge " << i << ": ";
	  (*this)[i]->ShowStatistics();
	}
    }
}

/**********************************************************/

template class EdgeInfoContainer<2>;
template class EdgeInfoContainer<3>;
}
