/**
*
* Copyright (C) 2004, 2011 by the Gascoigne 3D authors
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


#include  "hnstructureq12d.h"
#include  "gascoignemesh.h"


using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
HNStructureQ12d::HNStructureQ12d() : edges(NULL), wei(3)
{
  wei[0] =  0.5; 
  wei[1] =  0.5; 
  wei[2] =  0.; 

  lnoe[0][0]=0; lnoe[0][1]=1; lnoe[0][2]=3;
  lnoe[1][0]=1; lnoe[1][1]=3; lnoe[1][2]=2;
  lnoe[2][0]=3; lnoe[2][1]=2; lnoe[2][2]=0;
  lnoe[3][0]=2; lnoe[3][1]=0; lnoe[3][2]=1;

  lnop[0][0]=0; lnop[0][1]=2; lnop[0][2]=1;
  lnop[1][0]=2; lnop[1][1]=8; lnop[1][2]=5;
  lnop[2][0]=8; lnop[2][1]=6; lnop[2][2]=7;
  lnop[3][0]=6; lnop[3][1]=0; lnop[3][2]=3;
}

/*-----------------------------------------*/

void HNStructureQ12d::SparseStructureDiag(SparseStructure* S) const
{
  assert(edges);
  for(const_iterator p=edges->begin();p!=edges->end();p++)
    {
      int i = p->first;
      S->build_add(i,i);
    }
}

/*--------------------------------------------------------*/

void HNStructureQ12d::ReInit(const MeshInterface* M)
{
  const GascoigneMesh* GM = dynamic_cast<const GascoigneMesh*>(M);
  assert(GM);
  edges = GM->GetHangingIndexHandler().GetStructure();
}

/*-----------------------------------------*/

const fixarray<3,int>& HNStructureQ12d::regular_nodes(int i) const 
{
  map<int,EdgeVector>::const_iterator p = edges->find(i);

  assert(p!=edges->end());
//   if (p!=edges->end())
//     {
      return p->second; 
//     }
}

/*----------------------------------------------*/

void HNStructureQ12d::CondenseHanging(IntVector& indices) const
{
  assert(indices.size()==lnoe.size());
  for(int ii=0; ii<indices.size(); ii++)
    {
      fixarray<3,int> p = lnoe[ii];

      int& hang = indices[p[1]];

      if(!hanging(hang)) continue;

      const fixarray<3,int>& f = regular_nodes(hang);

      if      ( (indices[p[0]]==f[0]) || (indices[p[2]]==f[0]) )  hang = f[1];
      else if ( (indices[p[0]]==f[1]) || (indices[p[2]]==f[1]) )  hang = f[0];
      else  assert(0);
    }
}

/*----------------------------------------------*/

void HNStructureQ12d::CondenseHanging(EntryMatrix& E, IntVector& indices) const
{
  assert(indices.size()==lnoe.size());
  for(int ii=0; ii<indices.size(); ii++)
    {
      fixarray<3,int> p = lnoe[ii];

      int& hang = indices[p[1]];

      if(!hanging(hang)) continue;

      const fixarray<3,int>& f = regular_nodes(hang);

      if ((indices[p[2]]==f[0]) || (indices[p[2]]==f[1]) ) swap(p[0],p[2]);

      if      ( indices[p[0]]==f[0] )  hang = f[1];
      else if ( indices[p[0]]==f[1] )  hang = f[0];
      else assert(0);

      E.add_column     (p[0],p[1],weight(0));
      E.multiply_column(p[1],     weight(1));
      
      E.add_row        (p[0],p[1],weight(0));
      E.multiply_row   (p[1],     weight(1));
    }
}

/*----------------------------------------------*/

void HNStructureQ12d::CondenseHangingPatch(EntryMatrix& E, IntVector& indices) const
{
  assert(indices.size()==9);

  for(int ii=0; ii<4; ii++)
    {
      fixarray<3,int> p = lnop[ii];

      int i = indices[p[2]];
      if(!hanging(i)) continue;

      E.add_column     (p[0],p[2],weight(0));
      E.add_column     (p[1],p[2],weight(1));
      E.add_row        (p[0],p[2],weight(0));
      E.add_row        (p[1],p[2],weight(1));
      E.multiply_column(p[2],     0.);
      E.multiply_row   (p[2],     0.);
    }
}

/*----------------------------------------------*/

int HNStructureQ12d::hanging(int i) const 
{
  assert(edges);
  if (edges->find(i)!=edges->end()) return 1;
  return 0;
}

/*-----------------------------------------*/

void HNStructureQ12d::Zero(GlobalVector& u) const
{
  for(const_iterator p=edges->begin();p!=edges->end();p++)
    {
      u.zero_node(p->first);
    }
}

/*-----------------------------------------*/

bool HNStructureQ12d::ZeroCheck(const GlobalVector& u) const
{  
  for(const_iterator p=edges->begin();p!=edges->end();p++)
    {
      int i = p->first;
      for(int c=0; c<u.ncomp(); c++)
	{
	  if(u(i,c)!=0.) return 1;
	}
    }
  return 0;
}

/*-----------------------------------------*/

void HNStructureQ12d::Average(GlobalVector& u) const
{
  for(const_iterator p=edges->begin();p!=edges->end();p++)
    {
      const fixarray<3,int>& f = p->second;
      u.equ_node(p->first, wei[0], f[0], wei[1], f[1]);
    }
}

/*-----------------------------------------*/

void HNStructureQ12d::Distribute(GlobalVector& u) const
{
  for(const_iterator p=edges->begin();p!=edges->end();p++)
    {
      int i = p->first;

      u.add_node(p->second[0],wei[0],i);
      u.add_node(p->second[1],wei[1],i);

      u.zero_node(i);
    }
}

/*-----------------------------------------*/

void HNStructureQ12d::MatrixDiag(int ncomp, MatrixInterface& A) const
{
  nmatrix<double> M(ncomp);
  M.identity();
  for(const_iterator p=edges->begin();p!=edges->end();p++)
    {
      A.entry_diag(p->first,M);
    }
}
}
