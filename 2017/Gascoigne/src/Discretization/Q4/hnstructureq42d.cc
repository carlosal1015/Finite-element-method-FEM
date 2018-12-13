/**
*
* Copyright (C) 2006 by the Gascoigne 3D authors
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


#include "hnstructureq42d.h"
#include "gascoignemesh.h"

using namespace std;

namespace Gascoigne
{
/**********************************************************/

HNStructureQ42d::HNStructureQ42d() : HNStructureQ22d()
{
  w.resize(2);
  wq2.resize(2);

  w[0][0] =  0.2734375;  wq2[0][0] =  0.375;
  w[0][1] =  1.0937500;  wq2[0][1] =  0.75;
  w[0][2] = -0.5468750;  wq2[0][2] = -0.125,
  w[0][3] =  0.2187500;  wq2[0][3] =  0.;
  w[0][4] = -0.0390625;  wq2[0][4] =  0.;

  w[1][0] = -0.0390625;  wq2[1][0] = -0.125;
  w[1][1] =  0.4687500;  wq2[1][1] =  0.75;
  w[1][2] =  0.7031250;  wq2[1][2] =  0.375;
  w[1][3] = -0.1562500;  wq2[1][3] =  0.;
  w[1][4] =  0.0234375;  wq2[1][4] =  0.;

  M.resize(5,2);
  M(0,0) = w[0][0];   M(0,1) = w[1][0];
  M(1,0) = w[0][1];   M(1,1) = w[1][1];
  M(2,0) = w[0][2];   M(2,1) = w[1][2];
  M(3,0) = w[0][3]-1; M(3,1) = w[1][3];
  M(4,0) = w[0][4];   M(4,1) = w[1][4]-1;

  Mq2.resize(5,2);
  Mq2(0,0) = wq2[0][0];   Mq2(0,1) = wq2[1][0];
  Mq2(1,0) = wq2[0][1];   Mq2(1,1) = wq2[1][1];
  Mq2(2,0) = wq2[0][2];   Mq2(2,1) = wq2[1][2];
  Mq2(3,0) = wq2[0][3]-1; Mq2(3,1) = wq2[1][3];
  Mq2(4,0) = wq2[0][4];   Mq2(4,1) = wq2[1][4]-1;
}

/**********************************************************/

void HNStructureQ42d::add_column(EntryMatrix& A, const EntryMatrix&B, int j1, int j2, double s) const
{
  assert(A.Ndof()==B.Ndof());
  assert(A.Mdof()==B.Mdof());

  for(int i=0;i<A.Ndof();i++)
  {
    DoubleVector::iterator p1 = A.begin(i,j1);
    DoubleVector::const_iterator p2 = B.begin(i,j2);
    DoubleVector::const_iterator q1 = A.end(i,j1);
    while(p1!=q1) *p1++ += s * *p2++;
  }
}

/**********************************************************/

void HNStructureQ42d::add_row(EntryMatrix& A, const EntryMatrix&B, int i1, int i2, double s) const
{
  assert(A.Ndof()==B.Ndof());
  assert(A.Mdof()==B.Mdof());

  for(int j=0;j<A.Mdof();j++)
  {
    DoubleVector::iterator p1 = A.begin(i1,j);
    DoubleVector::const_iterator p2 = B.begin(i2,j);
    DoubleVector::const_iterator q1 = A.end(i1,j);
    while(p1!=q1) *p1++ += s * *p2++;
  }
}

/**********************************************************/

void HNStructureQ42d::GetHangingIndices(vector<int>& hang, const IntVector& indices) const
{
  for(int i=0; i<25; i++)
  {
    int x = i%5;
    int y = i/5;
    if((x%2==0)&&(y%2==0))
    {
      continue;
    }
    if(!hanging(indices[i]))
    {
      continue;
    }
    hang.push_back(i);
  }
}

/**********************************************************/

fixarray<5,int> HNStructureQ42d::local_nodes(int e,int n) const
{
  fixarray<5,int> R;

  int x = e%5;
  int y = e/5;
  int type = regular_nodes(n)[5];
  assert((type==0)||(type==1));

  assert(x%2+y%2==1);

  if(x%2==1)
  {
    for(int i=0; i<5; i++)
    {
      R[i]=5*y+(2*i)%5;
    }
  }
  else if(y%2==1)
  {
    for(int i=0; i<5; i++)
    {
      R[i]=(10*i)%25+x;
    }
  }
  else
  {
    abort();
  }
  if(((type==0)&&(R[4]==e))||((type==1)&&(R[3]==e)))
  {
    swap(R[0],R[2]);
    swap(R[3],R[4]);
  }
  assert(e==R[3+type]);

  return R;
}

/**********************************************************/

void HNStructureQ42d::modify_column_higher(EntryMatrix& E, const vector<int>& hang, const IntVector& indices) const
{
  EntryMatrix T=E;
  for(int i=0; i<hang.size(); i++)
  {
    int h = hang[i];
    int n = indices[h];
    int t = regular_nodes(n)[5];
    assert((t==0)||(t==1));

    fixarray<5,int> lnoe = local_nodes(h,n);

    for(int j=0; j<5; j++)
    {
      add_column(E,T,lnoe[j],h,M(j,t));
    }
  }
}

/**********************************************************/

void HNStructureQ42d::modify_column_lower(EntryMatrix& E, const vector<int>& hang, const IntVector& indices) const
{
  EntryMatrix T=E;
  for(int i=0; i<hang.size(); i++)
  {
    int h = hang[i];
    int n = indices[h];
    int t = regular_nodes(n)[5];
    assert((t==0)||(t==1));

    fixarray<5,int> lnoe = local_nodes(h,n);

    for(int j=0; j<5; j++)
    {
      add_column(E,T,lnoe[j],h,Mq2(j,t));
    }
  }
}

/**********************************************************/

void HNStructureQ42d::modify_row_higher(EntryMatrix& E, const vector<int>& hang, const IntVector& indices) const
{
  EntryMatrix T=E;
  for(int i=0; i<hang.size(); i++)
  {
    int h = hang[i];
    int n = indices[h];
    int t = regular_nodes(n)[5];
    assert((t==0)||(t==1));

    fixarray<5,int> lnoe = local_nodes(h,n);

    for(int j=0; j<5; j++)
    {
      add_row(E,T,lnoe[j],h,M(j,t));
    }
  }
}

/**********************************************************/

void HNStructureQ42d::modify_row_lower(EntryMatrix& E, const vector<int>& hang, const IntVector& indices) const
{
  EntryMatrix T=E;
  for(int i=0; i<hang.size(); i++)
  {
    int h = hang[i];
    int n = indices[h];
    int t = regular_nodes(n)[5];
    assert((t==0)||(t==1));

    fixarray<5,int> lnoe = local_nodes(h,n);

    for(int j=0; j<5; j++)
    {
      add_row(E,T,lnoe[j],h,Mq2(j,t));
    }
  }
}

/**********************************************************/

const fixarray<6,int>& HNStructureQ42d::regular_nodes(int i) const
{
  map<int,fixarray<6,int> >::const_iterator p = q4edges->find(i);
  if(p!=q4edges->end())
  {
    return p->second;
  }
  abort();
}

/**********************************************************/

void HNStructureQ42d::ReInit(const MeshInterface* m)
{
  HNStructureQ22d::ReInit(m);
  const GascoigneMesh* GM = dynamic_cast<const GascoigneMesh*>(m);
  assert(GM);
  q4edges = GM->GetHangingIndexHandler().GetQ4Structure();
}

/**********************************************************/

void HNStructureQ42d::CondenseHanging(IntVector& indices) const
{
  vector<int> hang(0);
  GetHangingIndices(hang,indices);

  for(int i=0; i<hang.size(); i++)
  {
    int h = hang[i];
    int n = indices[h];
    const fixarray<6,int>& regn = regular_nodes(n);

    int type = regn[5];
    indices[h] = regn[3+type];
  }
}

/**********************************************************/

void HNStructureQ42d::CondenseHanging(EntryMatrix& E, IntVector& indices) const
{
  vector<int> hang(0);
  GetHangingIndices(hang,indices);

  if(hang.size())
  {
    modify_row_higher(E,hang,indices);
    modify_column_higher(E,hang,indices);
    CondenseHanging(indices);
  }
}

/**********************************************************/

void HNStructureQ42d::CondenseHangingLowerHigher(EntryMatrix& E, IntVector& indices) const
{
  vector<int> hang(0);
  GetHangingIndices(hang,indices);

  if(hang.size())
  {
    modify_row_higher(E,hang,indices);
    modify_column_lower(E,hang,indices);
    CondenseHanging(indices);
  }
}

/**********************************************************/

void HNStructureQ42d::CondenseHangingHigherLower(EntryMatrix& E, IntVector& indices) const
{
  vector<int> hang(0);
  GetHangingIndices(hang,indices);

  if(hang.size())
  {
    modify_row_lower(E,hang,indices);
    modify_column_higher(E,hang,indices);
    CondenseHanging(indices);
  }
}

/**********************************************************/

void HNStructureQ42d::MatrixDiag(int ncomp, MatrixInterface& A) const
{
  nmatrix<double> M(ncomp);
  M.identity();
  for(const_iteratorq4 p=q4edges->begin();p!=q4edges->end();p++)
  {
    A.entry_diag(p->first,M);
  }
}

/**********************************************************/

void HNStructureQ42d::SparseStructureDiag(SparseStructure& S) const
{
  for(const_iteratorq4 p=q4edges->begin();p!=q4edges->end();p++)
  {
    int i = p->first;
    S.build_add(i,i);
  }
}

/**********************************************************/

void HNStructureQ42d::Zero(GlobalVector& u) const
{
  for(const_iteratorq4 p=q4edges->begin();p!=q4edges->end();p++)
  {
    u.zero_node(p->first);
  }
}

/**********************************************************/

void HNStructureQ42d::Average(GlobalVector& u) const
{
  for(const_iteratorq4 p=q4edges->begin();p!=q4edges->end();p++)
  {
    const fixarray<6,int>& f = p->second;
    int t = f[5];
    assert((t==0)||(t==1));
    u.equ_node(p->first, w[t][0], f[0], w[t][1], f[1], w[t][2], f[2]);
    u.add_node(p->first, w[t][3], f[3], w[t][4], f[4]);
  }
}

/**********************************************************/

void HNStructureQ42d::Distribute(GlobalVector& u) const
{
  for(const_iteratorq4 p=q4edges->begin();p!=q4edges->end();p++)
  {
    int i = p->first;
    int t = p->second[5];
    for(int j=0; j<5; j++)
    {
      u.add_node(p->second[j],w[t][j],i);
    }
    u.zero_node(i);
  }
}

/**********************************************************/

bool HNStructureQ42d::ZeroCheck(const GlobalVector& u) const
{
  for(const_iteratorq4 p=q4edges->begin();p!=q4edges->end();p++)
  {
    int i = p->first;
    for(int c=0; c<u.ncomp(); c++)
    {
      if(u(i,c)!=0.) return false;
    }
  }
  return true;
}

/**********************************************************/
}

