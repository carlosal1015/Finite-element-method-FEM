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


#include "hnstructureq43d.h"
#include "gascoignemesh.h"

using namespace std;

namespace Gascoigne
{
/**********************************************************/

HNStructureQ43d::HNStructureQ43d() : HNStructureQ23d()
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

  Medge.resize(5,2);
  Mq2edge.resize(5,2);
  Medge(0,0) = w[0][0];   Medge(0,1) = w[1][0];   Mq2edge(0,0) = wq2[0][0];   Mq2edge(0,1) = wq2[1][0];
  Medge(1,0) = w[0][3]-1; Medge(1,1) = w[1][3];   Mq2edge(1,0) = wq2[0][3]-1; Mq2edge(1,1) = wq2[1][3];
  Medge(2,0) = w[0][1];   Medge(2,1) = w[1][1];   Mq2edge(2,0) = wq2[0][1];   Mq2edge(2,1) = wq2[1][1];
  Medge(3,0) = w[0][4];   Medge(3,1) = w[1][4]-1; Mq2edge(3,0) = wq2[0][4];   Mq2edge(3,1) = wq2[1][4]-1;
  Medge(4,0) = w[0][2];   Medge(4,1) = w[1][2];   Mq2edge(4,0) = wq2[0][2];   Mq2edge(4,1) = wq2[1][2];

  Mface.resize(25,4);
  Mq2face.resize(25,4);
            // diese knoten bleiben erhalten
  for(int y=0; y<3; y++)
  {
    for(int x=0; x<3; x++)
    {
                // dof phi^n
      int n = y*10+2*x;
                // faces
      for(int i=0; i<4; i++)
      {
        Mface(n,i) = w[i%2][x]*w[i/2][y];
        Mq2face(n,i) = wq2[i%2][x]*wq2[i/2][y];
      }
    }
  }
            // die neuen
  for(int y=0; y<3; y++)
  {
    int n1 = 10*y+1;
    int n2 = 10*y+3;
    for(int i=0; i<4; i++)
    {
      Mface(n1,i) = w[i%2][3]*w[i/2][y];
      Mface(n2,i) = w[i%2][4]*w[i/2][y];
      Mq2face(n1,i) = wq2[i%2][3]*wq2[i/2][y];
      Mq2face(n2,i) = wq2[i%2][4]*wq2[i/2][y];
    }
  }
            // die neuen
  for(int x=0; x<3; x++)
  {
    int n1 = 2*x+5;
    int n2 = 2*x+15;
    for(int i=0; i<4; i++)
    {
      Mface(n1,i) = w[i%2][x]*w[i/2][3];
      Mface(n2,i) = w[i%2][x]*w[i/2][4];
      Mq2face(n1,i) = wq2[i%2][x]*wq2[i/2][3];
      Mq2face(n2,i) = wq2[i%2][x]*wq2[i/2][4];
    }
  }
  for(int x=0; x<2; x++)
  {
    for(int y=0; y<2; y++)
    {
      int n = 10*y+2*x+6;
      for(int i=0; i<4; i++)
      {
        Mface(n,i) = w[i%2][3+x]*w[i/2][3+y];
        Mq2face(n,i) = wq2[i%2][3+x]*wq2[i/2][3+y];
      }
    }
  }
  Mface(6,0)  -= 1; Mq2face(6,0)  -= 1;
  Mface(8,1)  -= 1; Mq2face(8,1)  -= 1;
  Mface(16,2) -= 1; Mq2face(16,2) -= 1;
  Mface(18,3) -= 1; Mq2face(18,3) -= 1;
}

/**********************************************************/

void HNStructureQ43d::add_column(EntryMatrix& A, const EntryMatrix&B, int j1, int j2, double s) const
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

void HNStructureQ43d::add_row(EntryMatrix& A, const EntryMatrix&B, int i1, int i2, double s) const
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

void HNStructureQ43d::GetHangingIndices(vector<int>& hang_e, vector<int>& hang_f, const IntVector& indices) const
{
  for(int i=0; i<125; i++)
  {
    int x = i%5;
    int y = (i/5)%5;
    int z = i/25;
    if((x%2==0)&&(y%2==0)&&(z%2==0))
    {
      continue;
    }
    int ht = hanging(indices[i]);
    if(!ht)
    {
      continue;
    }
    if(ht==2)
    {
      hang_e.push_back(i);
    }
    else if(ht==4)
    {
      hang_f.push_back(i);
    }
    else
    {
      abort();
    }
  }
}

/**********************************************************/

int HNStructureQ43d::hanging(int i) const
{
  if(q4edges->find(i)!=q4edges->end())
  {
    return 2;
  }
  if(q4faces->find(i)!=q4faces->end())
  {
    return 4;
  }
  return 0;
}

/**********************************************************/

fixarray<5,int> HNStructureQ43d::local_nodes_on_edge(int e, const IntVector& indices) const
{
  int node = indices[e];
  fixarray<5,int> R;
  int x = e%5;
  int y = (e/5)%5;
  int z = e/25;

  assert(x%2+y%2+z%2==1);

  const fixarray<6,int>& rege=regular_nodes_on_edge(node);
  int type = rege[5];
            // die Reihenfolge der local_nodes muss mit den regular
            // nodes uebereinstimmen.
  if(x%2==1)
  {
    for(int i=0;i<5;i++)
    {
      R[i]=25*z+5*y+i;
    }
  }
  else if(y%2==1)
  {
    for(int i=0;i<5;i++)
    {
      R[i]=25*z+5*i+x;
    }
  }
  else if(z%2==1)
  {
    for(int i=0;i<5;i++)
    {
      R[i]=25*i+5*y+x;
    }
  }
  else
  {
    abort();
  }
  if(((type==0)&&(R[3]==e))||((type==1)&&(R[1]==e)))
  {
    swap(R[0],R[4]);
    swap(R[1],R[3]);
  }
  assert(e==R[2*type+1]);
  assert(indices[R[0]]==rege[0]);
  assert(indices[R[4]]==rege[2]);

  return R;
}

/**********************************************************/

fixarray<25,int> HNStructureQ43d::local_nodes_on_face(int e, const IntVector& indices) const
{
  int node = indices[e];
  fixarray<25,int> R;
  int x = e%5;
  int y = (e/5)%5;
  int z = e/25;

  assert(x%2+y%2+z%2==2);

  const fixarray<26,int>& regf = regular_nodes_on_face(node);

            // die Reihenfolge der local_nodes muss mit den regular
            // nodes uebereinstimmen.
  if(z%4==0)
  {
    for(int i=0; i<25; i++)
    {
      R[i] = i+25*z;
    }
  }
  else if(x%4==0)
  {
    for(int i=0; i<25; i++)
    {
      R[i] = x + (25*(i%5)+5*(i/5));
    }
  }
  else if(y%4==0)
  {
    for(int i=0; i<25; i++)
    {
      R[i] = 5*y + (i%5+25*(i/5));
    }
  }
  assert((R[6]==e)||(R[8]==e)||(R[16]==e)||(R[18]==e));

  while(regf[0]!=indices[R[0]])
  {
    fixarray<25,int> R1=R;
    for(int y=0; y<5; y++)
    {
      for(int x=0; x<5; x++)
      {
        R[5*y+x] = R1[5*x+4-y];
      }
    }
  }

  assert(regf[0]==indices[R[0]]);
  if(regf[1]==indices[R[10]])
  {
    for(int y=1; y<5; y++)
    {
      for(int x=0; x<y; x++)
      {
        swap(R[5*y+x],R[5*x+y]);
      }
    }
  }
  assert(regf[1]==indices[R[2]]);
  assert(regf[5]==indices[R[10]]);

  return R;
}

/**********************************************************/

void HNStructureQ43d::modify_column_higher(EntryMatrix& E, const vector<int>& hang_e, const vector<int>& hang_f, const IntVector& indices) const
{
  EntryMatrix M=E;

            // kanten
  for(int i=0; i<hang_e.size(); i++)
  {
              // der lokale index des haengenden knoten
    int h    = hang_e[i];
    int node = indices[h];
    const fixarray<6,int>& rege = regular_nodes_on_edge(node);
    int type = rege[5];
    assert((type>=0)&&(type<2));
              // die 5 lokalen indices der regulars
    fixarray<5,int> lnoe = local_nodes_on_edge(h,indices);

    for(int j=0; j<5; j++)
    {
      add_column(E,M,lnoe[j],h,Medge(j,type));
    }
  }
            // flaechen
  for(int i=0; i<hang_f.size(); i++)
  {
              // der lokale index des haengenden knoten
    int h    = hang_f[i];
    int node = indices[h];
    const fixarray<26,int>& regf = regular_nodes_on_face(node);
    int type = regf[25];
    assert((type>=0)&&(type<4));
              // die 25 lokalen indices der regulars
    fixarray<25,int> lnof = local_nodes_on_face(h,indices);

    for(int j=0; j<25; j++)
    {
      add_column(E,M,lnof[j],h,Mface(j,type));
    }
  }
}

/**********************************************************/

void HNStructureQ43d::modify_column_lower(EntryMatrix& E, const vector<int>& hang_e, const vector<int>& hang_f, const IntVector& indices) const
{
  EntryMatrix M=E;

            // kanten
  for(int i=0; i<hang_e.size(); i++)
  {
              // der lokale index des haengenden knoten
    int h    = hang_e[i];
    int node = indices[h];
    const fixarray<6,int>& rege = regular_nodes_on_edge(node);
    int type = rege[5];
    assert((type>=0)&&(type<2));
              // die 5 lokalen indices der regulars
    fixarray<5,int> lnoe = local_nodes_on_edge(h,indices);

    for(int j=0; j<5; j++)
    {
      add_column(E,M,lnoe[j],h,Mq2edge(j,type));
    }
  }
            // flaechen
  for(int i=0; i<hang_f.size(); i++)
  {
              // der lokale index des haengenden knoten
    int h    = hang_f[i];
    int node = indices[h];
    const fixarray<26,int>& regf = regular_nodes_on_face(node);
    int type = regf[25];
    assert((type>=0)&&(type<4));
              // die 25 lokalen indices der regulars
    fixarray<25,int> lnof = local_nodes_on_face(h,indices);

    for(int j=0; j<25; j++)
    {
      add_column(E,M,lnof[j],h,Mq2face(j,type));
    }
  }
}

/**********************************************************/

void HNStructureQ43d::modify_row_higher(EntryMatrix& E, const vector<int>& hang_e, const vector<int>& hang_f, const IntVector& indices) const
{
  EntryMatrix M=E;

            // kanten
  for(int i=0; i<hang_e.size(); i++)
  {
              // der lokale index des haengenden knoten
    int h    = hang_e[i];
    int node = indices[h];
    const fixarray<6,int>& rege = regular_nodes_on_edge(node);
    int type = rege[5];
    assert((type>=0)&&(type<2));
              // die 5 lokalen indices der regulars
    fixarray<5,int> lnoe = local_nodes_on_edge(h,indices);

    for(int j=0; j<5; j++)
    {
      add_row(E,M,lnoe[j],h,Medge(j,type));
    }
  }
            // flaechen
  for(int i=0; i<hang_f.size(); i++)
  {
              // der lokale index des haengenden knoten
    int h    = hang_f[i];
    int node = indices[h];
    const fixarray<26,int>& regf = regular_nodes_on_face(node);
    int type = regf[25];
    assert((type>=0)&&(type<4));
              // die 25 lokalen indices der regulars
    fixarray<25,int> lnof = local_nodes_on_face(h,indices);

    for(int j=0; j<25; j++)
    {
      add_row(E,M,lnof[j],h,Mface(j,type));
    }
  }
}

/**********************************************************/

void HNStructureQ43d::modify_row_lower(EntryMatrix& E, const vector<int>& hang_e, const vector<int>& hang_f, const IntVector& indices) const
{
  EntryMatrix M=E;

            // kanten
  for(int i=0; i<hang_e.size(); i++)
  {
              // der lokale index des haengenden knoten
    int h    = hang_e[i];
    int node = indices[h];
    const fixarray<6,int>& rege = regular_nodes_on_edge(node);
    int type = rege[5];
    assert((type>=0)&&(type<2));
              // die 5 lokalen indices der regulars
    fixarray<5,int> lnoe = local_nodes_on_edge(h,indices);

    for(int j=0; j<5; j++)
    {
      add_row(E,M,lnoe[j],h,Mq2edge(j,type));
    }
  }
            // flaechen
  for(int i=0; i<hang_f.size(); i++)
  {
              // der lokale index des haengenden knoten
    int h    = hang_f[i];
    int node = indices[h];
    const fixarray<26,int>& regf = regular_nodes_on_face(node);
    int type = regf[25];
    assert((type>=0)&&(type<4));
              // die 25 lokalen indices der regulars
    fixarray<25,int> lnof = local_nodes_on_face(h,indices);

    for(int j=0; j<25; j++)
    {
      add_row(E,M,lnof[j],h,Mq2face(j,type));
    }
  }
}

/**********************************************************/

const fixarray<6,int>& HNStructureQ43d::regular_nodes_on_edge(int i) const
{
  map<int,EdgeVector>::const_iterator p = q4edges->find(i);
  if(p!=q4edges->end())
  {
    return p->second;
  }
  abort();
}

/**********************************************************/

const fixarray<26,int>& HNStructureQ43d::regular_nodes_on_face(int i) const
{
  map<int,FaceVector>::const_iterator p = q4faces->find(i);
  if(p!=q4faces->end())
  {
    return p->second;
  }
  abort();
}

/**********************************************************/

void HNStructureQ43d::ReInit(const MeshInterface* m)
{
  HNStructureQ23d::ReInit(m);
  const GascoigneMesh* GM = dynamic_cast<const GascoigneMesh*>(m);
  assert(GM);
  q4edges = GM->GetHangingIndexHandler().GetQ4Structure();
  q4faces = GM->GetHangingIndexHandler().GetQ4StructureFace();
}

/**********************************************************/

void HNStructureQ43d::CondenseHanging(IntVector& indices) const
{
  vector<int> hang_e,hang_f;
  GetHangingIndices(hang_e,hang_f,indices);

  for(int i=0; i<hang_e.size(); i++)
  {
    int h = hang_e[i];
    int n = indices[h];
    const fixarray<6,int>& regn = regular_nodes_on_edge(n);
    int type = regn[5];
    indices[h]=regn[3+type];
  }

  for(int i=0; i<hang_f.size(); i++)
  {
    int h = hang_f[i];
    int n = indices[h];
    const fixarray<26,int>& regn = regular_nodes_on_face(n);
    int type = regn[25];
    if(type==0)
    {
      indices[h]=regn[18];
    }
    else if(type==1)
    {
      indices[h]=regn[19];
    }
    else if(type==2)
    {
      indices[h]=regn[23];
    }
    else if(type==3)
    {
      indices[h]=regn[24];
    }
  }
}

/**********************************************************/

void HNStructureQ43d::CondenseHanging(EntryMatrix& E, IntVector& indices) const
{
  vector<int> hang_e(0),hang_f(0);
  GetHangingIndices(hang_e,hang_f,indices);

  if(hang_e.size()+hang_f.size())
  {
    modify_row_higher(E,hang_e,hang_f,indices);
    modify_column_higher(E,hang_e,hang_f,indices);
    CondenseHanging(indices);
  }
}

/**********************************************************/

void HNStructureQ43d::CondenseHangingLowerHigher(EntryMatrix& E, IntVector& indices) const
{
  vector<int> hang_e(0),hang_f(0);
  GetHangingIndices(hang_e,hang_f,indices);

  if(hang_e.size()+hang_f.size())
  {
    modify_row_higher(E,hang_e,hang_f,indices);
    modify_column_lower(E,hang_e,hang_f,indices);
    CondenseHanging(indices);
  }
}

/**********************************************************/

void HNStructureQ43d::CondenseHangingHigherLower(EntryMatrix& E, IntVector& indices) const
{
  vector<int> hang_e(0),hang_f(0);
  GetHangingIndices(hang_e,hang_f,indices);

  if(hang_e.size()+hang_f.size())
  {
    modify_row_lower(E,hang_e,hang_f,indices);
    modify_column_higher(E,hang_e,hang_f,indices);
    CondenseHanging(indices);
  }
}

/**********************************************************/

void HNStructureQ43d::MatrixDiag(int ncomp, MatrixInterface& A) const
{
  nmatrix<double> M(ncomp);
  M.identity();
  for(const_iteratorq4 p=q4edges->begin();p!=q4edges->end();p++)
  {
    A.entry_diag(p->first,M);
  }
  for(face_const_iteratorq4 p=q4faces->begin();p!=q4faces->end();p++)
  {
    A.entry_diag(p->first,M);
  }
}

/**********************************************************/

void HNStructureQ43d::SparseStructureDiag(SparseStructure& S) const
{
  for(const_iteratorq4 p=q4edges->begin();p!=q4edges->end();p++)
  {
    int i = p->first;
    S.build_add(i,i);
  }
  for(face_const_iteratorq4 p=q4faces->begin();p!=q4faces->end();p++)
  {
    int i = p->first;
    S.build_add(i,i);
  }
}

/**********************************************************/

void HNStructureQ43d::Zero(GlobalVector& u) const
{
  for(const_iteratorq4 p=q4edges->begin();p!=q4edges->end();p++)
  {
    u.zero_node(p->first);
  }
  for(face_const_iteratorq4 p=q4faces->begin();p!=q4faces->end();p++)
  {
    u.zero_node(p->first);
  }
}

/**********************************************************/

void HNStructureQ43d::Average(GlobalVector& u) const
{
  for(const_iteratorq4 p=q4edges->begin();p!=q4edges->end();p++)
  {
    const fixarray<6,int>& f = p->second;
    int t = f[5];
    assert((t==0)||(t==1));
    u.equ_node(p->first, w[t][0],f[0], w[t][1],f[1], w[t][2],f[2]);
    u.add_node(p->first, w[t][3],f[3], w[t][4],f[4]);
  }

  for(face_const_iteratorq4 p=q4faces->begin();p!=q4faces->end();p++)
  {
    const fixarray<26,int>& f = p->second;
    int t = f[25];
    assert((t>=0)&&(t<4));
    int tx = t%2;
    int ty = t/2;
    u.zero_node(p->first);
    for (int i=0;i<25;++i)
    {
      u.add_node(p->first, w[tx][i%5]*w[ty][i/5], f[i]);
    }
  }
}

/**********************************************************/

void HNStructureQ43d::Distribute(GlobalVector& u) const
{
  for(const_iteratorq4 p=q4edges->begin();p!=q4edges->end();p++)
  {
    int i = p->first;
    int t = p->second[5];
    for (int j=0; j<5; j++)
    {
      u.add_node(p->second[j],w[t][j],i);
    }
    u.zero_node(i);
  }
  for(face_const_iteratorq4 p=q4faces->begin();p!=q4faces->end();p++)
  {
    const fixarray<26,int>& f = p->second;
    int t = f[25];
    assert((t>=0)&&(t<4));
    int tx = t%2;
    int ty = t/2;
    for (int i=0;i<25;++i)
    {
      u.add_node(f[i],w[tx][i%5]*w[ty][i/5],p->first);
    }
    u.zero_node(p->first);
  }
}

/**********************************************************/

bool HNStructureQ43d::ZeroCheck(const GlobalVector& u) const
{
  for(const_iteratorq4 p=q4edges->begin();p!=q4edges->end();p++)
  {
    int i = p->first;
    for(int c=0; c<u.ncomp(); c++)
    {
      if(u(i,c)!=0.) return false;
    }
  }
  for(face_const_iteratorq4 p=q4faces->begin();p!=q4faces->end();p++)
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

