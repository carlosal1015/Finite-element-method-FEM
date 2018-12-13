/**
*
* Copyright (C) 2004, 2006 by the Gascoigne 3D authors
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


#include  "pi.h"

using namespace std;
namespace Gascoigne
{

/*-----------------------------------------*/

Pi::Pi() {}

/*-----------------------------------------*/

void Pi::vmult(CompVector<double>& y, const CompVector<double>& x, 
	       double s) const
{
  y.zero();
  int ncomp = x.ncomp();
  assert(ncomp==y.ncomp());
  {
    map<int,fixarray<2,int> >::const_iterator p;
    for(p=edge.begin();p!=edge.end();p++)
      {
	int i = p->first;
	const fixarray<2,int>& f = p->second;
	for(int c=0;c<ncomp;c++) y(i,c) = x(i,c) - 0.5*(x(f[0],c)+x(f[1],c));
      }
  }
  {
    map<int,fixarray<4,int> >::const_iterator p;
    for(p=face.begin();p!=face.end();p++)
      {
	int i = p->first;
	const fixarray<4,int>& f = p->second;
	for(int c=0;c<ncomp;c++) y(i,c) = x(i,c) - 0.25*(x(f[0],c)+x(f[1],c)+x(f[2],c)+x(f[3],c));
      }
  }
  {
    map<int,fixarray<8,int> >::const_iterator p;
    for(p=cell.begin();p!=cell.end();p++)
      {
	int i = p->first;
	const fixarray<8,int>& f = p->second;
	for(int c=0;c<ncomp;c++) y(i,c) = x(i,c) - 0.125*(x(f[0],c)+x(f[1],c)+x(f[2],c)+x(f[3],c)+x(f[4],c)+x(f[5],c)+x(f[6],c)+x(f[7],c));
      }
  }
}

/*-----------------------------------------*/

void Pi::Init(const MeshInterface* MP)
{
  edge.clear();
  face.clear();
  cell.clear();
  const GascoigneMesh3d* NMP = dynamic_cast<const GascoigneMesh3d*>(MP);
  if(NMP) 
    {
      Init3d(NMP);
    }
  else 
    {
      const GascoigneMesh2d* NMP2 = dynamic_cast<const GascoigneMesh2d*>(MP);
      assert(NMP2);
      Init2d(NMP2);
    }
}

/*-----------------------------------------*/

void Pi::Init2d(const GascoigneMesh2d* MP)
{
  assert(MP->HasPatch());
 
  for(int i=0;i<MP->npatches();i++)
    {
      const nvector<int>& ind = *MP->IndicesOfPatch(i);
      {
	fixarray<4,int> f;
	f[0] = ind[0];
	f[1] = ind[2];
	f[2] = ind[6];
	f[3] = ind[8];
	face.insert(make_pair(ind[4],f));
      }
      {
	fixarray<2,int> f;
	f[0] = ind[0];
	f[1] = ind[2];
	edge.insert(make_pair(ind[1],f));
	f[0] = ind[0];
	f[1] = ind[6];
	edge.insert(make_pair(ind[3],f));
	f[0] = ind[2];
	f[1] = ind[8];
	edge.insert(make_pair(ind[5],f));
	f[0] = ind[6];
	f[1] = ind[8];
	edge.insert(make_pair(ind[7],f));
      }
    }
}

/*-----------------------------------------*/

void Pi::Init3d(const GascoigneMesh3d* MP)
{
  assert(MP->HasPatch());
 
  for(int i=0;i<MP->npatches();i++)
    {
      const nvector<int>& ind = *MP->IndicesOfPatch(i);

      {
	fixarray<8,int> f;
	f[0] = ind[0];
	f[1] = ind[2];
	f[2] = ind[6];
	f[3] = ind[8];
	f[4] = ind[18];
	f[5] = ind[20];
	f[6] = ind[24];
	f[7] = ind[26];
	cell.insert(make_pair(ind[13],f));
      }
      {
	fixarray<4,int> f;
	f[0] = ind[0];
	f[1] = ind[2];
	f[2] = ind[6];
	f[3] = ind[8];
	face.insert(make_pair(ind[4],f));
	f[0] = ind[0];
	f[1] = ind[2];
	f[2] = ind[18];
	f[3] = ind[20];
	face.insert(make_pair(ind[10],f));
	f[0] = ind[6];
	f[1] = ind[8];
	f[2] = ind[24];
	f[3] = ind[26];
	face.insert(make_pair(ind[16],f));
	f[0] = ind[18];
	f[1] = ind[20];
	f[2] = ind[24];
	f[3] = ind[26];
	face.insert(make_pair(ind[22],f));

	f[0] = ind[2];
	f[1] = ind[20];
	f[2] = ind[8];
	f[3] = ind[26];
	face.insert(make_pair(ind[14],f));
	f[0] = ind[0];
	f[1] = ind[18];
	f[2] = ind[6];
	f[3] = ind[24];
	face.insert(make_pair(ind[12],f));
	
	

      }
      {
	fixarray<2,int> f;
	f[0] = ind[0];
	f[1] = ind[2];
	edge.insert(make_pair(ind[1],f));
	f[0] = ind[0];
	f[1] = ind[6];
	edge.insert(make_pair(ind[3],f));
	f[0] = ind[8];
	f[1] = ind[2];
	edge.insert(make_pair(ind[5],f));
	f[0] = ind[6];
	f[1] = ind[8];
	edge.insert(make_pair(ind[7],f));
	f[0] = ind[0];
	f[1] = ind[18];
	edge.insert(make_pair(ind[9],f));
	f[0] = ind[2];
	f[1] = ind[20];
	edge.insert(make_pair(ind[11],f));
	f[0] = ind[6];
	f[1] = ind[24];
	edge.insert(make_pair(ind[15],f));
	f[0] = ind[26];
	f[1] = ind[8];
	edge.insert(make_pair(ind[17],f));
	f[0] = ind[18];
	f[1] = ind[20];
	edge.insert(make_pair(ind[19],f));
	f[0] = ind[24];
	f[1] = ind[18];
	edge.insert(make_pair(ind[21],f));
	f[0] = ind[20];
	f[1] = ind[26];
	edge.insert(make_pair(ind[23],f));
	f[0] = ind[24];
	f[1] = ind[26];
	edge.insert(make_pair(ind[25],f));
      }
    }
}
}
