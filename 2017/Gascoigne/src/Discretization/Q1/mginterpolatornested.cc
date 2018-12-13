/**
*
* Copyright (C) 2004 by the Gascoigne 3D authors
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


#include  "mginterpolatornested.h"
#include  "gascoignemeshtransfer.h"


using namespace std;

/*--------------------------------------------------------*/

namespace Gascoigne
{
void MgInterpolatorNested::BasicInit(const MeshTransferInterface* MT)
{
  const GascoigneMeshTransfer* GT = dynamic_cast<const GascoigneMeshTransfer*>(MT);

  c2f.reservesize(GT->GetC2f().size());

  zweier = GT->GetZweier();
  vierer = GT->GetVierer();
  achter = GT->GetAchter();
  c2f    = GT->GetC2f();
}

/*-----------------------------------------*/
  
void MgInterpolatorNested::restrict_zero(GlobalVector& uL, const GlobalVector& ul) const
{
  for(int i=0;i<c2f.size();i++)  uL.equ_node(i,1.,c2f[i],ul);
  for(map<int,fixarray<2,int> >::const_iterator p=zweier.begin();
      p!=zweier.end();p++) {
    int il = p->first;
    fixarray<2,int> n2 = p->second;
    uL.add_node(n2[0],0.5,il,ul);
    uL.add_node(n2[1],0.5,il,ul);
  }
  for(map<int,fixarray<4,int> >::const_iterator p=vierer.begin();
      p!=vierer.end();p++) {
    int il = p->first;
    fixarray<4,int> n4 = p->second;
    uL.add_node(n4[0],0.25,il,ul);
    uL.add_node(n4[1],0.25,il,ul);
    uL.add_node(n4[2],0.25,il,ul);
    uL.add_node(n4[3],0.25,il,ul);
  }
  for(map<int,fixarray<8,int> >::const_iterator p=achter.begin();
      p!=achter.end();p++) {
    int il = p->first;
    fixarray<8,int> n8 = p->second;
    for (int i=0; i<8; i++)
      {
	uL.add_node(n8[i],0.125,il,ul);
      }
  }
}

/*-----------------------------------------*/

void MgInterpolatorNested::prolongate_add(GlobalVector& ul, const GlobalVector& uL) const
{
  for(int i=0;i<c2f.size();i++)  ul.add_node(c2f[i],1.,i,uL);
  for(map<int,fixarray<2,int> >::const_iterator p=zweier.begin();
      p!=zweier.end();p++) 
    {
      int il = p->first;
      fixarray<2,int> n2 = p->second;
      ul.add_node(il,0.5,n2[0],uL);
      ul.add_node(il,0.5,n2[1],uL);
    }
  for(map<int,fixarray<4,int> >::const_iterator p=vierer.begin();
      p!=vierer.end();p++) 
    {
      int il = p->first;
      fixarray<4,int> n4 = p->second;
      ul.add_node(il,0.25,n4[0],uL);
      ul.add_node(il,0.25,n4[1],uL);
      ul.add_node(il,0.25,n4[2],uL);
      ul.add_node(il,0.25,n4[3],uL);
    }
  for(map<int,fixarray<8,int> >::const_iterator p=achter.begin();
      p!=achter.end();p++) 
    {
      int il = p->first;
      fixarray<8,int> n8 = p->second;
      for (int i=0; i<8; i++)
	{
	  ul.add_node(il,0.125,n8[i],uL);
	}
    }
}

/*-----------------------------------------*/

void MgInterpolatorNested::SolutionTransfer(GlobalVector& uL, const GlobalVector& ul) const
{
  for(int i=0;i<c2f.size();i++) uL.equ_node(i,1.,c2f[i],ul);
}

/*-----------------------------------------*/

void MgInterpolatorNested::Pi(GlobalVector& u) const
{
  for(map<int,fixarray<2,int> >::const_iterator p=zweier.begin();
      p!=zweier.end();p++) {
    int il = p->first;
    fixarray<2,int> n2 = p->second;
    u.add_node(il,-0.5,c2f[n2[0]],u);
    u.add_node(il,-0.5,c2f[n2[1]],u);
  }
  for(map<int,fixarray<4,int> >::const_iterator p=vierer.begin();
      p!=vierer.end();p++) {
    int il = p->first;
    fixarray<4,int> n4 = p->second;
    u.add_node(il,-0.25,c2f[n4[0]],u);
    u.add_node(il,-0.25,c2f[n4[1]],u);
    u.add_node(il,-0.25,c2f[n4[2]],u);
    u.add_node(il,-0.25,c2f[n4[3]],u);
  }
  for(int i=0;i<c2f.size();i++)  u.node_zero(c2f[i]);
}
}
