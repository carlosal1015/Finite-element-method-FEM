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


#include  "hnstructureq22d.h"

using namespace std;

namespace Gascoigne
{
/*-----------------------------------------*/

HNStructureQ22d::HNStructureQ22d() : HNStructureQ12d(), q1wei(3)
{
  wei[0] =  0.375; q1wei[0] = 0.5;
  wei[1] =  0.75;  q1wei[1] = 0.5;
  wei[2] = -0.125; q1wei[2] = 0.;

  lnoe[0][0]=0; lnoe[0][1]=2; lnoe[0][2]=1;
  lnoe[1][0]=0; lnoe[1][1]=6; lnoe[1][2]=3;
  lnoe[2][0]=2; lnoe[2][1]=8; lnoe[2][2]=5;
  lnoe[3][0]=6; lnoe[3][1]=8; lnoe[3][2]=7;
}

/*-----------------------------------------*/

void HNStructureQ22d::Average(GlobalVector& u) const
{
  for(const_iterator p=edges->begin();p!=edges->end();p++)
    {
      const fixarray<3,int>& f = p->second;
      u.equ_node(p->first, wei[0], f[0], wei[1], f[1], wei[2], f[2]);
    }
}

/*-----------------------------------------*/

void HNStructureQ22d::Distribute(GlobalVector& u) const
{
  for(const_iterator p=edges->begin();p!=edges->end();p++)
    {
      int i = p->first;
      for (int j=0; j<3; j++)
	{
	  u.add_node(p->second[j],wei[j],i);
	}
      u.zero_node(i);
    }
}

/*-----------------------------------------*/

void HNStructureQ22d::CondenseHanging(EntryMatrix& E, IntVector& indices) const
{
  for(int ii=0; ii<4; ii++) // nur 4 kandiaten koennen haengen !!
    {
      int i = indices[2*ii+1];
      if(!hanging(i)) continue;

      const fixarray<3,int>& f = regular_nodes(i);

      fixarray<3,int> p = lnoe[ii];

      if ( (indices[p[0]]==f[1]) && (indices[p[1]]==f[0]) ) 
	{ 
	  swap(p[0],p[1]);
	} 

      indices[p[2]] = f[2];

      E.add_column     (p[0],p[2],weight(0));
      E.add_column     (p[1],p[2],weight(1));
      E.multiply_column(p[2],     weight(2));
      
      E.add_row        (p[0],p[2],weight(0));
      E.add_row        (p[1],p[2],weight(1));
      E.multiply_row   (p[2],     weight(2));
    }
}

/*-----------------------------------------*/

void HNStructureQ22d::CondenseHangingLowerHigher(EntryMatrix& E, IntVector& indices) const
{
  for(int ii=0; ii<4; ii++) // nur 4 kandiaten koennen haengen !!
    {
      int i = indices[2*ii+1];
      if(!hanging(i)) continue;

      const fixarray<3,int>& f = regular_nodes(i);

      fixarray<3,int> p = lnoe[ii];

      if ( (indices[p[0]]==f[1]) && (indices[p[1]]==f[0]) ) 
	{ 
	  swap(p[0],p[1]);
	} 

      indices[p[2]] = f[2];

      E.add_column     (p[0],p[2],q1wei[0]);
      E.add_column     (p[1],p[2],q1wei[1]);
      E.multiply_column(p[2],     q1wei[2]);
      
      E.add_row        (p[0],p[2],weight(0));
      E.add_row        (p[1],p[2],weight(1));
      E.multiply_row   (p[2],     weight(2));
    }
}

/*-----------------------------------------*/

void HNStructureQ22d::CondenseHangingHigherLower(EntryMatrix& E, IntVector& indices) const
{
  for(int ii=0; ii<4; ii++) // nur 4 kandiaten koennen haengen !!
    {
      int i = indices[2*ii+1];
      if(!hanging(i)) continue;

      const fixarray<3,int>& f = regular_nodes(i);

      fixarray<3,int> p = lnoe[ii];

      if ( (indices[p[0]]==f[1]) && (indices[p[1]]==f[0]) ) 
	{ 
	  swap(p[0],p[1]);
	} 

      indices[p[2]] = f[2];

      E.add_column     (p[0],p[2],weight(0));
      E.add_column     (p[1],p[2],weight(1));
      E.multiply_column(p[2],     weight(2));
      
      E.add_row        (p[0],p[2],q1wei[0]);
      E.add_row        (p[1],p[2],q1wei[1]);
      E.multiply_row   (p[2],     q1wei[2]);
    }
}

/*-----------------------------------------*/

void HNStructureQ22d::CondenseHangingMixed(EntryMatrix& E, IntVector& indices, int k) const
{
  assert(indices.size()==14);

  IntVector xx(9,-1);
  xx[1] = 9; xx[3] = 10; xx[5] = 12; xx[7] = 13;

  if      (k==0) { xx[0] = 0; xx[2] = 1; xx[6] = 3; xx[8] = 4;}
  else if (k==1) { xx[0] = 1; xx[2] = 2; xx[6] = 4; xx[8] = 5;}
  else if (k==2) { xx[0] = 3; xx[2] = 4; xx[6] = 6; xx[8] = 7;}
  else if (k==3) { xx[0] = 4; xx[2] = 5; xx[6] = 7; xx[8] = 8;}
  for(int ii=0; ii<4; ii++) // nur 4 kandiaten koennen haengen !!
    {
      int i = indices[xx[2*ii+1]];;

      if(!hanging(i)) continue;

      const fixarray<3,int>& f = regular_nodes(i);

      fixarray<3,int> p = lnoe[ii];

      if ( (indices[xx[p[0]]]==f[1]) && (indices[xx[p[1]]]==f[0]) ) 
	{ 
	  swap(p[0],p[1]);
	} 

      indices[xx[p[2]]] = f[2];

      E.add_column     (xx[p[0]],xx[p[2]],weight(0));
      E.add_column     (xx[p[1]],xx[p[2]],weight(1));
      E.multiply_column(xx[p[2]],         weight(2));
      
      E.add_row        (xx[p[0]],xx[p[2]],weight(0));
      E.add_row        (xx[p[1]],xx[p[2]],weight(1));
      E.multiply_row   (xx[p[2]],         weight(2));
    }
}

/*-----------------------------------------*/

/*void HNStructureQ22d::NewCondenseHanging(EntryMatrix& E, IntVector& indices1, IntVector& indices2) const
{
  for(int ii=0; ii<4; ii++) // nur 4 kandiaten koennen haengen !!
    {
      int i = indices1[2*ii+1];

      if(!hanging(i)) continue;

      const fixarray<3,int>& f = regular_nodes(i);

      fixarray<3,int> p = lnoe[ii];

      if ( (indices1[p[0]]==f[1]) && (indices1[p[1]]==f[0]) ) 
	{ 
	  swap(p[0],p[1]);
	} 

      indices1[p[2]] = f[2];

      E.add_row        (p[0],p[2],weight(0));
      E.add_row        (p[1],p[2],weight(1));
      E.multiply_row   (p[2],     weight(2));
    }
  for(int ii=0; ii<4; ii++) // nur 4 kandiaten koennen haengen !!
    {
      int i = indices2[2*ii+1];
      if(!hanging(i)) continue;

      const fixarray<3,int>& f = regular_nodes(i);

      fixarray<3,int> p = lnoe[ii];

      if ( (indices2[p[0]]==f[1]) && (indices2[p[1]]==f[0]) ) 
	{ 
	  swap(p[0],p[1]);
	} 

      indices2[p[2]] = f[2];

      E.add_column       (p[0],p[2],weight(0));
      E.add_column       (p[1],p[2],weight(1));
      E.multiply_column  (p[2],     weight(2));
    }
    }*/

/*-----------------------------------------*/

void HNStructureQ22d::CondenseHanging(IntVector& indices) const
{
  for(int ii=0; ii<4; ii++)
    {
      int j = lnoe[ii][2];
      int i = indices[j];

      if (hanging(i))
	indices[j] = regular_nodes(i)[2];
    }
}


}
