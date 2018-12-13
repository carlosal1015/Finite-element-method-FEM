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


#include  "hnstructureq23d.h"

using namespace std;

namespace Gascoigne
{
/*-----------------------------------------*/

HNStructureQ23d::HNStructureQ23d() : HNStructureQ13d(), q1wei(3)
{
  wei[0] =  0.375; q1wei[0] = 0.5;
  wei[1] =  0.75;  q1wei[1] = 0.5;
  wei[2] = -0.125; q1wei[2] = 0.;

  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    {
      fwei[3*i+j] = wei[i] * wei[j];
      fq1wei[3*i+j] = q1wei[i] * q1wei[j];
    }
  }

  lnoe[0][0]=0; lnoe[0][1]=1; lnoe[0][2]=2;
  lnoe[1][0]=0; lnoe[1][1]=3; lnoe[1][2]=6;
  lnoe[2][0]=2; lnoe[2][1]=5; lnoe[2][2]=8;
  lnoe[3][0]=6; lnoe[3][1]=7; lnoe[3][2]=8;
  lnoe[4][0]=0+18; lnoe[4][1]=1+18; lnoe[4][2]=2+18;
  lnoe[5][0]=0+18; lnoe[5][1]=3+18; lnoe[5][2]=6+18;
  lnoe[6][0]=2+18; lnoe[6][1]=5+18; lnoe[6][2]=8+18;
  lnoe[7][0]=6+18; lnoe[7][1]=7+18; lnoe[7][2]=8+18;
  lnoe[8][0]=0;  lnoe[8][1]=9;   lnoe[8][2]=18;
  lnoe[9][0]=2;  lnoe[9][1]=11;  lnoe[9][2]=20;
  lnoe[10][0]=6; lnoe[10][1]=15; lnoe[10][2]=24;
  lnoe[11][0]=8; lnoe[11][1]=17; lnoe[11][2]=26;

  lnop[0][0]=0; lnop[0][1]=2 ; lnop[0][2]=6; lnop[0][3]=8;  lnop[0][4]=4;
  lnop[1][0]=2; lnop[1][1]=20; lnop[1][2]=8; lnop[1][3]=26; lnop[1][4]=14;
  lnop[2][0]=6; lnop[2][1]=8;  lnop[2][2]=24;lnop[2][3]=26; lnop[2][4]=16;
  lnop[3][0]=0; lnop[3][1]=18; lnop[3][2]=6; lnop[3][3]=24; lnop[3][4]=12;
  lnop[4][0]=0; lnop[4][1]=2;  lnop[4][2]=18;lnop[4][3]=20; lnop[4][4]=10;
  lnop[5][0]=18;lnop[5][1]=20; lnop[5][2]=24;lnop[5][3]=26; lnop[5][4]=22;
}

/*-----------------------------------------*/

void HNStructureQ23d::Average(GlobalVector& u) const
{
  for(const_fiterator p=faces->begin(); p!=faces->end(); p++)
    {
      const fixarray<9,int>& f = p->second;
      u.equ_node(p->first, fwei[0], f[0], fwei[1], f[1], fwei[2], f[2]);
      u.add_node(p->first, fwei[3], f[3], fwei[4], f[4], fwei[5], f[5]);
      u.add_node(p->first, fwei[6], f[6], fwei[7], f[7], fwei[8], f[8]);
    }
  for(const_iterator p=edges->begin();p!=edges->end();p++)
    {
      const fixarray<3,int>& f = p->second;
      u.equ_node(p->first, wei[0], f[0], wei[1], f[1], wei[2], f[2]);
    }
}

/*-----------------------------------------*/

void HNStructureQ23d::Distribute(GlobalVector& u) const
{
  for(const_fiterator p=faces->begin(); p!=faces->end(); p++)
    {
      int i = p->first;
      const fixarray<9,int>& f = p->second;
      for (int j=0; j<9; j++)
	{
	  u.add_node(f[j],fwei[j],i);
	}
      u.zero_node(i);
    }
  for(const_iterator p=edges->begin();p!=edges->end();p++)
    {
      int i = p->first;
      const fixarray<3,int>& f = p->second;
      for (int j=0; j<3; j++)
	{
	  u.add_node(f[j],wei[j],i);
	}
      u.zero_node(i);
    }
}

/*-----------------------------------------*/

void HNStructureQ23d::CondenseHanging(EntryMatrix& E, IntVector& indices) const
{
  assert(indices.size()==27);

  CondenseHanging2er(E,indices);
  CondenseHanging4er(E,indices);
}

/*-----------------------------------------*/

void HNStructureQ23d::CondenseHangingLowerHigher(EntryMatrix& E, IntVector& indices) const
{
  assert(indices.size()==27);

  CondenseHanging2erLowerHigher(E,indices);
  CondenseHanging4erLowerHigher(E,indices);
}

/*-----------------------------------------*/

void HNStructureQ23d::CondenseHangingHigherLower(EntryMatrix& E, IntVector& indices) const
{
  assert(indices.size()==27);

  CondenseHanging2erHigherLower(E,indices);
  CondenseHanging4erHigherLower(E,indices);
}

/*----------------------------------------------*/

void HNStructureQ23d::CondenseHanging(IntVector& indices) const
{
  assert(indices.size()==27);

  for(int i=0; i<27; i++)
    {
      int h = indices[i];

      const_fiterator p=faces->find(h);
      if (p!=faces->end())
	{
	  indices[i] = p->second[8];
	}
      else
	{
	  const_iterator q=edges->find(h);
	  if (q!=edges->end())
	    indices[i] = q->second[2];	  
	}
    }
}

/*----------------------------------------------*/

void HNStructureQ23d::CondenseHanging2er
(EntryMatrix& E, nvector<int>& indices) const
{
  for(int i=0; i<12; i++)
    {
      fixarray<3,int>  p = lnoe[i];

      int elim = p[1];
      int h    = indices[elim];

      const_iterator q = edges->find(h);
      
      if (q==edges->end()) continue;

      const fixarray<3,int>& f = q->second;

      indices[elim] = f[2];

      if ( (indices[p[0]]==f[1]) && (indices[p[2]]==f[0]) ) 
	{ 
	  swap(p[0],p[2]);
	} 
      assert(indices[p[0]]==f[0]);
      assert(indices[p[2]]==f[1]);

      E.add_column     (p[0],elim,wei[0]);
      E.add_column     (p[2],elim,wei[1]);
      E.multiply_column(elim,     wei[2]);
      
      E.add_row        (p[0],elim,wei[0]);
      E.add_row        (p[2],elim,wei[1]);
      E.multiply_row   (elim,     wei[2]);
    }
}

/*----------------------------------------------*/

void HNStructureQ23d::CondenseHanging2erLowerHigher
(EntryMatrix& E, nvector<int>& indices) const
{
  for(int i=0; i<12; i++)
    {
      fixarray<3,int>  p = lnoe[i];

      int elim = p[1];
      int h    = indices[elim];

      const_iterator q = edges->find(h);
      
      if (q==edges->end()) continue;

      const fixarray<3,int>& f = q->second;

      indices[elim] = f[2];

      if ( (indices[p[0]]==f[1]) && (indices[p[2]]==f[0]) ) 
	{ 
	  swap(p[0],p[2]);
	} 
      assert(indices[p[0]]==f[0]);
      assert(indices[p[2]]==f[1]);

      E.add_column     (p[0],elim,q1wei[0]);
      E.add_column     (p[2],elim,q1wei[1]);
      E.multiply_column(elim,     q1wei[2]);
      
      E.add_row        (p[0],elim,wei[0]);
      E.add_row        (p[2],elim,wei[1]);
      E.multiply_row   (elim,     wei[2]);
    }
}

/*----------------------------------------------*/

void HNStructureQ23d::CondenseHanging2erHigherLower
(EntryMatrix& E, nvector<int>& indices) const
{
  for(int i=0; i<12; i++)
    {
      fixarray<3,int>  p = lnoe[i];

      int elim = p[1];
      int h    = indices[elim];

      const_iterator q = edges->find(h);
      
      if (q==edges->end()) continue;

      const fixarray<3,int>& f = q->second;

      indices[elim] = f[2];

      if ( (indices[p[0]]==f[1]) && (indices[p[2]]==f[0]) ) 
	{ 
	  swap(p[0],p[2]);
	} 
      assert(indices[p[0]]==f[0]);
      assert(indices[p[2]]==f[1]);

      E.add_column     (p[0],elim,wei[0]);
      E.add_column     (p[2],elim,wei[1]);
      E.multiply_column(elim,     wei[2]);
      
      E.add_row        (p[0],elim,q1wei[0]);
      E.add_row        (p[2],elim,q1wei[1]);
      E.multiply_row   (elim,     q1wei[2]);
    }
}

/*----------------------------------------------*/

void HNStructureQ23d::CondenseHanging4er
(EntryMatrix& E, nvector<int>& indices) const
{
  for(int i=0; i<6; i++)
    {
      fixarray<5,int> lf = lnop[i];

      int elim = lf[4];
      int h    = indices[elim];

      const_fiterator q = faces->find(h);
      
      if (q==faces->end()) continue;

      const fixarray<9,int>& gf = q->second;

      indices[elim] = gf[8];

      fixarray<8,int>  x = -1;

      for (int j=0; j<indices.size(); j++)
	{
	  int k = indices[j];
	  if      (k==gf[2]) x[2] = j;
	  else if (k==gf[5]) x[5] = j;
	  else if (k==gf[6]) x[6] = j;
	  else if (k==gf[7]) x[7] = j;
	}
      for (int j=0; j<4; j++)
	{
	  int k = indices[lf[j]];
	  if      (k==gf[0]) x[0] = lf[j];
	  else if (k==gf[1]) x[1] = lf[j];
	  else if (k==gf[3]) x[3] = lf[j];
	  else if (k==gf[4]) x[4] = lf[j];
	  else assert(0);
	}
      for (int j=0; j<8; j++) 
	{
	  assert(x[j]>=0);
	  E.add_column(x[j],elim,fwei[j]);
	  E.add_row   (x[j],elim,fwei[j]);
	}
      E.multiply_column(elim,fwei[8]);
      E.multiply_row   (elim,fwei[8]);
    }
}

/*----------------------------------------------*/

void HNStructureQ23d::CondenseHanging4erLowerHigher
(EntryMatrix& E, nvector<int>& indices) const
{
  for(int i=0; i<6; i++)
    {
      fixarray<5,int> lf = lnop[i];

      int elim = lf[4];
      int h    = indices[elim];

      const_fiterator q = faces->find(h);
      
      if (q==faces->end()) continue;

      const fixarray<9,int>& gf = q->second;

      indices[elim] = gf[8];

      fixarray<8,int>  x = -1;

      for (int j=0; j<indices.size(); j++)
	{
	  int k = indices[j];
	  if      (k==gf[2]) x[2] = j;
	  else if (k==gf[5]) x[5] = j;
	  else if (k==gf[6]) x[6] = j;
	  else if (k==gf[7]) x[7] = j;
	}
      for (int j=0; j<4; j++)
	{
	  int k = indices[lf[j]];
	  if      (k==gf[0]) x[0] = lf[j];
	  else if (k==gf[1]) x[1] = lf[j];
	  else if (k==gf[3]) x[3] = lf[j];
	  else if (k==gf[4]) x[4] = lf[j];
	  else assert(0);
	}
      for (int j=0; j<8; j++) 
	{
	  assert(x[j]>=0);
	  E.add_column(x[j],elim,fq1wei[j]);
	  E.add_row   (x[j],elim,fwei[j]);
	}
      E.multiply_column(elim,fq1wei[8]);
      E.multiply_row   (elim,fwei[8]);
    }
}

/*----------------------------------------------*/

void HNStructureQ23d::CondenseHanging4erHigherLower
(EntryMatrix& E, nvector<int>& indices) const
{
  for(int i=0; i<6; i++)
    {
      fixarray<5,int> lf = lnop[i];

      int elim = lf[4];
      int h    = indices[elim];

      const_fiterator q = faces->find(h);
      
      if (q==faces->end()) continue;

      const fixarray<9,int>& gf = q->second;

      indices[elim] = gf[8];

      fixarray<8,int>  x = -1;

      for (int j=0; j<indices.size(); j++)
	{
	  int k = indices[j];
	  if      (k==gf[2]) x[2] = j;
	  else if (k==gf[5]) x[5] = j;
	  else if (k==gf[6]) x[6] = j;
	  else if (k==gf[7]) x[7] = j;
	}
      for (int j=0; j<4; j++)
	{
	  int k = indices[lf[j]];
	  if      (k==gf[0]) x[0] = lf[j];
	  else if (k==gf[1]) x[1] = lf[j];
	  else if (k==gf[3]) x[3] = lf[j];
	  else if (k==gf[4]) x[4] = lf[j];
	  else assert(0);
	}
      for (int j=0; j<8; j++) 
	{
	  assert(x[j]>=0);
	  E.add_column(x[j],elim,fwei[j]);
	  E.add_row   (x[j],elim,fq1wei[j]);
	}
      E.multiply_column(elim,fwei[8]);
      E.multiply_row   (elim,fq1wei[8]);
    }
}

/*----------------------------------------------*/
}
