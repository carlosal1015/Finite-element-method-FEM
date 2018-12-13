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


#include  "mginterpolatormatrix.h"


using namespace std;

/*-----------------------------------------*/
  
namespace Gascoigne
{
void MgInterpolatorMatrix::restrict_zero(GlobalVector& uL, const GlobalVector& ul) const
{
  uL.zero();
  for(int i=0;i<ST.n();i++)
    {
      for(int pos=ST.start(i);pos<ST.stop(i);pos++)
	{
	  uL.add_node(ST.col(pos),val[pos],i,ul);
	}
    }
}

/*-----------------------------------------*/
  
void MgInterpolatorMatrix::prolongate_add(GlobalVector& ul, const GlobalVector& uL) const
{
  for(int i=0;i<ST.n();i++)
    {
      for(int pos=ST.start(i);pos<ST.stop(i);pos++)
	{
	  ul.add_node(i,val[pos],ST.col(pos),uL);
	}
    }
}

/*-----------------------------------------*/

void MgInterpolatorMatrix::SolutionTransfer(GlobalVector& uL, const GlobalVector& ul) const 
{
  uL.zero();
  for(int i=0;i<ST.n();i++)
    {
      for(int pos=ST.start(i);pos<ST.stop(i);pos++)
	{
	  if(val[pos]==1.)
	    uL.add_node(ST.col(pos),val[pos],i,ul);
	}
    }
}
}
