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


#include  "basicintegrator.h"


using namespace std;

/*---------------------------------------------------------*/

namespace Gascoigne
{
BasicIntegrator::BasicIntegrator() : IntegratorInterface() 
{
}

/*---------------------------------------------------------*/

void  Gascoigne::BasicIntegrator::universal_point(CellFunction& UCH, const LocalVector& UC, int i) const
{
  UCH.resize(UC.ncomp());
  for (int c=0; c<UC.ncomp(); c++)
    {
      UCH[c] = UC(i,c);
    }
}

/*---------------------------------------------------------*/

void BasicIntegrator::universal_point(const FemInterface& FEM, FemFunction& UH, const LocalVector& U) const
{
  UH.resize(U.ncomp());

  for (int c=0; c<UH.size(); c++)  UH[c].zero();

  for (int i=0; i<FEM.n(); i++)
    {
      FEM.init_test_functions(_NN,1.,i);
      for (int c=0; c<UH.size(); c++)
	{
	  UH[c].add(U(i,c),_NN);
	}
    }
}

/*---------------------------------------------------------*/

void BasicIntegrator::universal_point(FemFunction& UH, const LocalVector& U, const FemFunction& NN) const
{
  UH.resize(U.ncomp());
  for (int c=0; c<UH.size(); c++)  UH[c].zero();

  for (int i=0; i<NN.size(); i++)
    {
      for (int c=0; c<UH.size(); c++)
	{
	  UH[c].add(U(i,c),NN[i]);
	}
    }
}

/*---------------------------------------------------------*/

void BasicIntegrator::universal_point(const FemInterface& FEM, FemData& QH, const LocalData& Q) const
{
    QH.clear();
  LocalData::const_iterator p=Q.begin();
  for(; p!=Q.end(); p++)
    {
	universal_point(FEM, QH[p->first], p->second);
    }
}

/*---------------------------------------------------------*/

void  Gascoigne::BasicIntegrator::universal_point(CellData& QCH, const LocalData& QC, int i) const
{
  QCH.clear();
  LocalData::const_iterator p=QC.begin();
  for(; p!=QC.end(); p++)
    {
      universal_point(QCH[p->first], p->second, i);
    }
}
}
