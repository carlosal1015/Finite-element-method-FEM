/**
*
* Copyright (C) 2004, 2005, 2006, 2007, 2010 by the Gascoigne 3D authors
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


#include  "gascoignevisualization.h"
#include  "componentinformation.h"
#include  "compose_name.h"


/*-----------------------------------------*/

namespace Gascoigne
{

void GascoigneVisualization::AddVector(const ComponentInformation* CI, const GlobalVector* v) 
{
  assert(mesh);
  assert(v);
  _v=v;

  int ncomp = v->ncomp();
  assert(ncomp);

  VDI.Clear();

  VD .SetGlobalVector(v);

  //needed for when generating vectors
  CI->SetDimension(mesh->dimension());
  
  //VDI.AddScalars(ncomp);
  {
    int         ncomp2 = CI->GetNScalars();
    if(ncomp==ncomp2){
      std::string s_name;
      int i_notset=0;
      for(int i=0;i<ncomp;i++){
        s_name="notnamed";
        CI->GetScalarName(i,s_name);
        if ( s_name=="notnamed"){
          compose_name_without_dot(s_name,i_notset); 
          i_notset++;
        }
        VDI.AddScalar(i,s_name,i);
      }
    }else{
      std::string s_name;
      for(int i=0;i<ncomp;i++){
        s_name="u";
        compose_name_without_dot(s_name,i); 
        VDI.AddScalar(i,s_name,i);
      }
    }
  }  

  //add vectors
  {
    int             nvectors    = CI->GetNVectors();
    fixarray<3,int> fa_vectorindices;
    std::string     s_name;
    int i_notset=0;
    for(int i=0;i<nvectors;i++){
      s_name="notnamed";
      CI->GetVectorName(i,s_name);
      if ( s_name=="notnamed"){
        compose_name_without_dot(s_name,i_notset); 
        i_notset++;
      }
      CI->GetVectorIndices(i,fa_vectorindices);
      VDI.AddVector(i,s_name,fa_vectorindices);
    }
  }
}

/*-----------------------------------------*/

void GascoigneVisualization::AddVector(const GlobalVector* v) 
{
  assert(mesh);
  assert(v);
  _v=v;

  int ncomp = v->ncomp();
  assert(ncomp);

  VDI.Clear();

  VD .SetGlobalVector(v);
  VDI.AddScalars(ncomp);

  int vector_index=0;
  if (ncomp>=2) 
    {
      fixarray<3,int> ff;
      int dim = mesh->dimension();
      if (dim==3)
	{
	  ff[0] = 1; ff[1] = 2; ff[2] = 3;
	}
      else if ((dim==2) && (ncomp>2))
	{
	  ff[0] = 1; ff[1] = 2; ff[2] = -1;
	}
      else if ((dim==2) && (ncomp==2))
	{
	  ff[0] = 0; ff[1] = 1; ff[2] = -1;
	}
      VDI.AddVector(vector_index,"v",ff);
      vector_index++;
    }
}

/*-----------------------------------------*/

void GascoigneVisualization::AddPointVector(const ComponentInformation* CI, const GlobalVector* v) 
{
  AddVector(CI,v);

  SetPointData(&VD);
  SetPointDataInfo(&VDI);
}

/* -------------------------------------------------------*/

void GascoigneVisualization::AddPointVector(const GlobalVector* v) 
{
  AddVector(v);

  SetPointData(&VD);
  SetPointDataInfo(&VDI);
}

/* -------------------------------------------------------*/

void GascoigneVisualization::AddCellVector(const ComponentInformation* CI, const GlobalVector* v) 
{
  AddVector(CI,v);

  SetCellData(&VD);
  SetCellDataInfo(&VDI);
}

/* -------------------------------------------------------*/

void GascoigneVisualization::AddCellVector(const GlobalVector* v) 
{
  AddVector(v);

  SetCellData(&VD);
  SetCellDataInfo(&VDI);
}
}
