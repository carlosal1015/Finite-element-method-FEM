/**
*
* Copyright (C) 2004, 2009 by the Gascoigne 3D authors
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


#include  "q1lps3d.h"

namespace Gascoigne
{
/* ----------------------------------------- */
  
Q1Lps3d::~Q1Lps3d() 
{
  if (S) delete S;
  S=0;
}
  
/* ----------------------------------------- */

void Q1Lps3d::BasicInit(const ParamFile* paramfile)
{
  Q13d::BasicInit(paramfile);
  S = new Q1LpsStab3d;
  S  ->BasicInit(paramfile,HN);
}

/* ----------------------------------------- */

void Q1Lps3d::ReInit(const MeshInterface* M)
{
  Q13d::ReInit(M);
  S   ->ReInit(M);
}

/* ----------------------------------------- */

void Q1Lps3d::Structure(SparseStructureInterface* SI) const
{
  S->Structure(SI);
}

/* ----------------------------------------- */

void Q1Lps3d::StabForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  S->Form(f,u,EQ,d);
}

/* ----------------------------------------- */

void Q1Lps3d::Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  Q13d::Form(f,u,EQ,d);
  S   ->Form(f,u,EQ,d);
}

/* ----------------------------------------- */

void Q1Lps3d::Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double d) const
{
  Q13d::Matrix(A,u,EQ,d);
  S   ->Matrix(A,u,EQ,d);
}
}
