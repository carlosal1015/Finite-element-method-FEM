/**
*
* Copyright (C) 2004, 2005, 2006, 2008 by the Gascoigne 3D authors
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


#include  "galerkinglsintegratorq2.h"

using namespace std;

namespace Gascoigne
{
/*-----------------------------------------*/

template<int DIM>
void GalerkinGlsIntegratorQ2<DIM>::Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& U, 
    const LocalData& Q, const LocalData& QC) const
{
  GalerkinIntegratorQ2<DIM>::Form(EQ,F,FEM,U,Q,QC);
  Gls.                       Form(EQ,F,FEM,U,Q,QC);
}

/*-----------------------------------------*/

template<int DIM>
void GalerkinGlsIntegratorQ2<DIM>::Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, 
    const LocalData& Q, const LocalData& QC) const
{
  GalerkinIntegratorQ2<DIM>::Matrix(EQ,E,FEM,U,Q,QC);
  Gls.                       Matrix(EQ,E,FEM,U,Q,QC);
}

/*-----------------------------------------------------------*/

template class GalerkinGlsIntegratorQ2<2>;
template class GalerkinGlsIntegratorQ2<3>;
}
