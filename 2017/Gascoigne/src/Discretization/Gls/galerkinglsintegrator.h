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


#ifndef  __GalerkinGlsIntegrator_h
#define  __GalerkinGlsIntegrator_h

#include  "galerkinintegrator.h"
#include  "glsintegrator.h"

/*-----------------------------------------*/

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments GlsIntegrator

////
////
/////////////////////////////////////////////

template<int DIM>
class GalerkinGlsIntegrator : public GalerkinIntegrator<DIM>
{
 protected:

  GlsIntegrator<DIM>  Gls;

 public:


  GalerkinGlsIntegrator<DIM>() : GalerkinIntegrator<DIM>() {};
  ~GalerkinGlsIntegrator<DIM>() {}

  std::string GetName() const {return "GalerkinGlsIntegrator";}

  void Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector&U, const LocalData& Q, 
      const LocalData& QC) const;
  void Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, const LocalData& Q, 
      const LocalData& QC) const;
};
}

/*-----------------------------------------*/


#endif
