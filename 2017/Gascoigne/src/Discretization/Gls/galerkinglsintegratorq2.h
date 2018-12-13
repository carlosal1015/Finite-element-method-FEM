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


#ifndef  __GalerkinGlsIntegratorQ2_h
#define  __GalerkinGlsIntegratorQ2_h


#include  "galerkinintegratorq2.h"
#include  "glsintegratorq2.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments GalerkinGlsIntegratorQ2

////
////
/////////////////////////////////////////////

template<int DIM>
class GalerkinGlsIntegratorQ2 : public GalerkinIntegratorQ2<DIM>
{
protected:

  GlsIntegratorQ2<DIM>  Gls;

public:


//
////  Con(De)structor 
//

  GalerkinGlsIntegratorQ2<DIM>() : GalerkinIntegratorQ2<DIM>() {}
  ~GalerkinGlsIntegratorQ2<DIM>() {}

  std::string GetName() const {return "GalerkinGlsIntegratorQ2";}

  void Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector&U, 
      const LocalData& Q, const LocalData& QC) const;
  void Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, 
      const LocalData& Q, const LocalData& QC) const;

};
}

#endif
