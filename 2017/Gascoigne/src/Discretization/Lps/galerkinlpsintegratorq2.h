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


#ifndef  __GalerkinLpsIntegratorQ2_h
#define  __GalerkinLpsIntegratorQ2_h


/////////////////////////////////////////////
////
////@brief
////  ... comments GalerkinLpsIntegratorQ2

////
////
/////////////////////////////////////////////

#include  "galerkinintegratorq2.h"
#include  "lpsintegrator.h"

namespace Gascoigne
{
template<int DIM>
class GalerkinLpsIntegratorQ2 : virtual public GalerkinIntegratorQ2<DIM>
{
protected:

  LpsIntegratorQ2<DIM>  Lps;

public:

//
////  Con(De)structor 
//

  GalerkinLpsIntegratorQ2<DIM>() : GalerkinIntegratorQ2<DIM>() {}
  
  ~GalerkinLpsIntegratorQ2<DIM>() {}

  std::string GetName() const {return "GalerkinLpsIntegratorQ2";}

  void Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector&U, 
      const LocalData& Q, const LocalData& QC) const;
  void Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, 
      const LocalData& Q, const LocalData& QC) const;
};
}

#endif
