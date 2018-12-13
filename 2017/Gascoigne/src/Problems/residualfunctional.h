/**
*
* Copyright (C) 2004, 2005, 2007, 2009 by the Gascoigne 3D authors
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


#ifndef  __ResidualFunctional_h
#define  __ResidualFunctional_h

#include  "functional.h"
#include  <string>
#include  <set>
#include  "dirichletdata.h"
#include  "periodicdata.h"

/*-----------------------------------------*/


namespace Gascoigne
{
  class ResidualFunctional : public virtual Functional
    {
    protected:

      nvector<int>    __comps;
      nvector<double> __scales;
      std::set<int>   __cols;

      const DirichletData*   __DD;
      const PeriodicData*    __PD;
  
    public:
  
      ResidualFunctional();
      ~ResidualFunctional();
      ResidualFunctional(const ResidualFunctional& F) : Functional(F)
	{
	  __comps  = F.GetComps();
	  __scales = F.GetScales();
	  __cols   = F.GetColors();
	  __DD     = F.GetDirichletData();
          __PD     = F.GetPeriodicData();
	}

      std::string GetName() const {return "ResidualFunctional";}

      nvector<int>     GetComps()  const {return __comps;}
      nvector<int>&    GetComps()        {return __comps;}

      nvector<double>  GetScales() const { return __scales;}
      nvector<double>& GetScales()       { return __scales;}
  
      std::set<int>    GetColors() const {return __cols;}
      std::set<int>&   GetColors()       {return __cols;}
  

      const DirichletData* GetDirichletData()   const { return __DD;}
      const DirichletData*& GetDirichletDataPointer() { return __DD;}
      const PeriodicData* GetPeriodicData()     const { return __PD;}
      const PeriodicData*& GetPeriodicDataPointer()   { return __PD;}
  
    };
}

#endif
