/**
*
* Copyright (C) 2004, 2005, 2006, 2009 by the Gascoigne 3D authors
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


#ifndef  __BasicIntegrator_h
#define  __BasicIntegrator_h


#include  "gascoigne.h"
#include  "integratorinterface.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments BasicIntegrator

////
////
/////////////////////////////////////////////

class BasicIntegrator : public IntegratorInterface
{
 private:
  
  
 protected:
  
  mutable FemFunction   _NNN;
  mutable TestFunction  _NN;
  mutable FemFunction   _UH;
  mutable FemData       _QH;
  mutable CellData      _QCH;
  

 public:

  void  universal_point(const FemInterface& FEM, FemFunction& UH, const LocalVector& U) const;
  void  universal_point(CellFunction& UCH, const LocalVector& UC,int i=0) const;
  void  universal_point(FemFunction& UH, const LocalVector& U, const FemFunction& NN) const;
  
  void  universal_point(const FemInterface& FEM, FemData& QH, const LocalData& Q) const;
  void  universal_point(CellData& QCH, const LocalData& QC,int i=0) const;

  
  
  //
  ////  Con(De)structor 
  //

  BasicIntegrator();
  ~BasicIntegrator() {}
  
};
}

#endif
