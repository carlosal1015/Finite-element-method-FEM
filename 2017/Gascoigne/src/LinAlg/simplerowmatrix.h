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


#ifndef  __SimpleRowMatrix_h
#define  __SimpleRowMatrix_h

#include  "matrixinterface.h"
#include  "rowcolumnstencil.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments SimpleRowMatrix

////
////
/////////////////////////////////////////////

class SimpleRowMatrix : virtual public MatrixInterface
{
private:
  
  
protected:
  
  RowColumnStencil  ST;
  nvector<double> value; 
  
public:
  
  
  //
  ////  Con(De)structor 
  //
  
  SimpleRowMatrix() : MatrixInterface() {}
  ~SimpleRowMatrix() {}
  
  std::string GetName() const {return "SimpleRowMatrix";}
  
  std::ostream& Write(std::ostream& os) const;
  
  const StencilInterface* GetStencil() const { return &ST;}
  double& GetValue(int pos) {return value[pos];}
  const double& GetValue(int pos) const {return value[pos];}
  const double& GetValue(int i, int j) const {return value[ST.Find(i,j)];}
  
  void zero() {value.zero();}
  void ReInit(int n, int nentries);
  void ReInit(const SparseStructureInterface* S);
  void entry(niiterator start, niiterator stop, const EntryMatrix& M, double s=1.);
  void vmult(nvector<double>& y, const nvector<double>& x, double d=1.) const;
  
};
}

#endif
