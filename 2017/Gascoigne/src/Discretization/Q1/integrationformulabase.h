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


#ifndef  __IntegrationFormulaBase_h
#define  __IntegrationFormulaBase_h

#include  "integrationformulainterface.h"
#include  "gascoigne.h"

/*-----------------------------------------*/


namespace Gascoigne
{
template<int DIM>
class IntegrationFormulaBase : public IntegrationFormulaInterface
{
private:

  typedef Vertex<DIM>   VERTEX;

  int                 _in;
  DoubleVector     _iw;
  std::vector<VERTEX> _ic;

protected:

  void ReInit(int n) {
    _in = n;
    _iw.reserve(n);      
    _iw.resize (n);
    _ic.reserve(n);      
    _ic.resize (n);
  }

public:

  IntegrationFormulaBase<DIM>() : IntegrationFormulaInterface() {}
  IntegrationFormulaBase<DIM>(int n) {
    ReInit(n);
  } 
  IntegrationFormulaBase<DIM>(const IntegrationFormulaBase<DIM>& IF) 
    : _in(IF.n()), _iw(IF.w()), _ic(IF.c()) {}
 
  int    n()                 const { return _in;}
  double w(int k)            const { return _iw[k];}
  const VERTEX& c(int k) const { return _ic[k];}
  const DoubleVector& w() const { return _iw;}

  double& w(int k) { return _iw[k];}
  VERTEX& c(int k) { return _ic[k];}

  void xi(VERTEX& v, int k)  const { v = _ic[k];}
  const std::vector<VERTEX>& c()  const { return _ic;}

};
}

#endif
