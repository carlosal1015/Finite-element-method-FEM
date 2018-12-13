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


#ifndef  __Gascoigne_h
#define  __Gascoigne_h


/////////////////////////////////////////////
////
////@brief
////  ... comments Gascoigne

////
////
/////////////////////////////////////////////

#include  <set>
#include  <map>
#include  <string>
#include  "compvector.h"
#include  "nmatrix.h"
#include  "derivativevector.h"

namespace Gascoigne
{
  typedef CompVector<double>                                        GlobalVector;
  typedef CompVector<double>                                        LocalVector;
  typedef std::map<const std::string,const GlobalVector*>           GlobalData;
  typedef std::map<const std::string,LocalVector>                   LocalData;

  typedef nvector<double>                                           GlobalParameterVector;
  typedef nvector<double>                                           LocalParameterVector;
  typedef std::map<const std::string,const GlobalParameterVector*>  GlobalParameterData;
  typedef std::map<const std::string,LocalParameterVector>          LocalParameterData;

  typedef nvector<int>                                              IntVector;
  typedef nvector<double>                                           DoubleVector;
  typedef nmatrix<double>                                           DoubleMatrix;
  typedef std::set<int>                                             IntSet;
  typedef CompVector<double>::iterator                              VectorIterator;

  typedef nmatrix<double>                                           TimePattern;

  typedef DerivativeVector                                          TestFunction;
  typedef std::vector<TestFunction>                                 FemFunction;
  typedef std::map<const std::string,FemFunction>                   FemData;

  typedef nvector<double>                                           CellFunction;
  typedef std::map<const std::string,CellFunction>                  CellData;
}

#endif
