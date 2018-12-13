/**
*
* Copyright (C) 2004, 2005 by the Gascoigne 3D authors
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


#ifndef  __Q2Gls2d_h
#define  __Q2Gls2d_h

#include  "q22d.h"
#include  "q22dwithsecond.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments Q2Gls2d

////
////
/////////////////////////////////////////////

class Q2Gls2d : public Q22dWithSecond
{
public:

//
////  Con(De)structor 
//
  Q2Gls2d() : Q22dWithSecond() {}
  ~Q2Gls2d() {}

  std::string GetName() const {return "Q2Gls2d";}
  
  void BasicInit(const ParamFile* paramfile);
};

}
#endif
