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


#ifndef  __Q1Gls3d_h
#define  __Q1Gls3d_h

#include  "q13d.h"
#include  "glsintegrator.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments Q1Gls3d

////
////
/////////////////////////////////////////////

class Q1Gls3d : public Q13d
{
protected:

/*   GlsIntegrator<3> GlsInt;   */

public:

  //
  ////  Con(De)structor 
  //
  
  Q1Gls3d() : Q13d() {}
  ~Q1Gls3d() {}
  
  std::string GetName() const {return "Q1Gls3d";}
  
  void BasicInit(const ParamFile* pf);
};
}

#endif
