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


#ifndef  __Q2Lps3d_h
#define  __Q2Lps3d_h

#include  "q23d.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments Q2Lps3d

////
////
/////////////////////////////////////////////

class Q2Lps3d : public virtual Q23d
{
public:

//
////  Con(De)structor 
//

  Q2Lps3d() : Q23d() {}
  ~Q2Lps3d() {}

  std::string GetName() const {return "Q2Lps3d";}
  
  void BasicInit(const ParamFile* paramfile);
};

}
#endif
