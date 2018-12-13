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


#ifndef __dwrlps2d_h
#define __dwrlps2d_h

#include "q1lps2d.h"

namespace Gascoigne
{
/*-------------------------------------------------*/

class DwrLps2d : public Q1Lps2d
{
 protected:

  const PatchMesh* GetPatchMesh() const {
    const PatchMesh* MP = dynamic_cast<const PatchMesh*>(GetMesh());
    assert(MP);
    return MP;
  }

 public:

  DwrLps2d() : Q1Lps2d() {}
  ~DwrLps2d() {}

  std::string GetName() const {return "DwrLps2d";}
  void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
};
}
/*-------------------------------------------------*/

#endif
