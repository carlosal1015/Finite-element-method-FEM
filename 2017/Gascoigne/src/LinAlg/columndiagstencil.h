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


#ifndef  __ColumnDiagStencil_h
#define  __ColumnDiagStencil_h

#include  "columnstencil.h"


namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments ColumnDiagStencil

////
////
/////////////////////////////////////////////

class ColumnDiagStencil : public ColumnStencil
{
protected:

  IntVector   sdiag;

public:

//
////  Con(De)structor 
//
  ColumnDiagStencil() : ColumnStencil() {}
  ~ColumnDiagStencil() {}

  const IntVector&  diag() const { return sdiag; }
        IntVector&  diag()       { return sdiag; }
        int&           diag(int i)       { return sdiag[i]; } 
  const int&           diag(int i) const { return sdiag[i]; } 

  void memory(int n, int nt);
  void memory(const SparseStructureInterface*);
};
}

#endif
