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


#ifndef  __SimpleSparseStructureAdaptor_h
#define  __SimpleSparseStructureAdaptor_h

#include  "sparsestructureadaptor.h"


namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  ... comments SimpleSparseStructureAdaptor

///
///
/////////////////////////////////////////////

class SimpleSparseStructureAdaptor : public SparseStructureAdaptor
{
public:


private:


protected:


public:

  SimpleSparseStructureAdaptor() : SparseStructureAdaptor() {}

  std::string GetName() const {return "Simple";}

  int n() const {return n_base();} 
  int nentries() const {return nentries_base();} 
  void FillStencil(ColumnDiagStencil& S) const;

  int index(int i, int c) const {return i;}

  IntVector GetIndicesDirichlet(int inode, const std::vector<int>& cv) const{return IntVector(1,inode);}
};
}

#endif
