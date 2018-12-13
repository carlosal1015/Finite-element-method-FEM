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


#ifndef  __NodeSparseStructureAdaptor_h
#define  __NodeSparseStructureAdaptor_h

#include  "sparsestructureadaptor.h"


namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  ... comments NodeSparseStructureAdaptor

///
///
/////////////////////////////////////////////

class NodeSparseStructureAdaptor : public SparseStructureAdaptor
{
public:


private:


protected:

  int _ncomp;

public:


  NodeSparseStructureAdaptor(int ncomp) : SparseStructureAdaptor(), _ncomp(ncomp)  {}

  std::string GetName() const {return "Node";}

  int n() const { return _ncomp*n_base();}
  int nentries() const { return _ncomp*_ncomp*nentries_base();}

  int index(int i, int c) const {return i*_ncomp+c;}

  void FillStencil(ColumnDiagStencil& S) const;
  IntVector GetIndicesDirichlet(int inode, const std::vector<int>& cv) const;
};
}

#endif
