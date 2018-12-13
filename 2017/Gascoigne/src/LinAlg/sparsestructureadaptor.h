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


#ifndef  __SparseStructureAdaptor_h
#define  __SparseStructureAdaptor_h

#include  "columndiagstencil.h"
#include  "sparsestructure.h"
#include  <string>


namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  ... comments SparseStructureAdaptor

///
///
/////////////////////////////////////////////

class SparseStructureAdaptor
{
public:


private:

protected:

  int _nnodes;
  const SparseStructure* SSP;

  int n_base() const {assert(SSP);   return SSP->n();} 
  int nentries_base() const {assert(SSP); return SSP->ntotal();} 

public:


//
///  Constructor 
//
    SparseStructureAdaptor() : SSP(NULL) {}
    virtual ~SparseStructureAdaptor() {}

    virtual std::string GetName() const=0;

    void InitStructure(const SparseStructureInterface* SS) {
      SSP = dynamic_cast<const SparseStructure*>(SS);
      assert(SSP);
      _nnodes = SSP->n();
    }

    int nnodes() const {return _nnodes;}
    virtual int index(int i, int c) const=0;

    virtual int n() const=0; 
    virtual int nentries() const=0; 
    virtual void FillStencil(ColumnDiagStencil& S) const=0;
    virtual IntVector GetIndicesDirichlet(int inode, const std::vector<int>& cv)const=0;
};
}

#endif
