/**
*
* Copyright (C) 2004, 2005, 2006, 2009 by the Gascoigne 3D authors
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


#ifndef  __MeshAgentInterface_h
#define  __MeshAgentInterface_h


#include  "meshinterface.h"
#include  "multigridmeshinterface.h"
#include  "meshtransferinterface.h"
#include  "paramfile.h"
#include  "boundaryfunction.h"
#include  "stdperiodicmapping.h"
#include  <string>

namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ////
  ////@brief
  ////  ... comments MeshAgentInterface

  ////
  ////
  /////////////////////////////////////////////

  class MeshAgentInterface
  {
    private:

    protected:

    public:
      //
      ////  Con(De)structor 
      //
      MeshAgentInterface() {}
      virtual ~MeshAgentInterface() {}

      virtual int GetDimension() const=0;

      virtual void BasicInit(const ParamFile* pf)=0;
      virtual void BasicInit(const std::string& gridname, int dim, int patchdepth, int epatcher, bool goc2nc=false)=0;

      virtual int nnodes() const=0;
      virtual int ncells() const=0;
      virtual int nlevels() const=0;

      virtual void read_gup(const std::string& fname)=0;
      virtual void read_gip(const std::string& fname)=0;
      virtual void write_gup(const std::string& fname) const=0;
      virtual void write_gip(const std::string& fname) const=0;
      virtual void write_inp(const std::string& fname) const=0;
      virtual const MeshInterface* GetMesh(int l) const=0;

      virtual void global_patch_coarsen(int n)=0;
      virtual void global_refine(int n)=0;
      virtual void refine_nodes(IntVector& refnodes, IntVector& coarsenodes)=0;
      virtual void refine_nodes(IntVector& refnodes)=0;
      virtual void refine_cells(IntVector& refnodes)=0;
      virtual void random_patch_refine(double p, int n)=0;
      virtual void random_patch_coarsen(double p, int n)=0;
      virtual const MeshTransferInterface* GetTransfer(int l) const=0; 
      virtual const std::set<int> Cello2n(int i)const=0;
      virtual const int Cello2nFather(int i)const=0;
      virtual const bool Goc2nc()const=0;

      virtual void AddShape(int col, BoundaryFunction<2>* f) {std::cerr << "MeshAgentInterface::AddShape not written" << std::endl; abort();}
      virtual void AddShape(int col, BoundaryFunction<3>* f) {std::cerr << "MeshAgentInterface::AddShape not written" << std::endl; abort();}

      virtual void AddPeriodicMapping(int col, int col2, PeriodicMapping* map)
      {std::cerr << "MeshAgentInterface::AddPeriodicMapping not written" << std::endl; abort();}

  };
}

#endif
