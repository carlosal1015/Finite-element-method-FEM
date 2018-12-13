/**
*
* Copyright (C) 2008, 2010 by the Gascoigne 3D authors
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


#ifndef  __Algorithm_h
#define  __Algorithm_h

#include  "extrapolator.h"
#include  "solverinterface.h"
#include  "adaptordata.h"
#include  "meshagentinterface.h"
#include  "monitor.h"
#include  "visualization.h"
#include  "visudatacompvector.h"
#include  "visudatanvector.h"
#include  "numericinterface.h"

#include  "problemcontainer.h"
#include  "functionalcontainer.h"
#include  "stdiomanager.h"
#include  "stopwatch.h"
#include  "paramfile.h"
#include  "vectorinterface.h"
#include  "solverinfos.h"

/*-----------------------------------------*/

namespace Gascoigne
{

//////////////////////////////////////////////
//
///@brief
/// Implementation of diverse outer loops: mesh adaptation, time stepping...

///
///
//////////////////////////////////////////////

class Algorithm
{
private:

  MeshAgentInterface*        _MA;
  SolverInfos*               _SI;
  const NumericInterface*    _NI;

 protected:

  const ParamFile*           _paramfile;

  int                        _gmresmemsize;

        MeshAgentInterface*& GetMeshAgentPointer() { return _MA;}
  const NumericInterface*   GetNumeric()           { assert(_NI); return _NI;}
        SolverInfos*& GetSolverInfosPointer()      { return _SI;}
        SolverInfos* GetSolverInfos()              { assert(_SI); return _SI;}
  const SolverInfos* GetSolverInfos()      const   { assert(_SI); return _SI;}

  virtual const SolverInterface* GetSolver() const =0;
  virtual       SolverInterface* GetSolver()       =0;

  virtual void  ReInitVector(VectorInterface& u) const=0; 
  virtual void  DeleteVector(VectorInterface& u) const=0; 

  virtual void  AssembleMatrixAndIlu(VectorInterface& u) 
  { std::cerr << "Error: AssembleMatrixAndIlu" << std::endl; abort();}

  virtual void  LinearSolve(VectorInterface& du, const VectorInterface& y, CGInfo& cginfo) 
  { std::cerr << "Error: LinearSolve" << std::endl; abort();}

  void  Newton(VectorInterface& u, const VectorInterface& f, NLInfo& nlinfo);
  void  CopyVector(GlobalVector& dst, VectorInterface& src) const;
  void GmresSolve(VectorInterface& x, const VectorInterface& b, CGInfo& info);
  virtual void  Precondition(VectorInterface& x, VectorInterface& y);

public:

  Algorithm();
  virtual ~Algorithm();

  const MeshAgentInterface*  GetMeshAgent()  const { assert(_MA); return _MA;}
        MeshAgentInterface*  GetMeshAgent()        { assert(_MA); return _MA;}

  void PrintMeshInformation() const;

  virtual void BasicInit(const ParamFile* paramfile, const NumericInterface* NI);
};
}

/*-----------------------------------------*/

#endif
