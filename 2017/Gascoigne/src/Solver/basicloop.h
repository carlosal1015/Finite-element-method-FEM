/**
*
* Copyright (C) 2004, 2005, 2007 by the Gascoigne 3D authors
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


#ifndef  __BasicLoop_h
#define  __BasicLoop_h

#include  "extrapolator.h"
#include  "multilevelsolverinterface.h"
#include  "adaptordata.h"
#include  "meshagentinterface.h"
#include  "monitor.h"
#include  "visualization.h"
#include  "visudatacompvector.h"
#include  "visudatanvector.h"

#include  "problemdescriptorinterface.h"
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

class BasicLoop
{
private:

  MeshAgentInterface*        _MA;
  MultiLevelSolverInterface* _ML;
  SolverInfos*               _SI;

 protected:

  void WriteMeshAndSolution(const std::string& filename, const VectorInterface& u) const;
  void WriteSolution(const VectorInterface& u) const;
  void WriteMesh() const;
  void WriteMeshInp(const std::string& name) const;

  MeshAgentInterface*& GetMeshAgentPointer() { return _MA;}
  MultiLevelSolverInterface*& GetMultiLevelSolverPointer() { return _ML;}

  const MeshAgentInterface* GetMeshAgent() const { assert(_MA); return _MA;}
  const MultiLevelSolverInterface* GetMultiLevelSolver() const { assert(_ML); return _ML;}

  MeshAgentInterface* GetMeshAgent() { assert(_MA); return _MA;}
  MultiLevelSolverInterface* GetMultiLevelSolver() { assert(_ML); return _ML;}

        SolverInfos*& GetSolverInfosPointer() { return _SI;}
        SolverInfos* GetSolverInfos()         { assert(_SI); return _SI;}
  const SolverInfos* GetSolverInfos() const   { assert(_SI); return _SI;}

  mutable StopWatch   _clock_newmesh, _clock_solve, _clock_functionals, _clock_write, _clock_estimate;

  int    _niter, _iter;
  bool _writeVtk;
  bool _writeBupGup;
  bool _writeInp;

  std::string _reload, _initial;
  std::string _s_resultsdir;
  const ParamFile*  _paramfile;
  

  Monitor           Mon;
  StdIoManager      IOM;

  GlobalVector _GlobalErr;

  // new vectors

  virtual std::string Solve(VectorInterface& u, VectorInterface& f, std::string name);
  virtual std::string Solve(VectorInterface& u, VectorInterface& f){
     return Solve(u, f, _s_resultsdir+"/solve");
  }

  virtual void PrintMeshInformation(int outputlevel=0) const;

  virtual void Output(const VectorInterface& u, std::string name) const;
  virtual void Output(const VectorInterface& u){
    Output(u, _s_resultsdir+"/solve");
  }

  virtual void ComputeGlobalErrors(const VectorInterface& u);

  virtual void InitSolution(VectorInterface& u);
  virtual void CopyVector(GlobalVector& dst, VectorInterface& src);
  virtual void CopyVector(VectorInterface& dst, GlobalVector& src);

public:

  BasicLoop();
  virtual ~BasicLoop();

  virtual void BasicInit(const ParamFile* paramfile,
			 const ProblemContainer* PC,
			 const FunctionalContainer* FC=NULL);

  void run(const std::string& problemlabel);
  void ClockOutput() const;
};
}

/*-----------------------------------------*/

#endif
