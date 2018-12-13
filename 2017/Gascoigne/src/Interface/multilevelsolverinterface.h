/**
*
* Copyright (C) 2004, 2005, 2006, 2008 by the Gascoigne 3D authors
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


#ifndef  __MultiLevelSolverInterface_h
#define  __MultiLevelSolverInterface_h


#include  "meshagentinterface.h"
#include  "solverinterface.h"
#include  "monitor.h"
#include  "paramfile.h"
#include  "nlinfo.h"
#include  "vectorinterface.h"
#include  "problemcontainer.h"
#include  "functionalcontainer.h"

namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments MultiLevelSolverInterface

  ///
  ///
  /////////////////////////////////////////////

  class MultiLevelSolverInterface
  {
    private:

    protected:
      

    public:
      MultiLevelSolverInterface() {}
      virtual ~MultiLevelSolverInterface() {}

      virtual std::string GetName() const=0;
      virtual void BasicInit(const MeshAgentInterface* GMGM, const ParamFile* paramfile,
			     const ProblemContainer* PC, const FunctionalContainer* FC=NULL)=0;
      virtual void SetProblem(const std::string& problemlabel)=0;
      virtual void ReInit(const std::string& problemlabel)=0;
      virtual void SetMonitorPtr(Monitor* mon)=0;

      virtual const DoubleVector GetExactValues() const=0;
      virtual const DoubleVector ComputeFunctionals(VectorInterface& f, const VectorInterface& u)=0;
      virtual const DoubleVector ComputeFunctionals(VectorInterface& f, const VectorInterface& u,
						    FunctionalContainer* FC)=0;
      
      virtual void ReInitMatrix()=0;

      virtual int nlevels() const=0;

      virtual SolverInterface* GetSolver(int l)=0;
      virtual const SolverInterface* GetSolver(int l) const=0;
      virtual SolverInterface* GetSolver()=0;
      virtual const SolverInterface* GetSolver() const=0;
      virtual const ProblemContainer* GetProblemContainer() const=0;

//      virtual void SetState(const std::string& s)=0;
      virtual void AssembleMatrix(VectorInterface& u, NLInfo& nlinfo)=0;
      virtual void AssembleMatrix(VectorInterface& u)=0;
      virtual void ComputeIlu(VectorInterface& u)=0;
      virtual void ComputeIlu()=0;
      
      virtual void BoundaryInit(VectorInterface& u) const=0;

      //
      /// vector - manamgement
      //

      virtual void DeleteVector(VectorInterface& g)=0;
      virtual void RegisterVectors()=0;
      virtual void RegisterMatrix()=0;
      virtual void ReInitVector(VectorInterface& v)=0;
      virtual void ReInitVector(VectorInterface& v, int comp)=0;

      //
      /// vector 
      //

      virtual std::string LinearSolve(int level, VectorInterface& u, const VectorInterface& b, CGInfo& info)=0;
      virtual std::string LinearSolve(VectorInterface& u, const VectorInterface& b, CGInfo& info) {
        return LinearSolve(nlevels()-1,u,b,info);
      }

      virtual std::string Solve(int level, VectorInterface& x, const VectorInterface& b, NLInfo& nlinfo)=0;
      virtual std::string Solve(VectorInterface& x, const VectorInterface& b, NLInfo& nlinfo) {
        return Solve(nlevels()-1,x,b,nlinfo);
      }
      virtual void InterpolateSolution(VectorInterface& u, const GlobalVector& uold) const=0;
      virtual void InterpolateCellSolution(VectorInterface& u, const GlobalVector& uold) const=0;
      virtual double ComputeFunctional(VectorInterface& f, const VectorInterface& u, const std::string& label)=0;

      virtual void AssembleDualMatrix(VectorInterface& u)=0;
      virtual void vmulteq(VectorInterface& y, const VectorInterface&  x) const=0;

      virtual void Equ(VectorInterface& dst, double s, const VectorInterface& src)const=0;
      virtual void Zero(VectorInterface& dst)const=0;
      virtual void AddNodeVector(const std::string& name, VectorInterface& q)=0;
      virtual void DeleteNodeVector(const std::string& q)=0;

      virtual void SolutionTransfer(VectorInterface& u) const=0;
      virtual void Transfer(VectorInterface& u) const=0;

      virtual void newton(VectorInterface& u, const VectorInterface& f, VectorInterface& r, VectorInterface& w, NLInfo& info)=0;

  virtual std::vector<MgInterpolatorInterface*>& GetInterpolatorPointers()=0;
  virtual const std::vector<MgInterpolatorInterface*>& GetInterpolatorPointers() const=0;
  };
}

#endif
