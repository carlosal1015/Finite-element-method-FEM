/**
*
* Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009 by the Gascoigne 3D authors
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


#ifndef  __StdMultiLevelSolver_h
#define  __StdMultiLevelSolver_h

#include  "multilevelsolverinterface.h"
#include  "stdmultilevelsolverdata.h"
#include  "problemdescriptorinterface.h"
#include  "problemcontainer.h"
#include  "functionalcontainer.h"
#include  "monitor.h"
#include  "stopwatch.h"
#include  "mginterpolatorinterface.h"

/*-------------------------------------------------------------*/

namespace Gascoigne
{

//////////////////////////////////////////////
///
///@brief
/// Default nonlinear MultilevelSolver

/// - stores MultiGridMeshInterace
/// - stores array of MGInterpolator
/// - stores array of SolverInterface
///
//////////////////////////////////////////////

class StdMultiLevelSolver : public virtual MultiLevelSolverInterface
{
  private :

  std::vector<SolverInterface*>  _SP;
  const MeshAgentInterface* _MAP;
  std::vector<MgInterpolatorInterface*>   _Interpolator;

  protected :

  const MeshAgentInterface* GetMeshAgent() const {return _MAP;}

  std::vector<SolverInterface*>& GetSolverPointers() { return _SP; }
  const std::vector<SolverInterface*>& GetSolverPointers() const { return _SP; }
  
  mutable VectorInterface _cor, _res, _mg0, _mg1;

  mutable StopWatch   _clock_residual, _clock_solve;

  mutable int ComputeLevel;
  mutable int oldnlevels;

  const ParamFile*  _paramfile;

  Monitor*                          MON;
  StdMultiLevelSolverData*          DataP;
  const ProblemDescriptorInterface*      _PD;
  const ProblemContainer*                _PC;
  const FunctionalContainer*             _FC;
  
  virtual void NewSolvers();

  virtual SolverInterface* NewSolver(int solverlevel);
  virtual void NewMgInterpolator();
  virtual void SolverNewMesh();

  virtual const ProblemDescriptorInterface* GetProblemDescriptor() const { return _PD;}

  virtual const FunctionalContainer* GetFunctionalContainer()  const { return _FC; }
  virtual void SetFunctionalContainer(const FunctionalContainer* FC) { _FC=FC;     }

  const DoubleVector GetExactValues() const;
  const DoubleVector ComputeFunctionals(VectorInterface& f, const VectorInterface& u);
  const DoubleVector ComputeFunctionals(VectorInterface& f, const VectorInterface& u,
					FunctionalContainer* FC);
  
  virtual SolverInterface*& GetSolverPointer(int l) {assert(l<_SP.size()); return _SP[l];}
  virtual void SetComputeLevel(int level) {ComputeLevel=level;}

  virtual double NewtonNorm(const VectorInterface& u) const {
    return GetSolver(ComputeLevel)->NewtonNorm(u);
  }
  virtual void mgstep(std::vector<double>& res, std::vector<double>& rw, int l, int maxl, int minl, std::string& p0, std::string p, VectorInterface& u, VectorInterface& b, VectorInterface& v);

  virtual void Cg   (VectorInterface& x, const VectorInterface& f, CGInfo& info);
  virtual void Gmres(VectorInterface& x, const VectorInterface& f, CGInfo& info);

  virtual void ViewProtocoll() const;

  virtual void SolutionTransfer(int high, int low, VectorInterface& u) const;
  virtual void Transfer(int high, int low, VectorInterface& u) const;

 public:

  // Constructor

  StdMultiLevelSolver();
  ~StdMultiLevelSolver();

  std::string GetName() const {return "StdMultiLevelSolver";}

  std::vector<MgInterpolatorInterface*>& GetInterpolatorPointers() { return _Interpolator; }
  const std::vector<MgInterpolatorInterface*>& GetInterpolatorPointers() const { return _Interpolator; }

  void RegisterVectors();
  void RegisterMatrix();
  void ReInitMatrix();
  void ReInitVector(VectorInterface& v);
  void ReInitVector(VectorInterface& v, int comp);

  void BasicInit(const MeshAgentInterface* GMGM, const ParamFile* paramfile,
		 const ProblemContainer* PC,
		 const FunctionalContainer* FC=NULL);

  // Zugriff

//  virtual void SetState(const std::string& s) {
//    for(int l=0;l<_SP.size();l++) _SP[l]->SetState(s);
//  }

  int nlevels()                 const { assert(GetMeshAgent()); return GetMeshAgent()->nlevels();}
  virtual int FinestLevel  ()  const { return nlevels()-1;}
  virtual int CoarsestLevel()  const { return 0;}
  
  virtual const ProblemContainer* GetProblemContainer()        const { assert(_PC); return _PC; }
  virtual void SetProblemContainer(const ProblemContainer* PC)       { _PC=PC;     }

  SolverInterface* GetSolver(int l) {assert(l<_SP.size()); return _SP[l];}
  const SolverInterface* GetSolver(int l) const {assert(l<_SP.size()); return _SP[l];}
  SolverInterface* GetSolver() {assert(_SP.size()==nlevels()); return _SP[FinestLevel()];}
  const SolverInterface* GetSolver() const {assert(_SP.size()==nlevels()); return _SP[FinestLevel()];}

  void SetMonitorPtr(Monitor* mon) { MON = mon;}

  void ReInit(const std::string& problemlabel);
  void SetProblem(const std::string& label);

  // neue vektoren

  std::string LinearSolve(int level, VectorInterface& u, const VectorInterface& b, CGInfo& info);
  std::string Solve(int level, VectorInterface& x, const VectorInterface& b, NLInfo& nlinfo);
  void InterpolateSolution(VectorInterface& u, const GlobalVector& uold) const;
  void InterpolateCellSolution(VectorInterface& u, const GlobalVector& uold) const;

  virtual void NewtonVectorZero(VectorInterface& w) const;
  virtual double NewtonResidual(VectorInterface& y, const VectorInterface& x, const VectorInterface& b) const;
  virtual double NewtonUpdate(double& rr, VectorInterface& x, VectorInterface& dx, VectorInterface& r, const VectorInterface& f, NLInfo& nlinfo);
  virtual void NewtonLinearSolve(VectorInterface& x, const VectorInterface& b, CGInfo& info);
  virtual void NewtonMatrixControl(VectorInterface& u, NLInfo& nlinfo);
  virtual void NewtonOutput(NLInfo& nlinfo) const;
  virtual void NewtonPreProcess(VectorInterface& u, const VectorInterface& f,NLInfo& info) const;
  virtual void NewtonPostProcess(VectorInterface& u, const VectorInterface& f,NLInfo& info) const;


  void AssembleMatrix(VectorInterface& u, NLInfo& nlinfo);
  void AssembleMatrix(VectorInterface& u);
  /// not used in the library -- might be used in local
  void ComputeIlu(VectorInterface& u);
  void ComputeIlu();
  
  void BoundaryInit(VectorInterface& u) const;

  virtual void SolutionTransfer(VectorInterface& u) const;
  virtual void Transfer(VectorInterface& u) const;

  
  void vmulteq(VectorInterface& y, const VectorInterface&  x) const;
  
  virtual void LinearMg(int minlevel, int maxlevel, VectorInterface& u, const VectorInterface& f, CGInfo&);

  double ComputeFunctional(VectorInterface& f, const VectorInterface& u, const std::string& label);
  
  void AssembleDualMatrix(VectorInterface& u);

  // fuer gmres
  
  virtual void precondition(VectorInterface& x, VectorInterface& y);
  virtual void DeleteVector(VectorInterface& p);
  virtual void Equ(VectorInterface& dst, double s, const VectorInterface& src)const;
  void Zero(VectorInterface& dst)const;

  void AddNodeVector(const std::string& name, VectorInterface& q);
  void DeleteNodeVector(const std::string& q);

  void newton(VectorInterface& u, const VectorInterface& f, VectorInterface& r, VectorInterface& w, NLInfo& info);
};
}

#endif


