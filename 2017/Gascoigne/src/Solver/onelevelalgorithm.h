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


#ifndef  __OneLevelAlgorithm_h
#define  __OneLevelAlgorithm_h

#include  "algorithm.h"
#include  "stdsolver.h"

/*-----------------------------------------*/

namespace Gascoigne
{

//////////////////////////////////////////////
//
///@brief
///
///
//////////////////////////////////////////////

class OneLevelAlgorithm : public Algorithm
{
 private:

  SolverInterface*        _S;
  const ProblemContainer* _PC;

 protected:

  virtual const SolverInterface* GetSolver() const { return _S;}
  virtual       SolverInterface* GetSolver()       { return _S;}

  void JacobiSolver(VectorInterface& du, const VectorInterface& f, CGInfo& info);
  void IluSolver(VectorInterface& du, const VectorInterface& f, CGInfo& info);
  void CGSolver(VectorInterface& du, const VectorInterface& f, CGInfo& info, bool precondition);
  
  void  ReInitVector(VectorInterface& u) const { _S->ReInitVector(u);}
  void  DeleteVector(VectorInterface& u) const { _S->DeleteVector(u);}
  void  AssembleMatrixAndIlu(VectorInterface& u);
  void  LinearSolve(VectorInterface& du, const VectorInterface& y, CGInfo& cginfo);

public:

  OneLevelAlgorithm() :  _S(NULL) {}
  virtual ~OneLevelAlgorithm() { if (_S) delete _S;}

  virtual void BasicInit(const ParamFile *paramfile,
                         const NumericInterface *NI,
                         const ProblemContainer *PC);

  void Precondition(VectorInterface& x, VectorInterface& y);
  void RunLinear   (const std::string& problemlabel,
                    const std::string& solver = "ILU");
  void RunNonLinear(const std::string& problemlabel);
};
}

/*-----------------------------------------*/

#endif
