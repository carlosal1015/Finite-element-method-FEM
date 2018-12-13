/**
*
* Copyright (C) 2009 by the Gascoigne 3D authors
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


#ifndef  __SplittingMultiLevelSolver_h
#define  __SplittingMultiLevelSolver_h

#include  "multilevelsolver.h"

/*-----------------------------------------*/

namespace Gascoigne
{

//////////////////////////////////////////////
//
///@brief
///
///
//////////////////////////////////////////////

/*----------------------------------------------------------------------------*/

class SplittingMultiLevelSolver : public MultiLevelSolver
{
protected:

  std::vector<SolverInterface*>  _SP2;

  void NewSolvers()
  {
    int oldnlevels = _SP.size();
    
    if (oldnlevels>nlevels())
      {
	for (int l=oldnlevels-1; l>=nlevels(); l--)
	  {
	    delete _SP [l]; _SP [l] = NULL;
	    delete _SP2[l]; _SP2[l] = NULL;
	  }
      }
    _SP .resize(nlevels(),NULL);
    _SP2.resize(nlevels(),NULL);
    ComputeLevel = nlevels()-1;
    
    for (int level=0; level<nlevels(); ++level)  
      {
	int solverlevel = nlevels()-1-level;
	int dim         = GetMeshAgent()->GetDimension();
	
	// new Solvers
	if (_SP[solverlevel]==NULL) 
	  {
	    _SP[solverlevel] = _NI->NewSolver(solverlevel);
	    _SP[solverlevel]->BasicInit(_paramfile,dim);
	  }
	if (_SP2[solverlevel]==NULL) 
	  {
	    _SP2[solverlevel] = _NI->NewSolver(solverlevel);
	    _SP2[solverlevel]->BasicInit(_paramfile,dim);
	  }
      }
  }
  
public:
  
  SplittingMultiLevelSolver() : MultiLevelSolver(), _SP2(0) {}
  ~SplittingMultiLevelSolver() 
  {  
    for(int i=0; i<_SP.size(); i++)  { delete _SP2[i]; _SP2[i] = NULL;}
  }
  string GetName() const {  return "SplittingMultiLevelSolver";}

  void ReInit(const std::string& problem1, const std::string& problem2)
  {  
    NewMgInterpolator();
    NewSolvers();

    SetProblem(problem1);
    SolverNewMesh();
    RegisterMatrix();
    ReInitMatrix();

    SetProblem(problem2);
    SolverNewMesh();
    RegisterMatrix();
    ReInitMatrix();
  }
        SolverInterface* GetSolver(int l)
	  {
	    if (problem1) return _SP [l];
	    else                 _SP2[l];
	  }      
  const SolverInterface* GetSolver(int l) const
	  {
	    if (problem1) return _SP [l];
	    else                 _SP2[l];
	  }      
};

#endif
