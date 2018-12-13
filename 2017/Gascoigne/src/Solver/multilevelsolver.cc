/**
*
* Copyright (C) 2008 by the Gascoigne 3D authors
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


#include  "multilevelsolver.h"
#include  "mginterpolatornested.h"
#include  "stdtimesolver.h"

using namespace std;

/*-------------------------------------------------------------*/

namespace Gascoigne
{

/*-------------------------------------------------------------*/

MultiLevelSolver::MultiLevelSolver() : _SP(0), _Interpolator(0), _MAP(NULL), _paramfile(NULL),
				       _PC(NULL), _PD(NULL), _NI(NULL)
{
}

/*-------------------------------------------------------------*/

MultiLevelSolver::~MultiLevelSolver()
{
  for(int i=0;i<_SP.size();i++) 
    { 
      if (_SP[i]) 
        {
          delete _SP[i]; 
          _SP[i]=NULL; 
        }
    }
   for(int i=0; i<_Interpolator.size(); i++)  
    {
      if (_Interpolator[i]) 
        {
          delete _Interpolator[i]; 
          _Interpolator[i]=NULL;
        }
    }
}

/*-------------------------------------------------------------*/

void MultiLevelSolver::BasicInit(const NumericInterface* NI, const MeshAgentInterface* MAP, 
				 const ParamFile* paramfile, const ProblemContainer* PC)
{
  _NI        = NI;
  _MAP       = MAP;
  _paramfile = paramfile;
  _PC        = PC;
}

/*-------------------------------------------------------------*/

void MultiLevelSolver::ReInitVector(VectorInterface& v)
{
  for (int level=0; level<nlevels(); ++level)  
    {
      GetSolver(level)->ReInitVector(v);
    }
}

/*-------------------------------------------------------------*/

void MultiLevelSolver::DeleteVector(VectorInterface& v)
{
  for(int l=0; l<nlevels(); ++l)  
    {
      GetSolver(l)->DeleteVector(v);
    }
}

/*-------------------------------------------------------------*/

void MultiLevelSolver::NewSolvers()
{
  int oldnlevels = _SP.size();

  if (oldnlevels>nlevels())
    {
      for (int l=oldnlevels-1; l>=nlevels(); l--)
        {
          delete GetSolver(l);
          _SP[l] = NULL;
        }
    }
  _SP.resize(nlevels(),NULL);
  ComputeLevel = nlevels()-1;

  for (int level=0; level<nlevels(); ++level)  
    {
      int solverlevel = nlevels()-1-level;
      int dim         = GetMeshAgent()->GetDimension();

      // new Solvers
      if (GetSolver(solverlevel)==NULL) 
        {
          _SP[solverlevel] = _NI->NewSolver(solverlevel);
          GetSolver(solverlevel)->BasicInit(_paramfile,dim,_NI);
        }
    }
}

/*-------------------------------------------------------------*/

void MultiLevelSolver::RegisterMatrix()
{
  for(int level=0; level<nlevels(); ++level)  
    {
      GetSolver(level)->RegisterMatrix();
    }
}

/*-------------------------------------------------------------*/

void MultiLevelSolver::ReInitMatrix()
{
  for(int level=0; level<nlevels(); ++level)  
    {
      GetSolver(level)->ReInitMatrix();
    }
}

/*-------------------------------------------------------------*/

void MultiLevelSolver::SolverNewMesh()
{
  for(int level=0; level<nlevels(); ++level)  
    {
      const MeshInterface* MIP = GetMeshAgent()->GetMesh(level);
      assert(MIP);

      int solverlevel = nlevels()-1-level;
      GetSolver(solverlevel)->NewMesh(MIP);
    }
}

/*-------------------------------------------------------------*/

void MultiLevelSolver::NewMgInterpolator()
{
  for (int i=0;i<_Interpolator.size();++i)
    {
      assert(_Interpolator[i]!=NULL);
      delete _Interpolator[i];
      _Interpolator[i]=NULL;
    }
  _Interpolator.resize(nlevels()-1,NULL);

  for(int l=0; l<nlevels()-1; ++l)  
    {
      _Interpolator[l] = new MgInterpolatorNested;
    }
  //
  // Interpolator [l] :   interpoliert   (l+1)->l  (fein->grob)
  //
  for (int level=0;level<nlevels()-1;++level)
    {
      int sl = nlevels()-level-2;

      const MeshTransferInterface* MT = GetMeshAgent()->GetTransfer(sl);
      assert(MT);
      assert(_Interpolator[level]);
      GetSolver(level)->ConstructInterpolator(_Interpolator[level],MT);
    }
}

/*-------------------------------------------------------------*/

void MultiLevelSolver::SetProblem(const std::string& label)
{
  _PD = GetProblemContainer()->GetProblem(label);
  for(int level=0; level<nlevels(); ++level)  
    {
      int solverlevel = nlevels()-1-level;
      assert(GetSolver(solverlevel)) ;
      GetSolver(solverlevel)->SetProblem(*_PD);
    }
}

/*-------------------------------------------------------------*/

void MultiLevelSolver::ReInit(const std::string& problemlabel)
{  
  NewSolvers();
  SolverNewMesh();
  NewMgInterpolator();
  SetProblem(problemlabel);
  RegisterMatrix();
  ReInitMatrix();
}

/*-------------------------------------------------------------*/

void MultiLevelSolver::AssembleMatrix(VectorInterface& u)
{
  SolutionTransfer(u);
  for(int l=0; l<=ComputeLevel; l++)
    {
      GetSolver(l)->MatrixZero();
      GetSolver(l)->AssembleMatrix(u,1.);
    }
}

/*-------------------------------------------------------------*/

void MultiLevelSolver::AssembleDualMatrix(VectorInterface& u)
{
  for(int l=0; l<nlevels(); l++)
    {
      GetSolver(l)->AssembleDualMatrix(u,1.);
    }
}

/*-------------------------------------------------------------*/

void MultiLevelSolver::Transfer(int high, int low, VectorInterface& u) const
{
  for(int l=high;l>=low;l--)
    {
      GetSolver(l)->HNAverage(u);

      assert(_Interpolator[l-1]);
      _Interpolator[l-1]->SolutionTransfer(GetSolver(l-1)->GetGV(u),
					   GetSolver(l)  ->GetGV(u));
      
      GetSolver(l)->HNZero(u);
    }
}

/*-------------------------------------------------------------*/

void MultiLevelSolver::SolutionTransfer(VectorInterface& u) const
{
  Transfer(ComputeLevel,1,u);

  for(int l=ComputeLevel; l>=1; l--)
    GetSolver(l-1)->SetBoundaryVector(u);
}

/*-------------------------------------------------------------*/

void MultiLevelSolver::ComputeIlu()
{
  for(int l=0;l<=ComputeLevel;l++)
    {
      GetSolver(l)->ComputeIlu();
    }
}

/*-------------------------------------------------------------*/

void MultiLevelSolver::ComputeIlu(VectorInterface& u)
{
  SolutionTransfer(u);
  for(int l=0;l<=ComputeLevel;l++)
    {
      GetSolver(l)->ComputeIlu(u);
    }
}

/*-------------------------------------------------------------*/

}
