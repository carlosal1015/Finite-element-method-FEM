/**
*
* Copyright (C) 2004, 2011 by the Gascoigne 3D authors
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


#include  "q1.h"
#include  "gascoignemesh.h"
#include  "pressurefilter.h"

using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{
Q1::Q1() : CellDiscretization(), HN(NULL) 
{
}

/* ----------------------------------------- */

Q1::~Q1()
{
  if (HN) delete HN;
  HN = NULL;
}

/* ----------------------------------------- */

void Q1::ReInit(const MeshInterface* MP)
{
  CellDiscretization::ReInit(MP);
  HN->ReInit(MP);
}

/* ----------------------------------------- */

void Q1::LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const
{
  IntVector indices = GetLocalIndices(iq);
  HN->CondenseHanging(E,indices);
  IntVector::const_iterator  start = indices.begin();
  IntVector::const_iterator  stop  = indices.end();
  A.entry(start,stop,__E,s);
}

/* ----------------------------------------- */

void Q1::StrongDirichletMatrix(MatrixInterface& A, int col, const vector<int>& comp) const
{
  const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
  assert(GMP);
  const IntVector& bv = *GMP->VertexOnBoundary(col);
  for(int i=0;i<bv.size();i++)
    {
      A.dirichlet(bv[i], comp);
    }  
}

/* ----------------------------------------- */

void Q1::StrongDirichletMatrixOnlyRow(MatrixInterface& A, int col, const vector<int>& comp) const
{
  const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
  assert(GMP);
  const IntVector& bv = *GMP->VertexOnBoundary(col);
  for(int i=0;i<bv.size();i++)
    {
      A.dirichlet_only_row(bv[i], comp);
    }  
}

/* ----------------------------------------- */

void Q1::StrongDirichletVectorZero(GlobalVector& u, int col, const vector<int>& comp) const
{
  const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
  assert(GMP);
  const IntVector& bv = *GMP->VertexOnBoundary(col);

  for(int ii=0;ii<comp.size();ii++)
    {
      int c = comp[ii];
      if(c<0) {
	cerr << "negative component: " << c << endl;
        abort();
      } else if(c>=u.ncomp()){
	cerr << "unknown component: " << c << endl;
        abort();
      }
    }

  for(int i=0;i<bv.size();i++)
    {
      int index = bv[i];
      for(int ii=0;ii<comp.size();ii++)
	{
	  u( index,comp[ii] ) = 0.;
	}
    }  
}

/* ----------------------------------------- */

void Q1::InterpolateSolution(GlobalVector& u, const GlobalVector& uold) const
{
  InterpolateSolutionByPatches(u,uold);
  return;

  const IntVector& vo2n = *GetMesh()->Vertexo2n();
  assert(vo2n.size()==uold.n());

  DoubleVector habschon(GetMesh()->nnodes(),0.);  
  nvector<bool> oldnode(GetMesh()->nnodes(),0);

  u.zero();
  for(int i=0;i<vo2n.size();i++)
    {
      int in = vo2n[i];

      if(in>=0) 
	{
	  u.equ_node(in,1.,i,uold);
	  oldnode[in] = 1;
	}
    }

  nmatrix<double> w = GetLocalInterpolationWeights();

  for(int iq=0; iq<GetMesh()->ncells(); iq++)
    {
      IntVector v = GetMesh()->IndicesOfCell(iq);
      for(int iol=0; iol<v.size(); iol++)
	{
	  int io = v[iol];

	  if (oldnode[io])
	    {
	      for (int inl=0; inl<v.size(); inl++)
		{
		  if (iol==inl)        continue;
		  int in = v[inl];
		  if (oldnode[in])     continue;
		  if (habschon[in]>=1) continue;

		  double weight = w(iol,inl);

		  u.add_node(in,weight,io,uold);

		  habschon[in] += weight;
		}
	      
	    }
	}
    }
}

/*-----------------------------------------*/

void Q1::HNAverage(GlobalVector& x) const
{
  HN->Average(x);
}

/* ----------------------------------------- */

void Q1::HNDistribute(GlobalVector& x) const
{
  HN->Distribute(x);
}

/* ----------------------------------------- */

void Q1::HNZero(GlobalVector& x) const
{
  HN->Zero(x);
}

/* ----------------------------------------- */

bool Q1::HNZeroCheck(const GlobalVector& x) const
{
  return HN->ZeroCheck(x);
}

/* ----------------------------------------- */

void Q1::Structure(SparseStructureInterface* SI) const
{
  SparseStructure* S = dynamic_cast<SparseStructure*>(SI);
  assert(S);

  S->build_begin(n());
  for(int iq=0;iq<GetMesh()->ncells();iq++)
    {
      IntVector indices = GetLocalIndices(iq);
      HN->CondenseHanging(indices);
      S->build_add(indices.begin(), indices.end());
    }
  HN->SparseStructureDiag(S);
  S->build_end();
}

/* ----------------------------------------- */

void Q1::Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double d) const
{
  CellDiscretization::Matrix(A,u,EQ,d);

  HN->MatrixDiag(u.ncomp(),A);
}

/* ----------------------------------------- */

void Q1::MassMatrix(MatrixInterface& A) const
{
  CellDiscretization::MassMatrix(A);

  HN->MatrixDiag(1,A);  
}

/* ----------------------------------------- */

void Q1::InitFilter(DoubleVector& F) const
{
  PressureFilter* PF = static_cast<PressureFilter*>(&F);
  assert(PF);

  if (!PF->Active()) return;

  PF->ReInit(GetMesh()->nnodes(),GetMesh()->nhanging());
  nmatrix<double> T;
  for(int iq=0; iq<GetMesh()->ncells(); ++iq)
    {
      int nv = GetMesh()->nodes_per_cell(iq);
      EntryMatrix  E(nv,1);

      Transformation(T,iq);
      GetFem()->ReInit(T);

      double cellsize = GetIntegrator()->MassMatrix(E,*GetFem());
      PF->AddDomainPiece(cellsize);

      IntVector ind = GetMesh()->IndicesOfCell(iq);
      HN->CondenseHanging(E,ind);

      for(int i=0;i<ind.size();i++)
 	{
	  for(int j=0;j<ind.size();j++)
	    {
	      F[ind[j]] += E(i,j,0,0);
	    }
      	}
    }
}

/*----------------------------------------------*/

}
