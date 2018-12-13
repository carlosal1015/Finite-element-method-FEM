/**
*
* Copyright (C) 2004, 2005, 2006, 2009, 2011 by the Gascoigne 3D authors
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


#include  "q2.h"
#include  "gascoignemesh.h"
#include  "pressurefilter.h"

using namespace Gascoigne;
using namespace std;

/* ----------------------------------------- */

Q2::Q2() : PatchDiscretization(), HN(NULL)
{
}

/* ----------------------------------------- */

Q2::~Q2()
{
  assert(HN==NULL);
}

/* ----------------------------------------- */

void Q2::ReInit(const MeshInterface* MP)
{
  PatchDiscretization::ReInit(MP);
  HN->ReInit(MP);
}

/* ----------------------------------------- */

int Q2::n() const
{
  return GetMesh()->nnodes();
}

/* ----------------------------------------- */

int Q2::nc() const
{
  return GetMesh()->ncells();
}

/* ----------------------------------------- */

int Q2::n_withouthanging() const
{
  return GetMesh()->nnodes()-HN->nhnodes();
}

/* ----------------------------------------- */

void Q2::Structure(SparseStructureInterface* SI) const
{
  SparseStructure* S = dynamic_cast<SparseStructure*>(SI);
  assert(S);

  S->build_begin(n());
  for(int iq=0;iq<GetPatchMesh()->npatches();iq++)
    {
      nvector<int> indices = GetLocalIndices(iq);
      HN->CondenseHanging(indices);
      S->build_add(indices.begin(), indices.end());
    }
  HN->SparseStructureDiag(S);
  S->build_end();  
}

/* ----------------------------------------- */

void Q2::LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const
{
  nvector<int> indices = GetLocalIndices(iq);
  HN->CondenseHanging(E,indices);
  nvector<int>::const_iterator  start = indices.begin();
  nvector<int>::const_iterator  stop  = indices.end();
  A.entry(start,stop,E,s);
}

/* ----------------------------------------- */

void Q2::Interpolate(GlobalVector& u, const DomainInitialCondition& U) const
{
  std::cerr << "\"Q2::Interpolate\" not written!" << std::endl;
  abort();

//   if (&U==NULL) return;

//   for(int in=0; in<GetMesh()->nnodes(); ++in)
//     {
//       if (GetMesh()->dimension()==2)
// 	{
// 	  Vertex2d v = GetMesh()->vertex2d(in);
// 	  for(int c=0;c<u.ncomp();c++)
// 	    {
// 	      u(in,c) = U(c,v);
// 	    }
// 	}
//       else if (GetMesh()->dimension()==3)
// 	{
// 	  Vertex3d v = GetMesh()->vertex3d(in);
// 	  for(int c=0;c<u.ncomp();c++)
// 	    {
// 	      u(in,c) = U(c,v);
// 	    }
// 	}
//       else assert(0);
//     }
}

/* ----------------------------------------- */

void Q2::InitFilter(nvector<double>& F) const
{
  PressureFilter* PF = static_cast<PressureFilter*>(&F);
  assert(PF);

  if (!PF->Active()) return;

  PF->ReInit(GetMesh()->nnodes(),GetMesh()->nhanging());
  nmatrix<double> T;

  int nv = GetPatchMesh()->nodes_per_patch();
  EntryMatrix  E(nv,1);

  for(int iq=0;iq<GetPatchMesh()->npatches();++iq)
    {
      Transformation(T,iq);
      GetFem()->ReInit(T);

      double cellsize = GetIntegrator()->MassMatrix(E,*GetFem());
      PF->AddDomainPiece(cellsize);

      IntVector ind = *GetPatchMesh()->IndicesOfPatch(iq);
      HN->CondenseHanging(E,ind);

      for(int j=0; j<ind.size(); j++)
	{
	  int jj = ind[j];
	  for(int i=0; i<ind.size(); i++)
	    {
 	      F[jj] += E(i,j,0,0);
	    }
	}
    }
}

/* ----------------------------------------- */

void Q2::HNAverage(GlobalVector& x) const
{
  HN->Average(x);
}

/* ----------------------------------------- */

void Q2::HNDistribute(GlobalVector& x) const
{
  HN->Distribute(x);
}

/* ----------------------------------------- */

void Q2::HNZero(GlobalVector& x) const
{
  HN->Zero(x);
}

/* ----------------------------------------- */

bool Q2::HNZeroCheck(const GlobalVector& x) const
{
  return HN->ZeroCheck(x);
}

/* ----------------------------------------- */

void Q2::Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double d) const
{
  PatchDiscretization::Matrix(A,u,EQ,d);

  HN->MatrixDiag(u.ncomp(),A);
}

/* ----------------------------------------- */

void Q2::MassMatrix(MatrixInterface& A) const
{
  PatchDiscretization::MassMatrix(A);

  HN->MatrixDiag(1,A);  
}

/* ----------------------------------------- */

void Q2::StrongDirichletMatrix(MatrixInterface& A, int col, const vector<int>& comp) const
{
  const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
  const IntVector& bv = *GMP->VertexOnBoundary(col);
  for(int i=0;i<bv.size();i++)
    {
      A.dirichlet(bv[i], comp);
    }  
}

/* ----------------------------------------- */

void Q2::StrongDirichletMatrixOnlyRow(MatrixInterface& A, int col, const vector<int>& comp) const
{
  const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
  const IntVector& bv = *GMP->VertexOnBoundary(col);
  for(int i=0;i<bv.size();i++)
    {
      A.dirichlet_only_row(bv[i], comp);
    }  
}

/* ----------------------------------------- */

void Q2::StrongDirichletVector(GlobalVector& u, const DirichletData& BF, int col, const vector<int>& comp, double d) const
{
  const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
  nvector<double> ff(u.ncomp(),0.);
  const IntVector& bv = *GMP->VertexOnBoundary(col);

  FemData QH;

  GlobalToGlobalData();
  BF.SetParameterData(__QP);

  for(int i=0;i<bv.size();i++)
    {
      int index = bv[i];

      QH.clear();
      GlobalData::const_iterator p=GetDataContainer().GetNodeData().begin();
      for(; p!=GetDataContainer().GetNodeData().end(); p++)
      {
        QH[p->first].resize(p->second->ncomp());
        for(int c=0; c<p->second->ncomp(); c++)
        {
          QH[p->first][c].m() = p->second->operator()(index,c);
        }
      }

      assert(GetMesh()->dimension()==2 || GetMesh()->dimension()==3);

      if (GetMesh()->dimension()==2)
	{
	  const Vertex2d& v = GMP->vertex2d(index);
	  
	  BF(ff,v,col);
	  for(int iii=0;iii<comp.size();iii++)
	    {
	      int c = comp[iii];
	      u(index,c) = d * ff[c];
	    }
	}
      else if (GetMesh()->dimension()==3)
	{
	  const Vertex3d& v = GMP->vertex3d(index);
	  
	  BF(ff,v,col);
	  for(int iii=0;iii<comp.size();iii++)
	    {
	      int c = comp[iii];
	      u(index,c) = d * ff[c];
	    }
	}
    }
}

/* ----------------------------------------- */

void Q2::StrongDirichletVectorZero(GlobalVector& u, int col, const vector<int>& comp) const
{
  const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
  const IntVector& bv = *GMP->VertexOnBoundary(col);
  for(int i=0;i<bv.size();i++)
    {
      int index = bv[i];
      for(int iii=0;iii<comp.size();iii++)
	{
	  u( index,comp[iii] ) = 0.;
	}
    }  
}

/* ----------------------------------------- */

void Q2::StrongPeriodicVector(GlobalVector& u, const PeriodicData& BF, int col, const vector<int>& comp, double d) const
{
  const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
  nvector<double> ff(u.ncomp(),0.);
  const IntVector& bv = *GMP->VertexOnBoundary(col);

  FemData QH;

  GlobalToGlobalData();
  BF.SetParameterData(__QP);

  for(int i=0;i<bv.size();i++)
    {
      int index = bv[i];

      QH.clear();
      GlobalData::const_iterator p=GetDataContainer().GetNodeData().begin();
      for(; p!=GetDataContainer().GetNodeData().end(); p++)
      {
        QH[p->first].resize(p->second->ncomp());
        for(int c=0; c<p->second->ncomp(); c++)
        {
          QH[p->first][c].m() = p->second->operator()(index,c);
        }
      }

      if (GetMesh()->dimension()==2)
	{
	  const Vertex2d& v = GMP->vertex2d(index);
	  
	  BF(ff,v,col);
	  for(int iii=0;iii<comp.size();iii++)
	    {
	      int c = comp[iii];
	      u(index,c) = d * ff[c];
	    }
	}
      else if (GetMesh()->dimension()==3)
	{
	  const Vertex3d& v = GMP->vertex3d(index);
	  
	  BF(ff,v,col);
	  for(int iii=0;iii<comp.size();iii++)
	    {
	      int c = comp[iii];
	      u(index,c) = d * ff[c];
	    }
	}
      else abort();
      
    }
}

/* ----------------------------------------- */

void Q2::InterpolateSolution(GlobalVector& u, const GlobalVector& uold) const
{
  //
  // Das ist einfach von Q1 kopiert !!!!!
  //
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

  for(int iq=0; iq<GetMesh()->ncells(); iq++)
    {
      nmatrix<double> w = GetLocalInterpolationWeights(iq);

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
