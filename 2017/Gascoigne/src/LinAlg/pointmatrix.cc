/**
*
* Copyright (C) 2004, 2009, 2011 by the Gascoigne 3D authors
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


#include  "pointmatrix.h"
#include  "simplesparsestructureadaptor.h"
#include  "nodesparsestructureadaptor.h"
#include  "componentsparsestructureadaptor.h"
#include  "giota.h"


using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{
PointMatrix::PointMatrix(int ncomp, string type) : MatrixInterface(), _ncomp(ncomp)  
{
  if(type=="node")
    {
//       SSAP = new SimpleSparseStructureAdaptor;
      SSAP = new NodeSparseStructureAdaptor(_ncomp);
    }
  else if(type=="component")
    {
      SSAP = new ComponentSparseStructureAdaptor(_ncomp);
    }
  else
    {
      cerr << "PointMatrix::PointMatrix(): unknown type "<< type<<endl;
      abort();
    }

}

/* ----------------------------------------- */

PointMatrix::~PointMatrix() 
{
  if(SSAP) {delete SSAP; SSAP=NULL;}
}

/* ----------------------------------------- */

void PointMatrix::ReInit(const SparseStructureInterface* S)
{
  SSAP->InitStructure(S);
  SimpleMatrix::ReInit(SSAP->n(),SSAP->nentries());
  SSAP->FillStencil(ST);
}

/* ----------------------------------------- */

void PointMatrix::vmult(GlobalVector& y, const GlobalVector& x, double d) const
{
  assert(SSAP->GetName()=="Node");
  SimpleMatrix::vmult(y,x,d);
}

/* ----------------------------------------- */

void PointMatrix::vmult_transpose(GlobalVector& y, const GlobalVector& x, double d) const
{
  assert(SSAP->GetName()=="Node");
  SimpleMatrix::vmult_transpose(y,x,d);
}

/*-----------------------------------------*/

void PointMatrix::dirichlet(int inode, const vector<int>& cv)
{
  assert(SSAP);
  SimpleMatrix::dirichlet(SSAP->GetIndicesDirichlet(inode,cv));
}

/*-----------------------------------------*/

void PointMatrix::dirichlet_only_row(int inode, const vector<int>& cv)
{
  assert(SSAP);
  SimpleMatrix::dirichlet_only_row(SSAP->GetIndicesDirichlet(inode,cv));
}

/*-----------------------------------------*/

void PointMatrix::periodic(const map<int,int> &m_PeriodicPairs, const IntVector &iv_Components)
{
  const ColumnStencil S = *(dynamic_cast<const ColumnStencil*>(GetStencil()));

  for (int i = 0; i < iv_Components.size(); i++)
  {
    int comp = iv_Components[i];
    int first, second;

    for (map<int,int>::const_iterator p_pair = m_PeriodicPairs.begin(); p_pair!=m_PeriodicPairs.end(); p_pair++)
    {
      {
        // convert node and component to entry of matrix
        first  = p_pair->first  * _ncomp + comp;
        second = p_pair->second * _ncomp + comp;

        // normalize row "first" and row "second"
        for (int pos = S.start(first);pos<S.stop(first);pos++)
        {
          int j = S.col(pos);
          for (int pos2 = S.start(second);pos2<S.stop(second);pos2++)
          {
            if (S.col(pos2)==j)
            {
              value[pos2] = .5*value[pos2] + .5*value[pos];
              value[pos]  = value[pos2];
              break;
            }
          }

          // modify columns
          if (j != first && j != second)
          {
            for (int pos3 = S.start(j);pos3<S.stop(j);pos3++)
            {
              if (S.col(pos3) == first)
              {
                for (int pos4 = S.start(j);pos4<S.stop(j);pos4++)
                  if (S.col(pos4) == second)
                  {
                    value[pos4] = .5*value[pos4] + .5*value[pos3];
                    value[pos3] = value[pos4];
                    break;
                  }
                break;
              }
            }
          }
        }

        // finish modification of columns
        for (int pos = S.start(first);pos<S.stop(first);pos++)
        {
          if (S.col(pos) == first)
          {
            for (int pos2 = S.start(first);pos2<S.stop(first);pos2++)
              if (S.col(pos2) == second)
              {
                value[pos] += value[pos2];
                value[pos2] = 0.;
                break;
              }
            break;
          }
        }
        for (int pos = S.start(second);pos<S.stop(second);pos++)
        {
          if (S.col(pos) == first)
          {
            for (int pos2 = S.start(second);pos2<S.stop(second);pos2++)
              if (S.col(pos2) == second)
              {
                value[pos2] += value[pos];
                value[pos]   = 0.;
                break;
              }
            break;
          }
        }
      }
    }
  }
}

/*-----------------------------------------*/

void PointMatrix::entry_diag(int i, const nmatrix<double>& M)
{
  IntVector cv(_ncomp); iota(cv.begin(),cv.end(),0);
  dirichlet(i,cv);
}

/*-----------------------------------------*/

void PointMatrix::entry(niiterator start, niiterator stop, const EntryMatrix& M, double s)
{
  int n = stop-start;

  for(int ii=0;ii<n;ii++)
    {
      int i = *(start+ii);
      for(int c=0;c<_ncomp;c++)
        {
          int iglob = SSAP->index(i,c);
          for(int jj=0;jj<n;jj++)
            {
              int j = *(start+jj);
              for(int d=0;d<_ncomp;d++)
                {
                  int jglob = SSAP->index(j,d);
                  int pos = ST.Find(iglob,jglob);
                  value[pos] += s*M(ii,jj,c,d);
                }
            }
        }
    }
}

/*-----------------------------------------*/

void PointMatrix::AddMassWithDifferentStencil(const MatrixInterface* MP, const TimePattern& TP, double s)
{
  const SimpleMatrix* SM = dynamic_cast<const SimpleMatrix*>(MP);
  assert(SM);

  int n = SSAP->nnodes();

  const ColumnStencil*  SMS = dynamic_cast<const ColumnStencil*>(SM->GetStencil());
  const NodeSparseStructureAdaptor* NSMS = dynamic_cast<const NodeSparseStructureAdaptor*>(SSAP);
  assert(NSMS);

  assert(n==SMS->n());

  for(int i=0;i<n;i++)
    {
      for(int pos=SMS->start(i);pos<SMS->stop(i);pos++)
        {
          int j = SMS->col(pos);
          double m = s * SM->GetValue(pos);

          for(int c=0;c<_ncomp;c++)
            {
              int isystem = NSMS->index(i, c);
              for(int d=0;d<_ncomp;d++)
                {
                  int jsystem = NSMS->index(j, d);
                  bool found=0;
                  for(int pos2=ST.start(isystem);pos2<ST.stop(isystem);pos2++)
                    {
                      if(ST.col(pos2)==jsystem)
                        {
                          found = 1;
                          value[pos2] += m * TP(c,d);
                        }
                    }
                  if(!found) cerr << "not found ";
                }
            }
        }
    }


//   for(int i=0;i<n;i++)
//     {
//       for(int c=0;c<_ncomp;c++)
// 	{
// 	  intw iglob = SSAP->index(i,c);
// 	  for(int j=0;j<n;j++)
// 	    {
// 	      for(int d=0;d<_ncomp;d++)
// 		{
// 		  int jglob = SSAP->index(j,d);
// 		  int pos = ST.Find(iglob,jglob);
// 		  value[pos] += s * TP(c,d) * SM->GetValue(i,j);
// 		}
// 	    }
// 	}
//     }
}


/*-----------------------------------------*/

void PointMatrix::AddMassWithDifferentStencilJacobi(const MatrixInterface* MP, const TimePattern& TP, double s)
{
  const SimpleMatrix* SM = dynamic_cast<const SimpleMatrix*>(MP);
  assert(SM);

  int n = SSAP->nnodes();

  const ColumnStencil*  SMS = dynamic_cast<const ColumnStencil*>(SM->GetStencil());
  const NodeSparseStructureAdaptor* NSMS = dynamic_cast<const NodeSparseStructureAdaptor*>(SSAP);
  assert(NSMS);

  assert(n==SMS->n());

  vector<double> diag(n);
  for(int i=0; i<n; i++)
  {
    diag[i] = SM->GetValue(i,i);
  }

  for(int i=0;i<n;i++)
    {
      for(int pos=SMS->start(i);pos<SMS->stop(i);pos++)
        {
          int j = SMS->col(pos);
          double m = s * SM->GetValue(pos);

          m /= sqrt(diag[i] * diag[j]);

          for(int c=0;c<_ncomp;c++)
            {
              int isystem = NSMS->index(i, c);
              for(int d=0;d<_ncomp;d++)
                {
                  int jsystem = NSMS->index(j, d);
                  bool found=0;
                  for(int pos2=ST.start(isystem);pos2<ST.stop(isystem);pos2++)
                    {
                      if(ST.col(pos2)==jsystem)
                        {
                          found = 1;
                          value[pos2] += m * TP(c,d);
                        }
                    }
                  if(!found) cerr << "not found ";
                }
            }
        }
    }
}

/*-----------------------------------------*/

void PointMatrix::RestrictMatrix(const MgInterpolatorMatrix& I, const PointMatrix& Ah)
{
  std::cerr << "\"PointMatrix::RestrictMatrix\" not written!" << std::endl;
  abort();

//   const UnstructuredStencil& USH = GetStencil();
//   const UnstructuredStencil& USh = Ah.GetStencil();

//   const UnstructuredStencil& us = I.GetStencil();

//   for(int k=0;k<USh.n();k++)
//     {
//       for(int pos2=us.start(k);pos2<us.stop(k);pos2++)
// 	{	      
// 	  int i = us.col(pos2);
// 	  double alpha_ik = I.Alpha(pos2);

// 	  for(int pos=USh.start(k);pos<USh.stop(k);pos++)
// 	    {
// 	      int l = USh.col(pos);
	      
// 	      for(int pos3=us.start(l);pos3<us.stop(l);pos3++)
// 		{
// 		  int j = us.col(pos3);
// 		  double alpha_jl = I.Alpha(pos3);
// 		  int pos4 = USH.Find(i,j);
// 		  double d = alpha_ik*alpha_jl;

// 		  mat(pos4)->add(d, *Ah.mat(pos));
// 		}
// 	    }
// 	}
//     }
}
}
