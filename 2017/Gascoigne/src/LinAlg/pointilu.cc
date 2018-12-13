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


#include  "pointilu.h"
#include  "compareclass.h"
#include  "pointmatrix.h"
#include  "nodesparsestructureadaptor.h"
#include  "componentsparsestructureadaptor.h"
#include  "giota.h"


using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{
PointIlu::PointIlu(int ncomp, string type) : IluInterface(), SimpleIlu(), _ncomp(ncomp)  
{
  if(type=="node")
    {
      SSAP = new NodeSparseStructureAdaptor(_ncomp);
    }
  else if(type=="component")
    {
      SSAP = new ComponentSparseStructureAdaptor(_ncomp);
    }
  else
    {
      cerr << "PointIlu::PointIlu(): unknown type "<< type<<endl;
      abort();
    }
}

/* ----------------------------------------- */

PointIlu::~PointIlu() 
{
  if(SSAP) {delete SSAP; SSAP=NULL;}
}

/* ----------------------------------------- */

void PointIlu::ReInit(const SparseStructureInterface* S)
{
  SSAP->InitStructure(S);
  SimpleIlu::ReInit(SSAP->n(),SSAP->nentries());
}       

/* ----------------------------------------- */

void PointIlu::ConstructStructure(const IntVector& perm, const MatrixInterface& A)
{
  assert(p.size()==perm.size());
  assert(q.size()==perm.size());

  const ColumnDiagStencil* AS = dynamic_cast<const ColumnDiagStencil*>(A.GetStencil());
  assert(AS);

  ///////////////////////////////////////////////////////
  for(int i=0;i<AS->n();i++)
    {
      assert(AS->rowsize(i)>=0);
    }
  ///////////////////////////////////////////////////////
  int n    = AS->n();
  p = perm;
  for(int i=0;i<n;i++) q[p[i]] = i;

  int zmax = 1;
  for(int i=0;i<n;i++)
    {
      zmax = Gascoigne::max_int(zmax,AS->rowsize(i));
    }
  IntVector ppi(zmax), picol(zmax);

  ST.start(0) = 0;
  for(int i=0;i<n;i++)
    {
      int pi = p[i];
      int ni = AS->rowsize(pi);
      ST.stop(i) = ST.start(i) + ni;

      int count=0;
      for(int pos=AS->start(pi);pos<AS->stop(pi);pos++)
        {
          picol[count++] = q[AS->col(pos)];
        }
      iota(ppi.begin(),ppi.begin()+ni,0);
      sort(ppi.begin(),ppi.begin()+ni,CompareLess<IntVector >(picol));

      for(int ii=0;ii<ni;ii++)
        {
          ST.col(ST.start(i)+ ii) = picol[ppi[ii]];
        }
      for(int pos=ST.start(i);pos<ST.stop(i);pos++)
        {
          if(ST.col(pos)==i)
            {
              ST.diag(i) = pos;
              continue;
            }
        }
    }
} 
/*-------------------------------------------------------------*/

void PointIlu::modify(int c, double s)
{
  for(int i=0;i<ST.n();++i)
    {
      if( (i%_ncomp) == c )
        {
          double sum=0.;
          for(int pos=ST.start(i);pos<ST.stop(i);pos++)
            {
              sum += fabs(value[pos]);
            }
          sum -= fabs(value[ST.diag(i)]);
          value[ST.diag(i)] += s*sum;
        }
    }
}
}
