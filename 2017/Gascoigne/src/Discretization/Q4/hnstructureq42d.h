/**
*
* Copyright (C) 2006 by the Gascoigne 3D authors
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


#ifndef __HNStructureQ42d_h
#define __HNStructureQ42d_h

#include "hnstructureq22d.h"

namespace Gascoigne
{

/**********************************************************/

  class HNStructureQ42d : public HNStructureQ22d
  {
    protected:
      nmatrix<double>                       M,Mq2;
      nvector<fixarray<5,double> >          w,wq2;
      const std::map<int,fixarray<6,int> > *q4edges;

      typedef std::map<int,fixarray<6,int> >::iterator       iteratorq4;
      typedef std::map<int,fixarray<6,int> >::const_iterator const_iteratorq4;

      void add_column(EntryMatrix& A, const EntryMatrix&B, int j1, int j2, double s=1.) const;
      void add_row(EntryMatrix& A, const EntryMatrix&B, int i1, int i2, double s=1.) const;
      void GetHangingIndices(std::vector<int>& hang, const IntVector& indices) const;
      fixarray<5,int> local_nodes(int e,int n) const;
      void modify_column_higher(EntryMatrix& E, const std::vector<int>& hang, const IntVector& indices) const;
      void modify_column_lower(EntryMatrix& E, const std::vector<int>& hang, const IntVector& indices) const;
      void modify_row_higher(EntryMatrix& E, const std::vector<int>& hang, const IntVector& indices) const;
      void modify_row_lower(EntryMatrix& E, const std::vector<int>& hang, const IntVector& indices) const;
      const fixarray<6,int>& regular_nodes(int i) const;

    public:
      HNStructureQ42d();
      ~HNStructureQ42d() { }

      void ReInit(const MeshInterface* m);

      void CondenseHanging(IntVector& indices) const;
      void CondenseHanging(EntryMatrix& E, IntVector& indices) const;
      void CondenseHangingLowerHigher(EntryMatrix& E, IntVector& indices) const;
      void CondenseHangingHigherLower(EntryMatrix& E, IntVector& indices) const;
      void MatrixDiag(int ncomp, MatrixInterface& A) const;
      void SparseStructureDiag(SparseStructure& S) const;

      void Zero(GlobalVector& u) const;
      void Average(GlobalVector& u) const;
      void Distribute(GlobalVector& u) const;

      bool ZeroCheck(const GlobalVector& u) const;
  };

/**********************************************************/

}

#endif
