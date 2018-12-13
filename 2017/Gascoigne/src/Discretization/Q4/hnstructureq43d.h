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


#ifndef __HNStructureQ43d_h
#define __HNStructureQ43d_h

#include "hnstructureq23d.h"

namespace Gascoigne
{

/**********************************************************/

  class HNStructureQ43d : public HNStructureQ23d
  {
    protected:
      typedef fixarray<6,int>  EdgeVector;
      typedef fixarray<26,int> FaceVector;

      typedef std::map<int,EdgeVector>::iterator       iteratorq4;
      typedef std::map<int,EdgeVector>::const_iterator const_iteratorq4;
      typedef std::map<int,FaceVector>::iterator       face_iteratorq4;
      typedef std::map<int,FaceVector>::const_iterator face_const_iteratorq4;

      nmatrix<double>                 Medge,Mface,Mq2edge,Mq2face;
      nvector<fixarray<5,double> >    w,wq2;
      const std::map<int,EdgeVector> *q4edges;
      const std::map<int,FaceVector> *q4faces;

      void add_column(EntryMatrix& A, const EntryMatrix&B, int j1, int j2, double s=1.) const;
      void add_row(EntryMatrix& A, const EntryMatrix&B, int i1, int i2, double s=1.) const;
      void GetHangingIndices(std::vector<int>& hang_e, std::vector<int>& hang_f, const IntVector& indices) const;
      int hanging(int i) const;
      fixarray<5,int> local_nodes_on_edge(int e, const IntVector& indices) const;
      fixarray<25,int> local_nodes_on_face(int e, const IntVector& indices) const;
      void modify_column_higher(EntryMatrix& E, const std::vector<int>& hang_e, const std::vector<int>& hang_f, const IntVector& indices) const;
      void modify_column_lower(EntryMatrix& E, const std::vector<int>& hang_e, const std::vector<int>& hang_f, const IntVector& indices) const;
      void modify_row_higher(EntryMatrix& E, const std::vector<int>& hang_e, const std::vector<int>& hang_f, const IntVector& indices) const;
      void modify_row_lower(EntryMatrix& E, const std::vector<int>& hang_e, const std::vector<int>& hang_f, const IntVector& indices) const;
      const EdgeVector& regular_nodes_on_edge(int i) const;
      const FaceVector& regular_nodes_on_face(int i) const;

    public:
      HNStructureQ43d();
      ~HNStructureQ43d() { }

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
