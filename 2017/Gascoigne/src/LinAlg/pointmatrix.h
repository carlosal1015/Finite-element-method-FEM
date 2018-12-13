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


#ifndef  __PointMatrix_h
#define  __PointMatrix_h

#include  "matrixinterface.h"
#include  "simplematrix.h"
#include  "sparsestructureadaptor.h"
#include  "mginterpolatormatrix.h"

namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  ... comments PointMatrix

///
///
/////////////////////////////////////////////

class PointMatrix : public SimpleMatrix, virtual public MatrixInterface
{
protected:

  int _ncomp;
  SparseStructureAdaptor* SSAP;

public:

//
///  Constructor 
//
    PointMatrix(int ncomp, std::string type);
    virtual ~PointMatrix();

    std::string GetName() const {return "PointMatrix";}

    void zero() {
      SimpleMatrix::zero();
    }
    void vmult(GlobalVector& y, const GlobalVector& x, double d=1.) const;
    void vmult_transpose(GlobalVector& y, const GlobalVector& x, double d=1.) const;

    const StencilInterface* GetStencil() const { return SimpleMatrix::GetStencil();}
    void ReInit(const SparseStructureInterface* S);

    void entry(niiterator start, niiterator stop, const EntryMatrix& M, double s=1.);
    void entry_diag(int i, const nmatrix<double>& M);
    void dirichlet (int i, const std::vector<int>& cv);
    void dirichlet_only_row (int i, const std::vector<int>& cv);
    void periodic (const std::map<int,int> &m_PeriodicPairs, const IntVector &iv_Components);

    void transpose() {
      SimpleMatrix::transpose();
    }

    void AddMassWithDifferentStencil(const MatrixInterface* M, const TimePattern& TP, double s=1.);
    void AddMassWithDifferentStencilJacobi(const MatrixInterface* M, const TimePattern& TP, double s=1.);

    void RestrictMatrix(const MgInterpolatorMatrix& I, const PointMatrix& Ah);
};
}

#endif
