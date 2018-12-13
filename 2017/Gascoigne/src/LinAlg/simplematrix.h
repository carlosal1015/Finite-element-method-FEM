/**
*
* Copyright (C) 2004, 2005, 2008, 2011 by the Gascoigne 3D authors
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


#ifndef  __SimpleMatrix_h
#define  __SimpleMatrix_h

#include  "matrixinterface.h"
#include  "columndiagstencil.h"
#include  "sparsestructureadaptor.h"
#include  "compvector.h"


namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  ... comments SimpleMatrix

///
///
/////////////////////////////////////////////

class SimpleMatrix : virtual public MatrixInterface
{
protected:

  ColumnDiagStencil  ST;
  DoubleVector value; 
  DoubleVector _diag;

public:

//
///  Constructor 
//
    SimpleMatrix() : MatrixInterface() {}
    ~SimpleMatrix() {}

    std::string GetName() const {return "SimpleMatrix";}

    std::ostream& Write(std::ostream& os) const;

    const StencilInterface* GetStencil() const { return &ST;}
    double& GetValue(int pos) {return value[pos];}
    const double& GetValue(int pos) const {return value[pos];}
    const double& GetValue(int i, int j) const {return value[ST.Find(i,j)];}
    const DoubleVector& GetValues() const { return value; }

    void zero() {value.zero();}
    void ReInit(const SparseStructureInterface* S);
    void ReInit(int n, int nentries);
    void entry(niiterator start, niiterator stop, const EntryMatrix& M, double s=1.);
    void vmult(DoubleVector& y, const DoubleVector& x, double d=1.) const;
    void vmult_transpose(DoubleVector& y, const DoubleVector& x, double d=1.) const;
    void vmult_comp(int c, int d, GlobalVector& y, const GlobalVector& x, double s=1.) const;
    void vmult_comp_trans(int c, int d, GlobalVector& y, const GlobalVector& x, double s=1.) const;
    void vmult_time(GlobalVector& y, const GlobalVector& x, const TimePattern& TP, double s=1.) const;
    void dirichlet(const IntVector& indices);
    void dirichlet_only_row(const IntVector& indices);

    void transpose();
    void entry_diag(int i, const nmatrix<double>& M);
    
    void PrepareJacobi(double s=1.);
    void JacobiVector(GlobalVector &y) const;
    void Jacobi(GlobalVector &y) const;
    void vmult_time_Jacobi(GlobalVector& y, const GlobalVector& x, const TimePattern& TP, double s=1.) const;
    void copy_entries(const MatrixInterface&  A);

};
}

#endif
