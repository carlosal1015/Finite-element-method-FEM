/**
*
* Copyright (C) 2004, 2005, 2006, 2008, 2009, 2011 by the Gascoigne 3D authors
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


#ifndef __matrixinterface_h
#define __matrixinterface_h

#include  "entrymatrix.h"
#include  "compvector.h"
#include  "sparsestructureinterface.h"
#include  "stencilinterface.h"
#include  "periodicdata.h"
#include  <string>
#include  "matrixentrytype.h"

/*-------------------------------------------------------------*/

namespace Gascoigne
{
  class MatrixInterface
  {
    private:

    protected:

    public:
      MatrixInterface() { }
      virtual ~MatrixInterface() { }

      virtual std::string GetName() const=0;
  
      virtual const StencilInterface* GetStencil() const=0;
      virtual void ReInit(const SparseStructureInterface* S)=0;
      virtual void AddMassWithDifferentStencil(const MatrixInterface* M, const TimePattern& TP, double s=1.) {
        std::cerr << "\"MatrixInterface::AddMassWithDifferentStencil\" not written!" << std::endl;
        abort();
      }
      virtual void AddMassWithDifferentStencilJacobi(const MatrixInterface* M, const TimePattern& TP, double s=1.) {
        std::cerr << "\"MatrixInterface::AddMassWithDifferentStencil\" not written!" << std::endl;
        abort();
      }
      virtual void copy_entries(const MatrixInterface& S) {
        std::cerr << "\"MatrixInterface::copy_entries\" not written!" << std::endl;
        abort();
      }
      virtual void zero()=0;
      virtual void transpose()=0;

      virtual std::ostream& Write(std::ostream& os) const=0;

      //
      /// for matrix assembling
      //
      typedef IntVector::const_iterator niiterator;

      virtual void entry(nvector<int>::const_iterator start1, nvector<int>::const_iterator stop1,
		 nvector<int>::const_iterator start2, nvector<int>::const_iterator stop2,
		 const EntryMatrix& M, double s=1.) { 
        std::cerr << "\"MatrixInterface::entry\" not written!" << std::endl;
        abort();
      }
      virtual void entry(niiterator start, niiterator stop, const EntryMatrix& M, double s=1.)=0;
      virtual void entrydual(niiterator start, niiterator stop, const EntryMatrix& M, double s=1.) {
        std::cerr << "\"MatrixInterface::entrydual\" not written!" << std::endl;
        abort();
      }

      //
      /// for hanging nodes
      //
      virtual void entry_diag(int i, const nmatrix<double>& M)=0;
      
      //
      /// for boundary conditions
      //
      virtual void scale_diag(int i, const std::vector<int>& cv, double s) {
        std::cerr << "\"MatrixInterface::scale_diag\" not written!" << std::endl;
        abort();
      }
      virtual void dirichlet (int i, const std::vector<int>& cv) {
        std::cerr << "\"MatrixInterface::dirichlet\" not written!" << std::endl;
        abort();
      }
      virtual void dirichlet_only_row(int i, const std::vector<int>& indices) {
        std::cerr << "\"MatrixInterface::dirichlet_only_row\" not written!" << std::endl;
        abort();
      }
      virtual void dirichlet_only_column(int i, const std::vector<int>& indices) {
        std::cerr << "\"MatrixInterface::dirichlet_only_column\" not written!" << std::endl;
        abort();
      }
      virtual void dirichlet_only_row_no_diag(int i, const std::vector<int>& indices) {
        std::cerr << "\"MatrixInterface::dirichlet_only_row_no_diag\" not written!" << std::endl;
        abort();
      }
      virtual void periodic (const std::map<int,int> &m_PeriodicPairs, const IntVector &iv_Components) {
        std::cerr << "\"MatrixInterface::periodic\" not written!" << std::endl;
        abort();
      }
      virtual void vmult(GlobalVector& y, const GlobalVector& x, double s=1.) const {
        std::cerr << "\"MatrixInterface::vmult\" not written!" << std::endl;
        abort();
      }
      virtual void vmult_transpose(GlobalVector& y, const GlobalVector& x, double s=1.) const {
        std::cerr << "\"MatrixInterface::vmult_tranpose\" not written!" << std::endl;
        abort();
      }
      virtual void vmult_time(GlobalVector& y, const GlobalVector& x, const TimePattern& TP, double s=1.) const {
        std::cerr << "\"MatrixInterface::vmult_time\" not written!" << std::endl;
        abort();
      }

		/*-----------------------------------------------*/
      
      virtual void FillInterfaceList(const nvector<int>& elements,nvector<int>& start, nvector<MatrixEntryType>& values) const
	{ std::cerr << "\"MatrixInterface::FillInterfaceList\" not written!" << std::endl; abort(); }
      virtual void FurbishInterface (double d, const nvector<int>&   elements, const nvector<int>&   start, const nvector<MatrixEntryType>& values)
	{ std::cerr << "\"MatrixInterface::FurbishInterface\" not written!" << std::endl; abort(); }
      
      virtual void PrepareJacobi(double s)
	{ std::cerr << "\"MatrixInterface::PrepareJacobi\" not written!" << std::endl; abort(); }
      virtual void Jacobi           (GlobalVector& x) const
	{ std::cerr << "\"MatrixInterface::Jacobi\" not written!" << std::endl; abort(); }
      
		/*-----------------------------------------------*/
 };
}

/*-------------------------------------------------------------*/

#endif
