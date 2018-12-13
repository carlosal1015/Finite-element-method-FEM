/**
*
* Copyright (C) 2008, 2011 by the Gascoigne 3D authors
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


#ifndef __dynamicblockmatrix_h
#define __dynamicblockmatrix_h


#include  "dynamicstencil.h"
#include  "sparsestructure.h"
#include  "matrixinterface.h"
#include  "gascoigne.h"
#include  <list>

using namespace std;

/*-------------------------------------------------------------*/

namespace Gascoigne
{

    template<class B>
	class DynamicBlockMatrix : public MatrixInterface
    {
    protected:

	typedef typename list<int>::const_iterator    const_citerator;
	typedef typename list<int>::iterator          citerator;
	typedef typename list<B>::const_iterator      const_viterator;
	typedef typename list<B>::iterator            viterator;

	// it would be more easy to read if the stencil would be
	// embedded. However for ilu-sorting we need to access the stencil
	// form the solver... 
	// but this is also not nice, because the important functions 
	// for accessing the structure are not in the interface...
	DynamicStencil  DS;
    
	// the entries of the matrix, one list for every row!
	vector<list<B> >   smat;
	// number of components
	int              nc;
    
	void matrix_vector_trans(int p, double* yp, const double* xp, double s=1.) const;
  
	// number of entries
	int size() const {
          std::cerr << "\"DynamicBlockMatrix::size\" not written!" << std::endl;
          abort();
        }

    public:

	DynamicBlockMatrix<B>();
	DynamicBlockMatrix<B>(const DynamicBlockMatrix<B>& A);
	virtual ~DynamicBlockMatrix<B>() {}
  
	string GetName() const {return "DynamicBlockMatrix";}

	/////// Zugriff //////////////////////

	int   n()          const { assert(smat.size()==DS.n()); return DS.n(); }
	int   nentries()   const {
          std::cerr << "\"DynamicBlockMatrix::nentries\" not written!" << std::endl;
          abort();
        }
	int   ntotal()     const {
          std::cerr << "\"DynamicBlockMatrix::ntotal\" not written!" << std::endl;
          abort();
        }

	int  rowsize(int i)     const 
	{
	    assert(i<smat.size()); 
	    assert(i<DS.cols.size()); 
	    assert(smat[i].size()==DS.cols[i].size());
	    return smat[i].size(); 
	}


	//////////////////////////////
	//////////////////////////////
	//////////////////////////////
	// Working on the Stencil
	viterator       vfind(int i,int j)
	{
	    const_citerator cit = DS.cstart(i);
	    viterator       vit = vstart(i);
	    for (;(cit!=DS.cstop(i))&&(*cit!=j);++cit,++vit){}
	    assert(cit!=DS.cstop(i));
	    return vit;
	}
	const_viterator vfind(int i,int j) const
	{
	    const_citerator cit = DS.cstart(i);
	    const_viterator vit = smat[i].begin();
	    for (;(cit!=DS.cstop(i))&&(*cit!=j);++cit,++vit){}
	    assert(cit!=DS.cstop(i));
	    return vit;
	}
	citerator       cfind(int i,int j)
	{
	    citerator cit = DS.cstart(i);
	    for (;(cit!=DS.cstop(i))&&(*cit!=j);++cit){}
	    assert(cit!=DS.cstop(i));
	    return cit;
	}
	const_citerator cfind(int i,int j) const
	{
	    const_citerator cit = DS.cstart(i);
	    for (;(cit!=DS.cstop(i))&&(*cit!=j);++cit){}
	    assert(cit!=DS.cstop(i));
	    return cit;
	}
	viterator       vdiag(int i)        { return vfind(i,i); }
	const_viterator vdiag(int i) const  { return vfind(i,i); }
	citerator       cdiag(int i)        { return cfind(i,i); }
	const_citerator cdiag(int i) const  { return cfind(i,i); }
	const_citerator cstart(int i) const { return DS.cstart(i); }
	const_citerator cstop (int i) const { return DS.cstop(i);  }
	citerator       cstart(int i)       { return DS.cstart(i); }
	citerator       cstop (int i)       { return DS.cstop(i);  }
	const_viterator vstart(int i) const { assert(i<n()); return smat[i].begin(); }
	const_viterator vstop (int i) const { assert(i<n()); return smat[i].end(); }
	viterator       vstart(int i)       { assert(i<n()); return smat[i].begin(); }
	viterator       vstop (int i)       { assert(i<n()); return smat[i].end(); }

	// adds a coupling to the matrix, creates the entry in cols and smat (zero),
	// returns the iterator to the new value in smat
	viterator       add_coupling(int i,int j) 
	{
	    citerator cit = DS.cstart(i);
	    viterator vit = smat[i].begin();
	    for (;cit!=DS.cstop(i);++cit,++vit)
		if (*cit>=j)  break;
	    assert(*cit!=j);
	    DS.cols[i].insert(cit,j);
	    smat[i].insert(vit,B());
	    --vit;
	    return vit;
	}

	const StencilInterface* GetStencil() const { return &DS;}

	///// Methods //////////////////////

	//void copy_entries(const MatrixInterface& S);

	void AddMassWithDifferentStencil(const MatrixInterface* M, 
					 const TimePattern& TP, double s=1.);

	DynamicBlockMatrix& operator=(const DynamicBlockMatrix<B>& S); 
	void transpose();
	void ReInit   (const SparseStructureInterface*);
	void dirichlet(int i, const vector<int>& cv);
	void dirichlet_only_row(int i, const vector<int>& cv);

	void zero();
	void entry_diag(int i, const nmatrix<double>& M);
	void entry(nvector<int>::const_iterator start1, nvector<int>::const_iterator stop1,
		   nvector<int>::const_iterator start2, nvector<int>::const_iterator stop2,
		   const EntryMatrix& M, double s=1.){
          std::cerr << "\"DynamicBlockMatrix::entry\" not written!" << std::endl;
          abort();
        }
	void entry(nvector<int>::const_iterator start, nvector<int>::const_iterator stop, const EntryMatrix& M, double s=1.);
	void entrydual(nvector<int>::const_iterator start, nvector<int>::const_iterator stop, const EntryMatrix& M, double s=1.){
          std::cerr << "\"DynamicBlockMatrix::entrydual\" not written!" << std::endl;
          abort();
        }

	void vmult(GlobalVector& y, const GlobalVector& x, double s=1.) const;
	void vmult(GlobalVector& y, const GlobalVector& x, const TimePattern& TP, double s=1.)const;

	/*-----------------------------------------------*/
  
	void Jacobi(GlobalVector& x) const;
      
	/*-----------------------------------------------*/

	void FillInterfaceList(const nvector<int>& elements,nvector<int>& start, nvector<MatrixEntryType>& values) const
	{ 
          std::cerr << "\"DynamicBlockMatrix::FillInterfaceList\" not written!" << std::endl;
          abort();
        }
	void FurbishInterface (double d, const nvector<int>&   elements, const nvector<int>&   start, const nvector<MatrixEntryType>& values)
	{ 
          std::cerr << "\"DynamicBlockMatrix::FurbishInterface\" not written!" << std::endl;
          abort();
        }

	/*-----------------------------------------------*/

	ostream& Write(ostream &s) const
	{ 
          std::cerr << "\"DynamicBlockMatrix::Write\" not written!" << std::endl;
          abort();
        }
	friend   ostream& operator<<(ostream &s, const DynamicBlockMatrix<B>& A) 
	{ 
          std::cerr << "\"ostream& operator<<(ostream &s, const DynamicBlockMatrix<B>& A)\" not written!" << std::endl;
          abort();
        }
    };
}

#endif
