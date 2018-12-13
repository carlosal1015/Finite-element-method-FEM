/**
*
* Copyright (C) 2004, 2005, 2009 by the Gascoigne 3D authors
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


#ifndef __fmatrixblock_h
#define __fmatrixblock_h

#include  "entrymatrix.h"
#include  "nodematrix.h"
#include  "matrixentrytype.h"


/*****************************************************/

namespace Gascoigne
{
template<int N>
class FMatrixBlock : public NodeMatrix<N,MatrixEntryType>
{
  typedef typename NodeMatrix<N,MatrixEntryType>::iterator        iterator;
  typedef typename NodeMatrix<N,MatrixEntryType>::const_iterator  const_iterator;

  typedef nvector<double>::iterator        viterator;
  typedef nvector<double>::const_iterator  const_viterator;

public:

  int ncomp() const { return N;}

  inline  void   operator *= (const FMatrixBlock<N>&);
  void   operator *= (double s);

  void   transpose();
  void   transpose(FMatrixBlock<N>& A);
  void   copy_transpose(const FMatrixBlock<N>& A);

  void   zero_row(int);
  void   uno_diag(int);
  MatrixEntryType& diag(int i);
  void   getrow   (std::vector<double>& v, int i);
  void   getcolumn(std::vector<double>& v, int i);
  void   setrow   (std::vector<double>& v, int i);
  void   setcolumn(std::vector<double>& v, int i);

  void   DirichletRow (const std::vector<int>& cv);
  void   DirichletCol (const std::vector<int>& cv);
  void   DirichletDiag(const std::vector<int>& cv);
 


  void   entry     (const nmatrix<double>&);
  void   entry     (int i, int j, const EntryMatrix&, double s=1.);
  void   dual_entry(int i, int j, const EntryMatrix&, double s=1.);
  void   inverse ();
  inline void   vmult   (viterator) const;
  void   mult    (FMatrixBlock<N>&, const FMatrixBlock<N>&) const; 

  void   submult(const FMatrixBlock<N>& B, const FMatrixBlock<N>& C)
  {
    // this -= B*C
    nvector<MatrixEntryType>::iterator p(numfixarray<N*N,MatrixEntryType>::begin());
    for (char i=0; i<N; i++)
      {
	for (char j=0; j<N; j++)
	  {
	    nvector<MatrixEntryType>::const_iterator pC(C.begin()+j);
	    nvector<MatrixEntryType>::const_iterator pB(B.begin()+i*N);
	    nvector<MatrixEntryType>::const_iterator qB(pB+N);
	    //for (int k=0; k<N; k++)
	    for (; pB!=qB; pB++)
	      {
		//value(i,j) -= B(i,k) * C(k,j);
		*p -= *pB * *pC;
		pC += N;
	      }
	    p++;
        }      
      }
  }

  void add(double s, const FMatrixBlock<N>& A)
    {
      for (int i=0; i<N; i++)
	{
	  for (int j=0; j<N; j++)
	    {
	      NodeMatrix<N,MatrixEntryType>::value(i,j) += s*A(i,j);
	    }
	}
    }
  void adddiag(const nvector<double>& s, double l)
    {
      for (int i=0; i<N; i++)
	{
	  NodeMatrix<N,MatrixEntryType>::value(i,i) += s[i]*l;
	}
    }
  void add(double s, const TimePattern& TP);

  void cadd(double s, viterator p, const_viterator q0) const
    {
      const_iterator pm = numfixarray<N*N,MatrixEntryType>::begin();
      const_viterator pend = p+N;
      for ( ; p!=pend; p++)
	{
	  double sum = 0.;
	  const_viterator qend(q0+N);
	  for (const_viterator q=q0; q!=qend; q++)
	    {
	      sum +=  *pm++ * *q;
	    }
	  *p += s*sum;
	}
    }
  void caddtrans(double s, viterator p, const_viterator q0) const
    {
      const_iterator pm = numfixarray<N*N,MatrixEntryType>::begin();

      for (int k=0; k<N; k++)
	{
	  for (int h=0; h<N; h++)
	    {
	      //sum += M(k,h) * (*(q0+h));
	      *p++ += *pm++ * *q0;
	    }
	  p -= N;
	  q0++;
	}
    }
  void subtract(viterator p0, const_viterator q0) const
    {
      const_iterator pm = numfixarray<N*N,MatrixEntryType>::begin();

      for (viterator p(p0); p!=p0+N; p++)
	{
	  for (const_viterator q(q0); q!=q0+N; q++)
	    {
	      *p -= *pm++ * *q;
	    }
	}
    };
  std::ostream& print(std::ostream& s) const;

  // Zugriff auf Inhalt ueber ganzen Vektor, damits auch ohne
  // Struktur geht.
  void vector_get(nvector<MatrixEntryType>& v) const
    {
      v.resize(numfixarray<N*N,MatrixEntryType>::size());
      for (int i=0;i<numfixarray<N*N,MatrixEntryType>::size();++i)
        v[i]=NodeMatrix<N,MatrixEntryType>::operator[](i);
    }
  void vector_set(nvector<MatrixEntryType>& v)
    {
      assert(v.size()==this->size());
      for (int i=0;i<numfixarray<N*N,MatrixEntryType>::size();++i)
        NodeMatrix<N,MatrixEntryType>::operator[](i)=v[i];
    }
  void vector_add(double d, nvector<MatrixEntryType>& v)
    {
      assert(v.size()==N*N);
      for (int i=0;i<NodeMatrix<N,MatrixEntryType>::size();++i)
        NodeMatrix<N,MatrixEntryType>::operator[](i)+=d*v[i];
    }
};

/**********************************************************/

template<int N>
inline void FMatrixBlock<N>::operator *= (const FMatrixBlock<N>& B)
{
  numfixarray<N,double> vhelp;
  for (int i=0; i<N; i++)
    {
      for (int j=0; j<N; j++)
        {
          vhelp[j] = NodeMatrix<N,MatrixEntryType>::value(i,j);
        }
      for (int j=0; j<N; j++)
        {
          NodeMatrix<N,MatrixEntryType>::value(i,j) = 0.;
          for (int k=0; k<N; k++)
            {
              NodeMatrix<N,MatrixEntryType>::value(i,j) += vhelp[k] * B(k,j);
            }
        }      
    }
}

/**********************************************************/

template<int N>
inline void FMatrixBlock<N>::vmult(viterator p) const
{
  numfixarray<N,double> vhelp;
  // copy old entries of vector in vhelp
  //
  typename numfixarray<N,double>::iterator       a = vhelp.begin();
  typename numfixarray<N,double>::const_iterator b = vhelp.end();
  
  for ( ; a!=b; a++)
    {
      *a = *p++;
    }
  p -= N;
  a -= N;

  nvector<MatrixEntryType>::const_iterator q(NodeMatrix<N,MatrixEntryType>::begin());
  for (viterator c=p; c!=p+N; c++)
    {
      *c = 0.;
      for (; a!=b; a++)
	{
	  *c += *q++ * *a;
	} 
      a -= N;
    }
  
}


}

#endif


