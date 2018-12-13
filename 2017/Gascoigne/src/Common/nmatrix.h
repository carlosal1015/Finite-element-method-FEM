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


#ifndef __nmatrix_h
#define __nmatrix_h

#include  "nvector.h"
#include  "giota.h"

/**************************************************/

namespace Gascoigne
{
template<class T>
class nmatrix : public nvector<T>
{

  int nn,mm;
  
public:

  typedef typename  nvector<T>::const_iterator  const_iterator;
  typedef typename  nvector<T>::iterator        iterator;

  /*
    Constructeurs
   */

  nmatrix<T>()             : nvector<T>()   , nn(0), mm(0) {}
  nmatrix<T>(size_t n)        : nvector<T>(n*n), nn(n), mm(n) {}       
  nmatrix<T>(size_t n, size_t m) : nvector<T>(n*m), nn(n), mm(m) {}       
  nmatrix<T>(size_t n, size_t m, const T& d) : nvector<T>(n*m,d), nn(n), mm(m) {}       


  /*
    Fonctions memoire
   */
  //void reservesize(size_t n, size_t m, T s = static_cast<T>(0))
  void reservesize(size_t n, size_t m, T s)
    {
		nn = n;
		mm = m;
		nvector<T>::reservesize(n*m,s);
    }
  void reservesize(const nmatrix<T>& M)
    {
		nmatrix<T>::reservesize(M.n(),M.m(),0);
    }
  void memory(size_t n, size_t m)
    {
      //nmatrix<T>::reservesize(n,m);
      nmatrix<T>::reservesize(n,m,0);
    }
  void resize(size_t n)
    {
      nmatrix<T>::memory(n,n);
   }
   void resize(size_t n, size_t m)
    {
      nmatrix<T>::memory(n,m);
    }

  /*
    Fonctions d'acces
   */

  size_t      n()                     const  { return nn; }
  size_t      m()                     const  { return mm; }
  const T& operator()(int i,int j) const  { return (*this)[j+mm*i];   }
  T&       operator()(int i,int j)        { return (*this)[j+mm*i];   }
  const T& value     (int i,int j) const  { return (*this)[j+mm*i];   }
  T&       value     (int i,int j)        { return (*this)[j+mm*i];   }

  nmatrix<T>&   operator=(const nmatrix<T>& A)
    {
      if( (n()!=A.n()) || (m()!=A.m()) )
	{
	  std::cerr << "nmatrix<T>::operator=(const nmatrix<T>& A)\n";
	  std::cerr << "no = possible " << n()<<" "<<A.n()<<" "<<m()<<" "<<A.m()<<std::endl; abort();
	}
      for(int i=0;i<n();i++)
	{
	  for(int j=0;j<m();j++)
	    {
	      (*this)(i,j) = A(i,j);
	    }
	}
      return *this;
   }

  const_iterator rowstart(int i) const {return std::vector<T>::begin()+mm*i;}
  
  void identity()
    {
    nvector<T>::zero();
      for(int i=0;i<n();i++)
	{
	  (*this)(i,i) = 1.;
	}
    }

  /*
    Fonctions numeriques
   */
  T det() const
    {
      const nmatrix<T>& A = *this;
      T d;
      if( (n()==2) && (m()==2))
	{
	  d = A(0,0)*A(1,1)-A(0,1)*A(1,0);
	}
      else if( (n()==3) && (m()==3))
	{
	  d = A(0,0) * ( A(1,1)*A(2,2)-A(1,2)*A(2,1) ) 
	    - A(1,0) * ( A(0,1)*A(2,2)-A(0,2)*A(2,1) ) 
	    + A(2,0) * ( A(0,1)*A(1,2)-A(0,2)*A(1,1) );
	}
      else
	{
	  std::cerr << " cul de sac dans nmatrix::det() "<<n()<<" "<<m()<<std::endl;
	  exit(1);
	}
      return d;
    }

  void gram(const nmatrix<T>& B)
    {
      nvector<T>::zero();
      for(int i=0;i<n();i++)
	{
	  for(int j=0;j<m();j++)
	    {
	      for(int k=0;k<B.n();k++)
		{
		  (*this)(i,j) += B(k,i) * B(k,j);
		} 
	    }
	}
    }
  void transpose()
    {
      nmatrix<T> B(*this);
      resize(m(),n());
      for(int i=0;i<n();i++)
	{
	  for(int j=0;j<m();j++)
	    {
	      (*this)(i,j) = B(j,i);
	    }
	}
    }

  void mmult(nmatrix<T>& A, const nmatrix<T>& B) const
    {
      assert(this->m()==B.n());
      assert(A.n() == this->n());
      assert(A.m() == B.m());
      A.zero();
      for(int i=0;i<A.n();i++)
	{
	  for(int j=0;j<A.m();j++)
	    {
	      for(int k=0;k<m();k++)
		{
		  A(i,j) += (*this)(i,k) * B(k,j);
		} 
	    }
	}
    }

  void mmult_ad(nmatrix<T>& A, const nmatrix<T>& B) const
    {
      assert(this->n()==B.n());
      assert(A.n() == this->m());
      assert(A.m() == B.m());
      /* A = (*this)^T * B */
      A.zero();
      for(int i=0;i<A.n();i++)
	{
	  for(int j=0;j<A.m();j++)
	    {
	      for(int k=0;k<n();k++)
		{
		  A(i,j) += (*this)(k,i) * B(k,j);
		} 
	    }
	}
    }

  template<class VECTOR>
    void Mult(VECTOR& y, const VECTOR& x, double s=1.) const
    {
      const_iterator                    p  = std::vector<T>::begin();
      typename VECTOR::iterator         py = y.begin();
      typename VECTOR::const_iterator   px;
      
      while(p!=std::vector<T>::end())
	{
	  px = x.begin();
	  for(int j=0;j<m();j++)
	    {
	      *py += s * (*p++) * (*px++);
	    }
	  py++;
	}
    }
  template<class VECTOR>
    void mult(VECTOR& y, const VECTOR& x) const
    {
      const_iterator                    p  = std::vector<T>::begin();
      typename VECTOR::iterator         py = y.begin();
      typename VECTOR::const_iterator   px;
      
      while(p!=std::vector<T>::end())
	{
	  px = x.begin();
	  for(int j=0;j<m();j++)
	    {
	      *py += (*p++) * (*px++);
	    }
	  py++;
	}
    }
  template<class VECTOR>
    void multtrans(VECTOR& y, const VECTOR& x, double s=1.) const
    {
      assert(y.size()==this->n());
      const_iterator                              p  = std::vector<T>::begin();
      typename VECTOR::iterator         py = y.begin();
      typename VECTOR::const_iterator   px = x.begin();
      
      while(p!=std::vector<T>::end())
	{
	  py = y.begin();
	  for(int j=0;j<n();j++)
	    {
	      (*py) += s* (*p++) * (*px);
	      py++;
	    }
	  px++;
	}
    }

  template<class VECTOR>
    void multeq(VECTOR& y, const VECTOR& x, double s=1.) const
    {
      const_iterator                              p  = std::vector<T>::begin();
      typename VECTOR::iterator         py = y.begin();
      typename VECTOR::const_iterator   px;
      
      while(p!=std::vector<T>::end())
	{
	  px = x.begin();
	  *py = 0.;
	  for(int j=0;j<m();j++)
	    {
	      *py += s * (*p++) * (*px++);
	    }
	  py++;
	}
    }

  template<class VECTOR>
    void mult_ad(VECTOR& y, const VECTOR& x) const
    {
      // vmult with the adjoint matrix
      const_iterator                              p  = std::vector<T>::begin();
      typename VECTOR::iterator         py = y.begin();
      typename VECTOR::const_iterator   px = x.begin();
      
      while(p!=std::vector<T>::end())
	{
	  py = y.begin();
	  for(int j=0;j<n();j++)
	    {
	      (*py) += (*p++) * (*px);
	      py++;
	    }
	  px++;
	}
    }

  template<class ITER1, class ITER2>
    void mult_ad(ITER1 py, ITER2 px, double s=1.) const
    {
      // vmult with the adjoint matrix
      const_iterator    p  = std::vector<T>::begin();
      
      while(p!=std::vector<T>::end())
	{
	  for(int j=0;j<n();j++)
	    {
	      (*py) += s*(*p++) * (*px);
	      py++;
	    }
	  py -= m();
	  px++;
	}
    }

  template<class ITER1, class ITER2>
  void mult(ITER1 py, ITER2  px0, double s=1.) const
  {
    const_iterator    p  = std::vector<T>::begin();
    ITER2   px;
    
    while(p!=std::vector<T>::end())
      {
	px = px0;
	for(int j=0;j<m();j++)
	  {
	    *py += s* (*p++) * (*px++);
	  }
	py++;
      }
  }

  template<class VECTOR>
    void multeq_ad(VECTOR& y, const VECTOR& x) const
    {
      // vmulteq with the adjoint matrix
      const_iterator                              p  = std::vector<T>::begin();
      typename VECTOR::iterator         py = y.begin();
      typename VECTOR::const_iterator   px = x.begin();
      
      y.zero();
      while(p!=std::vector<T>::end())
	{
	  for(int j=0;j<n();j++)
	    {
	      (*py) += (*p++) * (*px);
	      py++;
	    }
	  py -= m();
	  px++;
	}
    }

  template<class ITER1, class ITER2>
    void multeq_ad(ITER1 py, ITER2 px) const
    {
      // vmult with the adjoint matrix
      const_iterator    p  = std::vector<T>::begin();
      
      for(int i=0;i<m();i++)
	{
	  (*py++) = 0.;
	}
      py -= m();
      while(p!=std::vector<T>::end())
	{
	  for(int j=0;j<n();j++)
	    {
	      (*py) += (*p++) * (*px);
	      py++;
	    }
	  py -= m();
	  px++;
	}
    }

  /**************************************************/
  
  void lu()
    {
      /* LU decomposition */
      
      for(int i=1;i<n();i++)
        {
          for(int k=0;k<i;k++)
            {
              value(i,k) /= value(k,k);
              for(int j=k+1;j<n();j++)
                {
                  value(i,j) -= value(i,k)*value(k,j);
                }
            }
        }
    }

  /**************************************************/
  
  void inverse2()
    {
      T a = value(0,0);
      T b = value(0,1);
      T c = value(1,0);
      T d = value(1,1);
      value(0,0) = d;
      value(0,1) = -b;
      value(1,0) = -c;
      value(1,1) = a;

      T idet = 1./(a*d-b*c);
      *this *= idet;
    }

  /**************************************************/
  
  void inverse3()
    {
      T idet = 1./det();

      T a = value(0,0);
      T b = value(0,1);
      T c = value(0,2);
      T d = value(1,0);
      T e = value(1,1);
      T f = value(1,2);
      T g = value(2,0);
      T h = value(2,1);
      T i = value(2,2);
      value(0,0) = e*i-f*h;
      value(0,1) = h*c-b*i;
      value(0,2) = b*f-e*c;
      value(1,0) = f*g-d*i;
      value(1,1) = a*i-g*c;
      value(1,2) = d*c-a*f;
      value(2,0) = d*h-g*e;
      value(2,1) = g*b-a*h;
      value(2,2) = a*e-b*d;

      *this *= idet;
    }

  /**************************************************/
  
  void inverse()
    {
      if (n()==2)
	{
	  inverse2();
	  return;
	}
      if (n()==3)
	{
	  inverse3();
	  return;
	}
      /* LU decomposition */
      
      for(int i=1;i<n();i++)
	{
	  for(int k=0;k<i;k++)
	    {
	      value(i,k) /= value(k,k);
	      for(int j=k+1;j<n();j++)
		{
		  value(i,j) -= value(i,k)*value(k,j);
		}
	    }
	}

      /* Inverse von L */
      
      for(int ncol=0;ncol<n()-1;ncol++)
	{
	  for(int i=ncol+1;i<n();i++)
	    {
	      value(i,ncol) = -value(i,ncol);
	      for(int k=ncol+1;k<i;k++)
		{
		  value(i,ncol) -= value(i,k)*value(k,ncol);
		}
	    }
	}


      /* Inverse von U */
      
      
      for(int nlin=0;nlin<n();nlin++)
	{
	  for(int j=nlin+1;j<n();j++)
	    {
	      value(nlin,j) /= -value(nlin,nlin);
	      for(int k=nlin+1;k<j;k++)
		{
		  value(nlin,j) -= value(nlin,k)*value(k,j);
		}
	      value(nlin,j) /= value(j,j);
	    }
	  value(nlin,nlin) = 1./value(nlin,nlin);
	}
      
      
      /* Inverse von A */
      
      for(int ncol=0;ncol<n();ncol++)
	{
	  for(int i=0;i<ncol+1;i++)
	    {
	      for(int k=ncol+1;k<n();k++)
		{
		  value(i,ncol) += value(i,k)*value(k,ncol);
		}
	    }
	  for(int i=ncol+1;i<n();i++)
	    {
	      value(i,ncol) *= value(i,i);
	      for(int k=i+1;k<n();k++)
		{
		  value(i,ncol) += value(i,k)*value(k,ncol);
		}
	    }
	}
    }

  void gauss_jordan()
    {
      nvector<int> p(n());
      iota(p.begin(),p.end(),0);
      
      for (int j=0;j<n();j++)
	{
	  double max = fabs(value(j,j));
	  int r = j;
	  for (int i=j+1;i<n();i++)
	    {
	      if (fabs(value(i,j)) > max)
		{
		  max = fabs(value(i,j));
		  r = i;
		}
	    }
	  if (r>j)
	    {
	      for (int k=0; k<n(); k++)
		{
		  //swap(value(j,k),value(r,k));
		  T h        = value(j,k);
		  value(j,k) = value(r,k);
		  value(r,k) = h;
		}
	      //swap(p[j],p[r]);
	      int h = p[j]; p[j] = p[r]; p[r] = h;
	    }
	  
	  double hr = 1./value(j,j);
	  value(j,j) = hr;
	  for (int k=0;k<n();k++)
	    {
	      if (k==j) continue;
	      for (int i=0;i<n();i++)
		{
		  if (i==j) continue;
		  value(i,k) -= value(i,j)*value(j,k)*hr;
		}
	    }
	  for (int i=0;i<n();i++)
	    {
	      value(i,j) *= hr;
	      value(j,i) *= -hr;
	    }
	  value(j,j) = hr;
	}
      nvector<double> hv(n());
      for (int i=0;i<n();i++)
	{
	  for (int k=0;k<n();k++) hv[p[k]] = value(i,k);
	  for (int k=0;k<n();k++) value(i,k) = hv[k];
	}
    }
};
}

#endif

