/**
*
* Copyright (C) 2004, 2005 by the Gascoigne 3D authors
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


#ifndef __nodematrix_h
#define __nodematrix_h

#include "numfixarray.h"


/**************************************************/

namespace Gascoigne
{
template<int N, class T>
class NodeMatrix : public numfixarray<N*N,T>
{
  
public:
  
  NodeMatrix<N,T>()           : numfixarray<N*N,T>() {}
  NodeMatrix<N,T>(const T& A) : numfixarray<N*N,T>(A) {}

  T&       operator()(int i,int j)       { return (*this)[j+i*N];}
  const T& operator()(int i,int j) const { return (*this)[j+i*N];}

  T&       value(int i,int j)       { return (*this)[j+i*N];}
  const T& value(int i,int j) const { return (*this)[j+i*N];}
  
  void reserve(int) const {};
  void resize(int) const {};

  int m() const { return N;}

  friend std::ostream& operator<<(std::ostream &s, const NodeMatrix<N,T>& A)
    {
      for(int c=0;c<N*N;c++)  s << A[c] << " ";
      s << "\t";
      return s;
    }

  void identity()
    {
    numfixarray<N*N,T>::zero();
      for(int d=0;d<N;d++)
        {
          value(d,d) = 1.;
        }
    }
  void zero_component(int c)
    {
      for(int d=0;d<N;d++)
        {
          value(c,d) = 0.;
        }
    }

/**************************************************/

  void addmult(double k, const NodeMatrix<N,T>& A, const NodeMatrix<N,T>& B)  
    // this += k*A*B
    {
      for(int i=0;i<N;i++)
        {
          for(int j=0;j<N;j++)
            {
              for(int k=0;k<N;k++)
                {
                  value(i,j) += k*A.value(i,k) * B.value(k,j); 
                }
            }
        }
    }

/**************************************************/

  void inverse() { inverse(*this);}

/**************************************************/

  void inverse(const NodeMatrix<N,T>& A)
    {
      /*for(int i=0;i<N;i++)
        {
          value(i,i) = 1./A.value(i,i);
        }
      return;*/

      *this = A;

      
      // LU decomposition

      for(int i=1;i<N;i++)
        {
          for(int k=0;k<i;k++)
            {
              value(i,k) /= value(k,k);
              for(int j=k+1;j<N;j++)
                {
                  value(i,j) -= value(i,k)*value(k,j);
                }
            }
        }

      // Inverse von L

      for(int ncol=0;ncol<N-1;ncol++)
	{
	  for(int i=ncol+1;i<N;i++)
	    {
	      value(i,ncol) = -value(i,ncol);
	      for(int k=ncol+1;k<i;k++)
		{
		  value(i,ncol) -= value(i,k)*value(k,ncol);
		}
	    }
	}
      

      // Inverse von U


      for(int nlin=0;nlin<N;nlin++)
	{
	  for(int j=nlin+1;j<N;j++)
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
      

      // Inverse von A

      for(int ncol=0;ncol<N;ncol++)
	{
	  for(int i=0;i<ncol+1;i++)
	    {
	      for(int k=ncol+1;k<N;k++)
		{
		  value(i,ncol) += value(i,k)*value(k,ncol);
		}
	    }
	  for(int i=ncol+1;i<N;i++)
	    {
	      value(i,ncol) *= value(i,i);
	      for(int k=i+1;k<N;k++)
		{
		  value(i,ncol) += value(i,k)*value(k,ncol);
		}
	    }
	}
      
    }

  void gauss_jordan()
    {
      std::vector<int> p(N);

      int i,j,k,r;
      double max, hr;
      
      for (i=0;i<N;i++) p[i] = i;
      
      for (j=0;j<N;j++)
	{
	  max = fabs(value(j,j));
	  r = j;
	  for (i=j+1;i<N;i++)
	    {
	      if (fabs(value(i,j)) > max)
		{
		  max = fabs(value(i,j));
		  r = i;
		}
	    }
	  if (r>j)
	    {
	      for (k=0;k<N;k++)
		{
		  hr = value(j,k) ; value(j,k) = value(r,k) ; value(r,k) = hr;
		}
	      i = p[j] ; p[j] = p[r] ; p[r] = i;
	    }
	  
	  hr = 1./value(j,j);
	  value(j,j) = hr;
	  for (k=0;k<N;k++)
	    {
	      if (k==j) continue;
	      for (i=0;i<N;i++)
		{
		  if (i==j) continue;
		  value(i,k) -= value(i,j)*value(j,k)*hr;
		}
	    }
	  for (i=0;i<N;i++)
	    {
	      value(i,j) *= hr;
	      value(j,i) *= -hr;
	    }
	  value(j,j) = hr;
	}
      std::vector<double> hv(N);
      for (i=0;i<N;i++)
	{
	  for (k=0;k<N;k++) hv[p[k]] = value(i,k);
	  for (k=0;k<N;k++) value(i,k) = hv[k];
	}
    }
  
};
}

#endif

