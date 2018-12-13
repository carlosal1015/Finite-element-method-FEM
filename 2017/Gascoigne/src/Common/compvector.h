/**
*
* Copyright (C) 2004, 2005, 2006, 2007 by the Gascoigne 3D authors
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


#ifndef  __compvector_h
#define  __compvector_h

#include  "nvector.h"
#include  "nmatrix.h"
#include  <cassert>

/*------------------------------------------------------*/

namespace Gascoigne
{
template<class T>
class CompVector : public nvector<T>
{
 protected:

  int N;

 public:

  typedef typename nvector<T>::const_iterator   const_iterator;
  typedef typename nvector<T>::iterator         iterator;

  int ncomp() const { return N; }
  int& ncomp() { return N; }

  ~CompVector() {}
  CompVector()                             : nvector<T>()      , N(0)  {}
  CompVector(int NN)                       : nvector<T>()      , N(NN) {}
  CompVector(int NN, size_t n)             : nvector<T>(NN*n)  , N(NN) {}
  CompVector(int NN, size_t n, const T& d) : nvector<T>(NN*n,d), N(NN) {}

  CompVector(const std::vector<T>& u)
    {
      N = 1;
      nvector<T>::reservesize(u.size());
      copy(u.begin(),u.end(),nvector<T>::begin());
    }

  CompVector(const CompVector& u)
    {
    ReInit(u);
    }

   CompVector& operator=(double d)
     {
       nvector<T>::operator=(d);
       return *this;
     }

  size_t  n() const { return nvector<T>::size()/N; } 

  const_iterator  start(int i) const { return std::vector<T>::begin() + i*N; }
  iterator        start(int i)       { return std::vector<T>::begin() + i*N; }
  const_iterator  stop (int i) const { return std::vector<T>::begin() + (i+1)*N; }
  iterator        stop (int i)       { return std::vector<T>::begin() + (i+1)*N; }

  const T& operator()(int i, int c) const {return *(start(i)+c);}
  T&       operator()(int i, int c)       {return *(start(i)+c);}

  void ReInit(size_t ncomp, size_t n) {assert(ncomp); N=ncomp; reservesize(n);}
  
  void ReInit(const CompVector& u) {
    N = u.ncomp();
    nvector<T>::reservesize(u.size());
    copy(u.begin(),u.end(),nvector<T>::begin());
  }
  
  void reservesize(size_t n, const T& s=0) { nvector<T>::reservesize(n*N,s); }
  void reservesize(const CompVector& u) 
    {
      ncomp() = u.ncomp();
      nvector<T>::reservesize(u.size()); 
    }
  void resize(size_t n, const T& s=0.) { nvector<T>::resize(n*N,s); }
  void total_reservesize(size_t n)    { nvector<T>::reservesize(n); }

  void equ_node(int i, double d0)
    {
      iterator       p = start(i);
      const_iterator q = start(i)+N;
      while(p!=q) *(p++) = d0;
   }
  void scale_comp(int c,double d0)                                                                                                                                    
    { 
      iterator       p = nvector<T>::begin()+c; 
      while(p<nvector<T>::end()) {*p *= d0; p+=N;} 
    } 
  void scale_node(int i, double d0)
    {
      iterator       p = start(i);
      const_iterator q = start(i)+N;
      while(p!=q) *(p++) *= d0;
   }

  void add_node(int i, double d0, const nvector<T>& u0)
    {
            iterator       p = start(i);
      const_iterator pp= p+N;
      const_iterator q = u0.begin();
      while(p!=pp) *(p++) += d0 * *(q++);
    }

  void equ_node(int i, int j, const CompVector<T>& u0)
    {
            iterator       p = start(i);
      const_iterator pp= p+N;
      const_iterator q = u0.start(j);
      while(p!=pp) *(p++) = *(q++);
    }

  void equ_node(int i, double d0, int i0, const CompVector<T>& u0)
    {
            iterator       p = start(i);
      const_iterator pp= p+N;
      const_iterator q = u0.start(i0);
      while(p!=pp) *(p++) = d0 * *(q++);
    }
  void equ_node(int i, double d0, int i0, double d1, int i1)
    {
            iterator       p = start(i);
      const_iterator pp= p+N;
      const_iterator q0 = start(i0);
      const_iterator q1 = start(i1);
      while(p!=pp) *(p++) = d0* *q0++ + d1* *q1++;
    }
  void equ_node(int i, double d0, int i0, double d1, int i1, double d2, int i2)
    {
            iterator       p = start(i);
      const_iterator pp= p+N;
      const_iterator q0 = start(i0);
      const_iterator q1 = start(i1);
      const_iterator q2 = start(i2);
      while(p!=pp) *(p++) = d0* *q0++ + d1* *q1++ + d2* *q2++;
    }
  void add_node(int i, double d0, int i0, double d1, int i1, double d2, int i2)
    {
            iterator       p = start(i);
      const_iterator pp= p+N;
      const_iterator q0 = start(i0);
      const_iterator q1 = start(i1);
      const_iterator q2 = start(i2);
      while(p!=pp) *(p++) += d0* *q0++ + d1* *q1++ + d2* *q2++;
    }
  void equ_node(int i, double d0, int i0, double d1, int i1, 
		       double d2, int i2, double d3, int i3)
    {
            iterator       p = start(i);
      const_iterator pp= p+N;
      const_iterator q0 = start(i0);
      const_iterator q1 = start(i1);
      const_iterator q2 = start(i2);
      const_iterator q3 = start(i3);
      while(p!=pp) *(p++) = d0* *q0++ + d1* *q1++ + d2* *q2++ + d3* *q3++;
    }
   void equ_node(int i, double d0, int i0, const CompVector<T>& u0,
		double d1, int i1, const CompVector<T>& u1)
    {
            iterator       p  = start(i);
      const_iterator q0 = u0.start(i0);
      const_iterator q1 = u1.start(i1);
      for(int c=0;c<N;c++)
        {
          *(p++) = d0* *(q0++) + d1* *(q1++) ;
        }
    }

  void zero_comp(int c)
    {
      iterator       p = std::vector<T>::begin()+c;
      while(p<nvector<T>::end()) {*p = 0.; p+=N;}
    }

  void zero_node(int i)
    {
            iterator       p = start(i);
      const_iterator q = p+N;;
      while(p<q) *(p++) = 0.;
    }

  void add_node(int i, double d0, int i0)
    {
            iterator       p = start(i);
      const_iterator q = start(i0);
      for(int c=0;c<N;c++) *(p++) += d0 * *(q++);
    }

  void add_node(int i, double d0, int i0, double d1, int i1)
    {
            iterator       p = start(i);
      const_iterator q0 = start(i0);
      const_iterator q1 = start(i1);
      for(int c=0;c<N;c++) *(p++) += d0 * *(q0++) + d1* *(q1++);
    }

  void add_node(int i, double d0, int i0, const CompVector<T>& u0)
    {
            iterator       p = start(i);
      const_iterator q = u0.start(i0);
      for(int c=0;c<N;c++) *(p++) += d0 * *(q++);
    }
  double CompScp(int c, const CompVector<T>& v) const
    {
      double d = 0.;
      const_iterator first  = std::vector<T>::begin()+c;
      const_iterator last   = std::vector<T>::end();
      const_iterator first2  = v.begin()+c;
      
      while( first < last)
        {
          d += *first * *first2;
          first += N;
          first2 += N;
        }
      return d;
    }
  double CompNormL8(int c) const
    {
      double d = 0;
      const_iterator first  = std::vector<T>::begin()+c;
      const_iterator last   = std::vector<T>::end();
      
      while( first < last)
        {
          d = Gascoigne::max( d, fabs(*first));
          first += N;
        }
      return d;
    }
  double CompMin(int c) const
    {
      double d = 1.e40;
      const_iterator first  = std::vector<T>::begin()+c;
      const_iterator last   = std::vector<T>::end();
      
      while( first < last)
        {
          d = Gascoigne::min( d, *first);
          first += N;
        }
      return d;
    }
  double CompMax(int c) const
    {
      double d = -1e14;
      const_iterator first  = std::vector<T>::begin()+c;
      const_iterator last   = std::vector<T>::end();
      
      while( first < last)
        {
          d = Gascoigne::max( d, *first);
          first += N;
        }
      return d;
    }
  //////////////////////////////////////
  void SetMax(int c, double val)
    {
      iterator first  = std::vector<T>::begin()+c;
      const_iterator last   = std::vector<T>::end();
      
      while( first < last)
	{
	  (*first) = Gascoigne::max(val,(*first));
	  first += N;
	}
    }
//////////////////////////////////////
  void FillLocal(int i, nvector<double>& uloc) const
    {
      assert(uloc.size()==N);
      const_iterator  first  = start(i);
      for(int ii=0;ii<N;++ii) uloc[ii] = *first++;
    }
  void node_zero(int i)
    {
      iterator  first  = start(i);
      const_iterator last   = stop(i);
      
      while( first != last)      {*first++ = 0.;}
    }
  void CompAdd(int c, double d)
    {
      iterator       first  = std::vector<T>::begin()+c;
      const_iterator last   = std::vector<T>::end();
      
      while( first < last)
        {
          *first += d;
          first += N;
        }
    }
  void CompAdd(int c1, double d, int c2, const CompVector& y)
    {
      iterator       first  = std::vector<T>::begin()+c1;
      const_iterator first2 = y.begin()+c2;
      const_iterator last   = std::vector<T>::end();
      int N2 = y.ncomp();
      
      while( first < last)
        {
          *first += d * *first2;
          first  += N;
          first2 += N2;
        }
    }
  void CompAdd(int c, double d, const CompVector& y)
    {
      CompAdd(c,d,c,y);
    }
  void CompEq(int c, double d)
    {
      iterator       first  = std::vector<T>::begin()+c;
      const_iterator last   = std::vector<T>::end();
      
      while( first < last)
        {
          *first = d;
          first += N;
        }
    }
  void CompEq(int c1, double d,  int c2, const CompVector& y)
    {
      iterator       first  = std::vector<T>::begin()+c1;
      const_iterator first2 = y.begin()+c2;
      const_iterator last   = std::vector<T>::end();
      int N2 = y.ncomp();
      
      while( first < last)
        {
          *first = d * *first2;
          first  += N;
          first2 += N2;
        }
    }
  void CompEq(int c, double d,  const CompVector& y)
    {
      CompEq(c,d,c,y);
    }
  double CompSum(int c) const
    {
      double d = 0.;
      const_iterator first  = std::vector<T>::begin()+c;
      const_iterator last   = std::vector<T>::end();
      
      while( first < last)
        {
          d += *first;
          first += N;
        }
      return d;
    }
  double CompNorm(int c) const
    {
      double d = 0.;
      const_iterator first  = std::vector<T>::begin()+c;
      const_iterator last   = std::vector<T>::end();
      
      while( first < last)
        {
          d += *first * *first;
          first += N;
        }
      return sqrt(d);
    }
  nvector<double> CompNorm() const
    {      
      nvector<double> d(N,0.);

      for(int c=0;c<N;c++)
        {
          const_iterator       first   = std::vector<T>::begin()  +c;
	  const_iterator last   = std::vector<T>::end();
      
	  while( first < last)
	    {
	      d[c] += *first * *first;
	      first += N;
	    }
	  d[c] = sqrt(d[c]);
	}
      return d;
    }
  void Add(const nvector<double>& scp, const CompVector<T>& y)
    {
      for(int c=0;c<N;c++)
        {
          iterator       first   = std::vector<T>::begin()  +c;
          const_iterator first2  = y.begin()+c;
          const_iterator last    = std::vector<T>::end();
          while( first < last) 
            {
              *first += scp[c] * *first2;
              first += N;   first2 += N;
            }
        }
    }
  void ScalarProductComp(nvector<double>& scp, const CompVector<T>& y) const
    {
      scp.resize(N);
      scp.zero();

      for(int c=0;c<N;c++)
        {
          const_iterator first   = std::vector<T>::begin()  +c;
          const_iterator first2  = y.begin()+c;
          const_iterator last    = std::vector<T>::end();
          while( first < last) 
            {
              scp[c] += *first * *first2;
              first += N;   first2 += N;
            }
        }
    }
  void ScalarProductCompMatrix(nmatrix<double>& scp, const CompVector<T>& y) const
    {
      scp.memory(N,N);
      scp.zero();

      for(int c=0;c<N;c++)
        {
          for(int d=0;d<N;d++)
            {
              if(d>c) continue;
              const_iterator first   = std::vector<T>::begin()  +c;
              const_iterator first2  = y.begin()+d;
              const_iterator last    = std::vector<T>::end();
              while( first < last) 
          {
            scp(c,d) += *first * *first2;
            first += N;   first2 += N;
          }
            }
        }
      for(int c=0;c<N;c++)
        {
          for(int d=0;d<N;d++)
            {
              if(d>c) scp(c,d) = scp(d,c);
            }
        }
    }
  void Read(std::istream& s, int c)
    {
      iterator first   = std::vector<T>::begin()+c;
      const_iterator last    = std::vector<T>::end();
      while( first < last) 
        {
          s >> *first;
          first += N;
        }
    }
  
  void BinWrite(std::ostream& out) const
    {
      out << ncomp() << " " << n() << std::endl << "[";
      
      int sizeT = sizeof(T);
      for(int i=0; i<nvector<T>::size(); i++)
        {
          out.write(reinterpret_cast<const char*>(&(nvector<T>::operator[](i))),sizeT);
        }
      out << "]"; 
    }

  void BinRead(std::istream& in)
    {
      char cc;
      int  c,n;
      in >> c >> n >> cc;
      ncomp() = c;
      resize(n);
      
      int sizeT = sizeof(T);
      for(int i=0; i<nvector<T>::size(); i++)
        {
          in.read(reinterpret_cast<char*>(&(nvector<T>::operator[](i))),sizeT);
        }
      in >> cc;
    }
};
}

#endif
