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


#ifndef __numfixarray_h
#define __numfixarray_h

#include  "fixarray.h"
#include  "gascoignemath.h"
#include  "nvector.h"

/*-------------------------------------------------*/

namespace Gascoigne
{
template<int N, class T>
class numfixarray  : public fixarray<N,T>
{
    typedef typename fixarray<N,T>::iterator         nvp;
    typedef typename fixarray<N,T>::const_iterator   cnvp;

  public:
    ~numfixarray()    {}
    numfixarray()                    : fixarray<N,T>(0.)    {}
/*   numfixarray(size_t n)               : fixarray<N,T>(n)   {} */
    numfixarray(const T& d)   : fixarray<N,T>(d) {}
    numfixarray(const numfixarray<N,T>& v) : fixarray<N,T>(v)   {}

    double operator*(const numfixarray& v) const;
    double operator*(const nvector<T>& v) const;
    numfixarray<N,T>&   operator=(const T&);
    numfixarray<N,T>&   operator=(const numfixarray&);
    numfixarray<N,T>&   operator*=(const numfixarray&);
    numfixarray<N,T>&   operator*=(double d);
    numfixarray<N,T>&   operator/=(double d);
    numfixarray<N,T>&   operator+=(double d) { add(d); return *this; }
    numfixarray<N,T>&   operator+=(const numfixarray& v) { add(1.,v); return *this; }
    numfixarray<N,T>&   operator-=(const numfixarray& v) { add(-1.,v); return *this; }
    numfixarray<N,T>    operator-(numfixarray v1) const
      { v1.add(-1.,*this);
	v1*=(-1);
	return v1;
      }
    
    int     n    () const { return fixarray<N,T>::size(); }
    void    zero ();
    void    equ  (const T&);
    void    equ  (double, const numfixarray&);
    void    equ  (double, const numfixarray&, double, const numfixarray&);
    void    equ  (double, const numfixarray&, double, const numfixarray&, 
		  double, const numfixarray&);
    void    equ  (double, const numfixarray&, double, const numfixarray&, 
		  double, const numfixarray&, double, const numfixarray&);
    void    equ  (double, const numfixarray&, double, const numfixarray&, 
		  double, const numfixarray&, double, const numfixarray&,
		  double, const numfixarray&, double, const numfixarray&);
    void    sequ (double, double, const numfixarray&);
    void    add  (double);
    void    add  (T, const numfixarray&);
    void    add  (double, const numfixarray&, double, const numfixarray&);
    void    add  (double, const numfixarray&, double, const numfixarray&, 
		  double, const numfixarray&);
    void    sadd (double, double, const numfixarray&);
    double  max  ()     const;
    double  norm ()     const;
    double  norm_l8 ()     const;
    double  norm_l1()     const;
    double  norm_l2()     const {return norm();}

    void normalise()
      {
	T d=0.;
	for(int i=0;i<n();i++)
	  {
	    d += (*this)[i]*(*this)[i];
	  }
	d = 1./sqrt(d);
	for(int i=0;i<n();i++)
	  {
	    (*this)[i] *= d;
	  }
      }

};


/**************************************************/

template<int N,class T>
inline double numfixarray<N,T>::norm_l8() const
{
  cnvp first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();

  double d=0.;
  while( first != last)
    {
      d = Gascoigne::max(d,fabs(*first));
      first++;
    }
  return d;
}

/**************************************************/

template<int N,class T>
inline double numfixarray<N,T>::norm() const
{
  cnvp first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();

  double n=0.;
  while( first != last)
    {
      n += ((*first)) * ((*first));
      first++;
    }
  return sqrt(n);
}

/**************************************************/

template<int N,class T>
inline double numfixarray<N,T>::norm_l1() const
{
  cnvp first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();

  double n=0.;
  while( first != last)
    {
      n += fabs((*first++));
    }
  return n;
}

/**************************************************/

template<int N,class T>
inline numfixarray<N,T>& numfixarray<N,T>::operator=(const T& d)
{
  nvp  first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();

  while( first != last)
    {
      *first++ = d;
    }
  return *this;
}

/**************************************************/

template<int N,class T>
inline numfixarray<N,T>& numfixarray<N,T>::operator=(const numfixarray<N,T>& v)
{
  nvp  first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();
  cnvp vfirst = v.fixarray<N,T>::begin();

  while( first != last)
    {
      *first++ = *vfirst++;
    }
  return *this;
}

/**************************************************/

template<int N,class T>
inline numfixarray<N,T>& numfixarray<N,T>::operator*=(const numfixarray<N,T>& d)
{
  nvp  first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();
  cnvp fd     = d.fixarray<N,T>::begin();

  while(first != last)
    {
      (*first++) *= (*fd++);
    }
  return *this;
}

/**************************************************/

template<int N,class T>
inline numfixarray<N,T>& numfixarray<N,T>::operator*=(double d)
{
  nvp  first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();

  while(first != last)
    {
      (*first++) *= d;
    }
  return *this;
}

/**************************************************/

template<int N,class T>
inline numfixarray<N,T>& numfixarray<N,T>::operator/=(double d)
{
  nvp  first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();

  while(first != last)
    {
      (*first++) /= d;
    }
  return *this;
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::zero()
{
  nvp  first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();

  while( first != last)
    {
      *first++ = 0.;
    }
}

/**************************************************/

template<int N,class T>
inline double numfixarray<N,T>::max() const
{
  double d = -1.;
  cnvp first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();

  while( first != last)
    {
      d = MAX( d, fabs((*first)));
      first++;
    }
  return d;
}

/**************************************************/

template<int N,class T>
inline double numfixarray<N,T>::operator* (const numfixarray<N,T>& v) const
{
  cnvp first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();
  cnvp first2 = v.fixarray<N,T>::begin();

  double d = 0.;
  while(first != last)
    {
      d += (*first++) * (*first2++);
    }
  return d;
}

/**************************************************/

template<int N,class T>
inline double numfixarray<N,T>::operator* (const nvector<T>& v) const
{
  cnvp first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();
  cnvp first2 = v.fixarray<N,T>::begin();

  double d = 0.;
  while(first != last)
    {
      d += (*first++) * (*first2++);
    }
  return d;
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::equ (const T& d)
{
  nvp  first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();

  while(first != last)
    {
      (*first++) = d;
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::equ (double d, const numfixarray<N,T>& v)
{
  nvp  first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();
  cnvp first2 = v.fixarray<N,T>::begin();

  while(first != last)
    {
      (*first++) = d*(*first2++);
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::equ (double d, const numfixarray<N,T>& v, 
				   double e, const numfixarray<N,T>& w)
{
  nvp  first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();
  cnvp first2 = v.fixarray<N,T>::begin();
  cnvp first3 = w.fixarray<N,T>::begin();

  while(first != last)
    {
      (*first++) = d*(*first2++) + e*(*first3++);
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::equ (double d, const numfixarray<N,T>& v, 
				   double e, const numfixarray<N,T>& w,
				   double f, const numfixarray<N,T>& x)
{
  nvp  first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();
  cnvp first2 = v.fixarray<N,T>::begin();
  cnvp first3 = w.fixarray<N,T>::begin();
  cnvp first4 = x.fixarray<N,T>::begin();

  while(first != last)
    {
      (*first++) = d*(*first2++) + e*(*first3++) + f*(*first4++);
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::equ (double d, const numfixarray<N,T>& v, 
				   double e, const numfixarray<N,T>& w,
				   double f, const numfixarray<N,T>& x,
				   double g, const numfixarray<N,T>& y	   )
{
  nvp  first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();
  cnvp first2 = v.fixarray<N,T>::begin();
  cnvp first3 = w.fixarray<N,T>::begin();
  cnvp first4 = x.fixarray<N,T>::begin();
  cnvp first5 = y.fixarray<N,T>::begin();

  while(first != last)
    {
      (*first++) = d*(*first2++) + e*(*first3++) + f*(*first4++) + g*(*first5++);
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::equ (double d, const numfixarray<N,T>& v, 
				   double e, const numfixarray<N,T>& w,
				   double f, const numfixarray<N,T>& x,
				   double g, const numfixarray<N,T>& y,
				   double h, const numfixarray<N,T>& z,
				   double i, const numfixarray<N,T>& zz)
{
  nvp  first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();
  cnvp first2 = v.fixarray<N,T>::begin();
  cnvp first3 = w.fixarray<N,T>::begin();
  cnvp first4 = x.fixarray<N,T>::begin();
  cnvp first5 = y.fixarray<N,T>::begin();
  cnvp first6 = z.fixarray<N,T>::begin();
  cnvp first7 = zz.fixarray<N,T>::begin();

  while(first != last)
    {
      (*first++) = d*(*first2++) + e*(*first3++) + f*(*first4++) + g*(*first5++) 
		   + h*(*first6++) + i*(*first7++);
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::sequ (double s, double d, const numfixarray<N,T>& v)
{
  equ(s,*this);
  add(d,v);
  return;

  nvp  first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();
  cnvp first2 = v.fixarray<N,T>::begin();

  while(first != last)
    {
      (*first) = s*(*first) + d*(*first2++);
      first++;
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::add (double d)
{
  nvp  first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();

  while(first != last)
    {
      (*first++) += d;
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::add (T d, const numfixarray<N,T>& v)
{
  nvp  first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();
  cnvp first2 = v.fixarray<N,T>::begin();

  while(first != last)
    {
      (*first++) += d*(*first2++);
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::add (double d, const numfixarray<N,T>& v,
				   double e, const numfixarray<N,T>& w)
{
  nvp  first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();
  cnvp first2 = v.fixarray<N,T>::begin();
  cnvp first3 = w.fixarray<N,T>::begin();

  while(first != last)
    {
      (*first++) += d*(*first2++) + e*(*first3++);
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::add (double d, const numfixarray<N,T>& v,
				   double e, const numfixarray<N,T>& w,
				   double f, const numfixarray<N,T>& x)
{
  nvp  first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();
  cnvp first2 = v.fixarray<N,T>::begin();
  cnvp first3 = w.fixarray<N,T>::begin();
  cnvp first4 = x.fixarray<N,T>::begin();

  while(first != last)
    {
      (*first++) += d*(*first2++) + e*(*first3++) + f*(*first4++);
    }
}

/**************************************************/

template<int N,class T>
inline void numfixarray<N,T>::sadd (double a, double d, const numfixarray<N,T>& v)
{
  nvp  first  = fixarray<N,T>::begin();
  cnvp last   = fixarray<N,T>::end();
  cnvp first2 = v.fixarray<N,T>::begin();

  while(first != last)
    {
      (*first) = a*(*first) + d*(*first2++);
      first++;
    }
}
}

#endif

