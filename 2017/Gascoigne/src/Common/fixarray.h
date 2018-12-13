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


#ifndef __fixarray_h
#define __fixarray_h

#include  <stdlib.h> 
#include  <iostream> 
#include  <iterator> 
#include  <algorithm> 

/*-------------------------------------------------*/

namespace Gascoigne
{
template<int N, class T>
class fixarray
{

public:

  typedef  T*        iterator;
  typedef  const T*  const_iterator;
  
 protected:
  
  T  val[N];
  
  void array_copy(const_iterator q)
    {
      iterator       p(begin());
      const_iterator pe(end());
      while(p!=pe)  *p++ = *q++;
    }
  
 public:
  
  fixarray<N,T>()      { BasicInit(T());}
  fixarray<N,T>(const T& d) { BasicInit(d);}
  fixarray<N,T>(const fixarray<N,T>& v)
    {
      BasicInit(T());
      array_copy(v.begin());
      //copy(v.begin(),v.end(),begin());
    }
  fixarray(const_iterator b)
    {
      BasicInit(T());
      array_copy(b);
    }
  
  virtual ~fixarray()
    {
//       Destroy(begin(),end());
    }
  
  void BasicInit(const T& d)
    {
      // Braucht man das wirklich ???
//       for(int i=0;i<N;i++)  construct(&(val[i]),d);
      for(int i=0;i<N;i++)  val[i]=d;
    }
  
  const T*  begin() const { return &(val[0]);}
  const T*  end  () const { return &(val[0])+N;}
  T*        begin()       { return &(val[0]);}
  T*        end  ()       { return &(val[0])+N;}
  
  size_t   size()            const { return N;}
  const T& operator[](int i) const { return val[i];}
  T&       operator[](int i)       { return val[i];}
  
  fixarray<N,T>& operator=(const T& d) 
    {
      iterator  p(end());
      while(p>begin()) *--p = d;
      return *this;
    } 
  
  fixarray<N,T>& operator=(const fixarray<N,T>& v) 
    {
      iterator        p(begin());
      const_iterator  q(v.begin());
      while(p<end()) *p++ = *q++;
      return *this;
    } 
  
  bool operator<(const fixarray<N,T>& v) const
    {
      const_iterator  p(  begin());
      const_iterator  q(v.begin());
      while(p<end())
	{
	  if (*p<*q) return 1;
	  if (*q<*p) return 0;
	  p++; q++;
	}
      return 0;
    }
  bool operator!=(const fixarray<N,T>& v) const
    {
      const_iterator  p(  begin());
      const_iterator  q(v.begin());
      while(p<end())
	{
	  if (*p!=*q) return 1;
	  p++; q++;
	}
      return 0;
    }
  
    
  std::ostream& put(std::ostream &s) const
    {
      copy(begin(),end(),std::ostream_iterator<T>(s," "));
      return s;
    }  
  
  std::istream& get(std::istream &s)
    {
    typename fixarray<N,T>::iterator p;
    for(p = begin();p!=end();p++) s >> *p;
      return s;
    }
  
  void read_data(std::istream& s)
    {
      size_t n;
      s >> n;
      if(size()!=n) 
	{
	  std::cerr << "read_data(): wrong size in fixarray" << N << " " << n << std::endl;
	  exit(1);
	}
      s >> *this;
    }
  
  void write_data(std::ostream& s) const
    {
      s << size() << std::endl;
      s << *this;
    }

  void BinWrite(std::ostream &s) const
  {
    int sizeT = sizeof(T);
    for (int i=0; i<N; i++)
    {
      s.write(reinterpret_cast<const char*>(&(operator[](i))),sizeT);
    }
  }

  void BinRead(std::istream &s)
  {
    int sizeT = sizeof(T);
    for (int i=0; i<N; i++)
    {
      s.read(reinterpret_cast<char*>(&(operator[](i))),sizeT);
    }
  }
};

/*-------------------------------------------------*/

template<int N, class T>
bool operator==(const fixarray<N,T>& x, const fixarray<N,T>& y) 
{
  return std::equal(x.begin(), x.end(), y.begin());
}

/*-------------------------------------------------*/

class fixarrayHash
{
 public:
  template<int N, class T>
    int operator()(const fixarray<N,T>& h) const { return static_cast<int>(h[0]);}
};

template<int N,class T>
std::ostream& operator<<(std::ostream &s, const fixarray<N,T>& A);
template<int N,class T>
std::istream& operator>>(std::istream &s, fixarray<N,T>& A);
}


#endif


