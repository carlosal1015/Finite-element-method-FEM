/**
*
* Copyright (C) 2004, 2008 by the Gascoigne 3D authors
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


#ifndef   __stlio_h
#define   __stlio_h

#include  "fixarray.h"
#include  <string>
#include  "compvector.h"
#include  "gascoigne.h"

#include  <vector>
#include  <set>
#include  <map>

#ifdef __OLDCOMPILER__
#include  <hash_map>
#define HASHMAP  hash_map
#else
#ifdef __NEWER_THAN_GCC_4_2__
#include <tr1/unordered_map>
#define HASHMAP   std::tr1::unordered_map
#else
#include  <ext/hash_map>
#define HASHMAP  __gnu_cxx::hash_map
#endif
#endif
/*---------------------------------------------*/

namespace Gascoigne
{

// std::ostream& operator<<(std::ostream &s, const std::map<int,fixarray<4,int> >& A);
// std::ostream& operator<<(std::ostream &s, const std::map<std::string,std::string>& A);
// std::ostream& operator<<(std::ostream &s, const std::map<std::string,int>& A);
// std::ostream& operator<<(std::ostream &s, const std::map<std::string,double>& A);
// std::ostream& operator<<(std::ostream &s, const std::map<std::pair<std::string,std::string>,int>& A);
// std::ostream& operator<<(std::ostream &s, const std::set<std::string>& A);
// std::ostream& operator<<(std::ostream &s, const std::set<int>& A);
// std::ostream& operator<<(std::ostream& s, const std::vector<std::pair<int,int> >& A);

// std::istream& operator>>(std::istream &s, std::set<std::string>& A);
// std::istream& operator>>(std::istream &s, std::set<int>& A);


/*----------------------------------------------*/

template<class T>
std::ostream& operator<<(std::ostream& s, const std::set<T>& A)
{
  std::ostream_iterator<T>  os(s, " ");
  copy(A.begin(),A.end(),os);
  return s;
}

template<class T>
std::istream& operator>>(std::istream& s, std::set<T>& A)
{
  std::ostream_iterator<T>  os(s, " ");
  copy(A.begin(),A.end(),os);
  return s;
}

/*---------------------------------------------*/

template<class T>
std::ostream& operator<<(std::ostream &s, const std::vector<T>& A)
{
  std::ostream_iterator<T>  os(s, " ");
  copy(A.begin(),A.end(),os);
  return s;
}

template<class T>
std::istream& operator>>(std::istream &s, std::vector<T>& A)
{
  typename std::vector<T>::iterator p = A.begin();
  while(p!=A.end())
    s >> *p++;
  return s;
}


/*---------------------------------------------*/

template<class T>
std::ostream& putvector(const std::vector<T>& v,std::ostream& s)
{
  for(typename std::vector<T>::const_iterator p = v.begin();p!=v.end();p++) p->put(s);
  return s;
}


/*------------------------------------------*/

template<class T, class S>
std::ostream& operator<<(std::ostream& os, const HASHMAP<T,S>& s)
{
  os  << s.size() << std::endl;
  for(typename HASHMAP<T,S>::const_iterator p = s.begin(); p!=s.end() ;p++)
    {
      os << p->first << "->" << p->second << " ";
    }
  os << std::endl;
  return os;
}

/*------------------------------------------*/

template<class T, class S>
std::ostream& operator<<(std::ostream& os, const std::map<T,S>& s)
{
  os  << s.size() << std::endl;
  for(typename std::map<T,S>::const_iterator p = s.begin(); p!=s.end() ;p++)
    {
      os << p->first << "->" << p->second << " ";
    }
  os << std::endl;
  return os;
}


/*---------------------------------------------*/


void write_data(const GlobalVector& v,std::ostream& s);
void read_data(       GlobalVector& v,std::istream& s);

void write_data(const int& v,std::ostream& s);
void read_data(int& v,std::istream& s)       ;

template<int N,class T>
void write_data(const fixarray<N,T>& v,std::ostream& s);
template<int N,class T>
void read_data(fixarray<N,T>& v,std::istream& s)       ;

void write_data(const double& v,std::ostream& s);
void read_data(double& v,std::istream& s)       ;

void write_data(const std::string& v,std::ostream& s);
void read_data(std::string& v,std::istream& s)       ;

/*---------------------------------------------*/

template<class T>
void write_data(const std::vector<T>& v,std::ostream& s);
template<class T>
void read_data(std::vector<T>& v,std::istream& s);

/*---------------------------------------------*/

template<class T>
void write_data(const std::set<T>& v,std::ostream& s);
template<class T>
void read_data(std::set<T>& v,std::istream& s);

/*---------------------------------------------*/

/*---------------------------------------------*/

template<class T,class S>
void write_data(const std::map<T,S>& v,std::ostream& s)
{
  s << v.size() << std::endl;
  for (typename std::map<T,S>::const_iterator it = v.begin();it!=v.end();++it)
    {
      write_data(it->first,s);  s << std::endl;
      write_data(it->second,s); s << std::endl;
    }
}

/*---------------------------------------------*/

template<class T,class S>
void read_data(std::map<T,S>& v,std::istream& s)
{
  size_t n;
  s >> n;
  for (int i=0;i<n;++i)
    {
      T fi;
      S se;
      read_data(fi,s);
      read_data(se,s);
      //Edit v[fi]=se;
      v.insert(std::make_pair<T,S>(fi,se));
    }
}

template<class T,class S>
void write_data(const HASHMAP<T,S>& v,std::ostream& s)
{
  s << v.size() << std::endl;
  for (typename HASHMAP<T,S>::const_iterator it = v.begin();it!=v.end();++it)
    {
      write_data(it->first,s);  s << std::endl;
      write_data(it->second,s); s << std::endl;
    }
}

/*---------------------------------------------*/

template<class T,class S>
void read_data(HASHMAP<T,S>& v,std::istream& s)
{
  size_t n;
  s >> n;
  for (int i=0;i<n;++i)
    {
      T fi;
      S se;
      read_data(fi,s);
      read_data(se,s);
      //Edit v[fi]=se;
      v.insert(std::make_pair<T,S>(fi,se));
    }
}
}

#undef HASHMAP

#endif


