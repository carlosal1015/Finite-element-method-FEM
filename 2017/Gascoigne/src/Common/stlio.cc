/**
*
* Copyright (C) 2004 by the Gascoigne 3D authors
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


#include  "stlio.h"
#include  "vertex.h"
#include  "compvector.h"

using namespace std;

namespace Gascoigne
{

// /*-------------------------------------------------------------------*/

// ostream& operator<<(ostream &s, const map<int,fixarray<4,int> >& A)
// {
//   for(map<int,fixarray<4,int> >::const_iterator p=A.begin();p!=A.end();p++)
//     {
//       s << "\t"<<p->first<<"\t"<<p->first<<"\t --> "<<p->second<<endl;
//     }
//   return s;
// }

// /*-------------------------------------------------------------------*/

// ostream& operator<<(ostream &s, const map<pair<string,string>,int>& A)
// {
//   for(map<pair<string,string>,int>::const_iterator p=A.begin();p!=A.end();p++)
//     {
//       s << "\t"<<p->first.first<<"\t"<<p->first.second<<"\t --> "<<p->second<<endl;
//     }
//   return s;
// }

// ostream& operator<<(ostream &s, const map<string,int>& A)
// {
//   for(map<string,int>::const_iterator p=A.begin();p!=A.end();p++)
//     {
//       s << "\t"<<p->first<<"\t --> "<<p->second<<endl;
//     }
//   return s;
// }

// ostream& operator<<(ostream &s, const map<string,double>& A)
// {
//   for(map<string,double>::const_iterator p=A.begin();p!=A.end();p++)
//     {
//       s << "\t"<<p->first<<"\t --> "<<p->second<<endl;
//     }
//   return s;
// }

// ostream& operator<<(ostream &s, const map<string,string>& A)
// {
//   for(map<string,string>::const_iterator p=A.begin();p!=A.end();p++)
//     {
//       s << "\t"<<p->first<<"\t --> "<<p->second<<endl;
//     }
//   return s;
// }

// /*-------------------------------------------------------------------*/

// ostream& operator<<(ostream &s, const set<string>& A)
// {
//   for(set<string>::const_iterator p=A.begin();p!=A.end();p++)
//     {
//       s << *p << " ";
//     }
//   return s;
// }

// ostream& operator<<(ostream &s, const set<int>& A)
// {
//   for(set<int>::const_iterator p=A.begin();p!=A.end();p++)
//     {
//       s << *p << " ";
//     }
//   return s;
// }

// /*-------------------------------------------------------------------*/

// ostream& operator<<(ostream& s, const vector<pair<int,int> >& A)
// {
//   typedef vector<pair<int,int> >::const_iterator it;
//   for(it p=A.begin();p!=A.end();p++)
//     {
//       s << "( "<<p->first <<" "<<p->second<<" )   ";
//     }
//   return s;
// }

void write_data(const GlobalVector& v,ostream& s)
{
  s << v.n() << endl
    << v.ncomp() << endl
    << v << endl;
}

void read_data(GlobalVector& v, istream& s)
{
  int n;
  int ncomp;
  s >> n;
  s >> ncomp;
  v.ncomp()=ncomp;
  v.resize(n,0);
  s >> v;
}

template<class T>
void write_data(const set<T>& v,ostream& s)
{
  s << v.size() << endl << v;
}

template<class T>
void read_data(set<T>& v, istream& s)
{
  size_t n;
  s >> n;
  for (int i=0;i<n;++i)
    {
      T a;
      s >> a;
      v.insert(a);
    }
}




void write_data(const int& v,ostream& s) { s << v; }
void read_data(int& v, istream& s)        { s >> v; }

template<int N,class T>
void write_data(const fixarray<N,T>& v,ostream& s) { s << v; }
template<int N,class T>
void read_data(fixarray<N,T>& v, istream& s)        { s >> v; }

void write_data(const double& v,ostream& s) { s << v; }
void read_data(double& v, istream& s)        { s >> v; }

void write_data(const string& v,ostream& s) { s << v; }
void read_data(string& v, istream& s)        { s >> v; }

template<class T>
void write_data(const vector<T>& v,ostream& s)
{
  s << v.size() << endl;
  for (int i=0;i<v.size();++i)
    {
      write_data(v[i],s);
      s << endl;
    }
  s << endl;
}


template<class T>
void read_data(vector<T>& v, istream& s)
{
  size_t n;
  s >> n;
  if(v.size()!=n) v.resize(n);
  for (int i=0;i<n;++i)
    read_data(v[i],s);
}

template void write_data(const fixarray<2,int>&,ostream& );
template void write_data(const fixarray<3,int>&,ostream& );
template void write_data(const fixarray<4,int>&,ostream& );
template void write_data(const fixarray<8,int>&,ostream& );
template void write_data(const fixarray<9,int>&,ostream& );

template void write_data(const vector<int>&,ostream& );
template void write_data(const vector<Vertex<2> >&,ostream& );
template void write_data(const vector<Vertex<3> >&,ostream& );
template void write_data(const vector<IntVector >&,ostream& );
template void write_data(const vector<fixarray<4,int> >&,ostream& );
template void write_data(const vector<fixarray<8,int> >&,ostream& );
template void write_data(const vector<set<int> >&,ostream& );

template void read_data(fixarray<2,int>&, istream& );
template void read_data(fixarray<3,int>&, istream& );
template void read_data(fixarray<4,int>&, istream& );
template void read_data(fixarray<8,int>&, istream& );
template void read_data(fixarray<9,int>&, istream& );

template void read_data(vector<int>&, istream& );
template void read_data(vector<Vertex<2> > &, istream& );
template void read_data(vector<Vertex<3> > &, istream& );
template void read_data(vector<IntVector >&, istream& );
template void read_data(vector<fixarray<4,int> >&, istream& );
template void read_data(vector<fixarray<8,int> >&, istream& );
template void read_data(vector<set<int> >&, istream& );
template void read_data(set<int> &, istream& );

}
