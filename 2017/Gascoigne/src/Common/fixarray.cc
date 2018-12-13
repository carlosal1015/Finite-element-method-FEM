/**
*
* Copyright (C) 2004, 2006 by the Gascoigne 3D authors
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


#include  "fixarray.h"

using namespace std;

/*------------------------------------------------*/
/*------------------------------------------------*/
/*------------------------------------------------*/

namespace Gascoigne
{
template<int N,class T>
ostream& operator<<(ostream &s, const fixarray<N,T>& A)
{
  copy(A.begin(),A.end(),ostream_iterator<T>(s," "));
  return s;
}


template<int N,class T> istream& operator>>(istream &s, fixarray<N,T>& A)
{
  typename fixarray<N,T>::iterator p = A.begin();
  while(p!=A.end())
    s >> *p++;
  return s;
}


template ostream& operator<<(ostream &s, const fixarray<1,float>& A);
template ostream& operator<<(ostream &s, const fixarray<1,double>& A);
template ostream& operator<<(ostream &s, const fixarray<2,double>& A);
template ostream& operator<<(ostream &s, const fixarray<3,double>& A);
template ostream& operator<<(ostream &s, const fixarray<4,double>& A);
template ostream& operator<<(ostream &s, const fixarray<6,double>& A);
template ostream& operator<<(ostream &s, const fixarray<8,double>& A);
template ostream& operator<<(ostream &s, const fixarray<9,float>& A);
template ostream& operator<<(ostream &s, const fixarray<16,float>& A);
template ostream& operator<<(ostream &s, const fixarray<25,float>& A);
template ostream& operator<<(ostream &s, const fixarray<27,float>& A);
template ostream& operator<<(ostream &s, const fixarray<400,float>& A);

template istream& operator>>(istream &s,  fixarray<1,float>& A);
template istream& operator>>(istream &s,  fixarray<1,double>& A);
template istream& operator>>(istream &s,  fixarray<2,double>& A);
template istream& operator>>(istream &s,  fixarray<3,double>& A);
template istream& operator>>(istream &s,  fixarray<4,double>& A);
template istream& operator>>(istream &s,  fixarray<6,double>& A);
template istream& operator>>(istream &s,  fixarray<8,double>& A);
template istream& operator>>(istream &s,  fixarray<9,float>& A);
template istream& operator>>(istream &s,  fixarray<16,float>& A);
template istream& operator>>(istream &s,  fixarray<25,float>& A);
template istream& operator>>(istream &s,  fixarray<400,float>& A);

template ostream& operator<<(ostream &, const fixarray<1,int>&);
template ostream& operator<<(ostream &, const fixarray<2,int>&);
template ostream& operator<<(ostream &, const fixarray<3,int>&);
template ostream& operator<<(ostream &, const fixarray<4,int>&);
template ostream& operator<<(ostream &, const fixarray<6,int>&);
template ostream& operator<<(ostream &, const fixarray<8,int>&);
template ostream& operator<<(ostream &, const fixarray<9,int>& );
template ostream& operator<<(ostream &, const fixarray<27,int>& );

template istream& operator>>(istream &, fixarray<1,int>& );
template istream& operator>>(istream &, fixarray<2,int>& );
template istream& operator>>(istream &, fixarray<3,int>& );
template istream& operator>>(istream &, fixarray<4,int>& );
template istream& operator>>(istream &, fixarray<6,int>& );
template istream& operator>>(istream &, fixarray<8,int>& );
template istream& operator>>(istream &, fixarray<9,int>& );
}
