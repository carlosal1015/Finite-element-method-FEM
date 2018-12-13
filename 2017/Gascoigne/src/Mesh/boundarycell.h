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


#ifndef __boundarycell_h
#define __boundarycell_h

#include  "cell.h"

#define NEDGES 1

namespace Gascoigne
{
template<int N>
class BoundaryCell : public Cell<N,NEDGES>
{
 protected:

  int mat, eiq, oq;

 public:

  BoundaryCell(int l = 0, int f = -1) : Cell<N,NEDGES>(l,f), mat(0), eiq(0), oq(0) {}
  BoundaryCell(const BoundaryCell& c) : Cell<N,NEDGES>(c)  , mat(c.material()), 
    eiq(c.edge_in_quad()), oq(c.of_quad()) {}

  int  nnchild ()     const { return N; }
  int  material()     const { return mat; }
  int& material()           { return mat; }
  int  edge_in_quad() const { return eiq; }
  int& edge_in_quad()       { return eiq; }
  int  of_quad()      const { return oq; }
  int& of_quad()            { return oq; }

  BoundaryCell<N>& operator=(const BoundaryCell<N>& c)
    {
      Cell<N,NEDGES>::operator=(c);
      mat = c.material();
      eiq = c.edge_in_quad();
      oq  = c.of_quad();
      return *this;
    }

  void BinWrite(std::ostream &s) const
  {
    Cell<N,NEDGES>::BinWrite(s);
    int sizeInt = sizeof(int);
    s.write(reinterpret_cast<const char*>(&(this->qlevel)),sizeInt);
    s.write(reinterpret_cast<const char*>(&(this->qfather)),sizeInt);
    s.write(reinterpret_cast<const char*>(&oq),sizeInt);
    s.write(reinterpret_cast<const char*>(&eiq),sizeInt);
    int nc = this->nchilds();
    s.write(reinterpret_cast<const char*>(&nc),sizeInt);
    for (int i=0; i<nc; i++)
    {
      s.write(reinterpret_cast<const char*>(&(this->qchilds[i])),sizeInt);
    }
  }

  void BinRead(std::istream &s)
  {
    Cell<N,NEDGES>::BinRead(s);
    int sizeInt = sizeof(int);
    s.read(reinterpret_cast<char*>(&(this->qlevel)),sizeInt);
    s.read(reinterpret_cast<char*>(&(this->qfather)),sizeInt);
    s.read(reinterpret_cast<char*>(&oq),sizeInt);
    s.read(reinterpret_cast<char*>(&eiq),sizeInt);
    int nc;
    s.read(reinterpret_cast<char*>(&nc),sizeInt);
    this->childs().resize(nc);
    for (int i=0; i<nc; i++)
    {
      s.read(reinterpret_cast<char*>(&(this->qchilds[i])),sizeInt);
    }
  }

  friend std::ostream& operator<<(std::ostream &s, const BoundaryCell& A)
    {
      s << " : ";
      s << A.vertex()  << " " << A.level()   << " ";
      s << A.father()  << " " << A.of_quad() << " ";
      s << A.edge_in_quad() << " @ ";
      s << A.nchilds() << " " << A.childs();
      
      
      /*s << "s-l-f-q  " << A.sleep() << " " << A.level();
	s << " " << A.father() << " " << A.of_quad() << std::endl;
	s << "Childs " << A.childs();
	s << std::endl;*/
      return s;
    }
  friend std::istream& operator>>(std::istream &s, BoundaryCell& A)
    {
      std::string symbol;
      int n;
      s >> symbol;
      if (symbol!=":")
	{
	  std::cout << "ERROR in BoundaryCell::operator>>" << std::endl;
	  exit(1);
	}
      s >> A.vertex() >> A.level();
      s >> A.father() >> A.of_quad() >> A.edge_in_quad();
      s >> symbol >> n;
      if (symbol!="@")
	{
	  std::cout << "ERROR in BoundaryCell::operator>>" << std::endl;
	  exit(1);
	}
      A.childs().resize(n);
      s >> A.childs();

      return s;
    }
};
}

#undef NEDGES

#endif
