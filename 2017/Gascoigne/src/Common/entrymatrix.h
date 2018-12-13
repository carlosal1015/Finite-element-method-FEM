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


#ifndef  __entrymatrix_h
#define  __entrymatrix_h

#include  "gascoigne.h"
#include  "gostream.h"

/*-------------------------------------------------------------*/

namespace Gascoigne
{
class EntryMatrix
{
  protected :

  typedef   DoubleVector::const_iterator  const_iterator;
  typedef   DoubleVector::iterator        iterator;
  typedef   std::vector<DoubleVector>     Vector;

  int              ndof, mdof, ncomp, mcomp, nmdof, nmcomp;
  int              ind;
  //  iterator         pind;
  DoubleVector  val;

  int dofindex (int i, int j) const {return nmcomp*(mdof*i+j);}
  int compindex(int c, int d) const {return mcomp*c+d;}
  int index(int i, int j,int c, int d) const {return dofindex(i,j)+compindex(c,d);}

  public :
    
    EntryMatrix() {}
  EntryMatrix(int nd, int nc)
    {
      SetDimensionDof(nd,nd);
      SetDimensionComp(nc,nc);
      resize();
    }
  EntryMatrix(int nd, int md, int nc, int mc)
    {
      SetDimensionDof(nd,md);
      SetDimensionComp(nc,mc);
      resize();
    }
  ~EntryMatrix()
    {
      val.clear();
      
    }

  const_iterator begin(int i, int j) const {return val.begin()+dofindex(i,j);}
  iterator       begin(int i, int j)       {return val.begin()+dofindex(i,j);}
  const_iterator end(int i, int j) const {return val.begin()+dofindex(i,j+1);}
  iterator       end(int i, int j)       {return val.begin()+dofindex(i,j+1);}

  void    SetDofIndex(int i, int j) { ind = dofindex(i,j);}
  double  operator()(int c, int d) const { return val[ind+compindex(c,d)];}
  double&  operator()(int c, int d) { return val[ind+compindex(c,d)];}

  void SetDimensionDof (int n, int m) {ndof =n; mdof =m; nmdof =n*m;}
  void SetDimensionComp(int n, int m) {ncomp=n; mcomp=m; nmcomp=n*m;}
  void resize() {val.reservesize(nmdof*nmcomp);}
  
  int Ncomp() const {return ncomp;}
  int Mcomp() const {return mcomp;}
  int Ndof() const {return ndof;}
  int Mdof() const {return mdof;}
  
  double   operator()(int i, int j, int c, int d) const { return val[index(i,j,c,d)];}
  double&  operator()(int i, int j, int c, int d)       { return val[index(i,j,c,d)];}
  double   operator()(int i, int j, int p) const { return val[dofindex(i,j)+p];}
  
  void zero() {val.zero();}
  
  friend std::ostream& operator<<(std::ostream& s, const EntryMatrix& A)
    {
    s << A.Ndof()  << "\t" << A.Mdof() << std::endl;
    s << A.Ncomp()  << "\t" << A.Mcomp() << std::endl;
    for(int i=0;i<A.Ndof();i++)
      {
      for(int j=0;j<A.Mdof();j++)
        {
        s << "\n[" << i << "," << j << "] ";
        for(int c=0;c<A.Ncomp();c++)
          {
          for(int d=0;d<A.Mcomp();d++)
            {
            s << A(i,j,c,d) << " ";
            }
          }
        }
      }
    return s;
    }
  
  
    void add(int il, int jl, int i, int j, const EntryMatrix& E)
    {
      iterator        p = begin(il,jl);
      const_iterator  q = end  (il,jl);
      const_iterator pE = E.begin(i,j);
      while(p!=q) *p++ += *pE++;
  }

    void zero_row(int i1)    
    {
      for(int j=0;j<mdof;j++)
	{
	  iterator       p1 = begin(i1,j);
	  const_iterator q1 = end(i1,j);
	  while(p1!=q1) *p1++ = 0.;
	}
    }
    void zero_column(int j1)
    {
      for(int i=0;i<ndof;i++)
	{
	  iterator p1 = begin(i,j1);
	  const_iterator q1 = end(i,j1);
	  while(p1!=q1) *p1++ = 0.;
	}
    }

    void add_row(int i1, int i2, double s=1.) 
    {
      for(int j=0;j<mdof;j++)
	{
	  iterator       p1 = begin(i1,j);
	  const_iterator p2 = begin(i2,j);
	  const_iterator q1 = end(i1,j);
	  while(p1!=q1) *p1++ += s * *p2++;
	}
    }

    void add_column(int j1, int j2, double s=1.)
    {
      for(int i=0;i<ndof;i++)
	{
	  iterator p1 = begin(i,j1);
	  const_iterator p2 = begin(i,j2);
	  const_iterator q1 = end(i,j1);
	  while(p1!=q1) *p1++ += s * *p2++;
	}
    }

    void add_column_row(int i1, int i2)    
    {
      add_column(i1,i2);
      add_row   (i1,i2);
    }

    void multiply_row(int i, double s=1.)    
      {
	for(int j=0;j<mdof;j++)
	  {
	    iterator p = begin(i,j);
	    const_iterator q = end(i,j);
	    while(p!=q) *p++ *= s;
	  }
      }
    void multiply_column(int j, double s=1.)    
    {
      for(int i=0;i<ndof;i++)
	{
	  iterator p = begin(i,j);
	  const_iterator q = end(i,j);
	  while(p!=q) *p++ *= s;
	}
    }

    void multiply_column_row(int j, double s=1.)    
    {
      multiply_column(j,s);
      multiply_row   (j,s);
    }

    void distribute_row(int s, int t1, int t2)    
    {
      add_column(t1, s, 0.5);
      add_column(t2, s, 0.5);
      multiply_column(s,0.);
    }

    void distribute_column(int s, int t1, int t2)    
    {
      add_row(t1, s, 0.5);
      add_row(t2, s, 0.5);
      multiply_row(s,0.);
    }
    void transpose(EntryMatrix& E)    
    {
      for(int i=0; i<ndof; i++)
	{
	  for(int j=0; j<mdof; j++)
	    {
	      SetDofIndex(i,j);
	      E.SetDofIndex(j,i);
	      for(int c=0; c<ncomp; c++)
		{
		  for(int d=0; d<mcomp; d++)
		    {
		      (*this)(c,d) = E(d,c);
		    }
		}
	    }
	}
    }
    void equ(double d)    
    {
      for(iterator p = val.begin(); p!=val.end(); p++)
	{
	  *p *= d;
	}
    }
    void mult(Vector& y, const Vector& x)
    {
      for(int i=0; i<ndof; i++)
	{
	  for(int j=0; j<mdof; j++)
	    {
	      SetDofIndex(i,j);
	      for(int c=0; c<ncomp; c++)
		{
		  for(int d=0; d<mcomp; d++)
		    {
		      y[i][c] += (*this)(c,d) * x[j][d];
		    }
		}
	    }
	}
    }
};
}
#endif
