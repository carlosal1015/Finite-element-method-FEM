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


#include  "sparsestructure.h"
#include  "gascoignemath.h"
#include  "stlio.h"

using namespace std;

/*----------------------------------------------*/

namespace Gascoigne
{
ostream& operator<<(ostream &s, const SparseStructure& A)
{
  A.statistics(s);
  s << endl;
  for(int i=0;i<A.n();i++)
    {
      s << i << " :: " << A.row(i)<<endl;
    }
  s << endl;

  return s;
}

/*----------------------------------------------*/

void SparseStructure::statistics(ostream &s) const
{
  s << n() << " " << ntotal() << "  " <<  ntotal()/static_cast<double>(n());
}

/*----------------------------------------------*/

SparseStructure& SparseStructure::operator=(const SparseStructure& A)
{
  sindices.reserve(A.n());
  sindices.resize (A.n());
  sntot = A.ntotal();
  sindices = A.indices();
  
  return *this;
}

/*----------------------------------------------*/

// SparseStructure::SparseStructure(const ColumnStencil& S) : US(&S)
// {
//   int n = S.n();
//   int ntot = S.nentries();
//   sindices.reserve(n);
//   sindices.resize (n);
//   sntot = ntot;

//   for(int i=0;i<n;i++) 
//     for(int pos=S.start(i);pos<S.stop(i);pos++)
//       sindices[i].insert(S.col(pos));
// }

/*----------------------------------------------*/

void SparseStructure::build_begin(int n)
{
  sindices.reserve(n);
  sindices.resize (n);
  for(int i=0;i<n;i++) sindices[i].clear();
}

/*----------------------------------------------*/

void SparseStructure::build_clear(int i)
{
  row(i).clear();
}

/*----------------------------------------------*/

void SparseStructure::hanging_node(int hi, int n1, int n2)
{
  // neu (eliminiert den hn in sich selbst und dann in den anderen), mit gascoigne nicht getestet!
  row(hi).erase(hi);


  for (set<int>::const_iterator p=rowbegin(hi); p!=rowend(hi); p++)
    {
      row(*p).erase(hi);
      row(*p).insert(n1);
      row(*p).insert(n2);
    }
  build_add(n1,rowbegin(hi),rowend(hi));
  build_add(n2,rowbegin(hi),rowend(hi));

  row(hi).clear();
  row(hi).insert(hi);


//   // alte version, lief mit gascoigne !
//   for (set<int>::const_iterator p=rowbegin(hi); p!=rowend(hi); p++)
//     {
//       if (*p!=hi)
// 	{
// 	  row(*p).erase(hi);
// 	  row(*p).insert(n1);
// 	  row(*p).insert(n2);
// 	}
//     }
//   build_add(n1,rowbegin(hi),rowend(hi));
//   build_add(n2,rowbegin(hi),rowend(hi));

//   row(hi).clear();
//   row(hi).insert(hi);
}

/*----------------------------------------------*/

void SparseStructure::build_end()
{
  sntot = 0;
  for(int i=0;i<n();i++)
    {
      sntot += row(i).size();
    }
}

/*----------------------------------------------*/

void SparseStructure::enlarge_lu()
{
  int maxbw = 0;
  for(int i=0;i<n();i++)
    {
      int imax = -1;
      int imin = n();
      for(set<int>::iterator p=rowbegin(i); p != rowend(i); ++p)
	{
	  int j = *p;
	  imin = Gascoigne::min_int(imin,j);
	  imax = Gascoigne::max_int(imax,j);
	}
      maxbw = Gascoigne::max_int(maxbw,abs(imax-i));
      maxbw = Gascoigne::max_int(maxbw,abs(imin-i));
    }
  //cerr << "bandbreite: " << maxbw << endl;
  for(int i=0;i<n();i++)
    {
      int imin = Gascoigne::max_int(0,i-maxbw);
      int imax = Gascoigne::min_int(i+maxbw,n()-1);
      for(int im=imin;im<=imax;im++)
	{
	  row(i).insert(im);
	}
    }
  build_end();
}

/*----------------------------------------------*/

void SparseStructure::enlarge(const SparseStructure& S)
{
  for(int i=0;i<n();i++)
    {
      for(set<int>::iterator p=S.rowbegin(i); p != S.rowend(i); ++p)
	{
	  int j = *p;
	  row(i).insert(S.rowbegin(j),S.rowend(j));
	}
    }
  build_end();
}

/*----------------------------------------------*/

void SparseStructure::enlarge_for_lu(const IntVector& p)
{
  assert (p.size()==n());
  vector<set<int> > transpose(n());

  IntVector q(n());
  for (int i=0;i<n();++i)
    q[p[i]]=i;
  
  for (int row=0;row<n();++row)
    for (set<int>::iterator col=rowbegin(row);col!=rowend(row);++col)
      transpose[*col].insert(row);
        
  for(int row=0; row<n(); row++)
    {
      int row1=p[row];
      
      for (set<int>::iterator col1=rowbegin(row1);col1!=rowend(row1);++col1)
	{
	  int col = q[*col1];
	  if (col>=row)
	    {
	      for(set<int>::iterator down1=transpose[row1].begin();down1!=transpose[row1].end();++down1)
		{
		  int down = q[*down1];
		  if (down>row)
		    {
		      this->row(*down1).insert(*col1);
		      transpose[*col1].insert(*down1);
		    }
		}
	    }
	}
      
    }
  build_end();
}
}
