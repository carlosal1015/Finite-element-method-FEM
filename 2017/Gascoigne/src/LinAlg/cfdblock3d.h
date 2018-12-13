/**
*
* Copyright (C) 2004, 2009, 2011 by the Gascoigne 3D authors
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


#ifndef __cfdblock3d_h
#define __cfdblock3d_h

#include "entrymatrix.h"
#include "matrixentrytype.h"

/*****************************************************/

namespace Gascoigne
{
class CFDBlock3d
{
  typedef nvector<double>::iterator        iterator;
  typedef nvector<double>::const_iterator  const_iterator;

 protected:

  MatrixEntryType lapx,dx,dy,dz,gx,gy,gz,s;
  MatrixEntryType lapy,lapz;
  
  void swap(MatrixEntryType a, MatrixEntryType b) const { MatrixEntryType c = a; a=b; b=c;}

public:

  int ncomp() const { return 4;}

  MatrixEntryType laplx  () const { return lapx;}
  MatrixEntryType laply  () const { return lapy;}
  MatrixEntryType laplz  () const { return lapz;}
  MatrixEntryType divx   () const { return dx;}
  MatrixEntryType divy   () const { return dy;}
  MatrixEntryType divz   () const { return dz;}
  MatrixEntryType stab   () const { return s;}
  MatrixEntryType gradx  () const { return gx;}
  MatrixEntryType grady  () const { return gy;}
  MatrixEntryType gradz  () const { return gz;}

  MatrixEntryType& laplx  ()  { return lapx;}
  MatrixEntryType& laply  ()  { return lapy;}
  MatrixEntryType& laplz  ()  { return lapz;}
  MatrixEntryType& divx   ()  { return dx;}
  MatrixEntryType& divy   ()  { return dy;}
  MatrixEntryType& divz   ()  { return dz;}
  MatrixEntryType& stab   ()  { return s;}
  MatrixEntryType& gradx  ()  { return gx;}
  MatrixEntryType& grady  ()  { return gy;}
  MatrixEntryType& gradz  ()  { return gz;}

  void operator  = (const CFDBlock3d&);
  void operator *= (double);
  void operator *= (const CFDBlock3d&);
  void operator += (const CFDBlock3d&);
  void operator -= (const CFDBlock3d&);

  MatrixEntryType  operator()(int i,int j) const;
  MatrixEntryType& diag(int i);

  void zero();
  void transpose();
  void transpose(CFDBlock3d& A);

  void copy_transpose(const CFDBlock3d& A) 
    { std::cerr << "CFDBlock3d::copy_transpose noch schreiben" << std::endl; }

  void   getrow   (std::vector<double>& v, int i) {abort();}
  void   getcolumn(std::vector<double>& v, int i) {abort();}
  void   setrow   (std::vector<double>& v, int i) {abort();}
  void   setcolumn(std::vector<double>& v, int i) {abort();}
  
  void   DirichletRow (const std::vector<int>& cv);
  void   DirichletCol (const std::vector<int>& cv);
  void   DirichletDiag(const std::vector<int>& cv);

  void   entry     (const nmatrix<double>&);
  void   entry     (int i, int j, const EntryMatrix&, double s=1.);
  void   dual_entry(int i, int j, const EntryMatrix&, double s=1.);
  void   inverse ();
  void   submult (const CFDBlock3d&, const CFDBlock3d&); 
  void   vmult   (iterator) const;
  void   mult    (CFDBlock3d&, const CFDBlock3d&) const; 

  void add(double s, const CFDBlock3d& A);
  void add(double s, const TimePattern& TP);
  void adddiag(const nvector<double>& s, double l);
 
  void cadd(double d, iterator p, const_iterator q0) const
    {
      // p and q are the vector offsets

      double a = s  * *(q0) + dx * *(q0+1) + dy * *(q0+2) + dz * *(q0+3);
      *p++ += d*a;

      a = gradx() * *q0  + laplx()* *(q0+1);
      *p++ += d*a;
      
      a = grady() * *q0  + laply()* *(q0+2);
      *p++ += d*a;

      a = gradz() * *q0  + laplz()* *(q0+3);
      *p += d*a;

      p -= 3;
    }

  void caddtrans(double s, iterator p, const_iterator q0) const
    { 
      std::cerr << "CFDBlock3d::caddtrans noch schreiben" << std::endl;
      abort();
    }
  void subtract(iterator p, const_iterator q0) const;

  std::ostream& print(std::ostream& s) const;

  // Zugriff auf Inhalt ueber ganzen Vektor, damits auch ohne
  // Struktur geht.
  void vector_get(nvector<MatrixEntryType>& v) const
    {
      std::cerr << "\"CFDBlock3d::vector_get\" not written!" << std::endl;
      abort();
    }
  void vector_set(nvector<MatrixEntryType>& v)
    {
      std::cerr << "\"CFDBlock3d::vector_set\" not written!" << std::endl;
      abort();
    }
  void vector_add(double d, nvector<MatrixEntryType>& v)
    {
      std::cerr << "\"CFDBlock3d::vector_add\" not written!" << std::endl;
      abort();
    }

  friend std::ostream& operator<<(std::ostream &s, const CFDBlock3d& A)
    {
      std::cerr << "\"std::ostream& operator<<(std::ostream &s, const CFDBlock3d& A)\" not written!" << std::endl;
      abort();
    }
};
}

#endif


