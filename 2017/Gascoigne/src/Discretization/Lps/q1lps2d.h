/**
*
* Copyright (C) 2004, 2005, 2006, 2009 by the Gascoigne 3D authors
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


#ifndef  __Q1Lps2d_h
#define  __Q1Lps2d_h

#include  "q12d.h"
#include  "q1lpsstab.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments Q1Lps2d

////
////
/////////////////////////////////////////////

/*----------------------------------------------*/

class Q1Lps2d : public virtual Q12d
{
 protected:

  Q1LpsStab2d*   S;

public:

//
////  Con(De)structor 
//

  Q1Lps2d() : Q12d() {}
  ~Q1Lps2d();

  std::string GetName() const {return "Q1Lps2d";}
  
  void BasicInit(const ParamFile* paramfile);
  void ReInit   (const MeshInterface* M);
  void Structure(SparseStructureInterface* SI) const;
  void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
  void Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double d) const;
  void StabForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
  void AddNodeVector(const std::string& name, const GlobalVector* q) const 
    {
      Q12d::AddNodeVector(name,q);
      S->   AddNodeVector(name,q);
    }
  void DeleteNodeVector(const std::string& name) const 
    {
      Q12d::DeleteNodeVector(name);
      S->   DeleteNodeVector(name);
    }
  void AddCellVector(const std::string& name, const GlobalVector* q) const 
    {
      Q12d::AddCellVector(name,q);
      S->   AddCellVector(name,q);
    }
  void DeleteCellVector(const std::string& name) const
    { 
      Q12d::DeleteCellVector(name);
      S->   DeleteCellVector(name);
    }
  void AddParameterVector(const std::string& name, const GlobalParameterVector* q) const 
    {
      Q12d::AddParameterVector(name,q);
      S->   AddParameterVector(name,q);
   }
  void DeleteParameterVector(const std::string& name) const
    {
      Q12d::DeleteParameterVector(name);
      S->   DeleteParameterVector(name);
    }
};

}

#endif
