/**
*
* Copyright (C) 2004, 2006, 2009 by the Gascoigne 3D authors
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


#ifndef  __boundarymanager_h
#define  __boundarymanager_h

#include  <map>
#include  <string>
#include  "stlio.h"
#include  "gascoigne.h"
#include  "paramfile.h"


namespace Gascoigne
{

//////////////////////////////////////////////
//
///@brief
/// Management of boundary colors

/// According to the input mesh, parts of the boundary are associated 
/// with a number (color). This number is used in different application 
/// classes to identify a certain part of the boundary.
/// The class BoundaryManager administrates a list of numbers for 
/// Dirichlet and Neumann colors.
//
//////////////////////////////////////////////

class BoundaryManager
{
 protected:

  IntSet                   _colsDirichlet, _colsRightHandSide, _colsEquation, _colsFunctional;
  IntVector                _colsPeriodic;
  std::map<int,IntVector>  _compsDirichlet;
  std::map<int,IntVector>  _compsPeriodic;


 public:

  BoundaryManager() {}
  virtual ~BoundaryManager() {}

  virtual void BasicInit(const ParamFile* pf);

  virtual std::string GetName() const {return "Std";}

  void AddDirichletData(int col, int c)    
    {
      _colsDirichlet.insert(col);
      _compsDirichlet[col].push_back(c);
    }
  void AddPeriodicData(int col, int c)
    {
      _compsPeriodic[col].push_back(c);
    }
  void AddBoundaryRightHandSide(int col)    
    {
      _colsRightHandSide.insert(col);
    }
  void AddBoundaryEquation(int col)    
    {
      _colsEquation.insert(col);
    }
  void AddBoundaryFunctional(int col)    
    {
      _colsFunctional.insert(col);
    }

  std::ostream& print(std::ostream& s) const;

  virtual const IntSet& GetBoundaryRightHandSideColors() const { return _colsRightHandSide; }
  virtual const IntSet& GetBoundaryEquationColors     () const { return _colsEquation; }
  virtual const IntSet& GetBoundaryFunctionalColors   () const { return _colsFunctional; }
  virtual const IntSet& GetDirichletDataColors        () const { return _colsDirichlet; }
  virtual const IntVector& GetPeriodicDataColors      () const { return _colsPeriodic; }

  virtual const IntVector& GetDirichletDataComponents(int c) const 
    {
      std::map<int,IntVector>::const_iterator p = _compsDirichlet.find(c);
      if(p==_compsDirichlet.end())
	{
	  std::cerr << "BoundaryManager::GetDirichletComponents()" << std::endl;
	  std::cerr << "No such color " << c <<std::endl;
	  std::cerr << "components = " << _compsDirichlet << std::endl;
	  abort();
	}
      return p->second;
    }

  virtual const IntVector& GetPeriodicDataComponents(int c) const
    {
      std::map<int,IntVector>::const_iterator p = _compsPeriodic.find(c);
      if(p==_compsPeriodic.end())
	{
	  std::cerr << "BoundaryManager::GetPeriodicComponents()" << std::endl;
	  std::cerr << "No such color " << c <<std::endl;
	  std::cerr << "components = " << _compsPeriodic << std::endl;
	  abort();
	}
      return p->second;
    }

};
}

#endif
