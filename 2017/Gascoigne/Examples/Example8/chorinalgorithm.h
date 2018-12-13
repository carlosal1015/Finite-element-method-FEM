/**
*
* Copyright (C) 2008, 2009 by the Gascoigne 3D authors
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


#ifndef  __ChorinAlgorithm_h
#define  __ChorinAlgorithm_h

#include  "splittingsolver.h"
#include  "multilevelsolver.h"
#include  "nonstationaryalgorithm.h"

/*-----------------------------------------*/

namespace Gascoigne
{

//////////////////////////////////////////////
//
///@brief
///
///
//////////////////////////////////////////////

class ChorinAlgorithm : public NonstationaryAlgorithm
{
 protected:

  int         niter;
  std::string initial, scheme, pproblem, vproblem;

  const SplittingSolver* GetSplittingSolver() const 
    { 
      const SplittingSolver* S = dynamic_cast<const SplittingSolver*>(GetSolver());
      assert(S);
      return S;
    }
       SplittingSolver* GetSplittingSolver()
    { 
      SplittingSolver* S = dynamic_cast<SplittingSolver*>(GetSolver());
      assert(S);
      return S;
    }
  void Chorin();
  void ChorinUzawa();
  void VanKan();
  void VelocityPredictor(VectorInterface& v, VectorInterface& fv, NLInfo& nlinfo, int iter);
  void PressurePoissonProblem(VectorInterface& q, VectorInterface& fp, CGInfo& cginfo);
  void VelocityProjection(VectorInterface& v, VectorInterface& q, 
			  VectorInterface& fv, CGInfo& cginfo, int iter);
  void PressureUpdate(VectorInterface& p, VectorInterface& q, int iter);

public:

  ChorinAlgorithm() : NonstationaryAlgorithm(), niter(0),initial(""),pproblem(""),vproblem("")  {}
  virtual ~ChorinAlgorithm() {}

  void Run(const std::string& problem1, const std::string& problem2);
};
}

/*-----------------------------------------*/

#endif
