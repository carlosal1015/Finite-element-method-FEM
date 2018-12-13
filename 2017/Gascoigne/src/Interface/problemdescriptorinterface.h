/**
*
* Copyright (C) 2004, 2005, 2006, 2007, 2009 by the Gascoigne 3D authors
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


#ifndef  __ProblemDescriptorInterface_h
#define  __ProblemDescriptorInterface_h

#include  "gostream.h"
#include  "stringutil.h"
#include  "filescanner.h"

#include  "boundarymanager.h"
#include  "equation.h"
#include  "faceequation.h"
#include  "dirichletdata.h"
#include  "periodicdata.h"
#include  "boundaryrighthandside.h"
#include  "boundaryequation.h"
#include  "exactsolution.h"
#include  "boundarymanager.h"
#include  "componentinformation.h"


namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments ProblemDescriptorInterface

  ///
  ///
  /////////////////////////////////////////////

  class ProblemDescriptorInterface
  {
    private:
      
    protected:

    public:
      ProblemDescriptorInterface() {}
      virtual ~ProblemDescriptorInterface() {}
  
      virtual void BasicInit(const ParamFile* pf) {}

      virtual std::string GetName() const=0;
      virtual std::ostream& OutputSettings(std::ostream& os) const=0;
      virtual void SetTime(double time, double dt) const=0;

      virtual const ParamFile* GetParamFile() const=0;

      virtual const Application*               GetRightHandSide           () const=0;
      virtual const BoundaryRightHandSide*     GetBoundaryRightHandSide   () const=0;
      virtual const Equation*                  GetEquation                () const=0;
      virtual const FaceEquation*              GetFaceEquation            () const=0;
      virtual const BoundaryEquation*          GetBoundaryEquation        () const=0;
      virtual const DirichletData*             GetDirichletData           () const=0;
      virtual const PeriodicData*              GetPeriodicData            () const=0;
      virtual const Application*               GetInitialCondition        () const=0;
      virtual const BoundaryInitialCondition*  GetBoundaryInitialCondition() const=0;
      virtual const ExactSolution*             GetExactSolution           () const=0;
      virtual const BoundaryManager*           GetBoundaryManager         () const=0;
      virtual const ComponentInformation*      GetComponentInformation    () const=0;
  };
}

#endif
