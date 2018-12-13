/**
*
* Copyright (C) 2006 by the Gascoigne 3D authors
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


#ifndef __DwrQ2Q4_h
#define __DwrQ2Q4_h

#include "solverinterface.h"

namespace Gascoigne
{

/**********************************************************/

  class DwrQ2Q4
  {
    private:
      SolverInterface                  &_S;
      const ProblemDescriptorInterface *_P;
      DiscretizationInterface          *_D;

    protected:
      DiscretizationInterface* CreateOtherDiscretization() const;

      double ScalarProduct(DoubleVector &eta, const GlobalVector &f, const GlobalVector &z) const;
      double ScalarProduct(DoubleVector &eta, const VectorInterface &gf, const VectorInterface &gz) const;
      double ScalarProductWithFluctuations(DoubleVector& eta, const VectorInterface &gf, const VectorInterface &gz) const;

      void PrimalResidualsHigher(VectorInterface &gf, const VectorInterface &gu);

      void DualResidualsHigher(VectorInterface &gf, const VectorInterface &gu, const VectorInterface &gz, const ProblemDescriptorInterface &PDI);

    public:
      DwrQ2Q4(SolverInterface &S);
      ~DwrQ2Q4() { }

      double Estimator(DoubleVector &eta, VectorInterface &gf, const VectorInterface &gu, const VectorInterface &gz, const ProblemDescriptorInterface &PDI);
  };

  /**********************************************************/

}

#endif
