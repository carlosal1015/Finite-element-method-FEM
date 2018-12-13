/**
*
* Copyright (C) 2008 by the Gascoigne 3D authors
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


#ifndef  __DwrAlgorithm_h
#define  __DwrAlgorithm_h

#include  "multilevelalgorithm.h"

/*-----------------------------------------*/

namespace Gascoigne
{

//////////////////////////////////////////////
//
///@brief
///
///
//////////////////////////////////////////////

class DwrAlgorithm : public MultiLevelAlgorithm
{
 protected:

  DiscretizationInterface* CreateOtherDiscretization() const;

  void PrimalResidualsHigher(VectorInterface& f, const VectorInterface& u);
  void DualResidualsHigher  (VectorInterface& f, const VectorInterface& u, const VectorInterface& z);

public:

  DwrAlgorithm() :  MultiLevelAlgorithm()  {}
  virtual ~DwrAlgorithm() {}
  void AdaptiveLoop(const std::string& problemlabel, const std::string& duallabel, Functional& J);
};

}

/*-----------------------------------------*/

#endif
