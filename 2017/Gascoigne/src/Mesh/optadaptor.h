/**
*
* Copyright (C) 2004 by the Gascoigne 3D authors
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


#ifndef __optadaptor_h
#define __optadaptor_h

#include "gascoigne.h"
#include "adaptordata.h"

/*********************************************************************/

namespace Gascoigne
{
class OptAdaptor
{
 protected:

  int d,p,p2,p4,n_aimed;
  int refined, double_refined, coarsened, marge, used;

  double co, dd, pp, factor;

  AdaptorData&             info;
  const DoubleVector&   vol;
        DoubleVector&   eta;

  void prepare();

public:

  OptAdaptor  (AdaptorData&, DoubleVector&, const DoubleVector&);

  void refine (IntVector&);
  void coarse (IntVector&);
  void RefineGnuplot (IntVector&);
};
}

/*********************************************************************/

#endif
