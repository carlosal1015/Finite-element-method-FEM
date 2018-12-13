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


#ifndef  __Q23d_h
#define  __Q23d_h

#include  "q2.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments Q23d

////
////
/////////////////////////////////////////////

class Q23d : public Q2
{
protected:

  int GetPatchNumber(const Vertex3d& p0, Vertex3d& p) const;
  void VertexTransformation(const Vertex3d& p0, Vertex3d& p, int iq) const;

public:

//
////  Con(De)structor 
//

  Q23d();
  ~Q23d();

  std::string GetName() const {return "Q23d";}
  
  void BasicInit(const ParamFile* paramfile);

  void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* GMT);

  nmatrix<double> GetLocalInterpolationWeights(int iq) const;
};
}

#endif
