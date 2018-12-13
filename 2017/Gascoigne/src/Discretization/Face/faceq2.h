/**
*
* Copyright (C) 2007 by the Gascoigne 3D authors
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


/*----------------------------   faceq2.h     ---------------------------*/
/*      $Id$                 */
#ifndef __faceq2_H
#define __faceq2_H
/*----------------------------   faceq2.h     ---------------------------*/


#include  "facediscretization.h"

namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments DiscretizationInterface

  ///
  ///
  /////////////////////////////////////////////

  template<int DIM>
    class FaceQ2 : public virtual FaceDiscretization
    {
    protected:

      void TransformationFace(FemInterface::Matrix& T1,FemInterface::Matrix& T2,int f) const;
	    
    public:
      FaceQ2() : FaceDiscretization() {}
      virtual ~FaceQ2() {}
      
      //
      //// Functions called from the Solver
      //
      virtual std::string GetName() const { return "Face Q2"; }
      
      virtual void BasicInit(const ParamFile* pf);
      virtual void ReInit   (const MeshInterface* M);
      virtual void build_faces();
    };
}



/*----------------------------   faceq2.h     ---------------------------*/
/* end of #ifndef __faceq2_H */
#endif
/*----------------------------   faceq2.h     ---------------------------*/
