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


#ifndef  __MgInterpolatorInterface_h
#define  __MgInterpolatorInterface_h

#include  "gascoigne.h"

/*--------------------------------------------------------*/

namespace Gascoigne
{
  class MgInterpolatorInterface
  {
    private:

    protected:

    public:
      MgInterpolatorInterface() {}
      virtual ~MgInterpolatorInterface() {}

      virtual void restrict_zero   (GlobalVector&, const GlobalVector&) const=0;
      virtual void prolongate_add  (GlobalVector&, const GlobalVector&) const=0;
      virtual void SolutionTransfer(GlobalVector&, const GlobalVector&) const=0;
      virtual void SolutionTransferUp(GlobalVector& ul, const GlobalVector& uL) const {
        prolongate_add(ul,uL); 
      }
      virtual void Pi(GlobalVector& u) const {
        std::cerr << "\"MgInterpolatorInterface::Pi\" not written!" << std::endl;
        abort();
      }
  };
}

#endif
