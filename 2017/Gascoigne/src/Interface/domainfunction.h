/**
*
* Copyright (C) 2005 by the Gascoigne 3D authors
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


#ifndef __DomainFunction_h
#define __DomainFunction_h

#include "application.h"
#include "vertex.h"

namespace Gascoigne
{

/**********************************************************/

class DomainFunction : public Application
{
  protected:

  public:
    DomainFunction() { }
    virtual ~DomainFunction() { }

    virtual int GetNcomp() const=0;
  
    virtual void F(DoubleVector& f, const Vertex2d &v) const
    {
      std::cerr << "DomainFunction::F not written\n";
      abort();
    }
    virtual void F(DoubleVector& f, const Vertex3d &v) const
    {
      std::cerr << "DomainFunction::F not written\n";
      abort(); 
    }
};

/**********************************************************/
}

#endif
