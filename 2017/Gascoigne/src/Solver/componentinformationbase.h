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


#ifndef  __ComponentInformationBase_h
#define  __ComponentInformationBase_h

#include  "gostream.h"
#include  "stringutil.h"
#include  "filescanner.h"

#include  "componentinformation.h"

namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments ComponentInformationBase

  ///
  ///
  /////////////////////////////////////////////

  class ComponentInformationBase : public ComponentInformation
  {
    private:
      
    protected:

    public:
      ComponentInformationBase():ComponentInformation() {}
      virtual ~ComponentInformationBase() {}
  
      virtual void BasicInit(const ParamFile* pf) {}

      virtual std::string GetName() const;


      virtual const int GetNScalars     () const;
      virtual void      GetScalarName   (int i, std::string& s_name) const;
      virtual const int GetNVectors     () const;
      virtual void      GetVectorName   (int i, std::string& s_name) const;
      virtual void      GetVectorIndices(int i, fixarray<3,int>& fa_vectorindices) const;
  };
}

#endif // __ComponentInformationBase_h
