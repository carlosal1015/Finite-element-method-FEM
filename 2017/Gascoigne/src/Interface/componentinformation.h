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


#ifndef  __ComponentInformation_h
#define  __ComponentInformation_h

#include  "gostream.h"
#include  "stringutil.h"
#include  "filescanner.h"
#include  "application.h"

namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments ComponentInformation

  ///
  ///
  /////////////////////////////////////////////

  class ProblemDescriptorInterface;
  class ComponentInformation : public virtual Application
  {
    private:
      
    protected:
      mutable int _i_dimension;
      ProblemDescriptorInterface* PDI;

    public:
      ComponentInformation() { PDI=NULL;}
      virtual ~ComponentInformation() {}
  
      virtual void BasicInit(const ParamFile* pf) {}

      virtual std::string GetName() const=0;
      virtual int         GetDimension()                const { return _i_dimension;      };
      virtual void        SetDimension(int i_dimension) const { _i_dimension=i_dimension; };
      ProblemDescriptorInterface*& GetProblemDescriptorInterface()       { return PDI;};
      ProblemDescriptorInterface*  GetProblemDescriptorInterface() const { return PDI;};

      const int ncomp   () const { return GetNScalars(); };
      const int GetNcomp() const { return GetNScalars(); };

      virtual const int GetNScalars     () const=0;
      virtual void      GetScalarName   (int i, std::string& s_name) const=0;
      virtual const int GetNVectors     () const=0;
      virtual void      GetVectorName   (int i, std::string& s_name) const=0;
      virtual void      GetVectorIndices(int i, fixarray<3,int>& fa_vectorindices) const=0;
  };
}

#endif // __ComponentInformation_h
