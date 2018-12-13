/**
*
* Copyright (C) 2005, 2006, 2008 by the Gascoigne 3D authors
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


#ifndef  __VectorInterface_h
#define  __VectorInterface_h

#include  <string>
#include  <iostream>
#include  <cassert>


namespace Gascoigne
{

  /////////////////////////////////////////////
  ////
  ////@brief
  ////  ... comments VectorInterface
  ////
  ////
  /////////////////////////////////////////////

  class VectorInterface : public std::string
  {
    private:
      std::string _type;

    protected:

    public:
      VectorInterface(const std::string& name) : std::string(name), _type("node")  { }
      VectorInterface(const std::string& name, const std::string& type) : std::string(name) {
        assert(type=="node" || type=="cell" || type=="parameter");
        GetType()=type;
      }
      VectorInterface(const VectorInterface& v) {
        SetName(v.GetName());
        SetType(v.GetType());
      }
      virtual ~VectorInterface() { }

      void SetName(const std::string& name) { GetName()=name; }
      std::string& GetName() { return *this; }
      const std::string& GetName() const { return *this; }

      void SetType(const std::string& type) {
        assert(type=="node" || type=="cell" || type=="parameter");
        GetType()=type;
      }
      std::string& GetType() { return _type; }
      const std::string& GetType() const { return _type; }

      friend std::ostream& operator<<(std::ostream& os, const VectorInterface& g) {
        os << "Name: '" << g.GetName() << "' ";
        os << "Type: '" << g.GetType() << "'" << std::endl;
        return os;
      }

  };
}

#endif
