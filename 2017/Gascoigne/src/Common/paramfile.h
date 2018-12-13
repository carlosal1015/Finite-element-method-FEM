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


#ifndef  __ParamFile_h
#define  __ParamFile_h


#include  <string>

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments ParamFile

////
////
/////////////////////////////////////////////

class ParamFile : public std::string
{
private:


protected:


public:


//
////  Con(De)structor 
//
  
  ParamFile() : std::string() {}
  ParamFile(const std::string& name) : std::string(name) {}
  std::string GetName() const {return *this;}
  void SetName(const std::string& name) { std::string::operator=(name);} 
};
}

#endif
