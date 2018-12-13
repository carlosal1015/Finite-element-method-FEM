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


#ifndef __backup_h
#define __backup_h

#include  "gascoigne.h"
#include  "compvector.h"
#include  <string>

/********************************************************************/

namespace Gascoigne
{
class WriteBackUp
{
 public:

  WriteBackUp(const GlobalVector&, const std::string&);
};

/********************************************************************/

class WriteBackUpBinary
{
 public:

  WriteBackUpBinary(const GlobalVector&, const std::string&);
};

/********************************************************************/

class ReadBackUp
{
 public:

  ReadBackUp(const std::string&, int&, int&);
  ReadBackUp(GlobalVector&, const std::string&);
};

/********************************************************************/

class ReadBackUpResize
{
 public:

  ReadBackUpResize(GlobalVector&, const std::string&);
};

/********************************************************************/

class ReadBackUpBinary
{
 public:

  ReadBackUpBinary(GlobalVector&, const std::string&);
};
}

/********************************************************************/

#endif
