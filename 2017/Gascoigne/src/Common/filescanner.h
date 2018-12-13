/**
*
* Copyright (C) 2004, 2007 by the Gascoigne 3D authors
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


#ifndef __filescanner_h
#define __filescanner_h

#include  "dataformathandler.h"
#include  "paramfile.h"

/***************************************************/

namespace Gascoigne
{
class FileScanner
{
  DataFormatHandler& DH;
  std::string        blocksymbol;
  bool               complain;


  void FormatToValue(const std::vector<std::string>& words);
  void print(const std::string& blockname) const;
  void _assert(bool b, const std::vector<std::string>& words) const;
  
public:
  
  int          _i_defaultvalues_level;
  int          _i_defaultvalues_save_all_to_file;
  std::string  _s_defaultvalues_save_filename;

  FileScanner(DataFormatHandler& D, const ParamFile* pf, const std::string& b="");
  FileScanner(DataFormatHandler& D);
  void readfile(const ParamFile* pf, const std::string& blockname);
  void NoComplain() { complain=0; }
};
}

/***************************************************/

#endif
