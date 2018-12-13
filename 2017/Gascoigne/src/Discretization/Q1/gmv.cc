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


#include  "visualization.h"
#include  "errormacros.h"

using namespace std;


/********************************************************************/

namespace Gascoigne
{
void Visualization::gmv(const string& bname) const
{
  string name = bname;
  name += ".gmv";

  ofstream file(name.c_str());
  FILE_ERROR(file,name);
  
  file << "gmvinput ascii" << endl;
  file << "nodes " << mesh->nnodes() << endl;
  
  output_vertexs_by_component(file);

  file << endl << "cells " << mesh->ncells() << endl;

  output_quads(file,"quad 4 ");
  output_hexs (file,"hex 8 ");

  if (PointData)
    {
      CheckPointData();
      file << "variable" << endl;
      for(VisuDataInfo::siterator p=PointDataInfo->sbegin();p!=PointDataInfo->send();++p)
	{
	  file << p->first <<" 1" << endl;
	  output_solution(file,p->second);
	  file << endl;
	}
     file << "endvars" << endl;
    }
  file << "endgmv" << endl;
  file.close();
 if (showoutput) cout << "[" << name << "]\n";

 if(compress)
 {
   string command = "gzip -f " + name;
   system(command.c_str());
 }
}
}
