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


#ifndef __visu_eps_h
#define __visu_eps_h

#include "patchmesh.h"
#include "string"
#include <map>
#include  "gascoigne.h"

/*-------------------------------------------------------------------------*/

namespace Gascoigne
{
class VisuEPS
{
  protected:  

  typedef std::pair<int,int> Line;
  typedef nvector<IntSet> Lines;

  const PatchMesh* M;
  int _i_compress;

  Lines lines;
  int   n_lines;
  Vertex2d offset;
  
  // EPS optionen
  std::map<int,int>    INTOPT;
  std::map<int,double> DOUBLEOPT;

  // sort line p (left,bottom) first
  void Lexiko(Line& p) const;

  // test if vertices a,b,c are aligned straightly
  bool InLine(int a,int b,int c) const;
  
  void CombineLines();
    
  public:

  VisuEPS(const ParamFile* paramfile);

  /**
   * Options for output:
   *
   * WRITE_PATCH:
   *   0 : write cells  (default)
   *   1 : write patchs
   *
   * LINEWIDTH:
   *   width of lines in pt. (0.1)
   *
   * WIDTH:
   *   horizontal size of output (int pt) (300)
   *
   * COMBINE_LINES
   *   1: straightly aligned lines are combined to one (def)
   *
   **/
  enum EPSOptions { WRITE_PATCH, LINEWIDTH, WIDTH, COMBINE_LINES };
  
  void SetOption(EPSOptions o, int v);
  void SetOption(EPSOptions o, double v);
  
  void SetMesh(const MeshInterface& PM)  { 
    const PatchMesh* PMP = dynamic_cast<const PatchMesh*>(&PM);
    assert(PMP);
    M = PMP;}
  void WriteGrid(std::string fname, int iter);
};
}

#endif
