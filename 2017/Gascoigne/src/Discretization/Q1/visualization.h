/**
*
* Copyright (C) 2004, 2005, 2006 by the Gascoigne 3D authors
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


#ifndef __visualization_h
#define __visualization_h

#include  "meshinterface.h"
#include  "visudata.h"
#include  "visudatainfo.h"
#include  "gnuplot.h"
#include  "paramfile.h"

#include  <stdio.h>
#include  <fstream>

/*--------------------------------------------------*/

namespace Gascoigne
{
class Visualization
{
 protected:

  /* Data */

  const MeshInterface*   mesh;
  const VisuData*        PointData;
  const VisuData*        CellData;
  const VisuDataInfo*    PointDataInfo;
  const VisuDataInfo*    CellDataInfo;

  std::string  filename;
  std::string  stepname;
  std::string  title;

  std::vector<GnuplotData> GP;

  bool    avsa, gmva, vigiea, vua, gnua, teca, vtka, compress;
  bool    b_rotatedvtk;
  int     i_rotatedvtk_slides;
  double  d_rotatedvtk_angle;
  int     pstep, showoutput;
  double  time, tstep, nexttime;

  /* Functions */

  void   BasicInit();

  void   avs        (const std::string&)  const;
  void   vu         (const std::string&)  const;
  void   gnuplot    (const std::string&)  const;
  void   gmv        (const std::string&)  const;
  void   vtk        (const std::string&)  const;
  void   rotatedvtk (const std::string&)  const;

  void   output_vertexs             (std::ofstream&) const;
  void   output_quads               (std::ofstream&, const std::string& s="") const;
  void   output_hexs                (std::ofstream&, const std::string& s="") const;
  void   output_vertexs_by_component(std::ofstream&) const;

  int    CheckPointData() const;
  int    CheckCellData() const;
  
  void output_solution(std::ofstream&, int) const;


  virtual void _vtk_pointdata(std::ofstream& out) const;
  virtual void _vtk_celldata(std::ofstream& out) const;
  virtual void _vtk_points(std::ofstream& out) const;
  virtual void _vtk_cells(std::ofstream& out) const;


  virtual void _rotatedvtk_pointdata(std::ofstream& out) const;
  virtual void _rotatedvtk_celldata(std::ofstream& out) const;
  virtual void _rotatedvtk_points(std::ofstream& out) const;
  virtual void _rotatedvtk_cells(std::ofstream& out) const;


 public:

  void SetPointData(const VisuData* VD) {
    assert(VD);
    PointData = VD;
  }
  void SetCellData(const VisuData* VD) {
    assert(VD);
    CellData = VD;
  }
  void SetPointDataInfo(const VisuDataInfo* VDI) {
    assert(VDI);
    PointDataInfo = VDI;
  }
  void SetCellDataInfo(const VisuDataInfo* VDI) {
    assert(VDI);
    CellDataInfo = VDI;
  }

  void NoOutput() {showoutput=0;}

  Visualization();
  virtual ~Visualization();
  Visualization(const Visualization&);
  Visualization& operator=(const Visualization&);

  void SetMesh(const MeshInterface* M) {mesh = M;}

  /* Functions */

  const std::string& get_name() const { return stepname;}
  void set_name       (const std::string&);

  void read_parameters(const ParamFile* pf);

  void set_time(double t)       
    {
      char ctitle[50];
      sprintf(ctitle,"Time: %6.4e",t); 
      title = ctitle; 
      time = t;
    }
  void set_pstep(int i)         { pstep = i; }
  void set_tstep(double t)      { pstep = -1; tstep = t; nexttime += t;}
  void set_gnuplotdata(const std::vector<GnuplotData>& gp) { GP=gp; }

  void  step   (int);
  int   active (int) const;
  void  format(const std::string&);

  void write();
};
}

#endif
