/**
*
* Copyright (C) 2004, 2005, 2006, 2011 by the Gascoigne 3D authors
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
#include  <set>
#include  "errormacros.h"
#include  "compose_name.h"
#include  "filescanner.h"


using namespace std;

/********************************************************************/

namespace Gascoigne
{
Visualization::Visualization() 
  : mesh(0), PointData(0), CellData(0),
    filename("none"), GP(0)
{
  BasicInit();
}

/********************************************************************/

Visualization::Visualization(const Visualization& W) : 
  mesh(W.mesh)
{
  *this = W;
}

/********************************************************************/

Visualization::~Visualization()
{}

/********************************************************************/

Visualization& Visualization::operator=(const Visualization& W)
{
  stepname       = W.stepname;
  title          = W.title;
  avsa           = W.avsa;
  gmva           = W.gmva;
  vua            = W.vua;
  vigiea         = W.vigiea;
  teca           = W.teca;
  vtka           = W.vtka;
  b_rotatedvtk         = W.b_rotatedvtk;
  i_rotatedvtk_slides  = W.i_rotatedvtk_slides;
  d_rotatedvtk_angle   = W.d_rotatedvtk_angle;
  //gnua         = W.gnua;
  pstep          = W.pstep;
  time           = W.time;
  tstep          = W.tstep;
  nexttime       = W.nexttime;
  showoutput     = W.showoutput;

  return *this;
}

/********************************************************************/

void Visualization::BasicInit()
{
  pstep = 1;
  avsa = gmva = vua = vigiea = gnua = teca = 0;
  time = 0.; tstep = 0.; nexttime = 0.;
  vtka = 1;
  b_rotatedvtk        = 0;
  i_rotatedvtk_slides = 10;
  d_rotatedvtk_angle  = 36;
  showoutput = 1;
  stepname = "solve";
  compress = false;
}

/********************************************************************/

void Visualization::read_parameters(const ParamFile* pf)
{
  double time;

  vector<string>  planes(0);
  DoubleVector gnupos(0);

  DataFormatHandler DH;
  DH.insert("step"       ,& pstep     ,1);
  DH.insert("tstep"      ,& time       ,0.);
  DH.insert("showoutput" ,& showoutput,1);
  DH.insert("gnuplot"    ,& planes);
  DH.insert("gnuposition",& gnupos);
  DH.insert("vtk"        ,& vtka, 1);
  DH.insert("rotatedvtk" ,& b_rotatedvtk, 0);
  DH.insert("rotatedvtk_slides" ,& i_rotatedvtk_slides, 10);
  DH.insert("rotatedvtk_angle" ,& d_rotatedvtk_angle, 36);
  DH.insert("gmv"        ,& gmva,0);
  DH.insert("vu"         ,& vua,0);
  DH.insert("vigie"      ,& vigiea,0);
  DH.insert("gnu"        ,& gnua,0);
  DH.insert("tecplot"    ,& teca,0);
  DH.insert("avs"        ,& avsa,0);
  DH.insert("compress"   ,& compress, 0);

  FileScanner FS(DH,pf,"Visualization");

  if (time>0.) set_tstep(time);

  // Gnuplot
  
  if (gnupos.size())
    {
      Vertex3d V(gnupos[0],gnupos[1],gnupos[2]);
      vector<GnuplotData> vgp;
      for (int i=0; i<planes.size(); i++)
	{
	  vgp.push_back(GnuplotData(planes[i],V));
	}
      if (vgp.size())
	{
	  set_gnuplotdata(vgp);
	}
    }
}

/********************************************************************/

void Visualization::set_name(const string& s) 
{ 
  stepname = s;
  filename = stepname;
}

/********************************************************************/

void Visualization::format(const string& s)
{
  if      (s=="vtk")            vtka         = 1;
  else if (s=="rotatedvtk")     b_rotatedvtk = 1;
  else if (s=="gmv")            gmva         = 1;
  else if (s=="vu")             vua          = 1;
  else if (s=="vigie")          vigiea       = 1;
  else if (s=="gnu")            gnua         = 1;
  else if (s=="tecplot")        teca         = 1;
  else if (s=="avs")            avsa         = 1;
  else {
    std::cerr << "Wrong format \"" <<s<<"\" in Visualization::format()" << std::endl;
    abort();
  }
}

/********************************************************************/

int Visualization::active(int i) const
{
  if (pstep<0)
    {
      if (time<nexttime) return 0;
    }
  else 
    {
      if ( (pstep==0) || (i%pstep) ) return 0;
    }
  return 1;
}

/********************************************************************/

void Visualization::step(int i)
{
  if ( (pstep==-1) ||  (!active(i)) )
    {
      if (pstep<0) nexttime = time+tstep;

      return;
    }
  filename = stepname;
  compose_name(filename,i);
}

/********************************************************************/

void Visualization::write()
{
  if(mesh==0)
    {
      cerr << "Visualization::write()\n";
      cerr << "mesh pointer not set\n";
      abort();
    }
  if(filename=="none")
    {
      cerr << "Visualization::write()\n";
      cerr << "no filename set [use \"step(i)\"]\n";
      abort();
    }
  if (avsa)         avs(filename);
  if (gmva)         gmv(filename);
  if (vua)          vu(filename);
  if (gnua)         gnuplot(filename);
  if (vtka)         vtk(filename);
  if (b_rotatedvtk) rotatedvtk(filename);
  if (showoutput) cout << "[" << filename << ".vtk]\n";
}

/********************************************************************/

int Visualization::CheckPointData() const
{
  if(PointData==NULL) return 0;

  set<int> comps;

  // scalars
  for(VisuDataInfo::siterator p=PointDataInfo->sbegin();p!=PointDataInfo->send();++p)
    {
      int q = p->second;

      assert(q>=0);
      assert(q<PointData->visucomp());

      comps.insert(q);
    }

  // vectors
  for(VisuDataInfo::viterator p=PointDataInfo->vbegin();p!=PointDataInfo->vend();++p)
    {
      for(int i=0;i<3;i++)
	{
	  int q = p->second[i];

	  assert(q>=-1);
	  assert(q<PointData->visucomp());

	  comps.insert(q);
	}
    }
  return comps.size();
}

/********************************************************************/

int Visualization::CheckCellData() const
{
  if(!CellData) return 0;

  set<int> comps;

  // scalars
  for(VisuDataInfo::siterator p=CellDataInfo->sbegin();p!=CellDataInfo->send();++p)
    {
      if( (p->second<0)||(p->second>=CellData->visucomp()))
	{
	  cerr << "Visualization::CheckCellData()\n";
	  cerr << "scalar does not exist "<<p->second<<endl;
	  exit(1);
	}
      comps.insert(p->second);
    }

  // vectors
  for(VisuDataInfo::viterator p=CellDataInfo->vbegin();p!=CellDataInfo->vend();++p)
    {
      for(int i=0;i<3;i++)
	{
	  if( (p->second[i]<-1)||(p->second[i]>=CellData->visucomp()))
	    {
	      cerr << "Visualization::CheckCellData()\n";
	      cerr << "vector component does not exist "<<p->second[i]<<endl;
	      exit(1);
	    }
	  comps.insert(p->second[i]);
	}
    }
  return comps.size();
}

/********************************************************************/

void Visualization::output_vertexs(ofstream& file) const
{
  if (mesh->dimension()==2)
    {
      for (int i=0; i<mesh->nnodes(); i++)
	{
	  file << mesh->vertex2d(i) << endl;
	}
    }
  else
    {
      for (int i=0; i<mesh->nnodes(); i++)
	{
	  file << mesh->vertex3d(i) << endl;
	}
    }
}

/********************************************************************/

void Visualization::output_vertexs_by_component(ofstream& file) const
{
  if (mesh->dimension()==2)
    {
      for (int i=0; i<mesh->nnodes(); i++)
	{
	  file << " " << mesh->vertex2d(i).x();
	}
      file << endl;
      for (int i=0; i<mesh->nnodes(); i++)
	{
	  file << " " << mesh->vertex2d(i).y();
	}
    }
  else
    {
      for (int i=0; i<mesh->nnodes(); i++)
	{
	  file << " " << mesh->vertex3d(i).x();
	}
      file << endl;
      for (int i=0; i<mesh->nnodes(); i++)
	{
	  file << " " << mesh->vertex3d(i).y();
	}
      file << endl;
      for (int i=0; i<mesh->nnodes(); i++)
	{
	  file << " " << mesh->vertex3d(i).z();
	}
    }
}

/********************************************************************/

void Visualization::output_quads(ofstream& file, const string& text) const
{
  if (mesh->dimension()==2)
    {
      for (int c=0; c<mesh->ncells(); c++)
	{
	  file << text;
	  for(int i=0;i<4;i++)
	    {
	      file << mesh->vertex_of_cell(c,i)+1 << " "; 
	    }
	  file<<endl;      
	}
    }
}

/********************************************************************/

void Visualization::output_hexs(ofstream& file, const string& text) const
{
  if (mesh->dimension()==3)
    {
      for (int c=0; c<mesh->ncells(); c++)
	{
	  file << text;
	  for(int i=0;i<8;i++)
	    {
	      file << mesh->vertex_of_cell(c,i)+1 << " "; 
	    }
	  file<<endl;      
	}
    }
}

/********************************************************************/

void Visualization::output_solution(ofstream& file, int c) const
{
  for (int ind=0; ind<PointData->visun(); ind++)
    {
      file << PointData->visudata(ind,c) << " ";
    }
}
}
