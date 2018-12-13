/**
*
* Copyright (C) 2004, 2005, 2007, 2011 by the Gascoigne 3D authors
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


#include "visu_eps.h"
#include "compose_name.h"
#include <fstream>
#include  "filescanner.h"

using namespace std;

/* -------------------------------------------------- */

namespace Gascoigne
{
VisuEPS::VisuEPS(const ParamFile* paramfile) : M(0)
{
  {
    DataFormatHandler DFH;
    DFH.insert("compress_eps" , &_i_compress, 0);
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(paramfile,"Visualization");
  }

  INTOPT[WRITE_PATCH]  = 0;
  INTOPT[COMBINE_LINES]  = 1;
  
  DOUBLEOPT[LINEWIDTH] = 0.1;
  DOUBLEOPT[WIDTH]     = 512;

  offset.x()=0;
  offset.y()=0;
}

/* -------------------------------------------------- */

void VisuEPS::SetOption(EPSOptions o, int v)
{
  if (INTOPT.find(o)==INTOPT.end())
    {
      if (DOUBLEOPT.find(o)!=DOUBLEOPT.end()) 
	{
	  SetOption(o,static_cast<double>(v));
	  return;
	}
      cerr << "void VisuEPS::SetOption(" << o << "," << v << ")!";
      abort();
    }
  INTOPT[o]=v;
}

/* -------------------------------------------------- */

void VisuEPS::SetOption(EPSOptions o, double v)
{
  assert(DOUBLEOPT.find(o)!=DOUBLEOPT.end());
  DOUBLEOPT[o]=v;
}

/* -------------------------------------------------- */

void VisuEPS::Lexiko(Line& p) const
{
  const Vertex2d& a = M->vertex2d(p.first);
  const Vertex2d& b = M->vertex2d(p.second);
  
  bool sm = (a.y()<b.y());
  if (a.y()==b.y()) sm = (a.x()<b.x());
  if (!sm)
    {
      int t = p.first;
      p.first = p.second;
      p.second = t;
    }
}

/* -------------------------------------------------- */

bool VisuEPS::InLine(int a,int b,int c) const
{
  double det = (M->vertex2d(a)-M->vertex2d(b))[0]*
	       (M->vertex2d(a)-M->vertex2d(c))[1]-
	       (M->vertex2d(a)-M->vertex2d(b))[1]*
	       (M->vertex2d(a)-M->vertex2d(c))[0];
  return (fabs(det)<1.e-13);
}

/* -------------------------------------------------- */

void VisuEPS::WriteGrid(string fname,int iter)
{
  assert(M);

  for (int i=0;i<lines.size();++i) lines[i].clear();
  lines.resize(M->nnodes());
  
  compose_name(fname,iter);
  fname += ".eps";


  n_lines = 0;
  
  double x_max,x_min,y_max,y_min;
  double scale = 1;
  
  x_min=x_max = M->vertex2d(0).x();
  y_min=y_max = M->vertex2d(0).y();
  
  for (int i=0;i<M->nnodes();++i)
    {
      x_min = Gascoigne::min(x_min,M->vertex2d(i).x());
      x_max = Gascoigne::max(x_max,M->vertex2d(i).x());
      y_min = Gascoigne::min(y_min,M->vertex2d(i).y());
      y_max = Gascoigne::max(y_max,M->vertex2d(i).y());
    }

  offset.x()=x_min;
  offset.y()=y_min;

  scale = DOUBLEOPT[WIDTH]/(x_max-x_min);

  if (!INTOPT[WRITE_PATCH])
    {
      for (int i=0;i<M->ncells();++i)
	for (int j=0;j<4;++j)
	  {
	    Line pu(M->vertex_of_cell(i,j),M->vertex_of_cell(i,(j+1)%4));
	    Lexiko(pu);
	    assert(pu.first!=pu.second);
	    assert(pu.first<lines.size());
	    if (lines[pu.first].find(pu.second)==lines[pu.first].end())
	      {
		++n_lines;
		lines[pu.first].insert(pu.second);
	      }
	}
    }
  else 
    {
      vector<int> bu(4);
      bu[0]=0;bu[1]=2;bu[2]=8;bu[3]=6;
      for (int i=0;i<M->npatches();++i)
	{
	  const IntVector& vop = *M->IndicesOfPatch(i);
	  
	  for (int j=0;j<4;++j)
	    {
	      pair<int,int>  pu(vop[bu[j]],vop[bu[(j+1)%4]]);
	      Lexiko(pu);
	      assert(pu.first!=pu.second);
	      assert(pu.first<lines.size());
	      if (lines[pu.first].find(pu.second)==lines[pu.first].end())
		{
		  ++n_lines;
		  lines[pu.first].insert(pu.second);	  
		}
	    }
	}
    }
  // Linien auf einer Gerade zusammenfassen
  if (INTOPT[COMBINE_LINES]) CombineLines();

  ofstream out(fname.c_str());
  cout << "[" << fname << "]" << endl;
  out << "%!PS-Adobe-2.0 EPSF-1.2" << endl
      << "%%Title: Gascoigne OutPut" << endl
      << "%%Creator: Gascoigne" << endl
      << "%%BoundingBox: "
	      // lower left corner
      << "-10 -10 "
	      // upper right corner
      << static_cast<unsigned int>( (x_max-x_min) * scale )+10
      << ' '
      << static_cast<unsigned int>( (y_max-y_min) * scale )+10
      << endl;
  
  // define some abbreviations to keep
  // the output small:
  // m=move turtle to
  // x=execute lineto
  out << "/m {moveto} bind def" << endl
      << "/x {lineto} bind def" << endl;
  
  //      << "/b {0 0 0 setrgbcolor} def" << endl
  //      << "/r {1 0 0 setrgbcolor} def" << endl;
  
  out << "%%EndProlog" << endl << endl;
  
  out << DOUBLEOPT[LINEWIDTH] << " setlinewidth" << endl;
  
  // Linienzuege malen
  for (int i=0;i<lines.size();++i)
    {
      int start = i;
      set<int>::iterator it;

      while (lines[i].size()>0)
	{
	  it = lines[i].begin();
	  int from = i;
	  int next = *it;
	  out << (M->vertex2d(start)-offset)[0]*scale
	      << " "  << (M->vertex2d(start)-offset)[1]*scale << " m";
	  while (lines[from].size()>0)
	    {
	      out << " " << (M->vertex2d(next)-offset)[0]*scale
		  << " " << (M->vertex2d(next)-offset)[1]*scale << " x";
	      it = lines[next].begin();
	      lines[from].erase(next);
	      from = next;
	      next = *it;
	    }
	  out << endl;
	}
    }
  out << "stroke" << endl;
  out.close();

  if(_i_compress) {
    string command = "gzip -f " + fname; 
    system(command.c_str()); 
  }

}

/* ----------------------------------------------------- */

void VisuEPS::CombineLines()
{
  int combine = 0;
  bool changed;
  int steps = 0;
  do 
    {
      ++steps;
      changed = false;
      for (int i=0;i<lines.size();++i)
	{
	  for (set<int>::iterator it = lines[i].begin();it!=lines[i].end();++it)
	    {
	      int next = *it;
	      bool found = false;
	      set<int>::iterator it1;
	      int nn;
	      for (it1 = lines[next].begin();((it1!=lines[next].end())&(!found));++it1)
		{
		  nn = *it1;
		  found = InLine(i,next,nn);
		}	
	      if (found)
		{
		  ++combine;
		  changed = true;
		      
		  lines[i].erase(next);
		  lines[i].insert(nn);
		  lines[next].erase(nn);
		  it = lines[i].begin();
		}
	    }
	}
    }
  while(changed);
//   cerr << "Combined: " << combine << " of " << n_lines << " lines in " << steps << " steps.\n";
}
}
