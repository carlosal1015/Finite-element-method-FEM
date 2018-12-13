/**
*
* Copyright (C) 2004, 2005, 2011 by the Gascoigne 3D authors
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
#include  "compvector.h"
#include  "compareclass.h"
#include  "gnuplot.h"

#include  <algorithm>
#include  "giota.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
GnuplotData::GnuplotData(const string& s, const Vertex3d& _pos) : 
  plane(s), pos(_pos)
{
  if ((plane!="x")  && (plane!="y")  && (plane!="z") &&
      (plane!="xz") && (plane!="yz") && (plane!="xy"))
    {
      cerr << "GnuplotData::GnuplotData() plane=" << plane << endl;
      abort();
    }
}
  
/*-----------------------------------------*/

GnuplotData::GnuplotData(const GnuplotData& GP) : 
  plane(GP.plane), 
  pos(GP.pos) 
{}

/*-----------------------------------------*/

void GnuplotData::SetName(string& filename) const
{
  if (plane=="x") filename += "_x";
  if (plane=="y") filename += "_y";
  if (plane=="z") filename += "_z";
  if (plane=="xz") filename += "_xz";
  if (plane=="yz") filename += "_yz";
  if (plane=="xy") filename += "_xy";
}

/*-----------------------------------------*/

bool GnuplotData::TestVertex(const Vertex2d& v) const
{
  if ((plane=="x") && (v.x()==pos.x())) return 1;
  if ((plane=="y") && (v.y()==pos.y())) return 1;
  return 0;
}

bool GnuplotData::TestVertex(const Vertex3d& v) const
{
  if ((plane=="xz") && (v.x()==pos.x()) && (v.z()==pos.z())) return 1;
  if ((plane=="yz") && (v.y()==pos.y()) && (v.z()==pos.z())) return 1;
  if ((plane=="xy") && (v.x()==pos.x()) && (v.y()==pos.y())) return 1;
  return 0;
}

/*-----------------------------------------*/

double GnuplotData::SetVertex(const Vertex2d& v) const
{
  assert(plane=="x" || plane=="y");

  if      (plane=="x") return v.y();
  else if (plane=="y") return v.x();
  else 
      return 0.;
}

double GnuplotData::SetVertex(const Vertex3d& v) const
{
  assert(plane=="xz" || plane=="yz" || plane=="xy");

  if      (plane=="xz") return v.y();
  else if (plane=="yz") return v.x();
  else if (plane=="xy") return v.z();
  else 
      return 0.;
}

/********************************************************************/

void Visualization::gnuplot(const string& name) const
{
  if (!PointData) return;

  for (int k=0; k<GP.size(); k++)
    {
      int i = 0;

      if (mesh->dimension()==3)
	{
	  for (int ind=0; ind<mesh->nnodes(); ind++)
	    {
	      const Vertex3d& V = mesh->vertex3d(ind);
	      if (GP[k].TestVertex(V)) i++;
	    }
	}
      else
	{
	  for (int ind=0; ind<mesh->nnodes(); ind++)
	    {
	      const Vertex2d& V = mesh->vertex2d(ind);
	      if (GP[k].TestVertex(V)) i++;
	    }
	}
      assert(i);
      if (!i) continue;

      DoubleVector x(i);
      int comp = PointData->visucomp();
      GlobalVector f(comp);
      f.resize(i);
      i = 0;
      if (mesh->dimension()==3)
	{
	  for (int ind=0; ind<mesh->nnodes(); ind++)
	    {
	      const Vertex3d& V = mesh->vertex3d(ind);
	      if (GP[k].TestVertex(V)) 
		{
		  x[i] = GP[k].SetVertex(V);
		  for (int c=0; c<comp; c++)
		    {
		      f(i,c) = PointData->visudata(ind,c);
		    }
		  i++;
		}
	    }
	}
      else
	{
	  for (int ind=0; ind<mesh->nnodes(); ind++)
	    {
	      const Vertex2d& V = mesh->vertex2d(ind);
	      if (GP[k].TestVertex(V)) 
		{
		  x[i] = GP[k].SetVertex(V);
		  for (int c=0; c<comp; c++)
		    {
		      f(i,c) = PointData->visudata(ind,c);
		    }
		  i++;
		}
	    }
	}


      IntVector C(x.size()); iota(C.begin(),C.end(),0);

      typedef CompareObject<DoubleVector >  CoC;
      
      sort(C.begin(),C.end(),CoC(x));
            
      string gnuname = name;
      
      GP[k].SetName(gnuname);
      gnuname += ".gpl";

      if (showoutput) cout << "[" << gnuname << "]\n";

      ofstream out(gnuname.c_str());
      FILE_ERROR(out,gnuname);
      
      for(int i=0;i<x.size();i++)
	{
	  out << x[C[i]];
	  for (int c=0; c<comp; c++)
	    {
	      out << "  " << f(C[i],c);
	    }
	  out <<  endl;
	}
      out.close();

      if(compress)
      {
        string command = "gzip -f " + gnuname;
        system(command.c_str());
      }
    }
}
}
