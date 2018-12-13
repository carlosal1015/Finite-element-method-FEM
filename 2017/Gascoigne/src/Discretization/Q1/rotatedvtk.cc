/**
*
* Copyright (C) 2006 by the Gascoigne 3D authors
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

/* ----------------------------------------- */

namespace Gascoigne
{
void Visualization::_rotatedvtk_pointdata(ofstream& out) const
{
  if (PointData) {
    int nn = mesh->nnodes()*i_rotatedvtk_slides;
 
    CheckPointData();
    //for(VisuDataInfo::siterator p=PointDataInfo->sbegin();p!=PointDataInfo->send();++p)
    for(int i=0; i< PointDataInfo->nscalars() ; i++) {
      VisuDataInfo::siterator p = (const_cast<VisuDataInfo*>(PointDataInfo))->GetSIterator(i);
      if(i==0) out << "POINT_DATA " << nn << endl;
      out << "SCALARS "<< p->first <<" DOUBLE "<< endl;
      out << "LOOKUP_TABLE default"<< endl;
      for (int i_slide=0; i_slide<i_rotatedvtk_slides; i_slide++) {
        for (int ind=0; ind<PointData->visun(); ind++) {
          if(mesh->dimension()==2) {
            out << PointData->visudata2(ind,p->second,mesh->vertex2d(ind)) << endl;
          } else {
            abort();
            out << PointData->visudata2(ind,p->second,mesh->vertex3d(ind)) << endl;
          }
        }
      }
      out << endl<< endl;
    }
    //for(VisuDataInfo::viterator p=PointDataInfo->vbegin();p!=PointDataInfo->vend();++p) {
    for(int i=0; i< PointDataInfo->nvectors() ; i++) {
      VisuDataInfo::viterator p = (const_cast<VisuDataInfo*>(PointDataInfo))->GetVIterator(i);
      out << "VECTORS "<< p->first <<" DOUBLE "<< endl;
      for (int i_slide=0; i_slide<i_rotatedvtk_slides; i_slide++) {
        for (int ind=0; ind<PointData->visun(); ind++) {
          for(int ii=0;ii<2;ii++) {
            if(mesh->dimension()==2) {
              abort();
              // out << PointData->visudata2(ind,p->second[ii],mesh->vertex2d(ind)) << " ";
            } else {
              abort();
            }
          }
          if(p->second[2]==-1) {
            out << 0. << " ";
          } else {
            out << PointData->visudata(ind,p->second[2]) << " ";
          }
          out << endl;
        }
      }
      out << endl<< endl;
    }
  }
}

/* ----------------------------------------- */

void Visualization::_rotatedvtk_celldata(ofstream& out) const
{
   if (CellData) {
     CheckCellData();

     // cout << "CellDataInfo->nscalars()" << CellDataInfo->nscalars() << endl;
     //
     // die reihenfolge der elemente per index ist wesentlich, es reicht nicht nur sie
     // per iterator aus CellDataInfo raus zu holen
     //for(VisuDataInfo::siterator p=CellDataInfo->sbegin();p!=CellDataInfo->send();++p){
     for(int i=0; i< CellDataInfo->nscalars() ; i++){
       VisuDataInfo::siterator p = (const_cast<VisuDataInfo*>(CellDataInfo))->GetSIterator(i);
       if(i==0) out << "CELL_DATA " << (i_rotatedvtk_slides-1)*mesh->ncells() << endl;
       //if(p==CellDataInfo->sbegin()) out << "CELL_DATA " << mesh->ncells() << endl;
       out << "SCALARS "<< p->first <<" DOUBLE "<< endl;
       out << "LOOKUP_TABLE default"<< endl;

       for (int i_slide=0; i_slide<i_rotatedvtk_slides; i_slide++) {
         for (int ind=0; ind<CellData->visun(); ind++) {
           out << CellData->visudata(ind,p->second) << endl;
         }
       }
       out << endl<< endl;
     }
     //cout << "CellDataInfo->nvectors()" << CellDataInfo->nvectors() << endl;
     //for(VisuDataInfo::viterator p=CellDataInfo->vbegin();p!=CellDataInfo->vend();++p){
     for(int i=0; i< CellDataInfo->nvectors() ; i++){
       VisuDataInfo::viterator p = (const_cast<VisuDataInfo*>(CellDataInfo))->GetVIterator(i);
       out << "VECTORS "<< p->first <<" DOUBLE "<< endl;
       for (int ind=0; ind<CellData->visun(); ind++) {
         abort();
         for(int ii=0;ii<2;ii++) {
           out << CellData->visudata(ind,p->second[ii]) << " ";
         }
         if(p->second[2]==-1) {
           out << 0. << " ";
         } else {
           out << CellData->visudata(ind,p->second[2]) << " ";
         }
         out << endl;
       }
       out << endl<< endl;
     }
   }
}

/* ----------------------------------------- */

void Visualization::_rotatedvtk_points(ofstream& out) const
{
  int nn = mesh->nnodes();
  out << "POINTS " << nn*i_rotatedvtk_slides << " DOUBLE " << endl;
  if(mesh->dimension()==2)
    { 
      for (int i_slide=0; i_slide<i_rotatedvtk_slides; i_slide++) {
        double d_phi = i_slide*d_rotatedvtk_angle;
        for (int i=0; i<nn; i++) {
          const Vertex2d & v = mesh->vertex2d(i);
          double r = v.x();
          double x = r*cos( 2.*M_PI *d_phi/360.);
          double z = r*sin( 2.*M_PI *d_phi/360.);
          double y = v.y();
          out<<  x << " " << y << " " << z << endl;
        }
      }  
    }
  else if(mesh->dimension()==3)
    { 
      abort();
    }
  else
    {
      abort();
    }
  out << endl;
}

/* ----------------------------------------- */

void Visualization::_rotatedvtk_cells(ofstream& out) const
{
  int ne = mesh->ncells();
  int nn = mesh->nnodes();

  int lenght=0;
  for (int c=0; c<ne; c++)
    {
      lenght += 2*mesh->nodes_per_cell(c)+1;
    }
 
  out << endl << "CELLS " << (i_rotatedvtk_slides-1)*ne <<" " << (i_rotatedvtk_slides-1)*lenght << endl;
  
  for (int i_slide=0; i_slide<i_rotatedvtk_slides-1; i_slide++) {
    for (int c=0; c<ne; c++) {
      int nle = mesh->nodes_per_cell(c);
      out << 2*nle << " ";
      for(int ii=0;ii<nle;ii++) {
        out << i_slide*nn + mesh->vertex_of_cell(c,ii) << " "; 
      }
      for(int ii=0;ii<nle;ii++) {
        out << (i_slide+1)*nn + mesh->vertex_of_cell(c,ii) << " "; 
      }
      out << endl;
    }
  }
  out << endl; 
  out << endl << "CELL_TYPES " << (i_rotatedvtk_slides-1)*ne << endl;
  for (int i_slide=0; i_slide<i_rotatedvtk_slides-1; i_slide++) {
    for (int c=0; c<ne; c++) {
      out << "12" << " "; // VTK Typ 12 = hexaeder
    }
  }  
  out << endl;
}

/* ----------------------------------------- */

void Visualization::rotatedvtk(const string& bname) const
{
  string name = bname;
  name += ".vtk";

  ofstream out(name.c_str());
  FILE_ERROR(out,name);

 //  Header

  out << "# vtk DataFile Version 4.2 "<<endl;
  out << "output from GascoigneStd, "<< title << endl;
  out << "ASCII" << endl;
  out << "DATASET UNSTRUCTURED_GRID" << endl;

 //  Mesh

  _rotatedvtk_points(out);
  _rotatedvtk_cells(out);

 //  Data

  _rotatedvtk_pointdata(out);
  _rotatedvtk_celldata(out);
 
 out.close();

 if(compress)
 {
   string command = "gzip -f " + name;
   system(command.c_str());
 }
}
}
