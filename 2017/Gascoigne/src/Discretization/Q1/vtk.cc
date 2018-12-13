/**
*
* Copyright (C) 2004, 2005, 2010, 2011 by the Gascoigne 3D authors
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
void Visualization::_vtk_pointdata(ofstream& out) const
{
  if (PointData) {
    int nn = mesh->nnodes();
 
    CheckPointData();
    //for(VisuDataInfo::siterator p=PointDataInfo->sbegin();p!=PointDataInfo->send();++p)
    for(int i=0; i< PointDataInfo->nscalars() ; i++) {
      VisuDataInfo::siterator p = (const_cast<VisuDataInfo*>(PointDataInfo))->GetSIterator(i);
      if(i==0) out << "POINT_DATA " << nn << endl;
      out << "SCALARS "<< p->first <<" DOUBLE "<< endl;
      out << "LOOKUP_TABLE default"<< endl;
      for (int ind=0; ind<PointData->visun(); ind++) {
        if(mesh->dimension()==2) {
          out << PointData->visudata2(ind,p->second,mesh->vertex2d(ind)) << endl;
        } else {
          out << PointData->visudata2(ind,p->second,mesh->vertex3d(ind)) << endl;
        }
      }
      out << endl<< endl;
    }
    //for(VisuDataInfo::viterator p=PointDataInfo->vbegin();p!=PointDataInfo->vend();++p) {
    for(int i=0; i< PointDataInfo->nvectors() ; i++) {
      VisuDataInfo::viterator p = (const_cast<VisuDataInfo*>(PointDataInfo))->GetVIterator(i);
      out << "VECTORS "<< p->first <<" DOUBLE "<< endl;
      for (int ind=0; ind<PointData->visun(); ind++) {
        for(int ii=0;ii<2;ii++) {
          if(mesh->dimension()==2) {
            out << PointData->visudata2(ind,p->second[ii],mesh->vertex2d(ind)) << " ";
          } else {
            out << PointData->visudata2(ind,p->second[ii],mesh->vertex3d(ind)) << endl;
          }
        }
        if(p->second[2]==-1) {
          out << 0. << " ";
        } else {
          out << PointData->visudata(ind,p->second[2]) << " ";
        }
        out << endl;
      }
      out << endl<< endl;
    }
  }
}

/* ----------------------------------------- */

void Visualization::_vtk_celldata(ofstream& out) const
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
       if(i==0) out << "CELL_DATA " << mesh->ncells() << endl;
       //if(p==CellDataInfo->sbegin()) out << "CELL_DATA " << mesh->ncells() << endl;
       out << "SCALARS "<< p->first <<" DOUBLE "<< endl;
       out << "LOOKUP_TABLE default"<< endl;

       for (int ind=0; ind<CellData->visun(); ind++) {
         out << CellData->visudata(ind,p->second) << endl;
       }
       out << endl<< endl;
     }
     //cout << "CellDataInfo->nvectors()" << CellDataInfo->nvectors() << endl;
     //for(VisuDataInfo::viterator p=CellDataInfo->vbegin();p!=CellDataInfo->vend();++p){
     for(int i=0; i< CellDataInfo->nvectors() ; i++){
       VisuDataInfo::viterator p = (const_cast<VisuDataInfo*>(CellDataInfo))->GetVIterator(i);
       out << "VECTORS "<< p->first <<" DOUBLE "<< endl;
       for (int ind=0; ind<CellData->visun(); ind++) {
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

void Visualization::_vtk_points(ofstream& out) const
{
  int nn = mesh->nnodes();
  out << "POINTS " << nn << " DOUBLE " << endl;

  assert(mesh->dimension()==2 || mesh->dimension()==3);

  if(mesh->dimension()==2)
    { 
      for (int i=0; i<nn; i++)
        {
          out<<  mesh->vertex2d(i) << " " << 0 << endl;
        }
    }
  else if(mesh->dimension()==3)
    { 
      for (int i=0; i<nn; i++)
        {
          out<<  mesh->vertex3d(i) << endl;
        }
    }
  out << endl;
}

/* ----------------------------------------- */

void Visualization::_vtk_cells(ofstream& out) const
{
  int ne = mesh->ncells();

  int lenght=0;
  for (int c=0; c<ne; c++)
    {
      lenght += mesh->nodes_per_cell(c)+1;
    }
 
  out << endl << "CELLS " << ne <<" " << lenght << endl;
  
  for (int c=0; c<ne; c++)
    {
      int nle = mesh->nodes_per_cell(c);
      out << nle << " ";
      for(int ii=0;ii<nle;ii++)
        {
          out << mesh->vertex_of_cell(c,ii) << " "; 
        }
      out << endl; 
    }     
  out << endl << "CELL_TYPES " << ne << endl;
  for (int c=0; c<ne; c++)
    {
      out << mesh->VtkType(c) << " ";
    }
  out << endl;
}

/* ----------------------------------------- */

void Visualization::vtk(const string& bname) const
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
  out << "FIELD FieldData 1" << endl;
    out << "TIME 1 1 double" << endl;
    out << time << endl;

 //  Mesh

  _vtk_points(out);
  _vtk_cells(out);

 //  Data

  _vtk_pointdata(out);
  _vtk_celldata(out);
 
 out.close();

 if(compress)
 {
   string command = "gzip -f " + name;
   system(command.c_str());
 }
}
}
