/**
*
* Copyright (C) 2004, 2011 by the Gascoigne 3D authors
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


#include  "gascoignemeshtransferconstructor.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
GascoigneMeshTransferConstructor2d::GascoigneMeshTransferConstructor2d
(const HierarchicalMesh2d* HM, GascoigneMeshTransfer* GMT,
 const LevelMesh2d* LMfine, const LevelMesh2d* LMcoarse)
{
  IntVector& c2f = GMT->GetC2f();
  map<int,fixarray<2,int> >& zweier = GMT->GetZweier();
  map<int,fixarray<4,int> >& vierer = GMT->GetVierer();
  map<int,int>&              CellEiner = GMT->GetCellEiner();
  map<int,fixarray<4,int> >& CellVierer = GMT->GetCellVierer();

  const QuadLawAndOrder& LaO = HM->QuadLawOrder();

  CellEiner.clear();
  CellVierer.clear();
  // Zellen
  for(int i=0;i<LMcoarse->ncells();i++)
    {
      int igq = LMcoarse->Quadl2g(i);
      if(LMfine->Quadg2l().find(igq)!=LMfine->Quadg2l().end())
	{
	  CellEiner[i] = LMfine->Quadg2l(igq);
	}
      else
	{
	  // kinder suchen
	  fixarray<4,int>  n4;
	  for(int ii=0;ii<4;ii++)
	    {
	      int ic = HM->quad(igq).child(ii);
	      n4[ii] = LMfine->Quadg2l(ic);
	    }
	  CellVierer[i] = n4;
	}
    }



  // 1er
  c2f.reservesize(LMcoarse->nnodes());
  for(int iL=0;iL<LMcoarse->nnodes();iL++)
    {
      int ig = LMcoarse->Vertexl2g(iL);
      int il = LMfine->Vertexg2l(ig);

      c2f[iL] = il;
    }

  // 2er 4er 8er (!)
  for(int i=0;i<LMcoarse->ncells();i++)
    {
      int igq = LMcoarse->Quadl2g(i);
      if(LMfine->Quadg2l().find(igq)!=LMfine->Quadg2l().end()) continue;
      const Quad& q = HM->quad(igq);

      // verfeinertes quad --> vierer

      fixarray<4,int>  n4;
      int igm = LaO.middle_vertex(q);
      int ilm = LMfine->Vertexg2l(igm);
      for(int ii=0;ii<4;ii++)
	{
	  int ig = q.vertex(ii);
	  int iL = LMcoarse->Vertexg2l(ig);
	  n4[ii] = iL;
	}
      vierer[ilm]=n4;


      // edges
      for(int ie=0;ie<4;ie++)
	{
	  int ige = LaO.edge_vertex(q,ie);
	  int ile = LMfine->Vertexg2l(ige);
	  if(LMcoarse->Vertexg2lCheck(ige)!=-2) continue;

	  fixarray<2,int> f;
	  LaO.globalvertices_of_edge(q,f,ie);


	  fixarray<2,int>  n2;
	  for(int ii=0;ii<2;ii++)
	    {
	      int iL = LMcoarse->Vertexg2l(f[ii]);
	      n2[ii] = iL;
	    }
	  zweier[ile]=n2;
	}
    }
//   ofstream file("MGI",ios::app);
//   cout << "Matrix:\n" << M << endl;
}

/*-----------------------------------------*/

GascoigneMeshTransferConstructor3d::GascoigneMeshTransferConstructor3d
(const HierarchicalMesh3d* HM, GascoigneMeshTransfer* GMT,
 const LevelMesh3d* LMfine, const LevelMesh3d* LMcoarse)
{
//   cerr << "GascoigneMeshTransferConstructor::Construct3d()\n";
//   cerr << "noch keine konstanten!\n";
//   abort();

  IntVector& c2f = GMT->GetC2f();
  map<int,fixarray<2,int> >& zweier = GMT->GetZweier();
  map<int,fixarray<4,int> >& vierer = GMT->GetVierer();
  map<int,fixarray<8,int> >& achter = GMT->GetAchter();

  const HexLawAndOrder& LaO = HM->HexLawOrder();

  // 1er
  c2f.reservesize(LMcoarse->nnodes());
  for(int iL=0;iL<LMcoarse->nnodes();iL++)
    {
      int ig = LMcoarse->Vertexl2g(iL);
      int il = LMfine->Vertexg2l(ig);
      c2f[iL] = il;
    }
  // 2er 4er 8er (!)
  for(int i=0;i<LMcoarse->ncells();i++)
    {
      int igq = LMcoarse->Hexl2g(i);
      if(LMfine->Hexg2l().find(igq)!=LMfine->Hexg2l().end()) continue;
      const Hex& q = HM->hex(igq);

      // verfeinertes hex

      fixarray<8,int>  n8;
      int igm = LaO.middle_vertex(q);
      int ilm = LMfine->Vertexg2l(igm);
      for(int ii=0;ii<8;ii++)
	{
	  int ig = q.vertex(ii);
	  // 8er
	  n8[ii] = LMcoarse->Vertexg2l(ig);
	}
      achter[ilm]=n8;
      // faces
      for(int ie=0;ie<6;ie++)
	{
	  int ige = LaO.face_vertex(q,ie);
	  int ile = LMfine->Vertexg2l(ige);
	  if(LMcoarse->Vertexg2lCheck(ige)!=-2) continue;
	  fixarray<4,int> f,n4;
	  LaO.globalvertices_of_face(q,f,ie);
	  for(int ii=0;ii<4;ii++)
	    {
	      // 4er
	      n4[ii] = LMcoarse->Vertexg2l(f[ii]);
	    }
	  vierer[ile]=n4;
	}
      //edges
      for(int ie=0; ie<12; ie++)
	{
	  int ige = LaO.edge_vertex(q,ie);
	  int ile = LMfine->Vertexg2l(ige);
	  if(LMcoarse->Vertexg2lCheck(ige)!=-2) continue;
	  fixarray<2,int> f,n2;
	  LaO.globalvertices_of_edge(q,f,ie);

	  n2[0] = LMcoarse->Vertexg2l(f[0]);
	  n2[1] = LMcoarse->Vertexg2l(f[1]);

	  zweier[ile]=n2;
	}
    }
}
}
