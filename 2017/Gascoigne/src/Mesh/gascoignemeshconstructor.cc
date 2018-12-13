/**
*
* Copyright (C) 2004, 2005, 2006, 2007 by the Gascoigne 3D authors
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


#include  "gascoignemeshconstructor.h"
#include  "gascoignemeshtransferconstructor.h"

using namespace std;

/*---------------------------------------------------*/

namespace Gascoigne
{
GascoigneMeshConstructor::GascoigneMeshConstructor
(const HierarchicalMesh* hm, GascoigneMultiGridMesh* gmg)
  : HM(hm), GMG(gmg), finestlevel(0)
{
}

/*---------------------------------------------------*/

void GascoigneMeshConstructor::BasicInit()
{
  if (HM->dimension()==2) Loop2d();
  else                    Loop3d();

  IntVector& v1 = *GMG->GetGascoigneMesh(0)->Vertexo2n();
  v1.reservesize(HM->Vertexo2n()->size());
  v1 = *HM->Vertexo2n();
}

/*---------------------------------------------------*/

void GascoigneMeshConstructor::Loop2d()
{
  finestlevel = 1;

  set<int>  newquads, oldquads;
  HM->GetAwakeCells(newquads);
 
  GascoigneMesh*  GM0 = GMG->GetGascoigneMesh(0);

  LevelMesh2d* FM=NULL;
  LevelMesh2d* CM=NULL;

  FM = LevelUpdate2d(GM0,newquads,oldquads);

  // Q4-PatchStruktur
  if(HM->patchdepth()>=2)
  {
    PatchIndexHandler &PIH = GM0->GetPatchIndexHandler();
    FM->ConstructCellIndOfPatch(_pl2g);
    PIH.GetHasQ4Patch() = 1;
    HangingIndexHandler &H = GM0->GetHangingIndexHandler();
    FM->ConstructHangingStructureQuartic(*(H.GetQ4Structure()));
  }

  //Wir brauchen local2global und global2local Vector aus FM
  _cl2g = FM->Quadl2g();
  _cg2l = FM->Quadg2l();
  const HierarchicalMesh2d* HMP = dynamic_cast<const HierarchicalMesh2d*>(HM);
  assert(HMP);

  for (int level=1; level<GMG->nlevels(); level++)
    {
      finestlevel = 0;

      GascoigneMesh*         GM = GMG->GetGascoigneMesh(level);
      GascoigneMeshTransfer* T  = GMG->GetTransfer(level-1);
      
      FM->construct_lists(newquads,oldquads);
      CM = LevelUpdate2d(GM,newquads,oldquads);
	  
      GascoigneMeshTransferConstructor2d(HMP,T,FM,CM);
      
      delete FM;
      FM = CM;
      CM = NULL;
    }
  delete FM;
  FM = NULL;
}

/*---------------------------------------------------*/

void GascoigneMeshConstructor::Loop3d()
{
  finestlevel = 1;

  set<int>  newquads, oldquads;
  HM->GetAwakeCells(newquads);
 
  GascoigneMesh*  GM0 = GMG->GetGascoigneMesh(0);

  LevelMesh3d* FM=NULL;
  LevelMesh3d* CM=NULL;

  FM = LevelUpdate3d(GM0,newquads,oldquads);

  // Q4-PatchStruktur
  if(HM->patchdepth()>=2)
  {
    PatchIndexHandler &PIH = GM0->GetPatchIndexHandler();
    FM->ConstructCellIndOfPatch(_pl2g);
    PIH.GetHasQ4Patch() = 1;
    HangingIndexHandler &H = GM0->GetHangingIndexHandler();
    FM->ConstructHangingStructureQuartic(*(H.GetQ4Structure()),*(H.GetQ4StructureFace()));
  }

  //Wir brauchen local2global und global2local Vector aus FM
  _cl2g =FM->Hexl2g();
  _cg2l = FM->Hexg2l();
  const HierarchicalMesh3d* HMP = dynamic_cast<const HierarchicalMesh3d*>(HM);
  assert(HMP);

  for (int level=1; level<GMG->nlevels(); level++)
    {
      finestlevel = 0;

      GascoigneMesh*         GM = GMG->GetGascoigneMesh(level);
      GascoigneMeshTransfer* T  = GMG->GetTransfer(level-1);
      
      FM->construct_lists(newquads,oldquads);
      //      cout << "*** construct_lists " << level << "\t" << newquads.size() << " " << oldquads.size() << endl;

      CM = LevelUpdate3d(GM,newquads,oldquads);
	  
      GascoigneMeshTransferConstructor3d(HMP,T,FM,CM);
      
      delete FM;
      FM = CM;
      CM = NULL;
    }
  delete FM;
  FM = NULL;
}

/*---------------------------------------------------*/

LevelMesh2d* GascoigneMeshConstructor::LevelUpdate2d(GascoigneMesh* GM, const IntSet& newquads, const IntSet& oldquads) const
{
  LevelMesh2d* CM = new LevelMesh2d(HM);
  CM->BasicInit(newquads,oldquads);
  Construct2d(GM,CM);
  return CM;
}

/*---------------------------------------------------*/

LevelMesh3d* GascoigneMeshConstructor::LevelUpdate3d(GascoigneMesh* GM, const IntSet& newquads, const IntSet& oldquads) const
{
  LevelMesh3d* CM = new LevelMesh3d(HM);
  CM->BasicInit(newquads,oldquads);
  Construct3d(GM,CM);
  return CM;
}

/*---------------------------------------------------*/

void GascoigneMeshConstructor::Construct2d
(GascoigneMesh* NNM, const LevelMesh2d* LM) const
{
  GascoigneMesh2d* NM = dynamic_cast<GascoigneMesh2d*>(NNM);

  assert(NM);

  IntVector& nc = NM->GetCellVector();
  vector<Vertex2d>& nx = NM->GetVertexVector();

  // zellen

  nc.reservesize(4*LM->ncells());
  for(int i=0;i<LM->ncells();i++)
    {
      for(int ii=0;ii<4;ii++) nc[4*i+ii] = LM->vertex_of_cell(i,ii);
    }

  // Koordinaten
  
  nx.reserve(LM->nnodes());
  nx.resize(LM->nnodes());
  for(int i=0;i<LM->nnodes();i++)
    {
      nx[i] = LM->vertex2d(i);
    }

  // PatchStructur

  PatchIndexHandler& PIH =NM->GetPatchIndexHandler();
  LM->ConstructIndOfPatch(PIH.GetIndex());
  PIH.GetDim() = 2;
  PIH.GetHasPatch() = 1;
  PatchToCell2d(PIH,LM);

  // BoundaryIndices

  LM->InitBoundaryHandler(NNM->GetBoundaryIndexHandler(),PIH);

  // Hanging nodes

  HangingIndexHandler& H = NNM->GetHangingIndexHandler();
  LM->ConstructHangingStructureQuadratic(*(H.GetStructure()));
}

/*-----------------------------------------*/

void GascoigneMeshConstructor::Construct3d
(GascoigneMesh* NNM, const LevelMesh3d* LM) const
{
  GascoigneMesh3d* NM = dynamic_cast<GascoigneMesh3d*>(NNM);

  assert(NM);

  IntVector& nc = NM->GetCellVector();
  vector<Vertex3d>& nx = NM->GetVertexVector();

  // zellen

  nc.reservesize(8*LM->ncells());
  for(int i=0;i<LM->ncells();i++)
    {
      for(int ii=0;ii<8;ii++) nc[8*i+ii] = LM->vertex_of_cell(i,ii);
    }

  // Koordinaten
  
  nx.reserve(LM->nnodes());
  nx.resize(LM->nnodes());
  for(int i=0;i<LM->nnodes();i++)
    {
      nx[i] = LM->vertex3d(i);
    }

  // PatchStructur

  PatchIndexHandler& PIH = NM->GetPatchIndexHandler();

  LM->ConstructIndOfPatch(PIH.GetIndex());

  PIH.GetDim() = 3;
  PIH.GetHasPatch() = 1;

  PatchToCell3d(PIH,LM);

  // BoundaryIndices

  LM->InitBoundaryHandler(NNM->GetBoundaryIndexHandler(),PIH);

  // Hanging nodes

  HangingIndexHandler& H = NNM->GetHangingIndexHandler();
  LM->ConstructHangingStructureQuadratic(*(H.GetStructure()),*(H.GetStructureFace()));
}

/*-----------------------------------------*/

void GascoigneMeshConstructor::PatchToCell2d
(PatchIndexHandler& PIH, const LevelMesh2d* LM) const
{
  IntVector ci;
  LM->ConstructCellIndOfPatch(ci);

  int np = ci.size();

  nvector<IntVector >& patch2cell=PIH.GetAllPatch2Cell();
  patch2cell.resize(np);

  const HierarchicalMesh2d* HM2d = dynamic_cast<const HierarchicalMesh2d*>(HM);
  assert(HM2d);

  for (int i=0;i<np;++i)
    {
      int j = ci[i];
      const Quad& q = HM2d->quad(j);
      patch2cell[i].resize(4);
      for (int p=0;p<4;++p)
	{
	  patch2cell[i][p]=LM->Quadg2l(q.child(p));
	  assert(patch2cell[i][p]>=0);
	}
    }
}

/*-----------------------------------------*/

void GascoigneMeshConstructor::PatchToCell3d
(PatchIndexHandler& PIH, const LevelMesh3d* LM) const
{
  IntVector ci;
  LM->ConstructCellIndOfPatch(ci);

  int np = ci.size();

  nvector<IntVector >& patch2cell = PIH.GetAllPatch2Cell();
  patch2cell.resize(np);

  const HierarchicalMesh3d* HM3d = dynamic_cast<const HierarchicalMesh3d*>(HM);
  assert(HM3d);

  for (int i=0;i<np;++i)
    {
      int j = ci[i];
      const Hex& q = HM3d->hex(j);
      patch2cell[i].resize(8);
      for (int p=0;p<8;++p)
	{
	  patch2cell[i][p]=LM->Hexg2l(q.child(p));
	  assert(patch2cell[i][p]>=0);
	}
    }
}
}

/*-----------------------------------------*/
