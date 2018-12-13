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


#include  "hnstructureq13d.h"
#include  "gascoignemesh.h"


using namespace std;

namespace Gascoigne
{
  

/*-----------------------------------------*/

  HNStructureQ13d::HNStructureQ13d() : HNStructureQ12d(), faces(NULL)
  {
    lnop[0][0]=0; lnop[0][1]=2 ; lnop[0][2]=6; lnop[0][3]=8;  lnop[0][4]=4;
    lnop[1][0]=2; lnop[1][1]=20; lnop[1][2]=8; lnop[1][3]=26; lnop[1][4]=14;
    lnop[2][0]=6; lnop[2][1]=8;  lnop[2][2]=24;lnop[2][3]=26; lnop[2][4]=16;
    lnop[3][0]=0; lnop[3][1]=18; lnop[3][2]=6; lnop[3][3]=24; lnop[3][4]=12;
    lnop[4][0]=0; lnop[4][1]=2;  lnop[4][2]=18;lnop[4][3]=20; lnop[4][4]=10;
    lnop[5][0]=18;lnop[5][1]=20; lnop[5][2]=24;lnop[5][3]=26; lnop[5][4]=22;

    lnoe[0][0]=0; lnoe[0][1]=1; lnoe[0][2]=2; 
    lnoe[1][0]=0; lnoe[1][1]=3; lnoe[1][2]=6; 
    lnoe[2][0]=2; lnoe[2][1]=5; lnoe[2][2]=8;
    lnoe[3][0]=6; lnoe[3][1]=7; lnoe[3][2]=8;
 
    lnoe[4][0]=18; lnoe[4][1]=19; lnoe[4][2]=20;
    lnoe[5][0]=18; lnoe[5][1]=21; lnoe[5][2]=24;
    lnoe[6][0]=20; lnoe[6][1]=23; lnoe[6][2]=26;
    lnoe[7][0]=24; lnoe[7][1]=25; lnoe[7][2]=26;
 
    lnoe[8][0]=0; lnoe[8][1]=9;  lnoe[8][2]=18;
    lnoe[9][0]=2; lnoe[9][1]=11; lnoe[9][2]=20;
    lnoe[10][0]=6; lnoe[10][1]=15; lnoe[10][2]=24;
    lnoe[11][0]=8; lnoe[11][1]=17; lnoe[11][2]=26;

  }

/*-----------------------------------------*/

  HNStructureQ13d::~HNStructureQ13d()
  {}

/*--------------------------------------------------------*/

  void HNStructureQ13d::ReInit(const MeshInterface* M)
  {
    const GascoigneMesh* GM = dynamic_cast<const GascoigneMesh*>(M);
    edges = GM->GetHangingIndexHandler().GetStructure();
    faces = GM->GetHangingIndexHandler().GetStructureFace();
    assert(edges);
    assert(faces);
  }

/*-----------------------------------------*/

  bool HNStructureQ13d::ZeroCheck(const GlobalVector& u) const
  {  
    bool r = HNStructureQ12d::ZeroCheck(u);
    if (r) return r;

    for(const_fiterator p=faces->begin(); p!=faces->end(); p++)
      {
	int i = p->first;
	for(int c=0;c<u.ncomp();c++) 
	  {
	    if(u(i,c)!=0.)  return 1;
	  }
      }
    return 0;
  }

/*-----------------------------------------*/

  void HNStructureQ13d::Zero(GlobalVector& u) const
  {
    HNStructureQ12d::Zero(u);

    for(const_fiterator p=faces->begin(); p!=faces->end(); p++)
      {
	int i = p->first;
	for(int c=0;c<u.ncomp();c++) u.zero_node(i);
      }
  }

/*-----------------------------------------*/

  void HNStructureQ13d::Average(GlobalVector& u) const
  {
    HNStructureQ12d::Average(u);

    for(const_fiterator p=faces->begin(); p!=faces->end(); p++)
      {
	int i = p->first;
	const fixarray<9,int>& f = p->second;
	u.equ_node(i, 0.25,f[0], 0.25,f[1], 0.25,f[3], 0.25,f[4]);
      }
  }

/*-----------------------------------------*/

  void HNStructureQ13d::Distribute(GlobalVector& u) const
  {
    HNStructureQ12d::Distribute(u);

    for(const_fiterator p=faces->begin(); p!=faces->end(); p++)
      {
	int i = p->first;
	const fixarray<9,int>& f = p->second;

	u.add_node(f[0],0.25,i);
	u.add_node(f[1],0.25,i);
	u.add_node(f[3],0.25,i);
	u.add_node(f[4],0.25,i);
	u.zero_node(i);
      }
  }

/*-----------------------------------------*/

  int HNStructureQ13d::hanging(int i) const 
  { 
    int r = HNStructureQ12d::hanging(i);

    if (r>0) return 2;
    if (faces->find(i)!=faces->end()) return 4;
    return 0;
  }

/*-----------------------------------------*/

  void HNStructureQ13d::MatrixDiag(int ncomp, MatrixInterface& A) const
  {
    HNStructureQ12d::MatrixDiag(ncomp,A);

    nmatrix<double> M(ncomp);
    M.identity();
    for(const_fiterator p=faces->begin(); p!=faces->end(); p++)
      {
	A.entry_diag(p->first,M);
      }
  }

/*-----------------------------------------*/

  void HNStructureQ13d::SparseStructureDiag(SparseStructure* S) const
  {
    HNStructureQ12d::SparseStructureDiag(S);

    for(const_fiterator p=faces->begin(); p!=faces->end(); p++)
      {
	int i = p->first;
	S->build_add(i,i);
      }
  }

/*----------------------------------------------*/

  fixarray<4,int> HNStructureQ13d::GetHangingFace(int i) const
  {
    const_fiterator p = faces->find(i);
    assert(p!=faces->end());
	  
    fixarray<4,int> Face;
    Face[0] = p->second[0];
    Face[1] = p->second[1];
    Face[2] = p->second[3];
    Face[3] = p->second[4];
  
    return Face;
  }

/*----------------------------------------------*/

  fixarray<2,int> HNStructureQ13d::GetHangingEdge(int i) const
  {
    map<int,EdgeVector>::const_iterator p = edges->find(i);
    assert(p!=edges->end());
	  
    fixarray<2,int> Edge;
    Edge[0] = p->second[0];
    Edge[1] = p->second[1];
  
    return Edge;
  }

/*----------------------------------------------*/

  void HNStructureQ13d::CondenseHanging(IntVector& indices) const
  {
    CondenseHanging2er(indices);
    CondenseHanging4er(indices);
  }

/*----------------------------------------------*/

  void HNStructureQ13d::CondenseHanging(EntryMatrix& E, IntVector& indices) const
  {
    CondenseHanging2er(E,indices);
    CondenseHanging4er(E,indices);
  }

/*----------------------------------------------*/

  void HNStructureQ13d::CondenseHanging4er(EntryMatrix& E, IntVector& indices) const
  {
    IntVector x(0), y(0);

    for(int ii=0; ii<8; ii++)
      {
	int j = indices[ii];

	if (hanging(j)==4) // 4er haengender Knoten
	  {
	    fixarray<4,int> Face = GetHangingFace(j);

	    x.push_back(ii);
	    for(int i=0; i<4; i++)
	      {
		int FaceIndex = Face[i];
			  //
			  // suche ob FaceIndex schon in indices sind
			  //
		int jj = 0;
		bool found = 0;
		while ( (jj<8) && !found)
		  {
		    found = (indices[jj] == FaceIndex);
		    jj++;
		  }
		jj--;
		if (found) y.push_back(jj); // merke Kopplung in Nachbar vertex
		else       indices[ii] = FaceIndex; // ersetze Kopplung
	      }
	  }
      }
    assert(y.size()==3*x.size());

    int counter = 0;
    for(int i=0;i<x.size();i++)
      {
	int i1 = x[i];           // new node !
      
	E.multiply_column_row(i1,0.25);
	int last = counter+3;
      
	assert (last<=y.size());

	for ( ; counter<last; counter++)
	  {
	    int i2 = y[counter];   // already there
	    E.add_column_row(i2,i1);	  
	  }
      }
  }

/*----------------------------------------------*/

  void HNStructureQ13d::Couplings(IntVector& indices) const
  {
	      // fuer Structure (nicht bigstencil)
	      //
	      // erst alle haengenden lines ersetzen
	      //
    int linecount = 0;
    int quadcount = 0;

    IntVector::const_iterator p0 = indices.begin();
    IntVector::const_iterator p1 = indices.end();

    for(int i=0; i<8; i++)
      {
	int& ind = indices[i];

	if (hanging(ind)!=2) continue;
      
	linecount++;

		  //const IntVector2& line = hang(ind);
	fixarray<2,int> line = GetHangingEdge(ind);
      
	for(int k=0; k<2; k++) 
	  {
		      // entweder gibt es newindex schon oder muss hinzugefuegt werden
		      //
	    int newindex = line[k];
	    IntVector::const_iterator p = find(p0,p1,newindex);
	    if (p==p1)
	      {
		ind = newindex;
		break;
	      }
	  }
      }
	      // jetzt duerften zu haengenden quads
	      // nur noch ein vertex nicht eingetragen sein
    for(int i=0; i<8; i++)
      {
	int& ind = indices[i];

	if (hanging(ind)!=4) continue;

	quadcount++;

	fixarray<4,int> face = GetHangingFace(ind);

	for(int k=0; k<4; k++) 
	  {
	    int newindex = face[k];
	    IntVector::const_iterator p = find(p0,p1,newindex);
	    if (p==p1)
	      {
		ind = newindex;
		break;
	      }
	  } 
      }
    assert (quadcount<=3);
    if (quadcount==1) assert (linecount>=2);
    if (quadcount>=2) assert (linecount==3);
  }

/*----------------------------------------------*/

  void HNStructureQ13d::CondenseHanging2er(EntryMatrix& E, IntVector& indices) const
  {
    IntVector x(0), y(0);

    for(int ii=0; ii<8; ii++)
      {
	int i = indices[ii];

	if (hanging(i)==2) // 2er haengender Knoten
	  {
	    fixarray<2,int> Edge = GetHangingEdge(i);

	    x.push_back(ii);

	    for(int iii=0; iii<2; iii++)
	      {
		int ir = Edge[iii];
		bool found = 0;
		for(int iiii=0; (iiii<8) && !found; iiii++)
		  {
		    if(ir==indices[iiii])
		      {
			found = 1;
			y.push_back(iiii);
		      }
		  }
		if(!found) indices[ii] = ir;
	      }
	  }
      }
    assert(x.size()==y.size());
    assert(x.size()<=3);

    for(int i=0;i<x.size();i++)
      {
	int i1 = x[i];        // new node !
	int i2 = y[i];        // already there
      
	E.multiply_column_row(i1,0.5);
	E.add_column_row (i2,i1);
      }
  }

/*----------------------------------------------*/

  void HNStructureQ13d::CondenseHanging2er(IntVector& indices) const
  {
    for(int ii=0; ii<8; ii++)
      {
	int i = indices[ii];

	if (hanging(i)==2) // 2er haengender Knoten
	  {
	    fixarray<2,int> Edge = GetHangingEdge(i);

	    for(int iii=0; iii<2; iii++)
	      {
		int ir = Edge[iii];
		bool found = 0;
		for(int iiii=0; (iiii<8) && !found; iiii++)
		  {
		    if(ir==indices[iiii])
		      {
			found = 1;
		      }
		  }
		if(!found) indices[ii] = ir;
	      }
	  }
      }
  }

/*----------------------------------------------*/

  void HNStructureQ13d::CondenseHanging4er(IntVector& indices) const
  {
    for(int ii=0; ii<8; ii++)
      {
	int j = indices[ii];

	if (hanging(j)==4) // 4er haengender Knoten
	  {
	    fixarray<4,int> Face = GetHangingFace(j);

	    for(int i=0; i<4; i++)
	      {
		int FaceIndex = Face[i];
			  //
			  // suche ob FaceIndex schon in indices sind
			  //
		int jj = 0;
		bool found = 0;
		while ( (jj<8) && !found)
		  {
		    found = (indices[jj] == FaceIndex);
		    jj++;
		  }
		jj--;
		if (!found) indices[ii] = FaceIndex; // ersetze Kopplung
	      }
	  }
      }
  }

/*----------------------------------------------*/

  void HNStructureQ13d::CondenseHangingPatch(EntryMatrix& E, IntVector& indices) const
  {
    assert(indices.size()==27);
  
    for(int ii=0; ii<6; ii++)
      {
	int i = indices[lnop[ii][4]];
      
	if(hanging(i)!=4) continue;
      
	const fixarray<5,int>& p = lnop[ii];
      
	int elim = p[4];
      
	E.add_column     (p[0],elim,0.25);
	E.add_column     (p[1],elim,0.25);
	E.add_column     (p[2],elim,0.25);
	E.add_column     (p[3],elim,0.25);
      
	E.add_row        (p[0],elim,0.25);
	E.add_row        (p[1],elim,0.25);
	E.add_row        (p[2],elim,0.25);
	E.add_row        (p[3],elim,0.25);
	E.multiply_column(elim,0.);
	E.multiply_row   (elim,0.);
      }
    for(int ii=0; ii<12; ii++)
      {
	int i = indices[lnoe[ii][1]];
      
	if(hanging(i)!=2) continue;
      
	const fixarray<3,int>& p = lnoe[ii];
      
	int elim = p[1];
      
	E.add_column     (p[0],elim,0.5);
	E.add_column     (p[2],elim,0.5);
	E.add_row        (p[0],elim,0.5);
	E.add_row        (p[2],elim,0.5);
	E.multiply_column(elim,0.);
	E.multiply_row   (elim,0.);
      }
  }
}

