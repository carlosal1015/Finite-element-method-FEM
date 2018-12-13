/**
*
* Copyright (C) 2004, 2005, 2006, 2007, 2009, 2010, 2011 by the Gascoigne 3D authors
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


#include  "q12d.h"
#include  "galerkinintegrator.h"
#include  "transformation2d.h"
#include  "finiteelement.h"
#include  "baseq12d.h"
#include  "sparsestructure.h"
#include  "gascoignemesh.h"
#include  "hnstructureq12d.h"
#include  "mginterpolatormatrix.h"
#include  "mginterpolatornested.h"
#include  "gascoignemeshtransfer.h"
#include  "hnstructureq12d.h"

using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{
Q12d::Q12d() : Q1()
{
}

/* ----------------------------------------- */

HNStructureInterface* Q12d::NewHNStructure()
{
  return new HNStructureQ12d;
}

/* ----------------------------------------- */

void Q12d::BasicInit(const ParamFile* pf)
{
  assert(HN==NULL);
  HN = NewHNStructure();
  assert(HN);

  if(!GetIntegratorPointer())
    GetIntegratorPointer() =  new GalerkinIntegrator<2>;
  assert(GetIntegrator());

  GetIntegratorPointer()->BasicInit();

  if(!GetFemPointer())
    {
      typedef Transformation2d<BaseQ12d>           TransQ1;
      typedef FiniteElement<2,1,TransQ1,BaseQ12d>  FiniteElement;
      CellDiscretization::GetFemPointer() =  new FiniteElement;
    }
  assert(GetFem());

  CellDiscretization::BasicInit(pf);
}

/* ----------------------------------------- */

nmatrix<double> Q12d::GetLocalInterpolationWeights() const
{
  // w(i,j) = interpolation weight of node i to node j
  int nn = 4;//GetMesh()->nodes_per_cell();
  nmatrix<double> w(nn,nn);
  w.zero();
  w(0,1) =  0.5  ; w(0,2) =  0.5  ; w(0,3) =  0.25;
  w(1,0) =  0.5  ; w(1,2) =  0.25 ; w(1,3) =  0.5 ;
  w(2,0) =  0.5  ; w(2,1) =  0.25 ; w(2,3) =  0.5 ;
  w(3,0) =  0.25 ; w(3,1) =  0.5  ; w(3,2) =  0.5 ;
  return w;
}

/* ----------------------------------------- */

void Q12d::StrongDirichletVector(GlobalVector& u, const DirichletData& BF, int col, const vector<int>& comp, double d) const
{
  const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
  assert(GMP);
  DoubleVector ff(u.ncomp(),0.);
  const IntVector& bv = *GMP->VertexOnBoundary(col);

  FemData QH;

  GlobalToGlobalData();
  BF.SetParameterData(__QP);

  for(int ii=0;ii<comp.size();ii++)
    {
      int c = comp[ii];
      if(c<0) {
        cerr << "negative component: " << c << endl;
        abort();
      } else if(c>=u.ncomp()){
        cerr << "unknown component: " << c << endl;
        abort();
      }
    }

  for(int i=0;i<bv.size();i++)
    {
      int index = bv[i];

      QH.clear();
      GlobalData::const_iterator p=GetDataContainer().GetNodeData().begin();
      for(; p!=GetDataContainer().GetNodeData().end(); p++)
      {
        QH[p->first].resize(p->second->ncomp());
        for(int c=0; c<p->second->ncomp(); c++)
        {
          QH[p->first][c].m() = p->second->operator()(index,c);
        }
      }

      BF.SetFemData(QH);

      const Vertex2d& v = GMP->vertex2d(index);
      
      BF(ff,v,col);
      for(int iii=0;iii<comp.size();iii++)
        {
          int c = comp[iii];
          u(index,c) = d * ff[c];
        }
    }
}

/* ----------------------------------------- */

void Q12d::StrongPeriodicVector(GlobalVector& u, const PeriodicData& BF, int col, const vector<int>& comp, double d) const
{
  const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(GetMesh());
  assert(GMP);
  DoubleVector ff(u.ncomp(),0.);
  const IntVector& bv = *GMP->VertexOnBoundary(col);

  FemData QH;
  GlobalToGlobalData();
  BF.SetParameterData(__QP);

  //for(int ii=0;ii<comp.size();ii++)
  //{
  //  int c = comp[ii];
  //  if(c<0) {
  //    cerr << "negative component: " << c << endl;
  //    abort();
  //  } else if(c>=u.ncomp()){
  //    cerr << "unknown component: " << c << endl;
  //    abort();
  //  }
  //}

  for(int i=0;i<bv.size();i++)
  {
    int index = bv[i];

    QH.clear();
    GlobalData::const_iterator p=GetDataContainer().GetNodeData().begin();
    for(; p!=GetDataContainer().GetNodeData().end(); p++)
    {
      QH[p->first].resize(p->second->ncomp());
      for(int c=0; c<p->second->ncomp(); c++)
      {
        QH[p->first][c].m() = p->second->operator()(index,c);
      }
    }

    BF.SetFemData(QH);

    const Vertex2d& v = GMP->vertex2d(index);

    BF(ff,v,col);
    for(int iii=0;iii<comp.size();iii++)
    {
      int c = comp[iii];
      u(index,c) = d * ff[c];
    }
  }
}

/* ----------------------------------------- */

void Q12d::Interpolate(GlobalVector& u, const DomainInitialCondition& U) const
{
  if (&U==NULL) return;

  for(int in=0; in<GetMesh()->nnodes(); ++in)
    {
      Vertex2d v = GetMesh()->vertex2d(in);
      for(int c=0;c<u.ncomp();c++)
        {
          u(in,c) = U(c,v);
        }
    }
}
/* ----------------------------------------- */

void Q12d::InterpolateDirac(GlobalVector& u, const GlobalVector& uold) const
{
  const IntVector& vo2n = *GetMesh()->Vertexo2n();

  assert(vo2n.size()==uold.n());
  assert(GetMesh()->nnodes()==u.n());
  assert(u.ncomp()==uold.ncomp());

  for(int i=0;i<vo2n.size();i++)
    {
      int in = vo2n[i];

      if(in>=0) 
        {
          u.equ_node(in,1.,i,uold);
        }
    }
}
/* ----------------------------------------- */

void Q12d::InterpolateSolutionByPatches(GlobalVector& u, const GlobalVector& uold) const
{
  const IntVector& vo2n = *GetMesh()->Vertexo2n();
  nvector<bool> habschon(GetMesh()->nnodes(),0);  

  assert(vo2n.size()==uold.n());
  assert(GetMesh()->nnodes()==u.n());
  assert(u.ncomp()==uold.ncomp());

  for(int i=0;i<vo2n.size();i++)
    {
      int in = vo2n[i];

      if(in>=0) 
        {
          u.equ_node(in,1.,i,uold);
          habschon[in] = 1;
        }
    }
  nvector<fixarray<3,int> > nodes(4);
  nodes[0][0] = 1; nodes[0][1] = 0;  nodes[0][2] = 2;
  nodes[1][0] = 3; nodes[1][1] = 0;  nodes[1][2] = 6;
  nodes[2][0] = 5; nodes[2][1] = 2;  nodes[2][2] = 8;
  nodes[3][0] = 7; nodes[3][1] = 6;  nodes[3][2] = 8;
 
  const PatchMesh* PM = dynamic_cast<const PatchMesh*>(GetMesh());
  assert(PM);

  for(int iq=0;iq<PM->npatches();++iq)
    {
      IntVector vi =* PM->IndicesOfPatch(iq);

      for(int j=0; j<nodes.size(); j++)
        {
          int v  = vi[nodes[j][0]];
          int v1 = vi[nodes[j][1]];
          int v2 = vi[nodes[j][2]];
          assert(habschon[v1]);
          assert(habschon[v2]);
          if (habschon[v]==0) 
            {
              u.equ_node(v,0.5,v1,uold);
              u.add_node(v,0.5,v2,uold);
              habschon[v] = 1;
            }
        }
      int v = vi[4];
      if (habschon[v]==0)
        {
          u.equ_node(v,0.25,vi[0],uold);
          u.add_node(v,0.25,vi[2],uold);	  
          u.add_node(v,0.25,vi[6],uold);	  
          u.add_node(v,0.25,vi[8],uold);	  
          habschon[v] = 1;
        }
    }  
}

/* ----------------------------------------- */

void Q12d::ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT)
{
  {
    MgInterpolatorNested* IP = dynamic_cast<MgInterpolatorNested*>(I);
    if(IP)
      {
        IP->BasicInit(MT);
        return;
      }
  }

  MgInterpolatorMatrix* IP = dynamic_cast<MgInterpolatorMatrix*>(I);
  assert(IP);
  const GascoigneMeshTransfer* GT = dynamic_cast<const GascoigneMeshTransfer*>(MT);
  assert(GT);

  const map<int,fixarray<2,int> >& zweier = GT->GetZweier();
  const map<int,fixarray<4,int> >& vierer = GT->GetVierer();
  const map<int,fixarray<8,int> >& achter = GT->GetAchter();
  const IntVector& c2f    = GT->GetC2f();

  int n  = c2f.size() +   zweier.size() +   vierer.size() +   achter.size();
  int nt = c2f.size() + 2*zweier.size() + 4*vierer.size() + 8*achter.size();

  ColumnStencil& ST = IP->GetStencil();
  DoubleVector& val = IP->GetAlpha();

  SparseStructure SS;

  SS.build_begin(n);
  for(int i=0;i<c2f.size();i++)
    {
      assert(c2f[i]>=0);

      SS.build_add(c2f[i],i);
    }
  for(map<int,fixarray<2,int> >::const_iterator p=zweier.begin();
      p!=zweier.end();p++) 
    {
      int il = p->first;
      fixarray<2,int> n2 = p->second;
      for(int ii=0;ii<2;ii++) SS.build_add(il,n2[ii]);
    }
  for(map<int,fixarray<4,int> >::const_iterator p=vierer.begin();
      p!=vierer.end();p++) {
    int il = p->first;
    fixarray<4,int> n4 = p->second;
    for(int ii=0;ii<4;ii++) SS.build_add(il,n4[ii]);
  }
  for(map<int,fixarray<8,int> >::const_iterator p=achter.begin();
      p!=achter.end();p++) {
    int il = p->first;
    fixarray<8,int> n8 = p->second;
    for(int ii=0;ii<8;ii++) SS.build_add(il,n8[ii]);
  }
  SS.build_end();

  assert(nt==SS.ntotal());

  ST.memory(&SS);

  val.reservesize(nt);

  for(int i=0;i<c2f.size();i++)
    {
      // ich weiss nicht, ob das richtig ist !!!!!
      int pos = ST.Find(c2f[i],i);
      assert(pos>=0);
      
      val[pos] = 1.;
    }
  for(map<int,fixarray<2,int> >::const_iterator p=zweier.begin();
      p!=zweier.end();p++) 
    {
      int il = p->first;
      fixarray<2,int> n2 = p->second;
      val[ST.Find(il,n2[0])] = 0.5;
      val[ST.Find(il,n2[1])] = 0.5;
    }
  for(map<int,fixarray<4,int> >::const_iterator p=vierer.begin();
      p!=vierer.end();p++) {
    int il = p->first;
    fixarray<4,int> n4 = p->second;
    val[ST.Find(il,n4[0])] = 0.25;
    val[ST.Find(il,n4[1])] = 0.25;
    val[ST.Find(il,n4[2])] = 0.25;
    val[ST.Find(il,n4[3])] = 0.25;
  }
  for(map<int,fixarray<8,int> >::const_iterator p=achter.begin();
      p!=achter.end();p++) {
    int il = p->first;
    fixarray<8,int> n8 = p->second;
    for (int i=0; i<8; i++)
      {
	val[ST.Find(il,n8[i])] = 0.125;
      }
  }
}

/* ----------------------------------------- */

void Q12d::EnergyEstimator(EdgeInfoContainerInterface& EIC, DoubleVector& eta, const GlobalVector& u, const Equation& EQ, const DomainRightHandSide* RHS, const std::string & s_energytype,double d_visc) const
{
  // RHS may be NULL
  //
  EnergyEstimatorIntegrator<2> EEI(s_energytype,d_visc);
  const HierarchicalMesh2d*    HM = dynamic_cast<const HierarchicalMesh2d*>(EIC.GetMesh());

  EdgeInfoContainer<2>& EICC = dynamic_cast<EdgeInfoContainer<2>&>(EIC);

  EEI.BasicInit();

  // Kanten initialisieren
  EEJumps(EICC,u,EEI,HM);
  
  // Kantenintegrale auswerten
  EEJumpNorm(EICC,eta,EEI,HM);

  // Residuenterme auswerten
  EEResidual(eta,u,EQ,RHS,EEI);
}

/* ----------------------------------------- */

void Q12d::EEJumps(EdgeInfoContainer<2>& EIC, const GlobalVector& u, const EnergyEstimatorIntegrator<2>& EEI, const HierarchicalMesh2d* HM) const
{
  fixarray<2,int> vertexes;
  nmatrix<double> T;

  for(int iq=0;iq<HM->ncells();++iq)
  {
    if (!(HM->sleep(iq)))
    {
      Transformation_HM(T,HM,iq);
      GetFem()->ReInit(T);

      GlobalToLocal_HM(__U,u,HM,iq);
      
      for (int ile=0; ile<4; ile++)
      {
        EEI.Jumps(__F,*GetFem(),__U,ile);

        int edgenumber = HM->edge_of_quad(iq,ile);

        if (EIC[edgenumber]==NULL)
        {
          const Edge& edge = HM->edge(edgenumber);
          HM->QuadLawOrder().globalvertices_of_edge(HM->quad(edge.master()),vertexes,edge.LocalMasterIndex());
          EIC[edgenumber] = new EdgeInfo<2>();
          EIC[edgenumber]->BasicInit(&edge,u.ncomp(),vertexes);
        }
        EIC[edgenumber]->AddNodes(__F);
      }
    }
  }
  EIC.ModifyHanging();
}

/* ----------------------------------------- */

void Q12d::EEJumpNorm(EdgeInfoContainer<2>& EIC, DoubleVector& eta, const EnergyEstimatorIntegrator<2>& EEI, const HierarchicalMesh2d* HM) const
{
  nmatrix<double> T;

  for (int iq=0; iq<HM->ncells(); iq++)
  {
    if (!(HM->sleep(iq)))
    {
      Transformation_HM(T,HM,iq);
      GetFem()->ReInit(T);

      double jump = 0.;
      for (int ile=0; ile<4; ile++)
      {
        int edgenumber = HM->edge_of_quad(iq,ile);
        if (EIC[edgenumber]->GetCount()==2)
        {
          jump += EEI.JumpNorm(*GetFem(),EIC[edgenumber]->GetNorm(),ile);
        }
      }
      double w =  0.25 * 0.5 * sqrt(jump);
      for (int in=0; in<4; in++)
      {
        int iv = HM->vertex_of_cell(iq,in);
        eta[iv] += w;
      }
    }
  }
}

/* ----------------------------------------- */

void Q12d::EEResidual(DoubleVector& eta, const GlobalVector& u, const Equation& EQ, const DomainRightHandSide* RHS, const EnergyEstimatorIntegrator<2>& EEI) const
{
  nmatrix<double> T;

  GlobalToGlobalData();
  if (RHS) RHS->SetParameterData(__QP);
  EQ.SetParameterData(__QP);

  for(int iq=0;iq<GetMesh()->ncells();++iq)
  {
    Transformation(T,iq);
    GetFem()->ReInit(T);
        
    //    GlobalToLocalData(iq);
    GlobalToLocal(__U,u,iq);
//    EQ.SetCellData(__QC);
//    if (RHS) EQ.SetCellData(__QC);

    // .cell is analogous to .point 
    // EQ.cell(GetMesh(),iq,__U,__QN); 
    double res = EEI.Residual(__U,*GetFem(),EQ,RHS,__QN);
    double w = 0.25 * sqrt(res);
    for (int in=0; in<4; in++)
    {
      eta[GetMesh()->vertex_of_cell(iq,in)] += w;
    }
  }
}

/* ----------------------------------------- */

int Q12d::GetCellNumber(const Vertex2d& p0, Vertex2d& p,int c0) const
{
  if(c0 != 0)
  {
    VertexTransformation(p0,p,c0,false);
    
    if((p[0]>0.-1.e-12)&&(p[0]<1.+1.e-12)&&(p[1]>0.-1.e-12)&&(p[1]<1.+1.e-12))
    {
      return c0;
    }
  }

  {
    int iq;
  
    for(iq=0; iq<GetMesh()->ncells(); ++iq)
    {
      bool found = true;
      
      for(int d=0; d<2; ++d)
      {
	double min=GetMesh()->vertex2d(GetMesh()->vertex_of_cell(iq,0))[d];
	double max=min;
	for(int j=1; j<4; ++j)
	{
	  double x = GetMesh()->vertex2d(GetMesh()->vertex_of_cell(iq,j))[d];
	  
	  min = Gascoigne::min(min,x);
	  max = Gascoigne::max(max,x);
	}
	if((p0[d]<min)||(p0[d]>max)) 
        {
	  found = false;
	  break;
	}
      }
      if(!found) continue;
    
      VertexTransformation(p0,p,iq);
    
      for(int d=0; d<2; ++d)
      {
	if((p[d]<0.-1.e-12)||(p[d]>1.+1.e-12))
	{
	  found = false;
	}
      }
      if(found) break;
    }

    if(iq<GetMesh()->ncells()) return iq;
    else                       return -1;
  }
}


/* ----------------------------------------- */

  void Q12d::VertexTransformation(const Vertex2d& p0, Vertex2d& p, int iq, bool abortiffail) const
{
  nmatrix<double> T;
  Transformation(T,iq);

  Transformation2d<BaseQ12d> Tr;
  Tr.init(T);

  Vertex2d res;
  
  p = 0.5;
  
  for(int niter=1; ;niter++)
  {
    Tr.point(p);
    
    res = p0;
    res.add(-1,Tr.x());

    if(res.norm()<1.e-13)
    {
      break;
    }
    if (niter>=10)
      {
	cerr << "void Q12d::VertexTransformation(const Vertex2d& p0, Vertex2d& p, int iq) const" << endl
	     << "\t did not converge " << res << "\t" << p << "\t" << p0 << endl;
	if(abortiffail)
	  abort();
	else 
	  return;
      }
    
    Tr.DTI().mult_ad(p,res);
  } 
}

/* ----------------------------------------- */

void Q12d::RhsCurve(GlobalVector &F, const Curve &C,int comp,int N) const
{
  double h = 1./N;        //Schrittweite (im Parameterbereich der Kurve)
  double t0=0.,t1;        //Zeitpunkte   (             "               )
  Vertex2d x0= C(t0),x1;  // zu t0, t1 korrespondierende Punkte auf der Kurve
//  cout << " x0 = " << x0 << endl;
  Vertex2d xr0,xr1;       //Zu x0, x1 gehoerende Punkte in der Referenzzelle
  int c0 = GetCellNumber(x0,xr0), c1;  //Nummern derjenigen Zellen, in denen x0, x1 liegen
//GetCellNumber -> bestimmt die Zelle, in der Punkt x0 liegt, und den korrespondierenden Punkt in der Referenzzelle xr0.

//  cout << " c0 = " << c0 << endl;

  for (int i=0;i<N;++i)
  {
    t1 = h*(i+1);
    x1 = C(t1);
    c1 = GetCellNumber(x1,xr1,c0);
    VertexTransformation(x1,xr1,c0);

    if (c1!=c0)
    {
      x1 = randpunkt(C,t0,t1,c1);
      --i;
      VertexTransformation(x1,xr1,c0);
    }
    cout.precision(16);

    nmatrix<double> T;
    Transformation(T,c0);
    GetFem()->ReInit(T);

    GetIntegrator()->RhsCurve(__F,*GetFem(),xr0,xr1,fabs(t1-t0),C.NormD(t0),C.NormD(t1),C.GetNcomp(),comp);
    BasicDiscretization::LocalToGlobal(F,__F,c0,1.);

    c0=c1;
    x0=x1;
    VertexTransformation(x0,xr0,c0);
    t0=t1;
  }
}

/* ----------------------------------------- */

Vertex2d Q12d::randpunkt(const Curve& C, double t0,double& t1,int& cc) const
{
  Vertex2d xr0,xr1,x_rand,abstand;
  int c0,c1,c_rand;

  Vertex2d x0 = C(t0);
  Vertex2d x1 = C(t1);

  c0 = GetCellNumber(x0,xr0);

  c1 = GetCellNumber(x1,xr1,c0);
  VertexTransformation(x1,xr1,c0);

  double norm;

  int niter = 0;
  do
    {
      double t = 0.5*(t0+t1);
      x1 = C(t);
      c_rand = GetCellNumber(x1,x_rand,c1);
      VertexTransformation(x1,x_rand,c0);
      if (c0==c_rand)
        {
          xr0 = x_rand;
          t0  = t;
        }
      else
        {
          xr1 = x_rand;
          t1 = t;
        }


      abstand = xr0;
      abstand -= xr1;
      norm = abstand.norm_l8();
      niter++;
      if(niter>56) abort();
    } while (fabs(norm)>1.e-12);

  cc = c_rand; 
  return C(t1);
}


/* ----------------------------------------- */

}
