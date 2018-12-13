/**
*
* Copyright (C) 2006, 2010 by the Gascoigne 3D authors
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


#include "integratorq2q4.h"
#include "patchintegrationformula.h"

using namespace std;

namespace Gascoigne
{

/**********************************************************/

template<>
double IntegratorQ2Q4<2>::Volume2MeshSize(double vol) const
{
  return sqrt(vol);
}

/**********************************************************/

template<>
double IntegratorQ2Q4<3>::Volume2MeshSize(double vol) const
{
  return cbrt(vol);
}

/**********************************************************/

template<int DIM>
int IntegratorQ2Q4<DIM>::PatchMeshNr2IntegratorNr(int in) const
{
  //
  // das ist nur von IntegratorQ1Q2 kopiert!!!
  //
  if(DIM==2)
  {
    return (1-(in>>1))*in+(in>>1)*(3-(in%2));
  }
  if(DIM==3)
  {
    int tmp= in-((in>>2)<<2);
    return ((in>>2)<<2)+(1-(tmp>>1))*tmp+(tmp>>1)*(3-(tmp%2));
  }
  std::cerr << "Integrator Q2/Q4: Dimension " << DIM << " nicht implementiert!" << std::endl;
  abort();
}

/**********************************************************/

template<int DIM>
void IntegratorQ2Q4<DIM>::Form(const Equation& EQ, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, const LocalVector& U, const LocalData& Q, const LocalData& QC) const
{
  assert(FemH.n()==FemL.n());

  F.ReInit(U.ncomp(),FemH.n());
  F.zero();

  IntegrationFormulaInterface* IF;
  if (DIM==2) IF = new PatchFormula2d<16,QuadGauss16>;
  else        IF = new PatchFormula3d<64,HexGauss64>;

  Vertex<DIM> x, xi;

  int numcells;
  if(DIM==2) numcells = 4;
  else numcells = 8;

  //Soviele Integr. Pkte je Zelle
  numcells = IF->n()/numcells;
  //Test, ist das sinnvollgewesen
  if(IF->n()%numcells != 0)
  {
    std::cerr << "Integrator Q2/Q4: Falsche Anzahl Zellen je Patch oder unzulaessige Integrationsformel!" << std::endl;
    abort();
  }

  for (int k=0; k<IF->n(); k++)
  {
    IF->xi(xi,k);
    FemH.point(xi);
    FemL.point(xi);
    double vol = FemL.J();
    double h  = Volume2MeshSize(vol);
    double weight  = IF->w(k) * vol;
    BasicIntegrator::universal_point(FemL,_UH,U);
    BasicIntegrator::universal_point(FemL,_QH,Q);
    FemL.x(x);

    EQ.SetFemData(_QH);
    //Die Zelldatensetzen
    if(k%numcells == 0)
    {
      //Wir sind am anfang einer neuen Zelle
      //Die CellData aus den LocalData hohlen.
      int IntegrCellNum = k/numcells;
      int PatchCellNum=PatchMeshNr2IntegratorNr(IntegrCellNum);

      universal_point(_QCH,QC,PatchCellNum);
      EQ.SetCellData(_QCH);
    }

    EQ.point(h,_UH,x);

    for (int i=0;i<FemH.n();i++)
    {
      FemH.init_test_functions(_NN,weight,i);
      EQ.Form(F.start(i),_UH,_NN);
    }
  }
  delete IF;
}

/**********************************************************/

template<int DIM>
void IntegratorQ2Q4<DIM>::AdjointForm(const Equation& EQ, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, const LocalVector& Z, const LocalData& Q, const LocalData& QC) const
{
  assert(FemH.n()==FemL.n());

  F.ReInit(Z.ncomp(),FemH.n());
  F.zero();

  _NNN.resize(FemH.n());
  EntryMatrix E;
  E.SetDimensionDof(FemH.n(),FemL.n());
  E.SetDimensionComp(Z.ncomp(),Z.ncomp());
  E.resize();
  E.zero();

  IntegrationFormulaInterface* IF;
  if (DIM==2) IF = new PatchFormula2d<16,QuadGauss16>;
  else        IF = new PatchFormula3d<64,HexGauss64>;

  Vertex<DIM> x, xi;
  TestFunction MM;
  for (int k=0; k<IF->n(); k++)
  {
    IF->xi(xi,k);
    FemH.point(xi);
    FemL.point(xi);
    double vol = FemL.J();
    double h  = Volume2MeshSize(vol);
    double weight  = IF->w(k) * vol;
    BasicIntegrator::universal_point(FemL,_UH,Z);
    BasicIntegrator::universal_point(FemL,_QH,Q);
    FemL.x(x);
    //EQ.pointmatrix(h,_QH["u"],_QH,x);
    EQ.pointmatrix(h,_QH["u"],x);
    for (int i=0;i<FemH.n();i++)
    {
      FemH.init_test_functions(_NNN[i],weight,i);
    }
    for (int j=0;j<FemL.n();j++)
    {
      FemL.init_test_functions(MM,1.,j);
      for (int i=0;i<FemH.n();i++)
      {
        E.SetDofIndex(j,i);
        EQ.Matrix(E,_QH["u"],_NNN[i],MM);
      }
    }
  }
  for (int i=0;i<FemH.n();i++)
  {
    for (int c=0; c<Z.ncomp(); c++)
    {
      double sum = 0.;
      for (int j=0; j<FemL.n(); j++)
      {
        for (int d=0; d<Z.ncomp(); d++)
        {
          sum += E(j,i,d,c)*Z(j,d);
        }
      }
      F(i,c) += sum;
    }
  }
  delete IF;
}

/**********************************************************/

template<int DIM>
void IntegratorQ2Q4<DIM>::BoundaryForm(const BoundaryEquation& BE, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, const LocalVector& U, int ile, int col, LocalData& Q, const LocalData& QC) const
{
  assert(FemH.n()==FemL.n());

  F.ReInit(BE.GetNcomp(),FemH.n());
  F.zero();

  IntegrationFormulaInterface* IF;
  if (DIM==2) IF = new PatchFormula1d<4,LineGauss4>;
  else        IF = new PatchFormula2d<16,QuadGauss16>;

  Vertex<DIM> x,n;
  Vertex<DIM-1> xi;

  for(int k=0; k<IF->n(); k++)
  {
    IF->xi(xi,k);
    FemH.point_boundary(ile,xi);
    FemL.point_boundary(ile,xi);
    BasicIntegrator::universal_point(FemL,_UH,U);
    BasicIntegrator::universal_point(FemL,_QH,Q);
    FemL.x(x);
    FemL.normal(n);
    double  h = FemL.G();
    double  weight = IF->w(k)*h;
    BE.SetFemData(_QH);
    BE.pointboundary(h,_UH,x,n);
    for(int i=0; i<FemH.n(); i++)
    {
      FemH.init_test_functions(_NN,weight,i);
      BE.Form(F.start(i),_UH,_NN,col);
    }
  }
  delete IF;
}

/**********************************************************/

template<int DIM>
void IntegratorQ2Q4<DIM>::Rhs(const DomainRightHandSide& RHS, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, const LocalData& Q, const LocalData& QC) const
{
  assert(FemH.n()==FemL.n());

  F.ReInit(RHS.GetNcomp(),FemH.n());
  F.zero();

  IntegrationFormulaInterface* IF;
  if (DIM==2) IF = new PatchFormula2d<16,QuadGauss16>;
  else        IF = new PatchFormula3d<64,HexGauss64>;

  Vertex<DIM> x, xi;
  for (int k=0; k<IF->n(); k++)
    {
      IF->xi(xi,k);
      FemH.point(xi);
      FemL.point(xi);
      double vol = FemL.J();
      double h  = Volume2MeshSize(vol);
      double weight  = IF->w(k) * vol;
      BasicIntegrator::universal_point(FemL,_QH,Q);
      RHS.SetFemData(_QH);
      RHS.SetCellSize(h);
      FemL.x(x);

      for (int i=0;i<FemH.n();i++)
        {
          FemH.init_test_functions(_NN,weight,i);
          RHS(F.start(i),_NN,x);
        }
    }
  delete IF;
}

/**********************************************************/

template<int DIM>
void IntegratorQ2Q4<DIM>::BoundaryRhs(const BoundaryRightHandSide& RHS, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, int ile, int col, const LocalData& Q, const LocalData& QC) const
{
  assert(FemH.n()==FemL.n());

  F.ReInit(RHS.GetNcomp(),FemH.n());
  F.zero();

  IntegrationFormulaInterface* IF;
  if (DIM==2) IF = new PatchFormula1d<4,LineGauss4>;
  else        IF = new PatchFormula2d<16,QuadGauss16>;

  Vertex<DIM> x, n;
  Vertex<DIM-1> xi;

  for(int k=0; k<IF->n(); k++)
  {
    IF->xi(xi,k);
    FemH.point_boundary(ile,xi);
    FemL.point_boundary(ile,xi);
    BasicIntegrator::universal_point(FemL,_QH,Q);
    RHS.SetFemData(_QH);
    FemL.x(x);
    FemL.normal(n);
    double  h = FemL.G();
    double  weight = IF->w(k)*h;
    RHS.SetCellSize(h);
    for(int i=0; i<FemH.n(); i++)
    {
      FemH.init_test_functions(_NN,weight,i);
      RHS(F.start(i),_NN,x,n,col);
    }
  }
  delete IF;
}

/**********************************************************/

template<int DIM>
void IntegratorQ2Q4<DIM>::DiracRhsPoint(LocalVector& b, const FemInterface& FemH, const FemInterface& FemL, const Vertex<DIM>& p, const DiracRightHandSide& DRHS, int j, const LocalData& Q, const LocalData& QC) const
{
  assert(FemH.n()==FemL.n());
  b.zero();

  IntegrationFormulaInterface* IF;
  if (DIM==2) IF = new PatchFormula2d<16,QuadGauss16>;
  else        IF = new PatchFormula3d<64,HexGauss64>;

  Vertex<DIM> x;
  FemH.point(p);
  FemL.point(p);
  FemL.x(x);
  BasicIntegrator::universal_point(FemL,_QH,Q);
  DRHS.SetFemData(_QH);

  for (int i=0; i<FemH.n(); i++)
  {
    FemH.init_test_functions(_NN,1.,i);
    DRHS.operator()(j,b.start(i),_NN,x);
  }
  delete IF;
}

/**********************************************************/

template<int DIM>
double IntegratorQ2Q4<DIM>::MassMatrix(EntryMatrix& E, const FemInterface& FemH, const FemInterface& FemL) const
{
  FemFunction NnnH, NnnL;

  NnnH.resize(FemH.n());
  NnnL.resize(FemL.n());

  E.SetDimensionDof(FemH.n(),FemL.n());
  int ncomp=1;
  E.SetDimensionComp(ncomp,ncomp);
  E.resize();
  E.zero();

  IntegrationFormulaInterface* IF;
  if (DIM==2) IF = new PatchFormula2d<16,QuadGauss16>;
  else        IF = new PatchFormula3d<64,HexGauss64>;

  Vertex<DIM> x, xi;
  double omega = 0.;
  for (int k=0; k<IF->n(); k++)
  {
    IF->xi(xi,k);
    FemH.point(xi);
    FemL.point(xi);
    double vol = FemL.J();
    double weight  = IF->w(k) * vol;
    omega += weight;
    for (int i=0;i<FemH.n();i++)
    {
      FemH.init_test_functions(NnnH[i],1.,i);
    }
    for (int i=0;i<FemL.n();i++)
    {
      FemL.init_test_functions(NnnL[i],1.,i);
    }
    for (int i=0;i<FemH.n();i++)
    {
      for (int j=0;j<FemL.n();j++)
      {
        E.SetDofIndex(i,j);
        E(0,0) += weight * NnnH[i].m() * NnnL[j].m();
      }
    }
  }
  delete IF;
  return omega;
}

/**********************************************************/

template<int DIM>
void IntegratorQ2Q4<DIM>::MassForm(const TimePattern& TP, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, const LocalVector& U) const
{
  assert(FemH.n()==FemL.n());

  F.ReInit(U.ncomp(),FemH.n());
  F.zero();

  IntegrationFormulaInterface* IF;
  if(DIM==2) IF = new PatchFormula2d<16,QuadGauss16>;
  else       IF = new PatchFormula3d<64,HexGauss64>;

  Vertex<DIM> xi;

  for(int k=0; k<IF->n(); k++)
  {
    IF->xi(xi,k);
    FemH.point(xi);
    FemL.point(xi);
    double vol = FemL.J();
    double weight = IF->w(k) * vol;
    BasicIntegrator::universal_point(FemL,_UH,U);
    for(int i=0;i<FemH.n();i++)
    {
      FemH.init_test_functions(_NN,weight,i);
      for(int m=0; m<TP.n(); m++)
      {
        for(int n=0; n<TP.n(); n++)
        {
          F(i,m) += TP(m,n) * _UH[n].m() * _NN.m();
        }
      }
    }
  }
  delete IF;
}

/**********************************************************/

template class IntegratorQ2Q4<2>;
template class IntegratorQ2Q4<3>;

}

