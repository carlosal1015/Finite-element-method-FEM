/**
*
* Copyright (C) 2004, 2005, 2006, 2008 by the Gascoigne 3D authors
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


#include "integratorwithsecond.h"

using namespace std;

namespace Gascoigne
{
/**********************************************************/

template<>
void IntegratorWithSecond<2>::point_hesse(const FemInterface& E, const Vertex<2>& v) const
{
  typedef FiniteElementWithSecond<2, 1, Transformation2d<BaseQ22dWithSecond>, BaseQ22dWithSecond>  FEWithSecond2d;
  
  const FEWithSecond2d* FE = dynamic_cast<const FEWithSecond2d*>(&E);

  FE->ComputeHesse(v);
}

/* ----------------------------------------- */
 
template<>
void IntegratorWithSecond<3>::point_hesse(const FemInterface& E, const Vertex<3>& v) const
{
  typedef FiniteElementWithSecond<3, 2, Transformation3d<BaseQ23dWithSecond>, BaseQ23dWithSecond> FEWithSecond3d;
  
  const FEWithSecond3d* FE = dynamic_cast<const FEWithSecond3d*>(&E);

  FE->ComputeHesse(v);
}

/* ----------------------------------------- */

template<>
void IntegratorWithSecond<2>::init_test_hesse(const FemInterface& E, TestFunction& N, double w, int i) const
{
  typedef FiniteElementWithSecond<2, 1, Transformation2d<BaseQ22dWithSecond>, BaseQ22dWithSecond>  FEWithSecond2d;
  
  const FEWithSecond2d* FE = dynamic_cast<const FEWithSecond2d*>(&E);
  
  FE->init_test_hesse(N,w,i);
}

/* ----------------------------------------- */

template<>
void IntegratorWithSecond<3>::init_test_hesse(const FemInterface& E,TestFunction& N, double w, int i) const
{
  typedef FiniteElementWithSecond<3, 2, Transformation3d<BaseQ23dWithSecond>, BaseQ23dWithSecond> FEWithSecond3d;
  
  const FEWithSecond3d* FE = dynamic_cast<const FEWithSecond3d*>(&E);
  
  FE->init_test_hesse(N,w,i);
}

/* ----------------------------------------- */

template<>
void IntegratorWithSecond<2>::hesse(const FemInterface& E, FemFunction& UH, const LocalVector& U) const
{
  UH.resize(U.ncomp());

  for (int c=0; c<UH.size(); c++)
  {
    UH[c].aux("xx") = 0;
    UH[c].aux("xy") = 0;
    UH[c].aux("yy") = 0;
  }

  for (int i=0; i<E.n(); i++)
    {
      init_test_hesse(E, _NN, 1., i);

      for (int c=0; c<UH.size(); c++)
	{
	  UH[c].aux("xx") += U(i,c) * _NN.aux("xx");
	  UH[c].aux("xy") += U(i,c) * _NN.aux("xy");
	  UH[c].aux("yy") += U(i,c) * _NN.aux("yy");
	}
    }
}
/* ----------------------------------------- */

template<>
void IntegratorWithSecond<3>::hesse(const FemInterface& E, FemFunction& UH, const LocalVector& U) const
{
  UH.resize(U.ncomp());

  for (int c=0; c<UH.size(); c++)
  {
    UH[c].aux("xx") = 0;
    UH[c].aux("xy") = 0;
    UH[c].aux("yy") = 0;
    UH[c].aux("xz") = 0;
    UH[c].aux("yz") = 0;
    UH[c].aux("zz") = 0;
  }

  for (int i=0; i<E.n(); i++)
    {
      init_test_hesse(E, _NN, 1., i);

      for (int c=0; c<UH.size(); c++)
	{
	  UH[c].aux("xx") += U(i,c) * _NN.aux("xx");
	  UH[c].aux("xy") += U(i,c) * _NN.aux("xy");
	  UH[c].aux("yy") += U(i,c) * _NN.aux("yy");
	  UH[c].aux("xz") += U(i,c) * _NN.aux("xz");
	  UH[c].aux("yz") += U(i,c) * _NN.aux("yz");
	  UH[c].aux("zz") += U(i,c) * _NN.aux("zz");	}
    }
}

/* ----------------------------------------- */

template<int DIM>
void IntegratorWithSecond<DIM>::hesse(const FemInterface& E, FemData& QH, const LocalData& Q) const
{
  LocalData::const_iterator p=Q.begin();
  for(; p!=Q.end(); p++)
    {
      hesse(E, QH[p->first], p->second);
    }
}

/* ----------------------------------------- */

template<int DIM>
double IntegratorWithSecond<DIM>::ComputeDomainFunctional(const DomainFunctional& F, const FemInterface& FEM, const LocalVector& U, 
    const LocalData& Q, const LocalData& QC) const
{
  const IntegrationFormulaInterface& IF = *GalerkinIntegratorQ2<DIM>::FormFormula();

  Vertex<DIM> x, xi;
  double j = 0.;
  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point(xi);
      point_hesse(FEM,xi);
      double vol = FEM.J();
      double weight  = IF.w(k) * vol;
      BasicIntegrator::universal_point(FEM,GalerkinIntegratorQ2<DIM>::_UH,U);
      BasicIntegrator::universal_point(FEM,GalerkinIntegratorQ2<DIM>::_QH,Q);

      hesse(FEM,GalerkinIntegratorQ2<DIM>::_UH,U);
      hesse(FEM,GalerkinIntegratorQ2<DIM>::_QH,Q);

      FEM.x(x);
      F.SetFemData(GalerkinIntegratorQ2<DIM>::_QH);
      j += weight * F.J(GalerkinIntegratorQ2<DIM>::_UH,x);
    }
  return j;
}

/* ----------------------------------------- */

template<int DIM>
void IntegratorWithSecond<DIM>::Rhs(const DomainRightHandSide& f, LocalVector& F, const FemInterface& FEM, 
    const LocalData& Q, const LocalData& QC) const
{
  F.ReInit(f.GetNcomp(),FEM.n());

  const IntegrationFormulaInterface& IF = *GalerkinIntegratorQ2<DIM>::FormFormula();

  F.zero();
  Vertex<DIM> x, xi;
  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point(xi);
      point_hesse(FEM,xi);
      double vol = FEM.J();
      double h  = GalerkinIntegratorQ2<DIM>::Volume2MeshSize(vol);
      double weight  = IF.w(k) * vol;
      GalerkinIntegratorQ2<DIM>::universal_point(FEM,GalerkinIntegratorQ2<DIM>::_QH,Q);
      hesse(FEM,GalerkinIntegratorQ2<DIM>::_QH,Q);

      f.SetFemData(GalerkinIntegratorQ2<DIM>::_QH);
      f.SetCellSize(h);
      FEM.x(x);
      for (int i=0;i<FEM.n();i++)
	{
        FEM.init_test_functions(GalerkinIntegratorQ2<DIM>::_NN,weight,i);
        init_test_hesse(FEM, GalerkinIntegratorQ2<DIM>::_NN, weight, i);
        f(F.start(i),GalerkinIntegratorQ2<DIM>::_NN,x);
	}
    }
}

template<int DIM>
void Gascoigne::IntegratorWithSecond<DIM>::
EstimateSecond(LocalVector& F, const FemInterface& FEM, const LocalVector& U) const
{
  F.ReInit(U.ncomp(),1);

  const IntegrationFormulaInterface& IF = *GalerkinIntegratorQ2<DIM>::FormFormula();

  F.zero();
  Vertex<DIM> xi;

  for (int k=0; k<IF.n(); k++)
  {
    IF.xi(xi,k);
    FEM.point(xi);
    point_hesse(FEM,xi);
    double vol = FEM.J();
    double weight  = IF.w(k) * vol;

    this->_UH.resize(U.ncomp());

    for (int c=0; c<this->_UH.size(); c++)  
    {
      this->_UH[c].zero();
      this->_UH[c].aux("xx")=0.;
      this->_UH[c].aux("xy")=0.;
      this->_UH[c].aux("yy")=0.;
    }


    for (int i=0;i<FEM.n();i++)
    {
      FEM.init_test_functions(GalerkinIntegratorQ2<DIM>::_NN,1.,i);
      init_test_hesse(FEM, GalerkinIntegratorQ2<DIM>::_NN,1., i);
      for (int c=0; c<this->_UH.size(); c++)
      {
        this->_UH[c].add(U(i,c),this->_NN);
        this->_UH[c].aux("xx") += U(i,c) * this->_NN.aux("xx");
        this->_UH[c].aux("xy") += U(i,c) * this->_NN.aux("xy");
        this->_UH[c].aux("yy") += U(i,c) * this->_NN.aux("yy");
      }
    }

    for (int c=0; c<this->_UH.size(); c++)
    {
      F(0,c) += vol*weight*(this->_UH[c].aux("xx")*this->_UH[c].aux("xx")
          +2*this->_UH[c].aux("xy")*this->_UH[c].aux("xy")
          +this->_UH[c].aux("yy")*this->_UH[c].aux("yy"));
    }
  }
}

/* ----------------------------------------- */

template class IntegratorWithSecond<2>;
template class IntegratorWithSecond<3>;

}
