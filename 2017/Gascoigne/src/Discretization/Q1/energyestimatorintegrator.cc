/**
*
* Copyright (C) 2004, 2005, 2006, 2007, 2010, 2011 by the Gascoigne 3D authors
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


#include "energyestimatorintegrator.h"

using namespace std;

/**********************************************************/

namespace Gascoigne
{

template<int DIM>
EnergyEstimatorIntegrator<DIM>::EnergyEstimatorIntegrator() : BasicIntegrator(), IF(NULL)
{ 
  _s_energytype = "energy_laplace";
  _d_visc       = 1;
}

template<int DIM>
EnergyEstimatorIntegrator<DIM>::EnergyEstimatorIntegrator(const std::string & s_energytype,double d_visc) : BasicIntegrator(), IF(NULL)
{ 
  _s_energytype = s_energytype;
  _d_visc       = d_visc;
}

/**********************************************************/

template<int DIM>
void EnergyEstimatorIntegrator<DIM>::BasicInit()
{
  if (DIM==2)
  {
    if (IF==NULL) IF = new QuadGauss4;
    _xi[0].x() = 0.;
    _xi[1].x() = 1.;
  }
  else
  {
    if (IF==NULL) IF = new HexGauss8;
    _xi[0].x() = 0.; _xi[0].y() = 0.;
    _xi[1].x() = 1.; _xi[1].y() = 0.;
    _xi[2].x() = 1.; _xi[2].y() = 1.;
    _xi[3].x() = 0.; _xi[3].y() = 1.;
  }
  assert(EnergyEstimatorIntegrator<DIM>::IF);
}

/**********************************************************/

template<int DIM>
EnergyEstimatorIntegrator<DIM>::~EnergyEstimatorIntegrator()
{
  if(IF)
  {
    delete IF;
    IF = NULL;
  }
}

/**********************************************************/

template<>
double EnergyEstimatorIntegrator<2>::Volume2MeshSize(double vol) const
{
  return sqrt(vol);
}

/**********************************************************/

template<>
double EnergyEstimatorIntegrator<3>::Volume2MeshSize(double vol) const
{
  return cbrt(vol);
}

/**********************************************************/

template<int DIM>
void EnergyEstimatorIntegrator<DIM>::Jumps(LocalVector& F, const FemInterface& FEM, const LocalVector& U, int ile) const
{
  Vertex<DIM> n;

  F.ReInit(U.ncomp(),2*DIM-2);
  F.zero();

  for (int i=0; i<2*DIM-2; i++)
  {
    FEM.point_boundary(ile,_xi[i]);
    FEM.normal(n);

    BasicIntegrator::universal_point(FEM,_UH,U);

    if( _s_energytype == "energy" || _s_energytype == "energy_laplace" ) {
      for (int c=0; c<U.ncomp(); c++) {
        F(i,c) += _d_visc * ( n.x() * _UH[c].x() + n.y() * _UH[c].y() );
        if (DIM==3) F(i,c) += _d_visc * ( n.z() * _UH[c].z() );
      }
    }else
    if( _s_energytype == "energy_stokes" ) {
      // jump of the pressure:
      F(i,0) += - n.x() * _UH[0].m() - n.y() * _UH[0].m();
      if (DIM==3) F(i,0) += - n.z() * _UH[0].m();
      
      // jump of the value visc * \del_n v :
      for (int c=1; c<U.ncomp(); c++) {
        F(i,c) += _d_visc * ( n.x() * _UH[c].x() + n.y() * _UH[c].y() );
        if (DIM==3) F(i,c) += _d_visc * ( n.z() * _UH[c].z() );
      }
    }else{
      std::cerr << " Bad EnergyEstimatorIntegrator type supplied: '" << _s_energytype  << "'.\n"; 
      std::cerr << " Please specify a known EnergyEstimatorIntegrator type: enery_laplace or enery_stokes.\n";
      abort();
    }
  }
}

/**********************************************************/

template<int DIM>
double EnergyEstimatorIntegrator<DIM>::JumpNorm(const FemInterface& FEM, fixarray<2*DIM-2,double> jumps, int ile) const
{
  double norm = 0.;
  for (int k=0; k<2*DIM-2; k++)
  {
    FEM.point_boundary(ile,_xi[k]);
    double h = Volume2MeshSize(FEM.J());
    double weight = (1.-DIM/4.)*h*FEM.G();
    norm += weight * jumps[k];
  }
  return norm;
}

/**********************************************************/

template<int DIM>
double EnergyEstimatorIntegrator<DIM>::Residual(const LocalVector& U, const FemInterface& FEM, const Equation& EQ, const DomainRightHandSide* RHS, const LocalData& Q) const
{
  double res = 0.;
  DoubleVector F(U.ncomp());
  
  F.zero();
  Vertex<DIM> x, xi;

  for (int k=0; k<GetFormula().n(); k++)
  {
    GetFormula().xi(xi,k);
    FEM.point(xi);
    double vol = FEM.J();
    double h  = Volume2MeshSize(vol);
    double weight  = GetFormula().w(k) * vol;
    BasicIntegrator::universal_point(FEM,_UH,U);
    BasicIntegrator::universal_point(FEM,_QH,Q);
    FEM.x(x);
    if (RHS) RHS->SetFemData(_QH);
    EQ.SetFemData(_QH);
    EQ.point(h,_UH,x);
    EQ.OperatorStrong(F,_UH);
    double value = 0.;
    if (RHS)
      {
        RHS->SetCellSize(h);
	for (int c=0; c<U.ncomp(); c++)
	  {
	    double rhs = (*RHS)(c,x);
	    value += (rhs-F[c])*(rhs-F[c]);
	  }
      }
    else
      {
	for (int c=0; c<U.ncomp(); c++)
	  {
	    value += F[c]*F[c];
	  }
      }
    res += h*h*weight * value;
  }
  return res;
}

/**********************************************************/

template class EnergyEstimatorIntegrator<2>;
template class EnergyEstimatorIntegrator<3>;
}
