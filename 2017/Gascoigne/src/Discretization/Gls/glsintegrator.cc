/**
*
* Copyright (C) 2004, 2005, 2006, 2010 by the Gascoigne 3D authors
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


#include  "glsintegrator.h"
#include  "glsequation.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
template<int DIM>
GlsIntegrator<DIM>::GlsIntegrator() : BasicIntegrator()
{
  if (DIM==2)
    IF = new QuadGauss4;
  else
    IF = new HexGauss8;
  assert(GlsIntegrator<DIM>::IF);
}

/* ----------------------------------------- */

template<>
double GlsIntegrator<2>::Volume2MeshSize(double vol) const
{
  return sqrt(vol);
}

/* ----------------------------------------- */

template<>
double GlsIntegrator<3>::Volume2MeshSize(double vol) const
{
  return cbrt(vol);
}

/*-----------------------------------------*/

template<int DIM>
void GlsIntegrator<DIM>::Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& U, 
    const LocalData& Q, const LocalData& QC) const
{
  assert(F.ncomp()==U.ncomp());

  const GlsEquation &GEQ = dynamic_cast<const GlsEquation &>(EQ);

  Vertex<DIM> x, xi;

  const IntegrationFormulaInterface& IFF = FormFormula();
  
  DoubleVector      Lu(U.ncomp());
  nmatrix<double>   A(F.ncomp(),F.ncomp());

  for (int k=0; k<IFF.n(); k++)
    {
      IFF.xi(xi,k);
      FEM.point(xi);

      double vol = FEM.J();
      double h  = Volume2MeshSize(vol);
      double weight  = IFF.w(k) * vol;
      BasicIntegrator::universal_point(FEM,_UH,U);
      BasicIntegrator::universal_point(FEM,_QH,Q);
      FEM.x(x);
      GEQ.SetFemData(_QH);
      GEQ.glspoint(h,_UH,x);
      Lu.zero();
      GEQ.L(Lu,_UH);
      for (int i=0;i<FEM.n();i++)
        {
	  FEM.init_test_functions(_NN,weight,i);
          A.zero();
	  GEQ.S(A,_UH,_NN);
	  
	  for (int c=0; c<F.ncomp(); c++)
	    {
	      for (int d=0; d<F.ncomp(); d++)
		{
		  F(i,c) += A(c,d) * Lu[d];
		}
	    }
	  //  vorher	  A.mult(F[i],Lu);
        }
    }
}

/*-----------------------------------------------------------*/

template<int DIM>
void GlsIntegrator<DIM>::Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, 
    const LocalData& Q, const LocalData& QC) const
{
  assert(E.Ndof()==FEM.n());
  assert(E.Mdof()==FEM.n());
  assert(E.Ncomp()==U.ncomp());

  const GlsEquation& GEQ = dynamic_cast<const GlsEquation&>(EQ);

  Vertex<DIM>     x, xi;
  FemFunction    NNN(FEM.n());

  const IntegrationFormulaInterface& IFF = FormFormula();

  int ncomp = U.ncomp();
  vector<nmatrix<double> >  LMat(FEM.n(),nmatrix<double>(ncomp,ncomp));
  nmatrix<double> SMat(ncomp,ncomp);

  for (int k=0; k<IFF.n(); k++)
    {
      IFF.xi(xi,k);
      FEM.point(xi);

      double vol = FEM.J();
      double h  = Volume2MeshSize(vol);
      double weight  = IFF.w(k) * vol;
      BasicIntegrator::universal_point(FEM,_UH,U);
      BasicIntegrator::universal_point(FEM,_QH,Q);
      FEM.x(x);
      GEQ.SetFemData(_QH);
      GEQ.glspointmatrix(h,_UH,x);
       
      for (int i=0; i<FEM.n(); i++)
        {
          FEM.init_test_functions(NNN[i],1.,i);
          LMat[i].zero();
          GEQ.LMatrix(LMat[i],_UH,NNN[i]);
        }
      for (int i=0; i<FEM.n(); i++)
        {
          SMat.zero();
          GEQ.S(SMat,_UH,NNN[i]);
	  SMat.equ(weight,SMat);

	  for(int d=0; d<ncomp; d++)
	    {
	      for (int j=0; j<FEM.n(); j++)
		{
		  E.SetDofIndex(i,j);
                  for(int c=0;c<ncomp;c++)  
                    {
		      for(int e=0; e<ncomp; e++)  
			{
			  E(c,d) += SMat(c,e) * LMat[j](e,d);
			}
		    }
		}
	    }
	}

      // derivative of S
      // does not enhance convergence for bench.param

      DoubleVector Lu(ncomp,0.);
      GEQ.L(Lu,_UH);
      Lu.equ(weight,Lu);
      
      vector<TestFunction>   NS(ncomp),MS(ncomp);
      for(int d=0; d<ncomp; d++)
	{
	  NS[d].zero();
	  MS[d].zero();
	}
      DoubleVector DS(ncomp,0.);
      for (int j=0; j<FEM.n(); j++)
	{
	  for(int d=0; d<ncomp; d++)
	    {
	      MS[d] = NNN[j];
	      for (int i=0; i<FEM.n(); i++)
		{
		  E.SetDofIndex(i,j);
		  for(int c=0; c<ncomp; c++)
		    {
		      NS[c] = NNN[i];
		      
		      GEQ.SMatrix(DS,_UH,MS,NS);
		      E(c,d) += Lu*DS;
		      
		      NS[c].zero();
		    }
		}      
	      MS[d].zero();
	    }
	}
    }
}

/*-----------------------------------------------------------*/

template class GlsIntegrator<2>;
template class GlsIntegrator<3>;
}
