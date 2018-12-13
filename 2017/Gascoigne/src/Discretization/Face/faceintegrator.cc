/**
*
* Copyright (C) 2007, 2011 by the Gascoigne 3D authors
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


#include "faceintegrator.h"

using namespace std;


namespace Gascoigne
{
  template<int DIM, int Q>
  FaceIntegrator<DIM, Q>::FaceIntegrator() : FaceIntegratorInterface()
  {
    if (!FaceFormulaPointer())
      {
        assert(2<=DIM && DIM<=3);
        assert(1<=Q && Q<=2);

	if ((DIM==2)&&(Q==1))
	  FaceFormulaPointer() = new LineGauss2;
	else if ((DIM==2)&&(Q==2))
	  FaceFormulaPointer() = new LineGauss3;
	else if ((DIM==3)&&(Q==1))
	  FaceFormulaPointer() = new QuadGauss4;
	else if ((DIM==3)&&(Q==2))
	  FaceFormulaPointer() = new QuadGauss9;
      }
  }

  //////////////////////////////////////////////////

  template<int DIM, int Q>
  void FaceIntegrator<DIM, Q>::universal_point_face(const FemInterface& FEM1,const FemInterface& FEM2, FemFunction& U1, FemFunction& U2, const LocalVector& U) const
  {
    U1.resize(U.ncomp());
    U2.resize(U.ncomp());
    
    for (int c=0; c<U1.size(); c++)  U1[c].zero();
    for (int c=0; c<U2.size(); c++)  U2[c].zero();
    
    assert(2<=DIM && DIM<=3);

    if (DIM==2)
      {
	assert(FEM1.n()==(Q+1)*(Q+1));
	for (int iy=0;iy<Q+1;++iy)
	  for (int ix=0;ix<Q+1;++ix)
	    {
	      FEM1.init_test_functions(__NN,1.,(Q+1)*iy+ix);
	      for (int c=0; c<U1.size(); c++)
		U1[c].add(U((2*Q+1)*iy+ix,c),__NN);
	      FEM2.init_test_functions(__NN,1.,(Q+1)*iy+ix);
	      for (int c=0; c<U1.size(); c++)
		U2[c].add(U((2*Q+1)*iy+ix+Q,c),__NN);
	    }
      }
    else if (DIM==3)
      {
	assert(FEM1.n()==(Q+1)*(Q+1)*(Q+1));
	for (int iz=0;iz<Q+1;++iz)
	  for (int iy=0;iy<Q+1;++iy)
	    for (int ix=0;ix<Q+1;++ix)
	      {
	      FEM1.init_test_functions(__NN,1.,(Q+1)*(Q+1)+iz+(Q+1)*iy+ix);
	      for (int c=0; c<U1.size(); c++)
		U1[c].add(U((Q+1)*(2*Q+1)*iz+(2*Q+1)*iy+ix,c),__NN);
	      FEM2.init_test_functions(__NN,1.,(Q+1)*(Q+1)*iz+(Q+1)*iy+ix);
	      for (int c=0; c<U1.size(); c++)
		U2[c].add(U((Q+1)*(2*Q+1)*iz+(2*Q+1)*iy+ix+Q,c),__NN);
	      }
      }
  }
  
  //////////////////////////////////////////////////
  
  template<int DIM, int Q>
  void FaceIntegrator<DIM, Q>::FaceForm(const FaceEquation& EQ, LocalVector& F, const FemInterface& FEM1, const FemInterface& FEM2, const LocalVector& U) const
  {
    int n_dof = U.n();
    
    F.ReInit(EQ.GetNcomp(),n_dof);

    const IntegrationFormulaInterface& IF = *FaceFormula();
  
    F.zero();
    Vertex<DIM> x1,x2, n;
    Vertex<DIM-1> xi;
  
    FemFunction  _U1,_U2;
    TestFunction _N1,_N2;

  
    for (int k=0; k<IF.n(); k++)
      {
	IF.xi(xi,k);
	FEM1.point_boundary(1,xi);
	FEM2.point_boundary(3,xi);
	universal_point_face(FEM1,FEM2,_U1,_U2,U);
	double h1  = FEM1.G();
	assert(fabs(h1-FEM2.G())<1.e-8);
	double weight  = IF.w(k) * h1;
	FEM1.x(x1);
	FEM2.x(x2);
	assert(fabs(x1[0]-x2[0])<1.e-8);
	assert(fabs(x1[1]-x2[1])<1.e-8);
	FEM1.normal(n);
	EQ.point_face(h1,_U1,_U2,x1,n);

        assert(2<=DIM && DIM<=3);

	if (DIM==2)
	  {
	    assert((Q+1)*(2*Q+1)==n_dof);
	    for (int iy=0;iy<Q+1;++iy)
	      for (int ix=0;ix<2*Q+1;++ix)
		{
		  _N1.zero();
		  _N2.zero();
		  if (ix<Q+1) FEM1.init_test_functions(_N1,weight,iy*(Q+1)+ix);
		  if (ix>Q-1) FEM2.init_test_functions(_N2,weight,iy*(Q+1)+ix-Q);
		  EQ.FaceForm(F.start((2*Q+1)*iy+ix),_U1,_U2,_N1,_N2);
		}
	  }
	else if (DIM==3)
	  {
	    assert((Q+1)*(Q+1)*(2*Q+1)==n_dof);
	    for (int iz=0;iz<Q+1;++iz)
	      for (int iy=0;iy<Q+1;++iy)
		for (int ix=0;ix<2*Q+1;++ix)
		{
		  _N1.zero();
		  _N2.zero();
		  if (ix<Q+1) FEM1.init_test_functions(_N1,weight,(Q+1)*(Q+1)*iz+iy*(Q+1)+ix);
		  if (ix>Q-1) FEM2.init_test_functions(_N2,weight,(Q+1)*(Q+1)*iz+iy*(Q+1)+ix-Q);
		  EQ.FaceForm(F.start((Q+1)*(2*Q+1)*iz+(2*Q+1)*iy+ix),_U1,_U2,_N1,_N2);
		}
	  }
      }
  }

  /* ----------------------------------------- */

  template<int DIM, int Q>
  void FaceIntegrator<DIM,Q>::FaceMatrix(const FaceEquation& EQ, EntryMatrix& E, const FemInterface& FEM1, const FemInterface& FEM2, const LocalVector& U) const
  {
    int n_dof = U.n();
    
    FemFunction NNN1(n_dof), NNN2(n_dof);
    FemFunction _U1,_U2;
    
    E.SetDimensionDof(n_dof,n_dof);
    E.SetDimensionComp(U.ncomp(),U.ncomp());
    E.resize();
    E.zero();

    const IntegrationFormulaInterface& IF = *FaceFormula();

    Vertex<DIM> x1,x2,n;
    Vertex<DIM-1> xi;
    for (int k=0; k<IF.n(); k++)
      {
	IF.xi(xi,k);
	FEM1.point_boundary(1,xi);
	FEM2.point_boundary(3,xi);
	universal_point_face(FEM1,FEM2,_U1,_U2,U);
	
	double h1  = FEM1.G();
	assert(fabs(h1-FEM2.G())<1.e-8);
	double weight  = IF.w(k) * h1;

	FEM1.x(x1);
	FEM2.x(x2);
	assert(fabs(x1[0]-x2[0])<1.e-8);
	assert(fabs(x1[1]-x2[1])<1.e-8);
	FEM1.normal(n);

	EQ.pointmatrix_face(h1,_U1,_U2,x1,n);

	double sw = sqrt(weight);

        assert(2<=DIM && DIM<=3);

	if (DIM==2)
	  {
	    assert(n_dof == (Q+1)*(2*Q+1));
	    
	    for (int iy=0;iy<(Q+1);++iy)
	      for (int ix=0;ix<2*Q+1;++ix)
		{
		  int i = (2*Q+1)*iy+ix;
		  NNN1[i].zero(); NNN2[i].zero();
		  if (ix<Q+1) FEM1.init_test_functions(NNN1[i],sw,iy*(Q+1)+ix);
		  if (ix>Q-1) FEM2.init_test_functions(NNN2[i],sw,iy*(Q+1)+ix-Q);
		}
	  }
	else if (DIM==3)
	  {
	    assert(n_dof == (Q+1)*(Q+1)*(2*Q+1));
	    
	    for (int iz=0;iz<(Q+1);++iz)
	      for (int iy=0;iy<(Q+1);++iy)
		for (int ix=0;ix<2*Q+1;++ix)
		  {
		    int i = (Q+1)*(2*Q+1)*iz+(2*Q+1)*iy+ix;
		    NNN1[i].zero();NNN2[i].zero();
		    if (ix<Q+1) FEM1.init_test_functions(NNN1[i],sw,iz*(Q+1)*(Q+1)+iy*(Q+1)+ix);
		    if (ix>Q-1) FEM2.init_test_functions(NNN2[i],sw,iz*(Q+1)*(Q+1)+iy*(Q+1)+ix-Q);
		  }
	  }

	for (int j=0; j<n_dof; j++)
	  {
	    for (int i=0; i<n_dof; i++)
	      {
		E.SetDofIndex(i,j);
		EQ.FaceMatrix(E,_U1,_U2,NNN1[j],NNN2[j],NNN1[i],NNN2[i]);
	      }
	  }
      }
  }




  template class FaceIntegrator<2,1>;
  template class FaceIntegrator<3,1>;
  template class FaceIntegrator<2,2>;
  template class FaceIntegrator<3,2>;

}
