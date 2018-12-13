/**
*
* Copyright (C) 2004, 2005 by the Gascoigne 3D authors
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


#ifndef __IntegrationFormulaSummed_h
#define __IntegrationFormulaSummed_h

#include  "integrationformula.h"

/*------------------------------------------------------------*/

namespace Gascoigne
{

template<class INT>
class IntegrationFormulaSummed1d : public IntegrationFormula1d
{
 protected:

  INT I;

 public:

  IntegrationFormulaSummed1d(int n) : IntegrationFormula1d(), I()
    {
      int    N = static_cast<int>(pow(2.,n));

      IntegrationFormula1d::ReInit(N*I.n());

      int nn = static_cast<int>(pow(2.,n));
      double d2 = pow(0.5,n);

      for(int i=0;i<N;i++)
	{
	  int ix = i%nn;
	  double dx = 1.*ix;

	  for(int ii=0;ii<I.n();ii++)
	    {
	      int index = i*I.n()+ii;
	      w(index) = d2 * I.w(ii);

	      double x = I.c(ii).x();

	      x += dx;
	      x *= d2;

	      c(index).x() = x;
	    }
	}
    }
};

/*------------------------------------------------------------*/

template<class INT>
class IntegrationFormulaSummed2d : public IntegrationFormula2d
{
 protected:

  INT I;

 public:

  IntegrationFormulaSummed2d(int n=4) : IntegrationFormula2d(), I()
    {
      int    N = static_cast<int>(pow(4.,n));

      IntegrationFormula2d::ReInit(N*I.n());

      int nn = static_cast<int>(pow(2.,n));
      double d2 = pow(0.5,n);
      double d4 = pow(0.25,n);

      for(int i=0;i<N;i++)
	{
	  int ix = i%nn;
	  int iy = i/nn;

	  double dx = 1.*ix;
	  double dy = 1.*iy;

	  for(int ii=0;ii<I.n();ii++)
	    {
	      int index = i*I.n()+ii;
	      w(index) = d4 * I.w(ii);

	      double x = I.c(ii).x();
	      double y = I.c(ii).y();

	      x += dx;
	      y += dy;

	      x *= d2;
	      y *= d2;

	      c(index).x() = x;
	      c(index).y() = y;
	    }
	}
    }
};

/*------------------------------------------------------------*/

template<class INT>
class IntegrationFormulaSummed3d : public IntegrationFormula3d
{
 protected:

  INT I;

 public:

  IntegrationFormulaSummed3d(int n=2) : IntegrationFormula3d(), I()
    {
      int    N = static_cast<int>(pow(8.,n));

      IntegrationFormula3d::ReInit(N*I.n());

      int nn = static_cast<int>(pow(2.,n));
      double d2 = pow(0.5,n);
      double d8 = pow(0.125,n);

      for(int i=0;i<N;i++)
	{
	  int ix = i%nn;
	  int j = (i-ix)/nn;

	  int iy = j%nn;
	  int iz = j/nn;

	  double dx = 1.*ix;
	  double dy = 1.*iy;
	  double dz = 1.*iz;

	  for(int ii=0;ii<I.n();ii++)
	    {
	      int index = i*I.n()+ii;
	      w(index) = d8 * I.w(ii);

	      double x = I.c(ii).x();
	      double y = I.c(ii).y();
	      double z = I.c(ii).z();

	      x += dx;
	      y += dy;
	      z += dz;

	      x *= d2;
	      y *= d2;
	      z *= d2;

	      c(index).x() = x;
	      c(index).y() = y;
	      c(index).z() = z;
	    }
	}
    }
};

}

#endif
