/**
*
* Copyright (C) 2004, 2008, 2011 by the Gascoigne 3D authors
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


#ifndef  __FemInterface_h
#define  __FemInterface_h

#include  "vertex.h"
#include  "nmatrix.h"
#include  "gascoigne.h"
#include  <string>

/*-----------------------------------------*/

namespace Gascoigne
{
  class FemInterface
  {
    private:

    protected:

    public:
      typedef nmatrix<double>    Matrix;

      FemInterface() {}
      virtual ~FemInterface() {}

      virtual std::string GetName() const=0;

      virtual int    n() const=0;
      virtual double J() const=0;
      virtual double G() const=0;

      virtual void x(Vertex2d& v) const {
        std::cerr << "\"FemInterface::x\" not written!" << std::endl;
        abort();
      }
      virtual void x(Vertex3d& v) const {
        std::cerr << "\"FemInterface::x\" not written!" << std::endl;
        abort();
      }

      virtual void normal(Vertex2d& v) const {
        std::cerr << "\"FemInterface::normal\" not written!" << std::endl;
        abort();
      }
      virtual void normal(Vertex3d& v) const {
        std::cerr << "\"FemInterface::normal\" not written!" << std::endl;
        abort();
      }

      virtual void point(const Vertex2d& v) const {
        std::cerr << "\"FemInterface::point\" not written!" << std::endl;
        abort();
      }
      virtual void point(const Vertex3d& v) const {
        std::cerr << "\"FemInterface::point\" not written!" << std::endl;
        abort();
      }

      virtual void  point_boundary(int ie, const Vertex1d& v) const {
        std::cerr << "\"FemInterface::point_boundary\" not written!" << std::endl;
        abort();
      }
      virtual void  point_boundary(int ie, const Vertex2d& v) const {
        std::cerr << "\"FemInterface::point_boundary\" not written!" << std::endl;
        abort();
      }

      virtual void ReInit(const Matrix& M) const=0;
      virtual void init_test_functions(TestFunction& Phi, double w, int i) const=0;
      virtual void Anisotropy(DoubleMatrix& A) const
      { 
        std::cerr << "\"FemInterface::Anisotropy\" not written!" << std::endl;
        abort();
      }

      virtual void GetCoordinates(DoubleMatrix& A) const { 
	std::cerr << "\"FemInterface::GetCoordinates\" not written!" << std::endl;
        abort();}
  };
}

#endif
