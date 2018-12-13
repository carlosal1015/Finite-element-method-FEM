/**
*
* Copyright (C) 2004 by the Gascoigne 3D authors
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


#ifndef  __ExactSolution_h
#define  __ExactSolution_h


#include  "vertex.h"
#include  <string>
#include  "application.h"

/*-----------------------------------------*/

namespace Gascoigne
{
  
  //////////////////////////////////////////////
  ///
  ///@brief
  /// Interface class for ExactSolution

  ///
  ///
  //////////////////////////////////////////////

  class ExactSolution : public virtual Application
  {
    private:

    protected:
      double eps;

    public:
      ExactSolution(): Application(), eps(1.e-6) {}
      virtual ~ExactSolution() {}

      ////////// 2d
      virtual double operator()(int c, const Vertex2d& v)const {
        std::cerr << "\"ExactSolution::operator()\" not written!" << std::endl;
        abort();
      } 

      virtual double x         (int c, const Vertex2d& v)const;
      virtual double y         (int c, const Vertex2d& v)const;
      virtual double z         (int c, const Vertex2d& v)const { abort(); }
      virtual double xx        (int c, const Vertex2d& v)const;
      virtual double yx        (int c, const Vertex2d& v)const;
      virtual double xy        (int c, const Vertex2d& v)const;
      virtual double yy        (int c, const Vertex2d& v)const;

      ////////// 3d
      virtual double operator()(int c, const Vertex3d& v)const {
        std::cerr << "\"ExactSolution::operator()\" not written!" << std::endl;
        abort();
      } 

      virtual double x         (int c, const Vertex3d& v)const;
      virtual double y         (int c, const Vertex3d& v)const;
      virtual double z         (int c, const Vertex3d& v)const;
      virtual double xx        (int c, const Vertex3d& v)const;
      virtual double yy        (int c, const Vertex3d& v)const;
      virtual double zz        (int c, const Vertex3d& v)const;
      virtual double xy        (int c, const Vertex3d& v)const;
      virtual double yz        (int c, const Vertex3d& v)const;
      virtual double xz        (int c, const Vertex3d& v)const;
  };
}

#endif
