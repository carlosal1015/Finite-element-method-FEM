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


#ifndef __iluinterface_h
#define __iluinterface_h

#include  "matrixinterface.h"

/*-------------------------------------------------------------*/

namespace Gascoigne
{
  class IluInterface
  {
    private:

    protected:

    public:
      IluInterface() {}
      virtual ~IluInterface() {};

      virtual int   n() const=0;
      virtual std::string GetName() const=0;
      virtual void ReInit(const SparseStructureInterface* A)=0;
      virtual void ConstructStructure(const IntVector& perm, const MatrixInterface& A)=0;

      virtual void modify(int c, double s) {
        std::cerr << "\"IluInterface::modify\" not written!" << std::endl;
        abort();
      }
      virtual void zero() {
        std::cerr << "\"IluInterface::zero\" not written!" << std::endl;
        abort();
      }
      virtual void compute_ilu () {
        std::cerr << "\"IluInterface::compute_ilu\" not written!" << std::endl;
        abort();
      }
      virtual void copy_entries(const MatrixInterface* A) {
        std::cerr << "\"IluInterface::copy_entries\" not written!" << std::endl;
        abort();
      }
      virtual void solve(GlobalVector& x) const {
        std::cerr << "\"IluInterface::solve\" not written!" << std::endl;
        abort();
      }
      virtual void solve_transpose(GlobalVector& x) const {
        std::cerr << "\"IluInterface::solve_transpose\" not written!" << std::endl;
        abort();
      }
      virtual std::ostream& Write(std::ostream &s) const {
        std::cerr << "\"IluInterface::Write\" not written!" << std::endl;
        abort();
      }
  };
}

#endif
