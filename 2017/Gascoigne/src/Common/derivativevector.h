/**
*
* Copyright (C) 2004, 2005, 2007 by the Gascoigne 3D authors
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


#ifndef  __derivativevector_h
#define  __derivativevector_h

#include  "numfixarray.h"

/*--------------------------------------------------------*/

namespace Gascoigne
{
class DerivativeVector : public numfixarray<6,double>
{
 private:
  std::map<std::string,double> _M;

 public:

  double m() const { return (*this)[0]; }
  double x() const { return (*this)[1]; }
  double y() const { return (*this)[2]; }
  double z() const { return (*this)[3]; }
  double n() const { std::cerr << "Normal derivative not written!" << std::endl;  abort(); }
  double D() const { return (*this)[5]; }  // fuer Laplace

  double& m() { return (*this)[0]; }
  double& x() { return (*this)[1]; }
  double& y() { return (*this)[2]; }
  double& z() { return (*this)[3]; }
  double& n() { std::cerr << "Normal derivative not written!" << std::endl;  abort(); }
  double& D() { return (*this)[5]; }

  double aux(const std::string &name) const
    {
      std::map<std::string,double>::const_iterator p = _M.find(name);
      if(p==_M.end())
        {
          std::cerr << name << " not found!" << std::endl;
          abort();
        }
      return p->second;
    }

  double &aux(const std::string &name) { return _M[name]; }

};
}

/*--------------------------------------------------------*/

#endif
