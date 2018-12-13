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


#ifndef __adaptordata_h
#define __adaptordata_h

#include  <string>
#include  <stdio.h>
#include  "paramfile.h"

namespace Gascoigne
{
class AdaptorData
{
  int     inc,inr,indr,inactual,idim,mnodes,iorder;
  double  irfactor,icfactor,idrfactor,iminf,imaxf,itol,ieta,_randomp;
  double  iglobal_conv,ilocal_conv;
  

  void    ErrorMessage(const std::string&, const std::string&) const;

  void    CheckEstimate  () const;
  void    CheckFunctional() const;
  void    CheckAdaptor   () const;

  std::string  _functional, _adaptor, _estimate;

public:

  /* Zugriff */

  const std::string& functional() const { return _functional;}
        std::string& functional()       { return _functional;}
  const std::string& adaptor   () const { return _adaptor;}
  const std::string& estimate  () const { return _estimate;}

  int      order()    const {return iorder;}
  int&     order()          {return iorder;}
  int      dim()       const {return idim;}
  int&     dim()             {return idim;}
  int      maxnodes()  const {return mnodes;}
  int      ncells()    const {return inactual;}
  int&     ncells()          {return inactual;}
  int      nc()        const {return inc;}
  int&     nc()              {return inc;}
  int      nr()        const {return inr;}
  int&     nr()              {return inr;}
  int      ndr()       const {return indr;}
  double   tol()       const {return itol;}
  double   drfactor()  const {return idrfactor;}
  double   rfactor()   const {return irfactor;}
  double&  rfactor()         {return irfactor;}
  double   global_conv()   const {return iglobal_conv;}
  double&  global_conv()         {return iglobal_conv;}
  double   local_conv()   const {return ilocal_conv;}
  double&  local_conv()         {return ilocal_conv;}
  double   cfactor()   const {return icfactor;}
  double   minf()      const {return iminf;}
  double&  minf()            {return iminf;}
  double   maxf()      const {return imaxf;}
  double&  maxf()            {return imaxf;}
  double   eta()       const { return ieta;}
  double&  eta()             { return ieta;}
  double   randomp()   const { return _randomp;}

  /*   Funktionen   */

  AdaptorData();

  void reset();
  void read (const ParamFile* pf);
};
}

#endif
