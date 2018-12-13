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


#ifndef  __Application_h
#define  __Application_h


#include  "gascoigne.h"


/*-------------------------------------------------------*/

namespace Gascoigne
{

  /////////////////////////////////////////////
  ////
  ////@brief
  ////  ... comments Application

  ////
  ////
  /////////////////////////////////////////////

  class Application
  {
    private:
      mutable double _dt, _time;

    protected:
      double GetTime() const {
        return _time;
      }
      double GetTimeStep() const {
        return _dt;
      }

    public:
      //
      ////  Con(De)structor 
      //

      Application() : _dt(0.),_time(0.) {}
      virtual ~Application() {}

      virtual std::string GetName() const=0;

      virtual void SetTime(double time, double dt) const {
        _time = time; 
        _dt = dt;
      }

      virtual void SetFemData(FemData& q) const {}
      virtual void SetCellData(CellData& q) const {}
      virtual void SetParameterData(LocalParameterData& q) const {}
 };
}

#endif

