/**
*
* Copyright (C) 2004, 2009 by the Gascoigne 3D authors
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


#ifndef  __SolverData_h
#define  __SolverData_h

#include  "gascoigne.h"
#include  "paramfile.h"
#include  <string>

/*----------------------------------------------*/

namespace Gascoigne
{
class SolverData
{
  protected:

    int           exact_lu, enlarge, iter_pre, iter_post, iter_exact;
    int           bicgstab_pstep;
    DoubleVector  _ilum;
    DoubleVector  vector_direction;
    IntVector     stream_direction, _pfilter;
    double        omega;
    std::string   ilusort, linear_smooth, bicgstab_residual;

    // gibt an wieviele Iterationen mindestens gemacht werden, und
    // ab wann im bicgstab ein Abbruchskriterium greift.
    // iter/miniter wird immer iteriert.
    double              bicgstab_miniter;
    
    int _cgMassMaxIter;
    double _cgMassTol, _cgMassGlobalTol;

  public:

    void BasicInit(const ParamFile* pf);

    const IntVector& GetPfilter()const { return _pfilter;}
    void SetPfilter(const IntVector& pfilter) {_pfilter=pfilter;}

    double GetIluModify(int c)      const { return _ilum[c];}
    const DoubleVector& GetIluModify() const { return _ilum;}
    void SetIluModify(const DoubleVector& ilum) {_ilum=ilum;}

    int    GetExactLu()             const { return exact_lu;}
    int    GetEnlarge()             const { return enlarge;}
    double GetOmega()           const { return omega;}
    int    GetIterPre ()        const { return iter_pre;}
    int    GetIterPost ()       const { return iter_post;}
    int    GetIterExact ()      const { return iter_exact;}
    const std::string& GetIluSort() const { return ilusort;}
    const IntVector& GetStreamDirection() const
      { return stream_direction; }
    const DoubleVector& GetVectorDirection() const
      { return vector_direction; }
    
    int    GetBiCGStabPStep()    const { return bicgstab_pstep;}
    std::string GetBiCGStabResidual()const { return bicgstab_residual;}
    std::string GetLinearSmooth()    const { return linear_smooth;}
    double GetBiCGStabMinIter()  const { return bicgstab_miniter;}

    int GetCgMassMaxIter() const { return _cgMassMaxIter; }
    double GetCgMassTol() const { return _cgMassTol; }
    double GetCgMassGlobalTol() const { return _cgMassGlobalTol; }
};
}

#endif
