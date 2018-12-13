/**
*
* Copyright (C) 2008, 2010 by the Gascoigne 3D authors
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


#ifndef  __TimeSolver_h
#define  __TimeSolver_h

#include  "stdsolver.h"
#include  "simplematrix.h"
#include  "cginfo.h"

/*-------------------------------------------------------------*/

namespace Gascoigne
{

//////////////////////////////////////////////
///
///@brief
/// Default nonlinear Solver for time-dependent Equations

///
///
//////////////////////////////////////////////


class TimeSolver : public StdSolver
{
private:
  
  MatrixInterface*  _MMP;

protected:

  double            theta, dt, time;
  TimePattern       _TP;
  MatrixInterface*& GetMassMatrixPointer() {return _MMP;}

public:

  TimeSolver() : StdSolver(), _MMP(NULL), theta(1.), dt(0.), time(0.) {}
  ~TimeSolver() { if (_MMP) { delete _MMP; _MMP=NULL;} };

  std::string GetName() const {  return "TimeSolver";}

  void SetTimeData(double _dt, double _theta, double _time);

  const MatrixInterface* GetMassMatrix() const {return _MMP;}
        MatrixInterface* GetMassMatrix()       {return _MMP;}

  void RegisterMatrix();
  void SetProblem(const ProblemDescriptorInterface& PDX);

  void ReInitTimePattern(const ProblemDescriptorInterface& PDX);
  void ReInitMatrix();

  MatrixInterface* NewMassMatrix(int ncomp, const std::string& matrixtype)
    {
      return new SimpleMatrix;
    }

  void AssembleMatrix(const VectorInterface& gu, double d);
  void Form(VectorInterface& gy, const VectorInterface& gx, double d) const;
  void MassMatrixVector(VectorInterface& gf, const VectorInterface& gu, double d) const;
  void InverseMassMatrix(VectorInterface& u, const VectorInterface& f, CGInfo& info);
  void precondition(VectorInterface& u, const VectorInterface& f);
  void cgvmult(VectorInterface& y, const VectorInterface& x, double d) const;
  void L2Projection(VectorInterface& Gu, VectorInterface& Gf);
};

/*-------------------------------------------------------------*/
}

#endif
