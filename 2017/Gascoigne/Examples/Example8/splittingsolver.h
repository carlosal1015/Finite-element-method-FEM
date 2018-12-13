/**
*
* Copyright (C) 2008 by the Gascoigne 3D authors
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


#ifndef  __SplittingSolver_h
#define  __SplittingSolver_h

#include  "timesolver.h"
#include  "filescanner.h"

/*-----------------------------------------*/

namespace Gascoigne
{

//////////////////////////////////////////////
//
///@brief
///
///
//////////////////////////////////////////////

class SplittingSolver : public TimeSolver
{
 protected:

  MatrixInterface*                 _MAP2;
  IluInterface*                    _MIP2;
  DiscretizationInterface*         _ZP2;
  const ProblemDescriptorInterface *PD1, *PD2;

  virtual DiscretizationInterface* CreateDiscretization_1() const=0;
  virtual DiscretizationInterface* CreateDiscretization_2() const=0;

 public:

  SplittingSolver() : TimeSolver(), _MAP2(NULL), _MIP2(NULL), _ZP2(NULL), PD1(NULL), PD2(NULL) {}
 ~SplittingSolver() 
   {
     if(_MAP2) delete _MAP2; _MAP2=NULL;
     if(_MIP2) delete _MIP2; _MIP2=NULL;
     if(_ZP2)  delete _ZP2;  _ZP2=NULL;
   }
 void SetProblem(const ProblemDescriptorInterface& PDX)
 {
   if      (PD1==NULL) { PD1 = &PDX;}
   else if (PD2==NULL) { PD2 = &PDX;}
   TimeSolver::SetProblem(PDX);
 }
 void NewMesh(const MeshInterface* mp)
 {
   TimeSolver::NewMesh(mp);
   _ZP2->ReInit(mp);
 }
 void RegisterMatrix()
 {
   if      (GetProblemDescriptor()==PD1) TimeSolver::RegisterMatrix();
   else if (GetProblemDescriptor()==PD2) 
     {
       const Equation*  EQ = GetProblemDescriptor()->GetEquation();
       assert(EQ);
       int ncomp = EQ->GetNcomp();
       if (_MAP2==NULL) _MAP2 = NewMatrix(ncomp, _matrixtype);
       if (_MIP2==NULL) _MIP2 = NewIlu   (ncomp, _matrixtype);
     }
   else abort();
 } 
 MatrixInterface* GetMatrix() const 
   { 
     if      (GetProblemDescriptor()==PD1) return TimeSolver::GetMatrix();
     else if (GetProblemDescriptor()==PD2) return _MAP2;
     abort();
   }
 IluInterface* GetIlu() const 
   { 
     if      (GetProblemDescriptor()==PD1) return TimeSolver::GetIlu();
     else if (GetProblemDescriptor()==PD2) return _MIP2;
     abort();
   }
  const DiscretizationInterface* GetDiscretization() const 
    {
      if      (GetProblemDescriptor()==PD1) return _ZP;
      else if (GetProblemDescriptor()==PD2) return _ZP2;
      abort();
    }
  DiscretizationInterface* GetDiscretization() 
    {
      if      (GetProblemDescriptor()==PD1) return _ZP;
      else if (GetProblemDescriptor()==PD2) return _ZP2;
      abort();
    }
  void BasicInit(const ParamFile* paramfile, const int dimension)
  {
    _paramfile = paramfile;
    _facediscname = "none";
    _useUMFPACK = 0;
    _matrixtype = "block";
    
    DataFormatHandler DFH;
    DFH.insert("ndirect"    , &_ndirect);
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(_paramfile,"Solver");

    GetDiscretizationPointer() = CreateDiscretization_1();
    _ZP2 = CreateDiscretization_2();

    GetDiscretization()->BasicInit(_paramfile);
    _ZP2->BasicInit(_paramfile);
  
    _Dat.BasicInit(_paramfile);
    _PF.SetComponents(_Dat.GetPfilter());
  }
};

}

#endif
