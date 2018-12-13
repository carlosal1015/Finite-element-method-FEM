/**
*
* Copyright (C) 2004, 2005, 2006, 2008, 2011 by the Gascoigne 3D authors
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


#include  "stdloop.h"
#include  "adaptordata.h"
#include  "diplomantenadaptor.h"
#include  "malteadaptor.h"
#include  "monitoring.h"
#include  "filescanner.h"
#include  "stdmultilevelsolver.h"
#include  "meshagent.h"
#include  "stdsolver.h"
#include  "energyestimator.h"
#include  "q22dwithsecond.h"


using namespace std;


/*-----------------------------------------*/

namespace Gascoigne
{
StdLoop::StdLoop() : BasicLoop()
{
  _estimator = _extrapolate = "none";
}

/*-----------------------------------------*/

StdLoop::~StdLoop() {}

/*-----------------------------------------*/

void StdLoop::ClockOutput() const
{
  BasicLoop::ClockOutput();
  cout << "  Functionals\t\t" << _clock_functionals.read() << endl;
  cout << "  Estimate\t\t" << _clock_estimate.read() << endl;
}

/*-----------------------------------------*/

void StdLoop::BasicInit(const ParamFile* paramfile, const ProblemContainer* PC, const FunctionalContainer* FC)
{
  BasicLoop::BasicInit(paramfile, PC, FC);

  DataFormatHandler DFH;

  DFH.insert("nmin",               &_nmin,              1000);
  DFH.insert("nmax",               &_nmax,              100000);
  DFH.insert("p",                  &_p,                 0.1);
  DFH.insert("random_coarsening",  &_random_coarsening, 0);
  DFH.insert("coarse",             &_coarse,            0);
  DFH.insert("refiner",            &_refiner,           "global");
  DFH.insert("estimator",          &_estimator,         "none");
  DFH.insert("extrapolate",        &_extrapolate,       "no");
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile,"Loop");
}

/*-------------------------------------------------------*/

DoubleVector StdLoop::ComputeFunctionals(VectorInterface& f, VectorInterface& u)
{
  DoubleVector j = GetMultiLevelSolver()->ComputeFunctionals(f,u);
  return j;
} 

/*-----------------------------------------*/

void StdLoop::EtaVisu(string name, int i, const DoubleVector& eta) const
{
  Visualization Visu;
  Visu.format("vtk");
  Visu.set_name(name);
  Visu.step(i);
  VisuDataInfo     VDI(1);
  VisuDataNVector  VD(eta);
  
  Visu.SetPointData(&VD);
  Visu.SetMesh(GetMeshAgent()->GetMesh(0));
  Visu.SetPointDataInfo(&VDI);
  Visu.write();
}

/*-----------------------------------------*/

void StdLoop::EtaCellVisu(string name, int i, const GlobalVector& eta) const
{
  Visualization Visu;
  Visu.format("vtk");
  Visu.set_name(name);
  Visu.step(i);
  VisuDataInfo     VDI(eta.ncomp());
  VisuDataCompVector  VD(eta);
  
  Visu.SetCellData(&VD);
  Visu.SetMesh(GetMeshAgent()->GetMesh(0));
  Visu.SetCellDataInfo(&VDI);
  Visu.write();
}

/*-------------------------------------------------*/

DoubleVector StdLoop::GetExactValues() const
{
  return GetMultiLevelSolver()->GetExactValues();
}

/*-------------------------------------------------*/

DoubleVector StdLoop::Functionals(VectorInterface& u, VectorInterface& f)
{
  DoubleVector J      = ComputeFunctionals(f,u);
  DoubleVector Jexact = GetExactValues();
  _JErr.resize(J.size());
  if (J.size())
    {
      cout << "\nvalue";
      cout.precision(16);
      for(int i=0; i<J.size(); i++) 
        {
          cout << "\t" << J[i];
        }
      cout << "\nerror";
      for(int i=0; i<J.size(); i++) 
        {
          _JErr[i] = Jexact[i] - J[i];
          cout << "\t" << _JErr[i] ;
        }
      cout << endl;
      if(_extrapolate=="yes")
        {
          Extra.NewValues(J);
          Extra.Print();
        }
      cout << endl;
    }
   return J;
}

/*-------------------------------------------------*/

double StdLoop::Estimator(DoubleVector& eta, VectorInterface& u, VectorInterface& f)
{
  double est = 0.;
  if (_estimator=="energy" || _estimator=="energy_laplace" || _estimator=="energy_stokes" )
    {
      StdSolver* S = dynamic_cast<StdSolver*>(GetMultiLevelSolver()->GetSolver());
      assert(S);
      
      MeshAgent* MA = dynamic_cast<MeshAgent*>(GetMeshAgent());
      assert(MA);
      S->GetHierarchicalMeshPointer() = MA -> GetHierarchicalMesh();

      EnergyEstimator E(*S);
      est = E.Estimator(eta,u,f);
      EtaVisu(_s_resultsdir+"/eta",_iter,eta);
    }
  else if(_estimator=="second")
    {
      int dim = GetMeshAgent()->GetMesh(0)->dimension();

      if(dim==2)
        {
        DiscretizationInterface *DIL;
        DIL = GetMultiLevelSolver()->GetSolver()->GetDiscretization();

        Q22dWithSecond DIH;
        DIH.BasicInit(_paramfile);

        GetMultiLevelSolver()->GetSolver()->SetDiscretization(DIH,true);
        GetMultiLevelSolver()->GetSolver()->HNAverage(u);

        eta.resize(GetMeshAgent()->nnodes());
        eta.zero();
        DIH.EstimateSecond(eta,GetMultiLevelSolver()->GetSolver()->GetGV(u));
        GetMultiLevelSolver()->GetSolver()->HNZero(u);

        GetMultiLevelSolver()->GetSolver()->SetDiscretization(*DIL);
        }
      else
        {
          cerr << "Estimator \"second\" not written for 3D!" << endl;
          abort();
        }
    }
  else
    {
      cerr << "Estimator type '" << _estimator << "' unknown. Use either energy_laplace or energy_stokes\n"; 
      abort();
    }
  return est;
}

/*-------------------------------------------------*/

void StdLoop::AdaptMesh(const DoubleVector& eta,string refine_or_coarsen_step)
{
  //das gleichzeitige vergroebern und verfeinern FUNKTIONIERT NICHT
  //wer das machen moechte, muss stattdessen in zwei getrennten laeufen 
  //das gitter verfeinern, reinit+interpolate und dann das gitter vergroebern
  if(refine_or_coarsen_step=="refine") ;
  else if(refine_or_coarsen_step=="coarsen") ;
  else {
    cerr<<"the variable refine_or_coarsen_step has to be set, either to 'refine' or 'coarsen'"<<endl;
    abort();
  }

  if (_refiner=="global") {
    if(refine_or_coarsen_step=="refine"){
      GetMeshAgent()->global_refine(1);
    }else{
      GetMeshAgent()->random_patch_refine(-0.1,0); 
    }
  }  
  else if(_refiner=="none")
    {
      GetMeshAgent()->global_refine(0); 
    }
  else if(_refiner=="random") 
    {
      if (GetMeshAgent()->nnodes()>_nmax) _p *= 0.5;
      if (GetMeshAgent()->nnodes()<_nmin) _p *= 1.1;

      if(refine_or_coarsen_step=="refine"){
        GetMeshAgent()->random_patch_refine(_p,0);
      }
      if(refine_or_coarsen_step=="coarsen"){
        if(_random_coarsening){
          GetMeshAgent()->random_patch_coarsen(_p,0);
        }else{
          GetMeshAgent()->random_patch_refine(-0.1,0); 
        }
      }
    }
  else if(_refiner=="random_refine") 
    {
      if(refine_or_coarsen_step=="refine"){
        if (GetMeshAgent()->nnodes()>_nmax) _p *= 0.5;
        if (GetMeshAgent()->nnodes()<_nmin) _p *= 1.1;
        GetMeshAgent()->random_patch_refine(_p,0);
      }
    }
  else if(_refiner=="random_coarsen") 
    {
      if(refine_or_coarsen_step=="coarsen"){
        if (GetMeshAgent()->nnodes()>_nmax) _p *= 0.5;
        if (GetMeshAgent()->nnodes()<_nmin) _p *= 1.1;
        GetMeshAgent()->random_patch_coarsen(_p,0);
      }
    }
  else if(_refiner=="eta") 
    {
      IntVector refnodes, coarsenodes,dummynodes;

      MalteAdaptor A(_paramfile,eta);
      A.refine(refnodes,coarsenodes);

      if(refine_or_coarsen_step=="coarsen") GetMeshAgent()->refine_nodes(dummynodes,coarsenodes);
      if(refine_or_coarsen_step=="refine")  GetMeshAgent()->refine_nodes(refnodes);
    }
  else if(_refiner=="dip") 
    {
      IntVector refnodes, coarsenodes;

      AdaptorData info;
      info.rfactor() = 1.; 
      DiplomandenAdaptor A(info,eta);
      A.refine(refnodes);
      if(refine_or_coarsen_step=="refine")  {
        GetMeshAgent()->refine_nodes(refnodes,coarsenodes);
      }else{
        GetMeshAgent()->random_patch_refine(-0.1,0); 
      }

    }
  else
    {
      cerr << "Unknown value \"" << _refiner << "\" for \"refiner\" in paramfile. " << endl;
      cerr << "Please use \"global\", \"none\", \"random\", \"random_refine\", \"random_coarsen\", \"eta\" or \"dip\" instead." << endl;
      abort();
    }
}

/*-------------------------------------------------*/

void StdLoop::AdaptMesh(const DoubleVector& eta)
{
  //das gleichzeitige vergroebern und verfeinern FUNKTIONIERT NICHT
  //wer das machen moechte, sollte stattdessen zwei getrennte laeufe durchfuehren:
  //das gitter vergroebern, reinit+interpolate und dann das gitter verfeinern
  //das entsprechend die methode AdaptMesh(eta,refine_or_coarsen_step) aufrufen
  if     (_refiner=="global") GetMeshAgent()->global_refine(1);
  else if(_refiner=="none") 
    {
      GetMeshAgent()->global_refine(0);
    }  
  else if(_refiner=="random") 
    {
      if (GetMeshAgent()->nnodes()>_nmax) _p *= 0.5;
      if (GetMeshAgent()->nnodes()<_nmin) _p *= 1.1;
      if(_random_coarsening){
        cerr<<"Das gleichzeitige Vergroebern und Verfeinern FUNKTIONIERT NICHT!"<<endl;
        cerr<<"Fuehre stattdessen zwei getrennte Laeufe durch: random_refine, random_coarsen"<<endl;
        cerr<<"und rufe dazu jewweils AdaptMesh(eta,refine_or_coarsen_step) auf."<<endl;
        abort();
      }
      GetMeshAgent()->random_patch_refine(_p,_random_coarsening);
    }
  else if(_refiner=="random_refine") 
    {
      if (GetMeshAgent()->nnodes()>_nmax) _p *= 0.5;
      if (GetMeshAgent()->nnodes()<_nmin) _p *= 1.1;
      GetMeshAgent()->random_patch_refine(_p,0);
    }
  else if(_refiner=="random_coarsen") 
    {
      if (GetMeshAgent()->nnodes()>_nmax) _p *= 0.5;
      if (GetMeshAgent()->nnodes()<_nmin) _p *= 1.1;
      GetMeshAgent()->random_patch_coarsen(_p,0);
    }
  else if(_refiner=="eta") 
    {
      IntVector refnodes, coarsenodes;

      MalteAdaptor A(_paramfile,eta);
      A.refine(refnodes,coarsenodes);

      if(refnodes.size()>0 && coarsenodes.size()>0){
        cerr<<"Das gleichzeitige Vergroebern und Verfeinern FUNKTIONIERT NICHT!"<<endl;
        cerr<<"Fuehre stattdessen zwei getrennte Laeufe durch, einmal vergroebern, einmal verfeinern"<<endl;
        cerr<<"und rufe dazu jewweils AdaptMesh(eta,refine_or_coarsen_step) auf."<<endl;
        abort();
      }

      GetMeshAgent()->refine_nodes(refnodes,coarsenodes);
    }
  else if(_refiner=="dip") 
    {
      IntVector refnodes, coarsenodes;

      AdaptorData info;
      info.rfactor() = 1.; 
      DiplomandenAdaptor A(info,eta);
      A.refine(refnodes);
      GetMeshAgent()->refine_nodes(refnodes,coarsenodes);
    }
  else
    {
      cerr << "Unknown value \"" << _refiner << "\" for \"refiner\" in paramfile. " << endl;
      cerr << "Please use \"global\", \"none\", \"random\", \"random_refine\", \"random_coarsen\", \"eta\" or \"dip\" instead." << endl;
      abort();
    }
}

/*-------------------------------------------------*/

void StdLoop::run(const std::string& problemlabel)
{
  VectorInterface u("u"), f("f");
  GlobalVector  ualt;

  Monitoring Moning;
  
  for (_iter=1; _iter<=_niter; _iter++)
    {
      cout << "\n================== " << _iter << " ================";
      PrintMeshInformation();
      Moning.SetMeshInformation(_iter,GetMeshAgent()->nnodes(),GetMeshAgent()->ncells());

      _clock_newmesh.start();

      GetMultiLevelSolver()->ReInit(problemlabel);
      GetMultiLevelSolver()->ReInitVector(u);
      GetMultiLevelSolver()->ReInitVector(f);
      GetMultiLevelSolver()->InterpolateSolution(u,ualt);
      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;

      _clock_newmesh.stop();

      if (_iter==1) 
        {
          GetMultiLevelSolver()->GetSolver()->OutputSettings();
          InitSolution(u);
          Moning.BasicInit(GetExactValues());
        }

      Solve(u,f);

      ComputeGlobalErrors(u);
      
      _clock_functionals.start();
      DoubleVector juh = Functionals(u,f);
      _clock_functionals.stop();

      DoubleVector eta;

      _clock_estimate.start();
      if (_estimator!="none")
        {
          double est = Estimator(eta,u,f);
          Moning.SetSolutionInformation(_iter,juh,est);
        }
      if (_iter<_niter) 
        {
          CopyVector(ualt,u);

          AdaptMesh(eta);
          //  wenn gleichzeitig verfeinert und vergroebert werden soll, dann sollte
          //  AdaptMesh zweimal folgendermassen aufgerufen werden :
          //     
          //    CopyVector(ualt,u);
          //    AdaptMesh(eta,"refine");
          //    GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
          //    GetMultiLevelSolver()->ReInit(problemlabel);
          //    GetMultiLevelSolver()->InterpolateSolution(u,ualt);
          //    CopyVector(ualt,u);
          //    AdaptMesh(eta,"coarsen");
          //   
        }
      _clock_estimate.stop();
     }

  ClockOutput();
}
}

/*-------------------------------------------------*/
