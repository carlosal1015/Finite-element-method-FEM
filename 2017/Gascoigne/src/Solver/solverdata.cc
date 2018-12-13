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


#include  "solverdata.h"
#include  "filescanner.h"


using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
  
  /*-----------------------------------------*/
  
  void SolverData::BasicInit(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    
    _pfilter.resize(0);

    DFH.insert("pfilter" , &_pfilter);
    DFH.insert("exact_lu", &exact_lu, 0);
    DFH.insert("enlarge" , &enlarge , 0);
    DFH.insert("iterpre" , &iter_pre , 4);
    DFH.insert("iterpost", &iter_post, 4);
    DFH.insert("iterexact",&iter_exact , 10);
    DFH.insert("omega"   , &omega, 1.);
    DFH.insert("ilum"    , &_ilum);
    
    //   DFH.insert("ilusort", &ilusort,"vectordirection");
    DFH.insert("ilusort", &ilusort,"cuthillmckee");
    DFH.insert("stream_direction",&stream_direction);
    DFH.insert("vector_direction",&vector_direction);
    
    DFH.insert("linear_smooth",     &linear_smooth,      "ilu");
    DFH.insert("bicgstab_residual" ,&bicgstab_residual , "approx");
    DFH.insert("bicgstab_pstep" ,   &bicgstab_pstep ,     0);
    DFH.insert("bicgstab_miniter",  &bicgstab_miniter,    1.);

    DFH.insert("cgmass_maxiter",    &_cgMassMaxIter,      100);
    DFH.insert("cgmass_tol",        &_cgMassTol,          1.e-8);
    DFH.insert("cgmass_globaltol",  &_cgMassGlobalTol,    1.e-14);
    
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(pf,"Solver");
    
    if ((ilusort=="streamdirection")&&(stream_direction.size()==0))
      {
        cerr << "Bei \n\tilusort\tstreamdiretion\nmuss" << endl
             << "\tstream_direction\n"
             << "mit Komponenten, nach denen sortiert wird, "
             << "angegeben werden.\n";
        abort();
      }
  }

  
}
