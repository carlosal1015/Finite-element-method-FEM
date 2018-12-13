/**
*
* Copyright (C) 2004, 2005, 2006 by the Gascoigne 3D authors
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


#ifndef  __BasicDiscretization_h
#define  __BasicDiscretization_h


#include  "discretizationinterface.h"
#include  "feminterface.h"
#include  "integratorinterface.h"

namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  ... comments BasicDiscretization

///
///
/////////////////////////////////////////////

class BasicDiscretization : public DiscretizationInterface
{
 private:
   const MeshInterface*  __MP;
   mutable DataContainer __q;
  
 protected:
   mutable EntryMatrix __E;
   mutable LocalVector __F;
   mutable LocalVector __U;

   mutable LocalData           __QN;
   mutable LocalData           __QC;
   mutable LocalParameterData  __QP;
   
   virtual const DataContainer& GetDataContainer() const {return __q;}
   virtual void SetDataContainer(const DataContainer& q) const {__q = q;}

   virtual const MeshInterface* GetMesh() const { assert(__MP); return __MP;}

   virtual void GlobalToGlobalData() const;
   virtual void GlobalToLocal(LocalVector& U, const GlobalVector& u, int iq) const {
     GlobalToLocalSingle(U,u,iq);
     GlobalToLocalData(iq);
   }
   virtual void GlobalToLocalData(int iq) const;
   virtual void GlobalToLocalSingle(LocalVector& U, const GlobalVector& u, int iq) const;
   virtual void GlobalToLocalCell(LocalVector& U, const GlobalVector& u, int iq) const;

   virtual void LocalToGlobal(GlobalVector& f, const LocalVector& F, int iq, double s) const;
   virtual void LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const;

   virtual IntVector GetLocalIndices(int iq) const=0;

 public:
   //
   ////  Constructor 
   //
   BasicDiscretization();
   ~BasicDiscretization();
   
   void BasicInit(const ParamFile* pf) {}
   void ReInit   (const MeshInterface* MP) {__MP=MP;}
   
   virtual void AddNodeVector(const std::string& name, const GlobalVector* q) const {
     __q.AddNodeVector(name,q);
   }
   virtual void DeleteNodeVector(const std::string& name) const {
     __q.DeleteNodeVector(name);
   }
   virtual void AddCellVector(const std::string& name, const GlobalVector* q) const {
     __q.AddCellVector(name,q);
   }
   virtual void DeleteCellVector(const std::string& name) const {
     __q.DeleteCellVector(name);
   }
   virtual void AddParameterVector(const std::string& name, const GlobalParameterVector* q) const {
     __q.AddParameterVector(name,q);
   }
   virtual void DeleteParameterVector(const std::string& name) const {
     __q.DeleteParameterVector(name);
   }

   void HNAverageData() const;
   void HNZeroData   () const;
};
}

#endif
