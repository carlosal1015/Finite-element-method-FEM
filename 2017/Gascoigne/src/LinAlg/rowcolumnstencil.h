/**
*
* Copyright (C) 2004 by the Gascoigne 3D authors
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


#ifndef  __RowColumnStencil_h
#define  __RowColumnStencil_h

#include  "columnstencil.h"


namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments RowColumnStencil

////
////
/////////////////////////////////////////////

class RowColumnStencil : virtual public StencilInterface, public ColumnStencil
{
private:


protected:

  nvector<int>   srow;

public:


//
////  Con(De)structor 
//

  RowColumnStencil() : ColumnStencil() {}
  ~RowColumnStencil() {}

  const nvector<int>&  row()    const { return srow; }
        nvector<int>&  row()          { return srow; }
        int&  row(int i)           { return srow[i]; } 
  const int&  row(int i)     const { return srow[i]; } 
  void memory(int n, int nt);
  void memory(const SparseStructureInterface*);
  
  virtual int Find(int i, int j) const
    {
      for(int ii=0;ii<n();ii++)
        {
          if(row(ii)==i)
            {
              return ColumnStencil::Find(ii,j);
            }
        }
      std::cerr << "ColumnStencil::Find()";
      std::cerr << "no such coupling: "<<i <<" "<<j<<std::endl;
      abort();
      return -1;
    }

  std::ostream& Write(std::ostream& os) const
  {
    os << n() << "\t" << nentries()<<std::endl<<std::endl;
    os << sstart<<std::endl<<std::endl;
    for(int i=0;i<n();i++)
      {
        for(int pos=start(i);pos<stop(i);pos++)
          {
            os << row(i) << ":" << col(pos)  << " ";
          }
        os << std::endl;
      }
    return os;
  }
};
}

#endif
