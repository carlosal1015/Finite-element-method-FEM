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


#ifndef  __ColumnStencil_h
#define  __ColumnStencil_h

#include  "stencilinterface.h"
#include  "sparsestructureinterface.h"
#include  "gascoigne.h"


namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments ColumnStencil

////
////
/////////////////////////////////////////////

class ColumnStencil : public virtual StencilInterface
{
protected:

  IntVector   scol, sstart;
  
public:

//
////  Con(De)structor 
//

  ColumnStencil() : StencilInterface() {}
  ~ColumnStencil() {}

  const IntVector&  col()    const { return scol; }
        IntVector&  col()          { return scol; }
  const IntVector&  start()  const { return sstart; }
        IntVector&  start()        { return sstart; }

  int  n()            const { return sstart.size()-1;}
  int  nentries()     const { return scol.size();}
  int  rowsize(int i) const { assert((i>=0)&&(i+1<sstart.size())); return sstart[i+1]-sstart[i];}

        int&  col(int pos)           { assert((pos>=0)&&(pos<scol.size())); return scol[pos]; } 
  const int&  col(int pos)     const { assert((pos>=0)&&(pos<scol.size())); return scol[pos]; } 
        int&  start(int i)           { assert((i>=0)&&(i<sstart.size()));   return sstart[i]; } 
  const int&  start(int i)     const { assert((i>=0)&&(i<sstart.size()));   return sstart[i]; } 
        int&  stop(int i)            { assert((i>=0)&&(i+1<sstart.size())); return sstart[i+1]; } 
  const int&  stop(int i)      const { assert((i>=0)&&(i+1<sstart.size())); return sstart[i+1]; } 

  void memory(const SparseStructureInterface*);
  void memory(int n, int nt);
  
  virtual int Find(int i, int j) const
    {
      for(int pos=start(i); pos<stop(i); pos++)
        {
          if (col(pos)==j) return pos;
        }
      std::cerr << "UnstructuredStencil::Find()";
      std::cerr << "no such coupling: "<< i <<" "<<j<<std::endl;
      abort();
      return -1;
    }

  std::ostream& Write(std::ostream& os) const;
};
}

#endif
