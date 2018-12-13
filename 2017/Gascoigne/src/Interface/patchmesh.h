/**
*
* Copyright (C) 2004, 2006 by the Gascoigne 3D authors
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


#ifndef  __PatchMesh_h
#define  __PatchMesh_h

#include  "meshinterface.h"

/*-----------------------------------------*/

namespace Gascoigne
{
  class PatchMesh : public virtual MeshInterface
  {
  
    private:

    protected:

    public:
      PatchMesh() {};
      virtual ~PatchMesh() {}
 
      virtual bool       HasPatch()                        const=0;
      virtual bool       HasQ4Patch()                      const=0;
      virtual int        npatches()                        const=0;
      virtual int        nq4patches()                      const=0;
      virtual int        nodes_per_patch()                 const=0;
      virtual int        nodes_per_q4patch()               const=0;
      virtual IntVector  Q2IndicesOfQ4Patch(int iq)        const=0;
      virtual IntVector  CoarseIndices(int iq)             const=0;
      virtual IntVector  CoarseIndicesQ4(int iq)           const=0;
      virtual const IntVector* IndicesOfPatch(int i)       const=0;
      virtual const IntVector* IndicesOfQ4Patch(int i)     const=0;
      virtual const IntVector* VertexOnBoundary(int color) const=0;
      virtual const IntVector* CellOnBoundary(int color)   const=0;
      virtual const IntVector* LocalOnBoundary(int color)  const=0;

      virtual bool CellIsCurved(int iq)                    const {
        return false;
      }

      // MPI
      virtual void send(int p) const {
        std::cerr << "\"PatchMesh::send\" not written!" << std::endl;
        abort();
      }
      virtual void recv(int p) {
        std::cerr << "\"PatchMesh::recv\" not written!" << std::endl;
        abort();
      }
  };
}

#endif
