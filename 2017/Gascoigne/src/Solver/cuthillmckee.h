/**
*
* Copyright (C) 2004, 2008, 2009 by the Gascoigne 3D authors
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


#ifndef __cuthillmckee_h
#define __cuthillmckee_h

#include "columnstencil.h"
#include "dynamicstencil.h"


/* CuthillMcKee
 *
 * macht nen CuthillMcKee fuer die ILU.
 * Bedeutet etwa, das immer Knoten mit
 * moeglichst wenig Nachbarn zuerst
 * genommen werden. Die Knoten einer Matrixzeile
 * stehen immer beieinander. Das soll die Bandbreite
 * klein halten.				
 *
 * So Gehts:
 *
 * CuthillMcKee cmc (UnstructuredStencilPointer);
 * cmc.Permutate (perm);   // reiner McKee
 * cmc.Permutate (perm,
 *                Vertex2d(1,0) ) // in x-Richtung
 * 
 */

namespace Gascoigne
{
class CuthillMcKee 
{
    const StencilInterface* S;
    const ColumnStencil*    CS;
    const DynamicStencil*   DS;
    
//     Vertex2d dir2d;
//     Vertex3d dir3d;
    int      dimension;
    std::vector<int> neighbors;

  public:
    CuthillMcKee (const StencilInterface *s);
    CuthillMcKee ();

    void Permutate    (IntVector &perm);
//     void Permutate    (IntVector &perm, const Vertex2d v);
//     void Permutate    (IntVector &perm, const Vertex3d v);
    
//     bool operator()(int i,int j) const;
    #ifdef __WITH_THREADS__
    void Permutate    (IntVector &perm, const IntVector &nodes_in_domain, const std::vector<std::vector<std::pair<int,int> > >& node2domain, int d);
    #endif
};
}

#endif
