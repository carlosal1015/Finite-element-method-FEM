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


#include  "nodematrix.h"
#include  "matrixentrytype.h"

/*-------------------------------------------------------*/
/*-------------------------------------------------------*/
/*-------------------------------------------------------*/

namespace Gascoigne
{
template class NodeMatrix<1 ,MatrixEntryType>;
template class NodeMatrix<2 ,MatrixEntryType>;
template class NodeMatrix<3 ,MatrixEntryType>;
template class NodeMatrix<4 ,MatrixEntryType>;
template class NodeMatrix<5 ,MatrixEntryType>;
template class NodeMatrix<6 ,MatrixEntryType>;
template class NodeMatrix<7 ,MatrixEntryType>;
template class NodeMatrix<8 ,MatrixEntryType>;
template class NodeMatrix<9 ,MatrixEntryType>;
template class NodeMatrix<16,MatrixEntryType>;
template class NodeMatrix<20,MatrixEntryType>;
template class NodeMatrix<25,MatrixEntryType>;
}
