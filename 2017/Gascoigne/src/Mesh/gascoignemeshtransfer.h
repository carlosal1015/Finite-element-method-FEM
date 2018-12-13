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


#ifndef  __GascoigneMeshTransfer_h
#define  __GascoigneMeshTransfer_h

#include  "meshtransferinterface.h"
#include  "fixarray.h"
#include  "gascoigne.h"
#include  <map>

/*-----------------------------------------*/

namespace Gascoigne
{
class GascoigneMeshTransfer : public MeshTransferInterface
{
protected:
  
  std::map<int,fixarray<2,int> >  zweier;
  std::map<int,fixarray<4,int> >  vierer;
  std::map<int,fixarray<8,int> >  achter;
  
  IntVector                c2f;
  std::map<int,int>             CellEiner;
  std::map<int,fixarray<4,int> >  CellVierer;
  std::map<int,fixarray<8,int> >  CellAchter;
  
public:
  
  const std::map<int,fixarray<2,int> >& GetZweier() const {return zweier;}
  const std::map<int,fixarray<4,int> >& GetVierer() const {return vierer;}
  const std::map<int,fixarray<8,int> >& GetAchter() const {return achter;}
  const IntVector&                   GetC2f()    const {return c2f;}
  
  std::map<int,fixarray<2,int> >& GetZweier() {return zweier;}
  std::map<int,fixarray<4,int> >& GetVierer() {return vierer;}
  std::map<int,fixarray<8,int> >& GetAchter() {return achter;}
  IntVector&                   GetC2f() {return c2f;}

  const std::map<int,int>             & GetCellEiner ()const  {return CellEiner;}
  const std::map<int,fixarray<4,int> >& GetCellVierer()const  {return CellVierer;}
  const std::map<int,fixarray<8,int> >& GetCellAchter()const  {return CellAchter;}

  std::map<int,int>             & GetCellEiner () {return CellEiner;}
  std::map<int,fixarray<4,int> >& GetCellVierer() {return CellVierer;}
  std::map<int,fixarray<8,int> >& GetCellAchter() {return CellAchter;}
  
  GascoigneMeshTransfer();
};
}

#endif
