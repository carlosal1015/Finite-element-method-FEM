/**
*
* Copyright (C) 2004, 2005, 2011 by the Gascoigne 3D authors
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


#ifndef __visudatainfo_h
#define __visudatainfo_h

#include  "fixarray.h"
#include  "visudata.h"

#include  <map>
#include  <string>

/*-------------------------------------------------------------------------*/

namespace Gascoigne
{
class VisuDataInfo
{
 protected:

  std::map<std::string,int>              scalars;
  std::map<std::string,fixarray<3,int> > vectors;
  std::map<std::string,int>              scalar_order;
  std::map<std::string,int>              vector_order;

 public:

  typedef std::map<std::string,int>::const_iterator                siterator;
  typedef std::map<std::string,fixarray<3,int> >::const_iterator   viterator;

  VisuDataInfo() {}
  VisuDataInfo(int ncomp) { AddScalars(ncomp);}
  VisuDataInfo(const VisuData& V, std::string def="U");
  VisuDataInfo(const VisuDataInfo& V) : scalars(V.Scalars()), vectors(V.Vectors()) {}
  VisuDataInfo& operator=(const VisuDataInfo& V);

  bool operator!=(const VisuDataInfo& V) const;

  void Clear() {
    scalar_order.clear();
    vector_order.clear();
    scalars.clear();
    vectors.clear();
  }

  siterator GetSIterator(int i) { 
    for(siterator p = sbegin() ; p!= send() ; p++){
      std::string s = p->first;
      if ( scalar_order[s]==i ) return p;
    }
    abort();
  }
  viterator GetVIterator(int i) {
    for(viterator p = vbegin() ; p!= vend() ; p++){
      std::string s = p->first;
      if ( vector_order[s]==i ) return p;
    }
    abort();
  }

  void AddScalar(int index,const std::string& name, int i)                    {scalar_order[name]=index;scalars[name]=i;}
  void AddVector(int index,const std::string& name, const fixarray<3,int>& i) {vector_order[name]=index;vectors[name]=i;}

  void AddScalars(int ncomp, std::string def="U");

  int nscalars() const {return scalars.size();}
  int nvectors() const {return vectors.size();}

  const std::map<std::string,int>&              Scalars() const {return scalars;}
  const std::map<std::string,fixarray<3,int> >& Vectors() const {return vectors;}

  siterator sbegin() const {return scalars.begin();}
  siterator send  () const {return scalars.end();}
  viterator vbegin() const {return vectors.begin();}
  viterator vend  () const {return vectors.end();}
};
}

#endif
