/**
*
* Copyright (C) 2006 by the Gascoigne 3D authors
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


#ifndef  __GlobalData_h
#define  __GlobalData_h


#include "gascoigne.h"

namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ////
  ////@brief
  ////  ... comments GlobalData

  ////
  ////
  /////////////////////////////////////////////

  class DataContainer
  {
    private:

    protected:

    public:
      GlobalData          _node;
      GlobalData          _cell;
      GlobalParameterData _parameter;

      //
      ////  Con(De)structor 
      //
      DataContainer() {}
      ~DataContainer() {}

      void AddNodeVector(const std::string& name, const GlobalVector* d) {
        assert(d!=NULL);
        if(!_node.insert(std::make_pair(name,d)).second)
        {
          std::cerr << "NodeVector \"" << name << "\" already added" << std::endl;
          abort();
        }
      }
      void AddCellVector(const std::string& name, const GlobalVector* d) {
        assert(d!=NULL);
        if(!_cell.insert(std::make_pair(name,d)).second)
        {
          std::cerr << "CellVector \"" << name << "\" already added" << std::endl;
          abort();
        }
      }
      void AddParameterVector(const std::string& name, const GlobalParameterVector* d) {
        assert(d!=NULL);
        if(!_parameter.insert(std::make_pair(name,d)).second)
        {
          std::cerr << "ParameterVector \"" << name << "\" already added" << std::endl;
          abort();
        }
      }
      
      void DeleteNodeVector(const std::string& name) {
        _node.erase(name);
      }
      void DeleteCellVector(const std::string& name) {
        _cell.erase(name);
      }
      void DeleteParameterVector(const std::string& name) {
        _parameter.erase(name);
      }

      const GlobalData& GetNodeData() const {
        return _node;
      }
      const GlobalData& GetCellData() const {
        return _cell;
      }
      const GlobalParameterData& GetParameterData() const {
        return _parameter;
      }
  };
}

#endif
