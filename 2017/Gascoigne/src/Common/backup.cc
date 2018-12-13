/**
*
* Copyright (C) 2004, 2011 by the Gascoigne 3D authors
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


#include "backup.h"
#include <fstream>

using namespace std;

namespace Gascoigne
{
/********************************************************************/

ReadBackUp::ReadBackUp(const string& name, int& size, int& comp)
{
  ifstream file;
  file.open(name.c_str());
  
  if(!file){
    cerr << "backup file '"<< name << "' not found" << endl;
    abort();
  }

  file >> size >> comp;
}

/********************************************************************/

ReadBackUp::ReadBackUp(GlobalVector& u, const string& name)
{
  ifstream file;
  file.open(name.c_str());
  
  if(!file){
    cerr << "backup file '"<< name << "' not found" << endl;
    abort();
  }

  int size, comp;

  file >> size;
  file >> comp;

  //  cout << "BackUp   : reading " << name << ", ";
  //  cout << comp<<" components, "<< size <<" nodes" <<endl;

  if (u.n()!=size)
    {
      cout << "Incompatibility u.n() size u.n()=" << u.n() << " file.n()=" << size << endl;
    }
  assert(u.n()==size);

  int v = max_int(u.ncomp(),comp);
  
  double d;
  for (int i=0; i<u.n(); i++)
    {
      for (int c=0; c<v; c++)  
        {
          double val = 0.;
          file >> val;
          u(i,c) += val;
        }
      for (int c=v; c<comp ;c++)  { file >> d;}
    }
  string test;
  file >> test;

  if(test!="BackUpEnd")
    {
      cout << "error in " <<__FILE__ << ":" << __LINE__ << " : error, test=='"<<test<<"' should be =='BackUpEnd'"<<endl;
      if( u.ncomp()!=comp)
        {
          cout << "probably, because: expected comp nr="<< u.ncomp() <<" NOT EQUAL the bup-file supplied comp nr="<<comp<<endl;
        }
      abort();
    }
}

/********************************************************************/

ReadBackUpResize::ReadBackUpResize(GlobalVector& u, const string& name)
{
  ifstream file;
  file.open(name.c_str());
  
  if(!file){
    cerr << "backup file '"<< name << "' not found" << endl;
    abort();
  }

  int size, comp;

  file >> size;
  file >> comp;

  cout << "BackUp   : reading " << name << ", ";
  cout << comp<<" components, "<< size <<" nodes" <<endl;


  u.ReInit(comp,size);

  if (u.n()!=size)
    {
      cout << "Incompatibility u.n() size " << u.n() << " " << size << endl;
      abort();
    }

  int v = max_int(u.ncomp(),comp);
  
  double d;
  for (int i=0; i<u.n(); i++)
    {
      for (int c=0; c<v; c++)  
	{
	  file >> u(i,c);
	}
      for (int c=v; c<comp ;c++)  { file >> d;}
    }
  string test;
  file >> test;

  if(test!="BackUpEnd")
    {
      cout << "error in " <<__FILE__ << ":" << __LINE__ << " : error, test=='"<<test<<"' should be =='BackUpEnd'"<<endl;
      if( u.ncomp()!=comp)
        {
          cout << "probably, because: expected comp nr="<< u.ncomp() <<" NOT EQUAL the bup-file supplied comp nr="<<comp<<endl;
        }
      abort();
    }
}

/********************************************************************/

WriteBackUp::WriteBackUp(const GlobalVector& u, const string& bname)
{
  string name = bname + ".bup";

  ofstream file;
  file.open(name.c_str());
  file.setf(ios::scientific,ios::floatfield);
  
  if(!file)
  {
    cerr << "BackUp: writing error" << endl;
    exit(10);
  }
  file << u.n() << " " << u.ncomp() << endl;

  file.precision(16);
  for (int i=0; i<u.n(); i++)
    {
      for (int c=0; c<u.ncomp(); c++)  
	{
	  file << u(i,c) << " ";
	}
      file << endl;
    }
  file << "BackUpEnd" << endl;
  file.close();
}

/********************************************************************/

WriteBackUpBinary::WriteBackUpBinary(const GlobalVector& u, const string& bname)
{
  string name(bname);
  name += ".bip";
  
  ofstream file;
  file.open(name.c_str());
  
  if(!file)
    {
      cerr << "BackUp: writing error" << endl;
      exit(10);
    }

  u.BinWrite(file);

  file << "BackUpEnd" << endl;
  file.close();
}

/********************************************************************/

ReadBackUpBinary::ReadBackUpBinary(GlobalVector& u, const string& bname)
{
  string name(bname);

  if(name.substr(name.size()-4,4) != ".bip")
    name += ".bip";

  ifstream file;
  file.open(name.c_str());
  
  if(!file);
  {
    cerr << "backup file '"<< name << "' not found" << endl;
    abort();
  }
  u.BinRead(file);

  cout << "BackUp   : reading " << name << ", ";
  cout << u.ncomp() <<" components, "<< u.n() <<" nodes " << endl;

  string test;
  file >> test;
  if(test!="BackUpEnd")
    {
      cout << "error in " <<__FILE__ << ":" << __LINE__ << " : error, test=='"<<test<<"' should be =='BackUpEnd'"<<endl;
      abort();
    }
}

/********************************************************************/
}
