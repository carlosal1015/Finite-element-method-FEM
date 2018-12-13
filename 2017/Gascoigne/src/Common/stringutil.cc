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


#include  "stringutil.h"
#include  "gostream.h"

#include  <iostream>
#include  <fstream>
#include  "stlio.h"

#ifdef __OLDCOMPILER__
#include  <strstream>
#include  <stdio.h>
#define STRINGSTREAM  strstream 
#define ISTRINGSTREAM istrstream 
#else
#include  <sstream>
#define STRINGSTREAM  stringstream
#define ISTRINGSTREAM istringstream
#endif

using namespace std;

/*--------------------------------------*/

namespace Gascoigne
{
string Int2String   (int a   )
{
#ifdef __OLDCOMPILER__
  int n = 1+(int) log(a);
  char* c;
  c = new char[n];
  sprintf(c,"%d",a);
  string r = c;
  delete c;
  return r;
#else
  STRINGSTREAM ss;
  ss << a;
  string r;
  r = ss.str();
  return r;
#endif
}

/*--------------------------------------*/

string Double2String(double a)
{
  STRINGSTREAM ss;
  ss << a;
  string r;
  r = ss.str();
  return r;
}

/*--------------------------------------*/

string GetBase(const char* buf, char sep)
{
  vector<string> all= StringSplit(buf,sep);
  assert(all.size());
  return all[0];
}

string GetTail(const char* buf, char sep)
{
  vector<string> all= StringSplit(buf,sep);
  assert(all.size());
  return all[all.size()-1];
}

/*--------------------------------------*/

vector<string> StringSplit(const char* buf, char sep)
{
  vector<string> words;
  ISTRINGSTREAM is(buf);
  while (!is.eof())
    {
      string t;
      getline(is,t,sep);
      if(t.size()==0) continue;
      if(t[0]!=sep) words.push_back(t);
    }
  return words;
}

vector<string> StringSplit(const char* buf, char sep1, char sep2)
{
  vector<string> words;
  ISTRINGSTREAM is(buf);
  while (!is.eof())
    {
      string t;
      getline(is,t,sep1);
      if(t.size()==0) continue;
      if(t[0]!=sep1) words.push_back(t);
    }
  vector<string> words2;
  for(int i=0;i<words.size();i++)
    {
      vector<string> s1 = StringSplit(words[i].c_str(),sep2);
      copy(s1.begin(),s1.end(),back_inserter(words2));
    }
      
  return words2;
}

/*-----------------------------------------*/

pair<string,vector<string> > SplitArgs(string s)
{
  vector<string> vs = StringSplit(s.c_str(),'_');
  if(!vs.size())
    {
      cerr << "SplitArgs()\t";
      cerr << s << " --> " <<  vs << endl;
      abort();
    }
  string name = vs[0];
  vector<string> args(vs.size()-1);
  for(int i=1;i<vs.size();i++)
    {
      args[i-1] = vs[i];
    }
  return make_pair(name,args);
}
}
