/**
*
* Copyright (C) 2004, 2005, 2008 by the Gascoigne 3D authors
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


#include  "monitor.h"
#include  <string>
#include  <fstream>
#include  <cstdio>
#include  <cstring>
#include  "filescanner.h"
#include  "stlio.h"


using namespace std;

/*****************************************************************/

namespace Gascoigne
{
Monitor::Monitor() : prec(6)
{
  aos = "Gascoigne";
  bos = "Adaptive";
  directory = ".";
  ps  = 1;
}

/*******************************************************************/

void Monitor::PrintInfoSummary(const CGInfo& info) const
{
  cout.precision(3);
  cout.setf(ios::fixed,ios::floatfield);
  cout << info.statistics().totaliter() << "\t";

  cout.setf(ios::scientific,ios::floatfield);
  cout.precision(2);
  cout << " { ";
  cout << info.control().residual() << " ";
  cout << info.control().firstresidual() << " (";
  cout << info.control().correction() << ") ";
  cout << "}\t";

  cout.precision(2);
  cout.setf(ios::fixed,ios::floatfield);
  cout << info.statistics().rate() << "\n";
}

/*******************************************************************/

void Monitor::PrintInfoSummary(const NLInfo& nlinfo) const
{
  const CGInfo& cginfo = nlinfo.GetLinearInfo();

  cout.precision(3);
  cout.setf(ios::fixed,ios::floatfield);
  cout << nlinfo.statistics().totaliter() << " [";
  cout << cginfo.statistics().totaliter() << "] M ";
  cout << nlinfo.statistics().totalmatrix() << " ";

  cout.setf(ios::scientific,ios::floatfield);
  cout.precision(2);
  cout << " { ";
  cout << nlinfo.control().residual() << " ";
  cout << nlinfo.control().firstresidual() << " (";
  cout << nlinfo.control().correction() << ") ";
  cout << "}\t";

  cout.precision(2);
  cout.setf(ios::fixed,ios::floatfield);
  cout << nlinfo.statistics().rate() << " ";
  cout << cginfo.statistics().rate() << "\n";
}

/*****************************************************************/

void Monitor::set_directory(const string& dir)
{
  //cout << "Monitor::set_directory()\tdirectory = " << dir << endl;
  directory = dir;
  texfile   = directory;
  texfile += "/TexResults";
  numfile   = directory;
  numfile += "/NumResults";
  protokoll = directory;
  protokoll += "/Protokoll";

  string command("mkdir -p ");
  command += dir;
  system(command.c_str());

//  ofstream tfile( texfile.c_str() );
//  if(!tfile)
//    {
//      error_io(texfile);
//    }
//  else
//    {
//      PrintHeader(tfile);
//    }
//  tfile.close();
  
//  ofstream pfile( protokoll.c_str() );
//  if(!pfile)
//    {
//      error_io(protokoll);
//    }
//  pfile.close();
  
//  ofstream nfile( numfile.c_str() );
//  if(!nfile)
//    {
//      error_io(numfile);
//    }
//  nfile.close();
  
  
  new_message.str("");
  

  //cout << "set_directory:" << protokoll <<  "#" << endl;
}

/*****************************************************************/

void Monitor::init(const ParamFile* pf, int c)
{
  control   = c;

  string dir;
  DataFormatHandler   DFH;
  DFH.insert("directory"  , &dir,  ".");
  DFH.insert("format"  , &format,  "latex");
  DFH.insert("header"  , &header);
  FileScanner FS(DFH, pf, "Monitor");

  if (header.size()>0) cout << "header:\n" << header << endl;

  if(dir!=".")  set_directory(dir);
}

/*****************************************************************/

void Monitor::error_io(const string& s) const
{
  cout << "Monitor::error_io()\n";
  cout << "Fehler beim Oeffnen von " << s << endl;
  cout << "directory = " << directory << endl;
  abort();
}

/*****************************************************************/

void Monitor::pre_monitor(char* s)
{
  if(control>2)
    {
      m_length[m_counter++] = strlen(s);
      cout << s;
    }
}

/*****************************************************************/

void Monitor::post_monitor()
{
  if(control>2)
    {
      int l = m_length[--m_counter];
      for(int i=0; i<l; i++) cout << "\b";
      for(int i=0; i<l; i++) cout << " " ;
      for(int i=0; i<l; i++) cout << "\b";
    }
}

/*****************************************************************/

void Monitor::print_message()
{
//  ofstream pfile( protokoll.c_str(),ios::app);

//  if(!pfile)
//  {
//    error_io(protokoll);
//  }

  if (control) cout << new_message.str() << endl;

  //  pfile << message << endl;

  new_message.str(""); 

//  pfile.close();
}

/*******************************************************************/

void Monitor::failed_step()
{
  print_message();


  new_message << " repeat";
  print_message();
}
  
/*******************************************************************/

void Monitor::mesh(int nc, int nl)
{
  new_message << " (N=" << nc << "[" << nl << "])";
  print_message();
}

/*******************************************************************/

void Monitor::post_nonlinear(const DoubleVector& Ju, double eta, int ncells,
			     int cg, int dcg)
{
  new_message.str(""); 
  new_message << "cg dnl: " << cg << " " << dcg 
	      << " [J(u)=" << Ju[0] << " eta=" << eta << "]";  
//   sprintf (message,"cg dnl: %d %d",cg,dcg);
//   sprintf (message,"%s [J(u)=%10.7e eta=%7.3e]",message,Ju[0],eta);
  
  print_message();
}

/*******************************************************************/

void Monitor::pre_nonlinear(int i)
{
  print_message();
  new_message.str("----------------------------------------");
  print_message();
  new_message.str("");
  new_message << aos.c_str() << "-" << bos.c_str() << ": ";
  new_message.width(4); new_message << i;

  //  sprintf(message,"%s-%s: %4d",aos.c_str(),bos.c_str(),i);
}

/*******************************************************************/

/*******************************************************************/

void Monitor::matrix_info(int i)
{
  if(!i)
    {
      for(int j=0; j<ps; j++) 
	new_message << " ";
    }
  else
    {
      if (_newmatrix) 
	new_message << "M";
      else            
	new_message << " ";
    }
}

/*******************************************************************/

void Monitor::nonlinear_step(const CGInfo& cginfo, const NLInfo& nlinfo)
{
  int i = nlinfo.control().iteration();

  matrix_info(i);

  if (!ps || (i%ps)) return;

  int  iter = cginfo.control().iteration(); 

  double r = nlinfo.control().residual();
  new_message.width(3);
  new_message << i << ": ";
  new_message.unsetf(ios::fixed);
  new_message.setf(ios::scientific); new_message.precision(2); 
  new_message << r;

  if (i && (control>1))
    {
      double c = nlinfo.control().correction();
      new_message << " " << c ;
    }
  if (i)
    {
      double rr = nlinfo.statistics().rate();
      double lr = nlinfo.statistics().lastrate();
      int    p  = nlinfo.control().relax();
      
      //sprintf(message,"%s [%4.2f %4.2f]",message,lr,rr);
      new_message << " [";
      new_message.unsetf(ios::scientific);
      new_message.setf(ios::fixed);
      new_message.width(3);
      new_message << lr << " " << rr << "]";
      if (p)
	{
	  new_message << " " << p;
	}
      else
	{
	  new_message << "  ";
	}
      
      rr = cginfo.statistics().rate();
      lr = cginfo.control().residual();
      double lc = cginfo.control().correction();
      
      new_message << " - ";
      new_message.setf(ios::scientific);
      new_message.unsetf(ios::fixed);
      new_message.width(4);
      new_message << lr;

      new_message.unsetf(ios::scientific);
      new_message.setf(ios::fixed);
      new_message.precision(3);	  
      if (control>1)
	{
	  new_message << " " << lc; 
	}
      new_message << " [" << rr << "] {" << iter;
      if (cginfo.control().status()=="running") 
	new_message << "*";
      new_message << "}";
    }
  if (cginfo.control().status()=="diverged")
    {
      new_message << " @";
    }
  print_message();
}

/*******************************************************************/

void Monitor::PrintResults(const string& s) const
{
//  ofstream   file(texfile.c_str(),ios::app);
//  if(!file)  error_io(texfile);
//  if(format=="latex")   PrintAscii(file,s);
  // ausser betrieb
//   if(format=="html")    PrintHtml (file,iv,dv,"\n");
//   if(format=="ascii")   PrintHtml (file,iv,dv,"\n");
//  file.close();
}

/*******************************************************************/

void Monitor::PrintResults(int i) const
{
//  ofstream   file(texfile.c_str(),ios::app);
//  if(!file)  error_io(texfile);
//  if(format=="latex")   PrintAscii(file,i);
//  file.close();
}

/*******************************************************************/

void Monitor::PrintResults(double d) const
{
//  ofstream   file(texfile.c_str(),ios::app);
//  if(!file)  error_io(texfile);
//  if(format=="latex")   PrintAscii(file,d);
//  file.close();
}

/*******************************************************************/

void Monitor::PrintResults(const IntVector& iv) const
{
//  ofstream   file(texfile.c_str(),ios::app);
//  if(!file)  error_io(texfile);
//  if(format=="latex")   PrintAscii(file,iv);
//  file.close();
}

/*******************************************************************/

void Monitor::PrintResults(const DoubleVector& dv) const
{
//  ofstream   file(texfile.c_str(),ios::app);
//  if(!file)  error_io(texfile);
//  if(format=="latex")   PrintAscii(file,dv);
//  file.close();
}

/*******************************************************************/

void Monitor::Print(const string& s, string se) const
{
//  ofstream   file(texfile.c_str(),ios::app);
//  if(!file)  error_io(texfile);

//  file << s << se;
//  cout << s << se;
//  file.close();
}

/*******************************************************************/

void Monitor::PrintHeader(ostream& os) const
{
  if(header.size()==0) return;
  for (int i=0; i<header.size()-1; i++) os << header[i] << " , "; 
  os << header[header.size()-1] << endl;
}


/*******************************************************************/

void Monitor::PrintAscii
(ostream& os, const IntVector& iv) const
{
  os.precision(prec);
  os << iv;
}

/*******************************************************************/

void Monitor::PrintAscii
(ostream& os, const DoubleVector& dv) const
{
  os.precision(prec);
  os.setf(ios::scientific,ios::floatfield);
  os << dv;
}

/*******************************************************************/

void Monitor::PrintAscii
(ostream& os, const string& s) const
{
  os << s;
}
}
