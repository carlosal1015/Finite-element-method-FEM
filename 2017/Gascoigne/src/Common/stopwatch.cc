/**
*
* Copyright (C) 2004, 2008 by the Gascoigne 3D authors
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


#include "stopwatch.h"

//#include  <strstream>
#include  <stdlib.h>
#include  <stdio.h>
#include  <sys/times.h>
#include  <sys/types.h>
#include  <sys/wait.h>
#include  <unistd.h>

using namespace std;

/***********************************************************
 *****   ueberarbeitete und fehlerbereinigte Version   ***** 
 *****   Michael Schmich                  15.01.2002   *****
 ***********************************************************/

/*----------------------------------------------------*/

namespace Gascoigne
{

double Time::GetTotalSeconds () const 
{
  return sec+60.*min+3600.*hour+24.*3600.*day;
}

/*----------------------------------------------------*/

void Time::add(double s)
{
  sec += s;
  int su = static_cast<int>(sec)/60;

  sec -= 60.*su;
  min += su;

  hour += min / 60;
  min   = min % 60;

  day += hour / 24;  
  hour = hour % 24;

}

/*----------------------------------------------------*/

ostream& operator<<(ostream &s, const Time& A)
{
  if(A.GetDays()>0)
    {
	s.precision(1);
	s << A.GetDays() << "-";
    }
  s.precision(2);
  if(A.GetHours()>0)
    {
	s << A.GetHours() << "-";
    }
  s << A.GetMinutes() << ":" << A.GetSeconds();
  return s;
}

/*----------------------------------------------------*/

StopWatch::StopWatch() : running(0), last_time(0), T() {}

/*----------------------------------------------------*/

void StopWatch::reset() 
{ 
  running = 0; last_time = 0; T.reset(); 
}

/*----------------------------------------------------*/

void StopWatch::start() 
{ 
  if (!running) { last_time = clock(); running = 1;}
}

/*----------------------------------------------------*/

double StopWatch::stop()  
{ 
  if (running) 
    {
	double s = static_cast<double>(clock() - last_time) / static_cast<double>(CLOCKS_PER_SEC);
	add(s);
	running = 0;
    }
  return T.GetTotalSeconds(); 
}


/*----------------------------------------------------*/

double StopWatch::read() const  
{
  if (running) return -1;
  return T.GetTotalSeconds(); 
} 


double StopWatch::read100() const  
{
  if (running) return -1;
  double X = static_cast<double> (0.01 * static_cast<int> (100 * read()));
  return X;

}

}
