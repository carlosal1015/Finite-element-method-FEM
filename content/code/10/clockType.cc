#include "clockType.hh"

void clockType::setTime(int hours, int minutes, int seconds)
{
    //Function to set the time.
    //The time is set according to the parameters.
    //Postcondition: hr = hours; min = minutes;
    //          sec = seconds;
    //          The function checks wheter the
    //          values of hours, minutes, and seconds
    //          default value 0 is assigned.
	if (0 <= hours && hours < 24)
	{
		hr = hours;
	}
	else
	{
		hr = 0;
	}
	if (0 <= minutes && minutes < 60)
	{
		min = minutes;
	}
	else
	{
		min = 0;
	}
	if (0 <=seconds && seconds < 60)
	{
		sec = seconds;
	}
	else
	{
		sec = 0;
	}
}