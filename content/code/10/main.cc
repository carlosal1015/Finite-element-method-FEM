//The program listing of the program that defines
//and uses the class clockType

#include <iostream>
#include "clockType.hh"

int main()
{
	clockType myClock;
	clockType yourClock;
	int hours;
	int minutes;
	int seconds;

	myClock.setTime(5, 4, 30);

	std::cout << "Line 9: myClock: ";
	// myClock.printTime();
	std::cout << std::endl;

	std::cout << "Line 12: yourClock: ";
	// yourClock.printTime();
	std::cout << std::endl;

	yourClock.setTime(5, 45, 16);

	std::cout << "Line 16: After setting, yourClock: ";
	// yourClock.printTime();
	std::cout << std::endl;

	// if (myClock.equalTime(yourClock))
	// {
	// 	std::cout << "Line 20: Both times are equal."
	// 			<< std::endl;
	// }
	// else
	// {
		// std::cout << "Line 22: The two times are not equal."
					// << std::endl;
	// }
	std::cout << "Line 23: Enter the hours, minutes, and "
				<< "seconds: ";
	// std::cin >> hours >> minutes >> seconds;
	std::cout << std::endl;

	myClock.setTime(hours, minutes, seconds);

	std::cout << "Line 27: New myClock: ";
	// myClock.printTime();
	std::cout << std::endl;

	// myClock.incrementHours();

	std::cout << "Line 31: After increamenting myClock by "
			<< "one second, myClock: ";
	// myClock.printTime();
	std::cout << std::endl;

	// myClock.getTime(hours, minutes, seconds);

	// std::cout << "Line 35: hours = " << hours
	// 		<< ", minutes = " << minutes
	// 		<< ", seconds = " << seconds << std::endl;

	return 0;	
}