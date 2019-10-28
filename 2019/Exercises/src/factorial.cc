#include <iostream>

int factorial(unsigned int number)
{
    if (number == 0)
      return 1;
    else
      return number*factorial(number - 1);
}


int main(int argc, char** argv)
{
  unsigned int number;
  std::cout << "Ingrese el número n: " ;
  std::cin >> number;

  std::cout << "n! es " << factorial(number) << "." << "\n";

  return 0;
}
