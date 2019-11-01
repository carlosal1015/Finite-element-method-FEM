#include <iostream>

int factorial(unsigned int number)
{
  if (number == 0)
    return 1;
  else
    return number * factorial(number - 1);
}


double euler(unsigned int n)
{
  float summation = 0;
  for(int k = 0; k < n; k++)
    summation += 1 / factorial(n);

  return summation;
}

int main(int argc, char** argv)
{
  unsigned int n;
  std::cout << "Por favor ingrese el valor de n: ";
  std::cin >> n;

  std::cout << "El valor de $e_n$ es " << euler(n) << std::endl;

  return 0;
}
//Turek
