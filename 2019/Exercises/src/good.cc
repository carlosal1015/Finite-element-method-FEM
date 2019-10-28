#include <iostream>

//namespace local

int answer()
{
  int answer = 42;
  return answer;
}

int main()
{
  //  int dummy = 0;

  double a = 1.0;
  double const b = 2.0;
  double c = a * b;

  std::cout << "The answer is: " << answer() << "\n";
  std::cout << "a = " << a << "\n";
  std::cout << "b = " << b << "\n";
  std::cout << "a * b = " << c << "\n";

  return 0;
} // main
