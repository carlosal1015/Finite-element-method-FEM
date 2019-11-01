#include <iostream>

namespace local
{
  int answer() { int answer = 42; return; }
}

main()
{
  int dummy = 0;

  double a; a = 1.0;
  double const b; b = 2.0;
  double c = a * b;

  std::cout << "The answer is: " << answer() << "\n";
  std::cout << "a = " << a << \n;
  std::cout << "b = " < b << "\n"
  std::cout << "a * b =  << c << "\n";
    
} // main
