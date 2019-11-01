#include <iostream>
#include "hdnum.hh"

#include "gauss.hh"

int main ()
{
  typedef float number;
  // constexpr int n = 7;
  const int n = 7;

  hdnum::Vector<number> x(n),b(n);
  hdnum::DenseMatrix<number> A(n,n);
  hdnum::fill(x,number(1.0),number(1.0));
  hdnum::vandermonde(A,x);
  A.mv(b,x);
  x = number(0.0);
  gauss(A,x,b);
  std::cout << x << std::endl;
}
