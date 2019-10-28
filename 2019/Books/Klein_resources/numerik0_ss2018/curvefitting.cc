#include <iostream>

#include "hdnum.hh"
#include "solveLeastSquares.hh"

int main ()
{
  // Setup matrix A
  typedef double NumberType;
  hdnum::DenseMatrix<NumberType> A;
  //hdnum::readMatrixFromFile("A6_linear.dat", A);
  hdnum::readMatrixFromFile("A6_quadratic.dat", A);
  A.scientific(false);
  std::cout << "A = " << A << std::endl;

  // Setup x
  std::size_t n = A.colsize();
  hdnum::Vector<NumberType> x(n);

  // Setup b
  hdnum::Vector<NumberType> b, b_orig;
  hdnum::readVectorFromFile("b6.dat", b);
  b.scientific(false);
  b_orig = b;
  std::cout << "b = " << b << std::endl;

  // Solve least squares problem
  solveLeastSquares(A, b, x);
  std::cout << "Solution x = " << x << std::endl;

  return 0;
}
