#include <iostream>
#include <fstream>

#include "hdnum.hh"
#include "solveLeastSquaresQR.hh"

template <typename T>
T polynomial(hdnum::Vector<T> coeff, const T point)
{
  T evaluation(0.0);
  T variable(1.0);
  for (std::size_t i=0; i<coeff.size(); ++i)
  {
    evaluation += coeff[i]*variable;
    variable *= point;
  }
  return evaluation;
}

int main ()
{
  // Setup matrix A
  typedef double NumberType;
  hdnum::DenseMatrix<NumberType> A;
  //hdnum::readMatrixFromFile("A6_linear.dat", A);
  hdnum::readMatrixFromFile("A6_quadratic.dat", A);
  A.scientific(false);
  // copy matrix for defect calculation afterwards because it gets overwritten
  hdnum::DenseMatrix<NumberType> A_orig(A);
  std::cout << "A = " << A << std::endl;

  // Setup x
  std::size_t n = A.colsize();
  hdnum::Vector<NumberType> x(n);

  // Setup b
  hdnum::Vector<NumberType> b, b_orig;
  hdnum::readVectorFromFile("b6.dat", b);
  b.scientific(false);
  // copy RHS because it gets overwritten
  b_orig = b;
  std::cout << "b = " << b << std::endl;

  // Solve least squares problem
  solveLeastSquares(A, b, x);
  std::cout << "Solution x = " << x << std::endl;

  // Defect
  b_orig *= -1.0;
  A_orig.umv(b_orig,x);
  std::cout << "quadratic error = ||Ax-b||_2 = " << b_orig.two_norm_2() << std::endl;

  return 0;
}
