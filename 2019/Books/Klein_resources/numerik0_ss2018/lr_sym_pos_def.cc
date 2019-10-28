#include <iostream>
#include "hdnum.hh"

#include "lr_no_pivot.hh"


// Funktion zum Aufstellen der Matrix
template<class NumberType>
void flussMatrix( hdnum::DenseMatrix<NumberType> &A )
{
  int M = A.rowsize();
  int N = A.colsize();
  if(M!=N)
    HDNUM_ERROR("Matrix muss quadratisch sein!");

  // Das gross N aus der Aufgabenstellung
  int N_orig = static_cast<int>(std::sqrt(M+1));

  // Die Matrix B aus der Vorlesung.
  hdnum::DenseMatrix<NumberType> B(2*N_orig*(N_orig-1), N);

  /*
    Wir nummerieren die Kanten und Knoten wiefolgt: (Beispiel für N=3)

    R---(0)---0---(1)---1
    |         |         |
   (6)       (7)       (8)
    |         |         |
    2---(2)---3---(3)---4
    |         |         |
   (9)       (10)      (11)
    |         |         |
    5---(4)---6---(5)---7

    wobei R der verwendete Refernzknoten ist. Alle Kanten sind so
    gerichtet, dass sie entweder nach rechts oder nach oben zeigen.

    Der Code stellt die Matrix B aus der Vorlesung auf und rechnet
    B^T*B aus.  Dazu werden erst all horizontalen Kanten, dann alle
    vertikalen Kanten betrachtet.  Man könnte das auch in einer
    Schleife machen und beim Matrixzugriff einen Indexshift
    einbauen...
  */
  for(int i=0; i < N_orig*(N_orig - 1); ++i)
  {
    int horizontal_v_out = i-1+i/(N_orig-1);
    int horizontal_v_in = horizontal_v_out + 1;

    if (horizontal_v_out>=0)
      B[i][horizontal_v_out] = -1.;
    B[i][horizontal_v_in] = +1.;
  }

  for(int i = N_orig*(N_orig - 1); i < 2* N_orig*(N_orig - 1); ++i)
  {
    int vertical_v_out = i-1-N_orig*(N_orig - 1);
    int vertical_v_in = vertical_v_out + N_orig;

    if (vertical_v_out>=0)
      B[i][vertical_v_out] = -1;
    B[i][vertical_v_in] = +1;
  }

  // Das Transponieren einer Matrix ist in hdnum nicht implementiert,
  // da es eine Operation ist die man im Scientific Computing üblicherweise
  // vermeidet (Kopie von großem Objekt erforderlich). Wir sind hier jetzt mal
  // nicht so und implementieren das von Hand.
  hdnum::DenseMatrix<NumberType> B_t(B.colsize(), B.rowsize());
  for(int i=0; i<B_t.rowsize(); ++i)
    for(int j=0; j<B_t.colsize(); ++j)
      B_t[i][j] = B[j][i];

  // Rechne B^T*B, das Ergebnis muss nicht zurückgegeben werden, da wir die Matrix
  // per Referenz übergeben bekommen haben.
  A = B_t*B;
}


int main()
{
  // Setup matrix
  // constexpr int N = 5;
  // constexpr int n = N*N-1;
  const int N = 5;
  const int n = N*N-1;
  typedef double number;
  hdnum::DenseMatrix<number> A(n,n);

  // Pretty printing
  A.scientific(false);
  A.width(15);

  // Fill matrix
  flussMatrix(A);
  // std::cout << A << std::endl;

  // Vector for storing pivot elements
  hdnum::Vector<number> pivotElements(n);

  // LU decomposition of A without pivoting
  lr_no_pivot(A, pivotElements);

  return 0;
}
