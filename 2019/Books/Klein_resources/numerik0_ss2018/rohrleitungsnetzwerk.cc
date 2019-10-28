// Kompilieren mit:
//
// g++ -I../hdnum/ -o rohrleitungsnetzwerk rohrleitungsnetzwerk.cc
//
// Das setzt voraus, dass ihr Programm in einem Ordner parallel zur
// entpackten HDNum Bibliothek liegt.

#include <iostream>
#include <cstdlib>
#include "hdnum.hh"


// Funktion zum Aufstellen der Matrix
template<class NumberType>
void flussMatrix( hdnum::DenseMatrix<NumberType> &A )
{
  int M( A.rowsize() );
  int N( A.colsize() );
  if(M!=N)
    HDNUM_ERROR("Matrix muss quadratisch sein!");

  // TODO Implementieren Sie hier das Aufstellen der Matrix des
  // Rohleitungsgleichungssystems
}


// Funktion zur Berechnung der Frobenius-Norm einer Matrix
template<class NumberType>
NumberType frobeniusNorm(const hdnum::DenseMatrix<NumberType> &A)
{
  // Error checking
  int M(A.rowsize());
  int N(A.colsize());
  if(M!=N)
    HDNUM_ERROR("Matrix muss quadratisch sein!");

  NumberType result=0.0;

  // TODO Implementieren Sie hier die Frobeniusnorm

  return result;
}


// Funktion zur Berechnung des betragsgrößten Eigenwertes mit Potenzmethode
template<class NumberType>
NumberType maxEigenwert(const hdnum::DenseMatrix<NumberType> &A)
{
  // Error checking
  int M(A.rowsize());
  int N(A.colsize());
  if(M!=N)
    HDNUM_ERROR("Matrix muss quadratisch sein!");

  // TODO Implementieren Sie hier die Potenzmethode
}

// Hauptprogramm
int main(int argc, char ** argv)
{

  // Anzahl der Knoten
  const int N(3);
  std::cout << "Knotenanzahl N: " << N << std::endl;

  // Größe der Matrix
  const int n(N*N-1);

  // Datentyp für die Matrix
  typedef double REAL;

  // Matrix initialisieren
  hdnum::DenseMatrix<REAL> A(n,n);

  // Pretty-printing einmal setzen für alle Matrizen
  A.scientific(false);
  A.width(15);

  flussMatrix(A);
  if (N<=4)
    std::cout << A << std::endl;

  // Bei Schwierigkeiten mit Teilaufgabe a) können Sie Teilaufgaben b)
  // und c) mit folgender Matrix testen
  int size_b = 4;
  hdnum::DenseMatrix<REAL> B(size_b,size_b);
  for (std::size_t i=0; i<size_b; ++i)
  {
    for (std::size_t j=0; j<size_b; ++j)
    {
      B[i][j] = i+j;
    }
  }

  // TODO Geben Sie hier alle verfügbaren Matrixnormen von A (bzw. B)
  // sowie den maximalen Eigenwert aus

  return 0;
}
