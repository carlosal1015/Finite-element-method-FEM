#include <iostream>  
#include <cmath>
#include "hdnum.hh" // hdnum header

/*
 * The class that contains all important parts of solving the Poisson problem.
 * 
 */
class Poisson{
public:
  Poisson(double h, double alpha,double a, int p = 1)
  :h(h),
   alpha(alpha),
   a(a),
   p(p)
  {}
  
  void run(); 
private:
  void make_grid();
  void setup_system();
  void assemble_system();
  void solve();
  void output_results() const;
  
  int    p;       //polynomial degree
  int    n_dofs;  //number of degrees of freedom 
  int    n_elems; //number of finite elements (intervals)
  
  double h;       //grid size
  double alpha;   //model parameter
  double a;       //rhs factor
  
  hdnum::Vector<double> rhs;      // right hand side (f)
  hdnum::Vector<double> solution; // solution vector (y)
  hdnum::Vector<double> grid;     // vector of grid points (x)
  
  hdnum::DenseMatrix<double> system_matrix; //system matrix   (A)
  
};

/*
 * This function calculates the subdivision of [0,1] 
 * based on the given grid size.
 * The positions of the grid points x_i are calculated
 * and saved in grid.
 */
void Poisson::make_grid()
{
   n_elems = 1.0/h; 
   std::cout << "==================== Making Grid ====================\n";
   if ( p == 1 )
   {
     n_dofs = n_elems+1;
     grid = hdnum::Vector<double>(n_dofs);
     for ( int i = 0; i < n_dofs ; i++ )
     {
       grid[i] = i*h;
     }
   }
   else{
     /*
      * TODO: 
      * 1) Change the #DoFs to match quadratic elements
      * 2) Change grid to also contain the x-values of the interval-midpoints
      * 3) Remove the existing two lines since everything is now implemented
      */
     std::cout << "Grid for p=2 not yet implemented" << std::endl;
     exit(0);
   }
   std::cout << "\t #DoFs:     " << n_dofs << "\n"
             << "\t #Elems:    " << n_elems << "\n"
             << "\t Grid size: " << h << "\n"
             << "Done." << std::endl;
}
/*
 * This functions initializes 
 * the system matrix and 
 * right hand side (rhs) to the right sizes i.e. #DoFs
 */
void Poisson::setup_system()
{
  std::cout << "==================== Setup System ===================\n";
  
  rhs = hdnum::Vector<double>(n_dofs,0.0); 
  solution = hdnum::Vector<double>(n_dofs,0.0); 
  system_matrix = hdnum::DenseMatrix<double>(n_dofs,n_dofs,0.0);
  
  std::cout << "Done." << std::endl;
}
/*
 * This is the "heart" of the Problem
 * Here the Matrix A and right hand side f are assembled according to
 * 1) The Problem
 * 2) The polynomial degree of the finite elements
 * 
 */
void Poisson::assemble_system()
{
  std::cout << "==================== Assemble System ================\n";
  std::cout << "Assembly" << std::endl;
  double h_inv = 1.0/h;
  
  //First assemble the matrix and right hand side by going over all elements (subintervals)
  for ( int j = 0 ; j < n_elems ; j++ )
  {
    if ( p == 1 )
    {
      //Index of the grid point on the left end of the interval/element
      int i = j; 
      
      //Add the local system matrix of the current (linear) element to the global system matrix 
      system_matrix[i][i]     += h_inv;
      system_matrix[i+1][i]   += -h_inv;
      system_matrix[i][i+1]   += -h_inv;
      system_matrix[i+1][i+1] += h_inv;
      //Add the local rhs to the global rhs
      rhs[i]   += -a*h/2.0;
      rhs[i+1] += -a*h/2.0;
    }
    else{
      /*
       * TODO:
       * 1) The index for the grid point on the left end of each element is different for p=2
       * 2) The quadratic element matrix has to be added instead of the linear one
       * 3) Same for the right hand side
       * 4) Remove the existing two lines since everything is now implemented
       */
      std::cout << "Assembly for p=2 not yet implemented" << std::endl;
      exit(0);
    }
  } 
  // Insert the model parameter into the system matrix
  // Be careful: this only works in this easy fashion
  // because alpha is constant. In the other, more
  // general case, you need to integrate alpha locally
  // in the element matrix.
  system_matrix *= alpha;
  std::cout << "Applying boundary conditions" << std::endl;
  //Now apply zero Dirichlet boundary conditions
  rhs[0] = 0.0;
  rhs[n_dofs-1] = 0.0;
  system_matrix[0][1] = 0.0;
  system_matrix[1][0] = 0.0;
  system_matrix[n_dofs-1][n_dofs-2] = 0.0;
  system_matrix[n_dofs-2][n_dofs-1] = 0.0;
  if ( p == 2 )
  {
    /*
     * TODO:
     * 1) Remove further entries coupling the left and right boundary to the other grid points
     * 2) Remove the existing two lines since everything is now implemented
     */
    std::cout << "Additional boundary conditions for p=2 not yet implemented" << std::endl;
    exit(0);
  }
  
  std::cout << "Done." << std::endl;
}
/*
 * This functions solves the linear equation system
 * A*y=f 
 * using LU-decomposition
 * 
 */
void Poisson::solve()
{
  std::cout << "==================== Solving System =================\n";
  hdnum::Vector<long unsigned int> p(n_dofs);
  hdnum::Vector<long unsigned int> q(n_dofs);
  hdnum::Vector<double> s(n_dofs);
  hdnum::row_equilibrate(system_matrix,s);
  hdnum::lr_fullpivot(system_matrix,p,q);
  hdnum::apply_equilibrate(s,rhs);
  hdnum::permute_forward(p,rhs);
  hdnum::solveL(system_matrix,rhs,rhs);
  hdnum::solveR(system_matrix,solution,rhs);
  hdnum::permute_backward(q,solution);
  std::cout << "Done." << std::endl;
}
/*
 * This function writes the grid and solution
 * in a file that can be plotted using gnuplot
 * All parameters are added to the filename
 * 
 */
void Poisson::output_results() const
{
//   std::cout << solution << std::endl;
  std::cout << "==================== Output of Solution ==============\n";
  std::ofstream myfile;
  std::ostringstream filename;
  filename << "solution_h_" << h << "_alpha_"<< alpha << "_a_" << a << "_p_" << p << ".txt";
  myfile.open(filename.str());
  myfile << "plot '-' u 1:2 with lines" << std::endl;
  for(int i = 0 ; i < n_dofs ; i++ )
  {
    myfile << grid[i] << "\t" << solution[i] << std::endl; 
  }
  std::cout << "Output written to " << filename.str() << std::endl;
  std::cout << "Plot this using:\n" 
            << "gnuplot -p < " << filename.str() << std::endl;
  myfile.close();
}
/*
 * This function calls all internal functions in the proper order
 */
void Poisson::run()
{
  make_grid();
  setup_system();
  assemble_system();
  solve();
  output_results();
}

/*
 * Main function
 * 1) Reads in the parameters from the command line call
 * 2) Makes sure the parameters are in the correct range
 * 3) Creates an instance of the Poisson class and calls run()
 */
int main (int argc, char** argv)
{
  if( argc < 4 )
  {
    std::cout << "Not enough parameters provided!" << std::endl;
    std::cout << "Usage:\n" << argv[0] << " <h> " << " <alpha> "<< " <a> " << " [<p>] " << std::endl;
    std::cout << "\t <h> is the grid size\n" 
              << "\t <alpha> is the model parameter\n"
              << "\t <a> is the factor for the right hand side f = -a\n"
	      << "\t <p> (optional) is the polynomial degree (by default 1)"<< std::endl;
    return 1;
  }
  double h = 0;
  double alpha = 0.0;
  double a = 0.0;
  int    p = 1;
  h     = std::atof(argv[1]);
  alpha = std::atof(argv[2]);
  a     = std::atof(argv[3]);
  if( argc == 5 )
  {
     p = std::atoi(argv[4]);
  }
  if ( h > 1 )
  {
    std::cout << "Grid size must be smaller than 1! " << std::endl;
    return 1;
  }
  if ( p > 2 || p < 0 )
  {
    std::cout << "Polynomial degree not in the supported range, choose either p=1 or p=2" << std::endl; 
    return 1;
  }
  
  Poisson problem(h,alpha,a,p);
  problem.run();
}
