#include <iostream>
#include <vector>
#include "hdnum.hh"

using namespace hdnum;

/** @brief Example class for a differential equation model

    The model is

    u_0'(t) = -u_1(t), t>=t_0, u_0(t_0) = U_0
    u_1'(t) =  u_0(t), t>=t_0, u_1(t_0) = U_1

    \tparam T a type representing time values
    \tparam N a type representing components of states and f-values
*/
template<class T, class N=T>
class LinearOscillatorProblem
{
public:
  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export time_type */
  typedef T time_type;

  /** \brief export number_type */
  typedef N number_type;

  //! constructor
  LinearOscillatorProblem ()
  {}

  //! return number of components for the model
  std::size_t size () const
  {
    return 2;
  }

  //! set initial state including time value
  void initialize (T& t0, hdnum::Vector<N>& u) const
  {
    t0 = 0;
    u[0] = 1.0;
    u[1] = 0.0;
  }

  //! model evaluation
  void f (const T& t, const hdnum::Vector<N>& u, hdnum::Vector<N>& result) const
  {
    result[0] = -u[1];
    result[1] =  u[0];
  }

  //! jacobian evaluation needed for implicit solvers
  void f_x (const T& t, const hdnum::Vector<N>& x, hdnum::DenseMatrix<N>& result) const
  {
    result[0][0] = 0.0; result[0][1] = -1.0;
    result[1][0] = 1.0; result[1][1] =  0.0;
  }
};

template<class Solver, class Number>
Number run_scheme (Solver solver, std::string filename, Number T)
{
  std::vector<Number> times;           // store time values here
  std::vector<Vector<Number> > states; // store states here
  times.push_back(solver.get_time());  // initial time
  states.push_back(solver.get_state()); // initial state

  int steps=0;
  while (solver.get_time()<T-1e-8) // the time loop
    {
      solver.step();                  // advance model by one time step
      times.push_back(solver.get_time()); // save time
      states.push_back(solver.get_state()); // and state
      steps++;
    }

  auto t = solver.get_time();
  auto u0 = std::cos(t);
  auto u1 = std::sin(t);
  auto u = solver.get_state();
  Number error = std::sqrt((u0-u[0])*(u0-u[0])+(u1-u[1])*(u1-u[1]));
  //std::cout << "file=" << filename << " steps=" << steps << " error=" << error << std::endl;
  //gnuplot(filename,times,states); // output model result
  return error;
}

int main (int argc, char** argv)
{
  typedef double Number;               // define a number type
  int k=40;
  Number T = 10*2*M_PI;

  // instantiate model
  typedef LinearOscillatorProblem<Number> Model; // Model type
  Model model;                         // instantiate model

  // do a refinement study
  Number errors[11][9];
  int steps[11];
  for (int i=0; i<11; i++)
    {
      steps[i] = k;
      const Number dt = T/steps[i];

      // instantiate explicit solvers
      EE<Model> solver0(model);
      Heun2<Model> solver1(model);
      Heun3<Model> solver2(model);
      RungeKutta4<Model> solver3(model);
      solver0.set_dt(dt);
      solver1.set_dt(dt);
      solver2.set_dt(dt);
      solver3.set_dt(dt);

      // instantiate implicit solvers
      Newton newton;
      newton.set_maxit(20);
      newton.set_verbosity(0);    
      newton.set_reduction(1e-12);
      newton.set_abslimit(1e-15);
      newton.set_linesearchsteps(3);  
      DIRK<Model,Newton> solver4(model,newton,"Implicit Euler");               
      DIRK<Model,Newton> solver5(model,newton,"Midpoint Rule");               
      DIRK<Model,Newton> solver6(model,newton,"Alexander");               
      DIRK<Model,Newton> solver7(model,newton,"Crouzieux");               
      solver4.set_dt(dt);
      solver5.set_dt(dt);
      solver6.set_dt(dt);
      solver7.set_dt(dt);

       // Gauß with s=3
      DenseMatrix<double> A(3,3,0.0);
      A[0][0] = 5.0/36.0;
      A[0][1] = (10.0-3.0*sqrt(15.0))/45.0;
      A[0][2] = (25.0-6.0*sqrt(15.0))/180.0;
      A[1][0] = (10.0+3.0*sqrt(15.0))/72.0;
      A[1][1] = 2.0/9.0;
      A[1][2] = (10.0-3.0*sqrt(15.0))/72.0;
      A[2][0] = (25.0+6.0*sqrt(15.0))/180.0;
      A[2][1] = (10.0+3.0*sqrt(15.0))/45.0;
      A[2][2] = 5.0/36.0;

      Vector<double> B(3,0.0);
      B[0] = 5.0/18.0;
      B[1] = 4.0/9.0;
      B[2] = 5.0/18.0;

      Vector<double> C(3, 0.0);
      C[0] = (5.0-sqrt(15.0))/10.0;
      C[1] = 0.5;
      C[2] = (5.0+sqrt(15.0))/10.0;

      RungeKutta<Model> solver8(model,A,B,C);
      solver8.set_dt(dt);

      // run the solvers
      errors[i][0] = run_scheme(solver0,"output0.dat",T);
      errors[i][1] = run_scheme(solver1,"output1.dat",T);
      errors[i][2] = run_scheme(solver2,"output2.dat",T);
      errors[i][3] = run_scheme(solver3,"output3.dat",T);
      errors[i][4] = run_scheme(solver4,"output4.dat",T);
      errors[i][5] = run_scheme(solver5,"output5.dat",T);
      errors[i][6] = run_scheme(solver6,"output6.dat",T);
      errors[i][7] = run_scheme(solver7,"output7.dat",T);
      errors[i][8] = run_scheme(solver8,"output8.dat",T);

      // double number of steps
      k = 2*k;
    }

  // report convergence rates
  for (int s=0; s<9; s++)
    {
      std::cout << "scheme " << s << std::endl;
      for (int i=1; i<11; i++)
      {
	auto rate = std::log(errors[i][s]/errors[i-1][s])/std::log(0.5);
	std::cout << "steps=" << steps[i] << " error=" << errors[i][s] << " rate=" << rate << std::endl;
      }
    }
	
  return 0;
}
