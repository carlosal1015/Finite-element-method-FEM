#include <iostream>
#include <vector>
#include "hdnum.hh"

using namespace hdnum;

template<class T, class N=T>
class VanDerPolProblem
{
public:
  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export time_type */
  typedef T time_type;

  /** \brief export number_type */
  typedef N number_type;

  //! constructor stores parameter lambda
  VanDerPolProblem () : eps(1.0E-3), ctr(0)
  {
  }

  VanDerPolProblem (N eps_) : eps(eps_), ctr(0)
  {
  }

  //! return number of componentes for the model
  std::size_t size () const
  {
    return 2;
  }

  //! set initial state including time value
  void initialize (T& t0, Vector<N>& x0) const
  {
    t0 = 0;
    x0[0] = 1.0;
    x0[1] = 2.0;
  }

  //! model evaluation
  void f (const T& t, const Vector<N>& x, Vector<N>& result) const
  {
    result[0] = -x[1];
    result[1] = (x[0]-x[1]*x[1]*x[1]/N(3.0)+x[1])/eps;
    ctr++;
  }

  //! model evaluation
  void f_x (const T& t, const Vector<N>& x, DenseMatrix<N>& result) const
  {
    result[0][0] = 0.0;      result[0][1] = -1.0;
    result[1][0] = 1.0/eps;  result[1][1] = (1.0-x[1]*x[1])/eps;
  }

  size_type get_count () const
  {
    return ctr;
  }

private:
  N eps;
  mutable size_type ctr;
};



int main ()
{
  typedef double Number;               // define a number type

  typedef VanDerPolProblem<Number> Model; // Model type
  Model model(1e-3);                         // instantiate model

  Newton newton;
  newton.set_maxit(20);
  newton.set_verbosity(0);    
  newton.set_reduction(1e-12);
  newton.set_abslimit(1e-15);
  newton.set_linesearchsteps(3);  

  /* Runge-Kutta Fehlberg solver */
  // typedef RKF45<Model> Solver;         
  // Solver solver(model);                
  // solver.set_TOL(0.2); 

  /* DIRK solver 

     allows for  "Alexander"
                 "Midpoint Rule"
                 "Crouzieux"
                 "Implicit Euler"
                 "Fractional Step Theta"
  */
  typedef DIRK<Model,Newton> Solver;         
  Solver solver(model,newton,"Implicit Euler");               
  //Solver solver(model,newton,"Alexander");               

  const Number dt = 1.0/16.0;
  solver.set_dt(dt);             // set initial time step
  //solver.set_verbosity(1);

  std::vector<Number> times;           // store time values here
  std::vector<Vector<Number> > states; // store states here
  std::vector<Number> dts;             // store delta t
  times.push_back(solver.get_time());  // initial time
  states.push_back(solver.get_state()); // initial state
  dts.push_back(solver.get_dt());      // initial dt

  Number T = 10;
  int steps = 0;
  Number min_dt = dt;
  while (solver.get_time()<T-1e-8) // the time loop
    {
      solver.step();                  // advance model by one time step
      times.push_back(solver.get_time()); // save time
      states.push_back(solver.get_state()); // and state
      dts.push_back(solver.get_dt());      // used dt
      min_dt = std::min(dts.back(),min_dt);
      steps++;
    }

  std::cout << "Solution at time " << solver.get_time() 
            << " : " << solver.get_state() << std::endl;

  solver.get_info();

  std::cout << " number of time steps: " << steps
	    << " number of f evaluations: " << model.get_count() 
            << " minimum step width: " << min_dt
	    << std::endl;
  //gnuplot("vanderpol_rk45_loose.dat",times,states,dts); // output model result
  gnuplot("vanderpol_IE_loose.dat",times,states,dts); // output model result

  return 0;
}
