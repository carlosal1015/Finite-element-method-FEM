#include <iostream>
#include <vector>
#include "hdnum.hh"
#include <cmath>

using namespace hdnum;

int main ()
{
  typedef double Number;    // define a type for numbers
  typedef Vector<Number> V; // define a type for states u

  Number lambda = -1.0; // model problem u' = lambda * u
  Number dt = 0.1;      // time step \Delta t
  Number time = 0.0;    // initial time t_0
  V u(2);               // a state with one component (d=1)
  u[0] = 1.0;           // assign initial value

  std::vector<Number> times; // store all time values here
  std::vector<V> states;     // store all states here
  times.push_back(time);     // initial time
  states.push_back(u);       // initial state
  std::cout << "Enter the value of Delta t: ";
  std::cin >> dt;
  while (time<4.0-1e-6) // the time loop
    {
      u[0] = u[0] + dt*u[1];
      u[1] = u[1] + dt*lambda*u[0]*std::sin(u[0]); // compute a new state
      time += dt;                   // advance time
      times.push_back(time);        // save new time
      states.push_back(u);          // and new state
    }

  gnuplot("eemp.dat",times, states); // output model result

  return 0;
}
