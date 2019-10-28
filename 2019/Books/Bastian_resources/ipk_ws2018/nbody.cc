#include<iostream>
#include<vector>
#include<map>
#include<chrono>

// for SDL display
#include "sdlwrapper.hh"

// for random positions
#include"io.hh"

/* data structures */

// structure for 2D points
struct Point
{
  // coordinates
  double x;
  double y;

  // addition operators
  Point& operator += (const Point& other)
  {
    x += other.x;
    y += other.y;
    return *this;
  }

  Point operator + (const Point& other) const
  {
    Point output = *this;
    output += other;
    return output;
  }

  // subtraction operators
  Point& operator -= (const Point& other)
  {
    x -= other.x;
    y -= other.y;
    return *this;
  }

  Point operator - (const Point& other) const
  {
    Point output = *this;
    output -= other;
    return output;
  }

  // scaling operators
  Point& operator *= (double alpha)
  {
    x *= alpha;
    y *= alpha;
    return *this;
  }

  Point operator * (double alpha) const
  {
    Point output = *this;
    output *= alpha;
    return output;
  }
};

// scaling operator with scalar in front
Point operator * (double alpha, const Point& point)
{
  Point output = point;
  output *= alpha;
  return output;
}

// structure representing celestial body
struct Body
{
  Point position;
  Point velocity;
  double mass;
  std::array<int,3> color;
};

/* display and transform functions */

// draw current positions on canvas
  template<typename Body>
void displayBodies(SDLCanvas& canvas, const std::vector<Body>& bodies)
{
  for (unsigned int i = 0; i < bodies.size(); i++)
    canvas.drawPixel(512 + bodies[i].position.x, 384 - bodies[i].position.y,bodies[i].color);

  canvas.display();
}

// switch to center of mass system
  template<typename Body>
void centerMass(std::vector<Body>& bodies)
{
  Point center    = {0.,0.};
  Point centerVel = {0.,0.};
  double totalMass = 0.;

  // calculate weighted center and total mass
  for (unsigned int i = 0; i < bodies.size(); i++)
  {
    center    += bodies[i].mass * bodies[i].position;
    centerVel += bodies[i].mass * bodies[i].velocity;
    totalMass += bodies[i].mass;
  }

  center    *= 1./totalMass;
  centerVel *= 1./totalMass;

  // correct position and velocity
  for (unsigned int i = 0; i < bodies.size(); i++)
  {
    bodies[i].position -= center;
    bodies[i].velocity -= centerVel;
  }
}

// merge bodies that are sufficiently close
template<typename Body>
void mergeMass(std::vector<Body>& bodies)
{
  //check until nothing is merged
  bool dropped = true;
  while(dropped)
  {
    dropped = false;

    // loop over all pairs
    for (unsigned int i = 0; i < bodies.size(); i++)
      for (unsigned int j = i+1; j < bodies.size(); j++)
        // merge if distance is small enough
        if (std::pow(bodies[j].position.x-bodies[i].position.x,2)
                  + std::pow(bodies[j].position.y-bodies[i].position.y,2) < 10.)
        {
          // weighted mean of position, velocity, mass
          bodies[i].position = 1./(bodies[i].mass + bodies[j].mass)
            * (bodies[i].mass * bodies[i].position
                + bodies[j].mass * bodies[j].position);

          bodies[i].velocity = 1./(bodies[i].mass + bodies[j].mass)
            * (bodies[i].mass * bodies[i].velocity
                + bodies[j].mass * bodies[j].velocity);

          bodies[i].mass += bodies[j].mass;

          // mean of colors
          for (unsigned int k = 0; k < 3; k++)
            bodies[i].color[k] = 1./2. * (bodies[i].color[k] + bodies[j].color[k]);

          // swap to end and then drop last entry
          std::swap(bodies[j],bodies[bodies.size()-1]);
          bodies.pop_back();
          dropped = true;
      }
  }
}

/* time stepping functions */

// explicit Euler method
// u(t + dt) = u(t) + dt * u'(t)
  template<typename Force, typename Body>
void explEulerStep (const Force& force, std::vector<Body>& bodies,
    double t, double dt)
{
  // calculate acceleration
  std::vector<Point> accel;
  for (unsigned int i = 0; i < bodies.size(); i++)
    accel.push_back(force(bodies, i, t));

  // update position and velocity
  for (unsigned int i = 0; i < bodies.size(); i++)
  {
    bodies[i].position += dt * bodies[i].velocity;
    bodies[i].velocity += dt * accel[i];
  }
}

// implicit Euler method
// u(t + dt) = u(t) + dt * u'(t + dt)
  template<typename Force, typename Body>
void implEulerStep (const Force& force, std::vector<Body>& bodies,
    double t, double dt)
{
  std::vector<Body> oldBodies = bodies;
  std::vector<Body> prevBodies;
  double residual = 1.;
  int count = 0;
  // repeat while residual is to large
  while (residual > 1e-10 && count < 20)
  {
    prevBodies = bodies;
    bodies     = oldBodies;

    // calculate approximate acceleration at t+dt
    std::vector<Point> accel;
    for (unsigned int i = 0; i < bodies.size(); i++)
      accel.push_back(force(prevBodies, i, t+dt));

    // update position and velocity
    for (unsigned int i = 0; i < bodies.size(); i++)
    {
      bodies[i].position += dt * prevBodies[i].velocity;
      bodies[i].velocity += dt * accel[i];
    }

    // calculate residual (== mismatch)
    residual = 0.;
    for (unsigned int i = 0; i < bodies.size(); i++)
    {

      residual += std::pow(bodies[i].position.x - prevBodies[i].position.x,2)
        + std::pow(bodies[i].position.y - prevBodies[i].position.y,2)
        + std::pow(bodies[i].velocity.x - prevBodies[i].velocity.x,2)
        + std::pow(bodies[i].velocity.y - prevBodies[i].velocity.y,2);
    }
    residual = std::sqrt(residual);
    count++;
  }

}

// leapfrog method
// update position, then velocity, then position...
  template<typename Force, typename Body>
void leapfrogStep (const Force& force, std::vector<Body>& bodies,
    double t, double dt)
{
  // calculate acceleration
  std::vector<Point> accel;
  for (unsigned int i = 0; i < bodies.size(); i++)
    accel.push_back(force(bodies, i, t));

  // update half of velocity, position, other half of velocity
  for (unsigned int i = 0; i < bodies.size(); i++)
  {
    bodies[i].velocity += accel[i] * (dt/2.);
    bodies[i].position += bodies[i].velocity * dt;
    bodies[i].position += accel[i] * (dt*dt/2.);
    bodies[i].velocity += accel[i] * (dt/2.);
  }
}

// Crank-Nicolson method
// u(t + dt) = u(t) + 1/2 * [ dt * u'(t) + dt * u'(t + dt) ]
  template<typename Force, typename Body>
void crankNicolsonStep (const Force& force, std::vector<Body>& bodies,
    double t, double dt)
{
  std::vector<Body> oldBodies = bodies;
  std::vector<Body> prevBodies;

  // calculate explicit Euler acceleration
  std::vector<Point> oldAccel;
  for (unsigned int i = 0; i < bodies.size(); i++)
    oldAccel.push_back(force(oldBodies,i,t));

  // repeat implicit Euler part until converged
  double residual = 1.;
  int count = 0;
  while (residual > 1e-10 && count < 20)
  {
    prevBodies = bodies;
    bodies     = oldBodies;

    std::vector<Point> prevAccel;
    for (unsigned int i = 0; i < bodies.size(); i++)
      prevAccel.push_back(force(prevBodies,i,t+dt));

    for (unsigned int i = 0; i < bodies.size(); i++)
    {
      bodies[i].position += 1./2. * (oldBodies[i].velocity + prevBodies[i].velocity) * dt;
      bodies[i].velocity += 1./2. * (oldAccel[i] + prevAccel[i]) * dt;
    }

    residual = 0.;
    for (unsigned int i = 0; i < bodies.size(); i++)
    {

      residual += std::pow(bodies[i].position.x - prevBodies[i].position.x,2)
        + std::pow(bodies[i].position.y - prevBodies[i].position.y,2)
        + std::pow(bodies[i].velocity.x - prevBodies[i].velocity.x,2)
        + std::pow(bodies[i].velocity.y - prevBodies[i].velocity.y,2);
    }
    residual = std::sqrt(residual);
    count++;
  }

}

// Runge-Kutta (RK4) method
template<typename Force, typename Body>
void rk4Step (const Force& force, std::vector<Body>& bodies,
    double t, double dt)
{
  std::vector<std::vector<Body>> k(4,bodies);

  // perform explicit Euler steps of various lengths
  for (unsigned int stage = 0; stage < 4; stage++)
  {
    double dtStage = dt;
    if (stage == 1 or stage == 2)
    {
      dtStage *= 1./2.;
      for (unsigned int i = 0; i < bodies.size(); i++)
      {
        k[stage][i].position += 1./2. * k[stage-1][i].position;
        k[stage][i].velocity += 1./2. * k[stage-1][i].velocity;
      }
    }
    else if (stage == 3)
    {
      for (unsigned int i = 0; i < bodies.size(); i++)
      {
        k[stage][i].position += k[stage-1][i].position;
        k[stage][i].velocity += k[stage-1][i].velocity;
      }
    }

    explEulerStep(force,k[stage],t,dtStage);

    for (unsigned int i = 0; i < bodies.size(); i++)
    {
      k[stage][i].position -= bodies[i].position;
      k[stage][i].velocity -= bodies[i].velocity;
    }
  }

  // update position and velocity based on stages
  for (unsigned int i = 0; i < bodies.size(); i++)
  {
    bodies[i].position += 1./6. * (k[0][i].position + 2. * k[1][i].position
        + 2. * k[2][i].position + k[3][i].position);
    bodies[i].velocity += 1./6. * (k[0][i].velocity + 2. * k[1][i].velocity
        + 2. * k[2][i].velocity + k[3][i].velocity);
  }
}


/* simulation and main function */

// choose method and simulate system
template<typename Force>
void simulateNBody(std::vector<Body>& bodies, Force force)
{
  // choose numerical method
  std::cout << "Choose scheme:\n"
    << "\t1) explicit Euler\n"
    << "\t2) implicit Euler\n"
    << "\t3) leapfrog\n"
    << "\t4) Crank-Nicolson\n"
    << "\t5) Runge-Kutta (RK4)\n"
    << std::endl;

  int scheme;
  std::cin >> scheme;

  // choose postprocessing
  std::cout << "Choose post-processing:\n"
    << "\t1) none\n"
    << "\t2) center of mass\n"
    << "\t3) merge close points\n"
    << "\t4) both\n"
    << std::endl;

  int processing;
  std::cin >> processing;

  // create canvas
  SDLCanvas canvas("screen", 1024, 768);

  // save current time
  const auto start = std::chrono::steady_clock::now();
  unsigned int count = 0;
  double t = 0.;
  double dt = 1.;
  // repeat until stop requested
  while (true)
  {
    // perform step
    switch (scheme)
    {
      case 1: explEulerStep    (force,bodies,t,dt); break;
      case 2: implEulerStep    (force,bodies,t,dt); break;
      case 3: leapfrogStep     (force,bodies,t,dt); break;
      case 4: crankNicolsonStep(force,bodies,t,dt); break;
      case 5: rk4Step          (force,bodies,t,dt); break;

      default: std::cout << "input error" << std::endl; return;
    }
    // perform post-processing
    switch (processing)
    {
      case 1: /* (nothing) */                         break;
      case 2: centerMass(bodies);                     break;
      case 3: mergeMass (bodies);                     break;
      case 4: mergeMass (bodies); centerMass(bodies); break;

      default: std::cout << "input error" << std::endl; return;
    }

    // update time
    t += dt;

    // show results
    displayBodies(canvas,bodies);

    if (++count % 1000 == 0)
    {
      // get current time, print frames per second
      const auto now = std::chrono::steady_clock::now();
      const auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(now-start).count();
      std::cout << "Elapsed time: " << diff/1000. << " s, fps: " << count*1000./diff << std::endl;

      // calculate angular momentum to check its conservation
      double angularMomentum = 0.;
      for (const Body& body : bodies)
        angularMomentum += body.mass
          * (body.position.x * body.velocity.y
              - body.position.y * body.velocity.x);
      std::cout << "Angular momentum: " << angularMomentum << std::endl;
    }

    // stop if window was closed
    if (canvas.windowClosed())
      break;
  }
}

// choose bodies and force, then simulate
int main()
{
  // define different choices of bodies
  std::vector<Body> solarBodies = {
    {{   0.,  0.}, {   0., 0.}, 1e3, {255,  0,  0}},
    {{ 100.,  0.}, {   0.,0.3}, 10., {  0,255,  0}},
    {{-200.,  0.}, {   0.,0.2}, 10., {  0,  0,255}},
    {{   0.,250.}, {-0.25, 0.}, 10., {255,255,  0}},
  };
  std::vector<Body> symmetricBodies = {
    {{  150.,    0.}, {  0.,-0.2}, 1e3, {255,  0,  0}},
    {{ -150.,    0.}, {  0., 0.2}, 1e3, {  0,255,  0}},
    {{    0.,  150.}, { 0.2,  0.}, 1e3, {  0,  0,255}},
    {{    0., -150.}, {-0.2,  0.}, 1e3, {255,255,  0}},
  };
  std::vector<Body> randomPosBodies;
  for (unsigned int i = 0; i < 20; i++)
   {
     std::vector<double> vals = normal_distribution (random_seed(),3,0.,1.);
     std::vector<double> cols = uniform_distribution(random_seed(),3,0.,256.);
     randomPosBodies.push_back({
         {100.*vals[0],100.*vals[1]},
         {0.,0.},
         100.+10.*vals[3],
         {int(cols[0]),int(cols[1]),int(cols[2])}
         });
   }
  std::vector<Body> randomVelBodies;
  for (unsigned int i = 0; i < 20; i++)
   {
     std::vector<double> vals = normal_distribution (random_seed(),5,0.,1.);
     std::vector<double> cols = uniform_distribution(random_seed(),3,0.,256.);
     randomVelBodies.push_back({
         {100.*vals[0],100.*vals[1]},
         {0.1*vals[2],0.1*vals[3]},
         100.+10.*vals[4],
         {int(cols[0]),int(cols[1]),int(cols[2])}
         });
   }

  // main loop
  bool finished = false;
  while (not finished)
  {
    // choose set of bodies for simulation
    std::vector<Body> bodies;
    std::cout << "Choose setup: \n"
      << "\t1) solar system\n"
      << "\t2) symmetric\n"
      << "\t3) random positions\n"
      << "\t4) random pos. + vel.\n"
      << "\t5) (quit)\n"
      << std::endl;

    int setup;
    std::cin >> setup;
    switch (setup)
    {
      case 1: bodies = solarBodies;     break;
      case 2: bodies = symmetricBodies; break;
      case 3: bodies = randomPosBodies; break;
      case 4: bodies = randomVelBodies; break;
      case 5: finished = true;          break;

      default: std::cout << "input error" << std::endl; continue;
    }

    // stop if "quit" chosen
    if (finished)
      break;


    // forces

    // gravity
    auto gravity = [](auto& bodies, int i, double t)
      {
        Point output = {0.,0.};

        for (unsigned int j = 0; j < bodies.size(); j++)
          if (i != j)
          {
            const double distFactor = 0.01 * bodies[j].mass
              / std::pow(std::pow(bodies[j].position.x-bodies[i].position.x,2)
                  + std::pow(bodies[j].position.y-bodies[i].position.y,2),1.5);

            output += distFactor * (bodies[j].position - bodies[i].position);
          }

        return output;
      };

    // elastic force (springs)
    auto elastic = [](auto& bodies, int i, double t)
      {
        Point output = {0.,0.};

        for (unsigned int j = 0; j < bodies.size(); j++)
          if (i != j)
          {
            const double distFactor = 0.0001 / bodies[i].mass;

            output += distFactor * (bodies[j].position - bodies[i].position);
          }

        return output;
      };

    // magnetic (Lorenz) force
    auto magnetic = [](auto& bodies, int i, double t)
      {
        Point output = {0.,0.};

        ///*
        for (unsigned int j = 0; j < bodies.size(); j++)
          if (i != j)
          {
            const double distFactor = - 0.01 * bodies[j].mass
              / std::pow(std::pow(bodies[j].position.x-bodies[i].position.x,2)
                  + std::pow(bodies[j].position.y-bodies[i].position.y,2),1.5);;

            output += distFactor * (bodies[j].position - bodies[i].position);
          }
          //*/

        output.x += - 0.005 * bodies[i].velocity.y;
        output.y +=   0.005 * bodies[i].velocity.x;

        return output;
      };

    // class for force resulting from Rosenbrock potential
    auto potential = [](double x, double y)
      {
        return 0.1 * (std::pow(1.-x/100.,2) + 100. * std::pow(y/100. - std::pow(x/100.,2),2));
      };

    auto rosenbrock_force = [potential](auto& bodies, int i, double t)
      {
        const double eps = 1e-6;
        const Point p = bodies[i].position;
        const double dVx = (potential(p.x+eps,p.y) - potential(p.x-eps,p.y))/(2.*eps);
        const double dVy = (potential(p.x,p.y+eps) - potential(p.x,p.y-eps))/(2.*eps);
        const Point output = {-dVx/bodies[i].mass,-dVy/bodies[i].mass};

        return output;
      };

    // choose acting force and run simulation
    std::cout << "Choose force: \n"
      << "\t1) gravity\n"
      << "\t2) elastic\n"
      << "\t3) electromagnetic\n"
      << "\t4) gravity + elastic\n"
      << "\t5) force field\n"
      << std::endl;

    int force;
    std::cin >> force;
    switch (force)
    {
      case 1: simulateNBody(bodies,gravity);                 break;
      case 2: simulateNBody(bodies,elastic);                 break;
      case 3: simulateNBody(bodies,magnetic);                break;
      case 4: simulateNBody(
        bodies,
        [=](auto& bodies, int i, double t)
        {
          return gravity(bodies,i,t) + elastic(bodies,i,t);
        });
                                                             break;
      case 5: simulateNBody(bodies,rosenbrock_force);        break;

      default: std::cout << "input error" << std::endl; continue;
    }
  }

  return 0;
}
