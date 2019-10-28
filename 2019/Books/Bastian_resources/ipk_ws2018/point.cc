#include "point.hh"

Point::Point()
  : _x(0.0)
  , _y(0.0)
{}

Point::Point(double x, double y)
  : _x(x)
  , _y(y)
{}

double Point::x() const
{
  return _x;
}

double Point::y() const
{
  return _y;
}
