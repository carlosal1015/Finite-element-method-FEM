#ifndef POINT_HH
#define POINT_HH

// Interface und Implementierung sind hier der Einfachheit
// halber in einer Datei

class Point
{

public:

  // default constructor
  Point();

  // constructor für übergebene Koordinaten
  Point(double x, double y);

  // Accessors

  double x() const;

  double y() const;

private:

  double _x;
  double _y;

};

#endif // POINT_HH
