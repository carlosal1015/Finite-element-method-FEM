#ifndef CANVAS_HH
#define CANVAS_HH

#include <vector>

#include "point.hh"
#include "pgm.hh"

// example implementation that combines declaration
// and definition
//
// benefit:
// code in one place, easy to get overview and see
// interconnection of code parts
//
// drawback:
// interface linked to implementation, no longer
// possible to provide more than one implementation

// simple canvas class for image generation
class Canvas
{

  // member variables are private
private:

  const Point _center;
  const double _width;
  const double _height;
  const int _horPixels;
  const int _vertPixels;
  std::vector<std::vector<int> > _pixels;

  // Constructor and method definitions
public:

  Canvas(const Point& center, double width, double height,
         int horPixels, int vertPixels)
    : _center(center)
    , _width(width)
    , _height(height)
    , _horPixels(horPixels)
    , _vertPixels(vertPixels)
    , _pixels(horPixels,std::vector<int>(vertPixels))
  {}

  // Accessors for member variables

  Point center() const
  {
    return _center;
  }

  double width() const
  {
    return _width;
  }

  double height() const
  {
    return _height;
  }

  int horPixels() const
  {
    return _horPixels;
  }

  int vertPixels() const
  {
    return _vertPixels;
  }

  // returns coordinates of given pixel
  Point coord(int i, int j) const
  {
    double x = (i - _horPixels/2.)*_width / _horPixels;
    double y = (j - _vertPixels/2.)*_height / _vertPixels;
    return {x + _center.x(), y + _center.y()};
  }

  // read-only access to pixel
  // to use it with a variable called mycanvas,
  // write mycanvas(i_coord,j_coord)
  int operator()(int i, int j) const
  {
    return _pixels[i][j];
  }

  // read/write access to pixel
  // to use it with a variable called mycanvas,
  // write mycanvas(i_coord,j_coord) = value;
  int& operator()(int i, int j)
  {
    return _pixels[i][j];
  }

  // create image file from canvas content
  void write(const std::string& filename) const
  {
    write_pgm(_pixels,filename);
  }
};

#endif //CANVAS_HH
