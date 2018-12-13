#ifndef  __ZeroPeriodicData_h
#define  __ZeroPeriodicData_h

#include  "periodicdata.h"

/*-----------------------------------------*/


namespace Gascoigne
{
class ZeroPeriodicData : public PeriodicData
{
protected:

public:

  ZeroPeriodicData() : PeriodicData() {}
  std::string GetName() const {return "Zero";}
  void operator()(DoubleVector& b, const Vertex2d& v, int col) const {}
  void operator()(DoubleVector& b, const Vertex3d& v, int col) const {}
};
}

#endif
