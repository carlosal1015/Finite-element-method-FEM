#include <iostream>
#include <array>
#include <list>
#include <vector>
#include <algorithm>
#include <iterator>

#include "io.hh"
#include "statistics.hh"


template<typename T>
void print_statistics(const T& t)
{
  std::cout << "mean: " << mean(t) << std::endl;
  std::cout << "median: " << median(t) << std::endl;
  std::cout << "2nd moment: " << moment(t,2) << std::endl;
  std::cout << "std_dev: " << standard_deviation(t) << std::endl;
}


int main()
{
  auto vec = uniform_distribution(0,10000,0.0,1.0);
  print_statistics(vec);

  std::array<int,6> int_array = {1,2,3,4,5,6};
  print_statistics(int_array);

  // only enable the following block for task c)
  /*
  std::list<float> float_list;
  std::copy(vec.begin(),vec.end(),std::back_inserter(float_list));
  print_statistics(float_list);
  */

  return 0;
}
