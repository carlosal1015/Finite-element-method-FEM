#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <iterator>

std::vector<double> reversed(const std::vector<double>& v)
{
  std::vector<double> r;
  std::copy(v.rbegin(),v.rend(),std::back_inserter(r));
  return r;
}


int main()
{
  auto v = std::vector<double>{ 1,2,3,4 };
  auto vr = reversed(v);
  
  for(auto entry : vr)
    std::cout << entry << std::endl;

  std::map<std::string,int> kurse;

  kurse["ipk"] = 170;
  kurse["ipi"] = 150;

  for (std::pair<const std::string,int> kurs : kurse)
    std::cout << kurs.first << ": " << kurs.second << std::endl;

  for (auto& kurs : kurse)
    kurs.second -= 10;

  for (std::pair<const std::string,int> e : kurse)
    std::cout << e.first << ": " << e.second << std::endl;

  return 0;
}
