#include <fstream>
#include<cmath>

#include "pgm.hh"

std::vector<std::vector<int> > read_pgm(const std::string& filename)
{
  std::ifstream file(filename);
  
  std::string token;
  file >> token;
  if (token != "P2")
    throw std::runtime_error("wrong file format");
  
  int xSize, ySize, maxVal;
  file >> xSize >> ySize >> maxVal;

  std::vector<std::vector<int> > data(xSize,std::vector<int>(ySize));

  for (unsigned int j = 0; j < data[0].size(); j++)
    for (unsigned int i = 0; i < data.size(); i++)
      file >> data[i][data[0].size()-1-j];

  return data;
}

void write_pgm(const std::vector<std::vector<int> >& data, const std::string& filename)
{
  // check that data is nonempty
  if (data.empty() || data[0].empty())
    throw std::invalid_argument("empty pixel array");

  // check that data array is rectangular
  for (const std::vector<int>& row : data)
    if (row.size() != data[0].size())
      throw std::invalid_argument("row sizes inconsistent");


  // write file
  std::ofstream file(filename);

  file << "P2\n"
    << data.size() << " " << data[0].size() << "\n";
  
  // renormalize if necessary
  int maxVal = 0;
  for (const std::vector<int>& row : data)
    for (const int& entry : row)
      maxVal = std::max(entry,maxVal);
  
  // renormalizing output
  if (maxVal > 65535)
  {
    file << 65535 << "\n";
    for (unsigned int j = 0; j < data[0].size(); j++)
    {
      for (unsigned int i = 0; i < data.size(); i++)
      {
        if (data[i][j] < 0)
          throw std::range_error("negative entry at indices ("
              + std::to_string(i) + ", " + std::to_string(j) + ")");
      
        file << std::floor(data[i][j] * (65535./maxVal)) << " ";
      }
      file << "\n";
    }
  }
  // normal output
  else
  {
    file << maxVal << "\n";
    for (unsigned int j = 0; j < data[0].size(); j++)
    {
      for (unsigned int i = 0; i < data.size(); i++)
      {
        if (data[i][j] < 0)
          throw std::range_error("negative entry at indices ("
              + std::to_string(i) + ", " + std::to_string(j) + ")");
      
        file << data[i][data[0].size()-1-j] << " ";
      }
      file << "\n";
    }
  }
}
