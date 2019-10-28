#ifndef IO_HH
#define IO_HH

#include <vector>

// reads numbers from stdin, one number per line until it detects an empty line
// returns all read numbers in a vector
std::vector<double> read_vector();

// generates N uniformly distributed numbers between min and max.
// seed is an arbitrary number that initializes the random number generator. Using the
// same number always creates the same output. If you want to have different results each
// time, use random_seed() to generate seed.
std::vector<double> uniform_distribution(int seed, int N, double min, double max);

// generates N normally distributed numbers with given average and standard deviation.
// seed is an arbitrary number that initializes the random number generator. Using the
// same number always creates the same output. If you want to have different results each
// time, use random_seed() to generate seed.
std::vector<double> normal_distribution(int seed, int N, double avg, double std_dev);

// returns a random seed (different on each call)
int random_seed();

#endif // IO_HH
