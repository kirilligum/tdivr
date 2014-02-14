#include "boltzmann_path_integral.hpp"
#include <vector>
#include <functional>
#include <complex>

std::complex<double> boltzmann_path_integral(double beta, std::function<double(std::vector<double>)> pot, int sample_size, double mass, double x_min, double x_max) { 
  return std::complex<double>(1.0,0.0);
}
