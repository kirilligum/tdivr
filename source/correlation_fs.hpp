#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <functional>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <limits>
#include <tuple>
#include <iterator>
#include <complex>
//#include <mpi.h>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/range/numeric.hpp>
#include <Eigen/Eigen>
#include <Eigen/StdVector>

#include "potential/harmonic_oscillator.hpp"

struct correlation_fs {
  double gamma = 1.0;
  double hb = 1.0;
  std::complex<double> ii = {0.0,1.0};
  std::complex<double> ir = {1.0,0.0};
  std::complex<double> operator()( 
      std::vector<double> traj1,
      std::vector<double> traj2,
      std::vector< std::vector<double>> mm1,
      std::vector< std::vector<double>> mm2,
      double s1,
      double s2,
      std::complex<double> a,
      std::complex<double> b
      ) {
    std::complex<double> corr;
    auto preb = preexp(mm1)*preexp(mm2)*exp(ii*(s1-s2)/hb)*a;
    if((harmonic_oscillator{}.q(traj1)+harmonic_oscillator{}.q(traj2))/2>0) ///> dirty approximation for slow varying 
      corr = preb*b;
    else
      corr = preb*0.0;
    return corr;
  }
  std::complex<double> preexp( const std::vector< std::vector<double>> & mm) {
    using namespace Eigen;
    Map<const MatrixXd> 
      mm0(mm[0].data(),2,2),
      mm1(mm[1].data(),2,2),
      mm2(mm[2].data(),2,2),
      mm3(mm[3].data(),2,2);
    return sqrt((0.5*(mm0*ir+mm1*ir-mm2*ii*gamma*hb+mm3*ii/gamma/hb)).determinant());
  }
};
