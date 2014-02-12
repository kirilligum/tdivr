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
#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/density.hpp>
#include <boost/accumulators/statistics/stats.hpp>
//#include <boost/math/constants/constants.hpp>
#include <ro/ro.h>

#include "potential/harmonic_oscillator.hpp"
#include "ode_step.hpp"
#include "correlation_fs.hpp"
#include "ohmic_bath_param.hpp"
#include "make_bath.hpp"
#include "make_initial_state.hpp"


int main(int argc, char const *argv[]) {
  using namespace std;
  using namespace std::placeholders;
  //using namespace boost;
  using namespace boost::accumulators;
  namespace b = boost;
  typedef accumulator_set<double, features<tag::density> > acc;
  typedef b::iterator_range<std::vector<std::pair<double, double> >::iterator > histogram_type;

  cout << " hi \n";
  size_t world_rank=1, seed = 1;
  std::default_random_engine gen(world_rank+seed);
  ohmic_bath_param  bp; std::vector<double> initial_bath,bw,bc;
  tie(initial_bath,bw,bc) = make_bath(bp,gen);
  ofstream ibo("initial_bath.txt");
  for (size_t i = 0; i < bw.size(); ++i) {
    ibo << bw[i] << setw(15) << bc[i] << setw(15) << initial_bath[2*i] << setw(15) << initial_bath[2*i+1] << endl;
  }
  acc myAccumulatorp( tag::density::num_bins = 20, tag::density::cache_size = 10);
  acc myAccumulatorq( tag::density::num_bins = 20, tag::density::cache_size = 10);
  for (size_t i = 0; i < initial_bath.size(); ++i) {
    std::ref(myAccumulatorp)(initial_bath[2*i]);
    std::ref(myAccumulatorq)(initial_bath[2*i+1]);
  }
  //std::for_each(lines.begin(), lines.end(),std::ref(myAccumulator));
  histogram_type histp = density(myAccumulatorp);
  histogram_type histq = density(myAccumulatorq);
  double totalp = 0.0;
  for( int i = 0; i < histp.size(); i++ )
  {
    std::cout << "Bin lower bound: " << histp[i].first << ", Value: " << histp[i].second << std::endl;
    totalp += histp[i].second;
  }
  std::cout << "totalp: " << totalp << std::endl; //should be 1 (and it is)
  double totalq = 0.0;
  for( int i = 0; i < histp.size(); i++ )
  {
    std::cout << "Bin lower bound: " << histp[i].first << ", Value: " << histp[i].second << std::endl;
    totalq += histp[i].second;
  }
  std::cout << "totalq: " << totalq << std::endl; //should be 1 (and it is)


  //pcet_param sp;
  //auto initial_state = make_initial_state(sp,gen);
  //initial_state.insert(end(initial_state),begin(initial_bath),end(initial_bath));
  //std::copy(std::begin(initial_state), std::end(initial_state),std::ostream_iterator<double>(cout," ")); cout << endl;
  //std::vector<double> initial_state = { 0.0, 1.0, 0.0};///> time, q , p
  //std::copy(std::begin(initial_state), std::end(initial_state),std::ostream_iterator<double>(std::cout," ")); std::cout << std::endl;
  //std::cout << harmonic_oscillator{}(initial_state) << "\n";
  //size_t seed = 1;
  //size_t time_steps = 4e1;
  //double end_time =1e1;
  //double dt = end_time/time_steps;
  //std::vector<std::vector<double>> traj(time_steps);
  //traj[0]=initial_state;
  //boost::partial_sum(traj, traj.begin(),ode_step(dt,harmonic_oscillator{}));
  //ofstream ot("traj.txt"); for(auto i:traj) { std::copy(std::begin(i), std::end(i),std::ostream_iterator<double>(ot," ")); ot << std::endl;}
  //std::vector<double> traj1(initial_state),traj2(initial_state);
  //std::vector<double> identity_matrix_asv = {1.0,0.0,0.0,1.0}; ///> Identity matrix as a vector
  //std::vector< std::vector<double>> mm1(4,identity_matrix_asv),mm2(4,identity_matrix_asv);///> had to write matrices as vector so that eigen can use them without copying and may be better for odeint
  //double s1,s2;
  //std::complex<double> a(1.0,0.0),b(1.0,0.0);
  //auto corr = correlation_fs{}(traj1,traj2,mm1,mm2,s1,s2,a,b);
  //std::cout << corr << "  <-- corr \n";
  return 0;
}
