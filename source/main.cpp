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
#include <ro/ro.h>
//#include <boost/math/constants/constants.hpp>

#include "potential/harmonic_oscillator.hpp"
#include "ode_step.hpp"
#include "correlation_fs.hpp"
#include "ohmic_bath_param.hpp"
#include "make_bath.hpp"
#include "make_initial_state.hpp"


int main(int argc, char const *argv[]) {
  using namespace std;
  namespace b = boost;
  cout << " hi \n";
  size_t world_rank=1, seed = 1;
  std::default_random_engine gen(world_rank+seed);
  ohmic_bath_param  bp; std::vector<double> initial_bath,bw,bc;
  tie(initial_bath,bw,bc) = make_bath(bp,gen);
  //copy(bw,ostream_iterator<double>(cout," "));
  ofstream ibo("initial_bath.txt");
  //for (size_t i = 0; i < bw.size(); ++i) {
    //ibo << bw[i] << setw(15) << bc[i] << setw(15) << initial_bath[2*i] << setw(15) << initial_bath[2*i+1] << endl;
  //}
  std::for_each(
      b::make_zip_iterator(b::make_tuple(
          std::begin(bw),std::begin(bc)
          )),
      b::make_zip_iterator(b::make_tuple(
          std::end(bw),std::end(bc)
          )),
        [&ibo](b::tuple<const double &,const double &> t){
        ibo << t.get<0>() << "  " << t.get<1>() << "\n";}
      );

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
