#include "histogram_bath_pq.hpp"

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

#include <vector>
#include <string>
void histogram_bath_pq(const std::vector<double> &bath_pq, int nbin, int cache_size) {
  using namespace std;
  using namespace std::placeholders;
  //using namespace boost;
  using namespace boost::accumulators;
  namespace b = boost;
  typedef accumulator_set<double, features<tag::density> > acc;
  typedef b::iterator_range<std::vector<std::pair<double, double> >::iterator > histogram_type;

  acc myAccumulatorp( tag::density::num_bins = 20, tag::density::cache_size = 10);
  acc myAccumulatorq( tag::density::num_bins = 20, tag::density::cache_size = 10);
  for (size_t i = 0; i < bath_pq.size(); ++i) {
    std::ref(myAccumulatorp)(bath_pq[2*i]);
    std::ref(myAccumulatorq)(bath_pq[2*i+1]);
  }
  //std::for_each(lines.begin(), lines.end(),std::ref(myAccumulator));
  histogram_type histp = density(myAccumulatorp);
  histogram_type histq = density(myAccumulatorq);
  ofstream ohbp("bath_p_histogram.txt");
  ofstream ohbq("bath_q_histogram.txt");
  for( int i = 0; i < histp.size(); i++ ) {
    ohbp << histp[i].first << "  " << histp[i].second << std::endl;
  }
  for( int i = 0; i < histq.size(); i++ ) {
    ohbq << histq[i].first << "  " << histq[i].second << std::endl;
  }
}
