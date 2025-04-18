/* mixture_time.cpp */

#include <random>
#include <numeric>
#include <cpp11.hpp>
#include <cpp11/function.hpp>

#include "NEXTNetR/NEXTNetR_types.h"
#include "nextnet/random.h"

[[cpp11::linking_to("BH")]]
[[cpp11::linking_to("NEXTNetR")]]

using namespace cpp11;

struct mixture_time_impl : public virtual transmission_time {
  virtual interval_t sample(rng_t &rng, interval_t t, double m) const override {
    if ((t == 0.0) && (m == 1.0))
      return times.at(pick(rng)).get()->sample(rng, 0.0, 1.0);
    return transmission_time::sample(rng, t, m);
  }
  
  virtual double survivalprobability(interval_t tau) const override {
    double p = 0.0;
    for(std::size_t i=0; i < times.size(); ++i)
      p += weights.at(i) * times.at(i).get()->survivalprobability(tau);
    return p;
  }

  virtual double density(interval_t tau) const override {
    double p = 0.0;
    for(std::size_t i=0; i < times.size(); ++i)
      p += weights.at(i) * times.at(i).get()->density(tau);
    return p;
  }
  std::vector<transmission_time_R> times;
  std::vector<double> weights;
  mutable std::discrete_distribution<std::size_t> pick;
};

[[cpp11::register]]
SEXP mixture_time(list times, doubles weights) {
  // Validate parameters
  if (times.size() != weights.size())
    throw std::runtime_error("number of distributions and weights must agree");

  // Create object
  auto r = std::make_unique<mixture_time_impl>();
  
  // Set individual distributions and weights
  r->times.reserve(times.size());
  r->weights.reserve(times.size());
  const double ws = std::reduce(weights.begin(), weights.end());
  for(R_xlen_t i=0; i < times.size(); ++i) {
    r->times.push_back((transmission_time_R)times[i]);
    r->weights.push_back(weights[i] / ws);
  }
  r->pick = std::discrete_distribution<std::size_t>(weights.begin(), weights.end());

  // Return created object
  return transmission_time_R(r.release());
}

