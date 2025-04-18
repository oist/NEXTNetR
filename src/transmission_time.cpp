#include <optional>
#include <memory>
#include <utility>

#include <cpp11.hpp>
#include <cpp11/function.hpp>
#include <cpp11/external_pointer.hpp>

#include "NEXTNetR_types.h"
#include "rng.h"

#include "nextnet/types.h"
#include "nextnet/random.h"

using namespace nextnetR;

using namespace cpp11;
namespace writable = cpp11::writable;

[[cpp11::register]]
doubles nextnetR_time_sample(int n, const transmission_time_R& ttr, double t, double m) {
    if (!ttr) throw std::runtime_error("time distribution cannot be NULL"); 

    RNG_SCOPE_IF_NECESSARY;
    transmission_time& tt = *(ttr.get());
    writable::doubles r;
    r.reserve(n);
    for(int i=0; i < n; ++i)
        r.push_back(tt.sample(rng_engine(), t, m));
    return r;
}

[[cpp11::register]]
doubles nextnetR_time_density(const transmission_time_R& ttr, doubles taus, doubles ts, doubles ms) {
    if (!ttr)
      throw std::runtime_error("time distribution cannot be NULL");
    if ((taus.size() != ts.size()) || (taus.size() != ms.size()))
      throw std::runtime_error("taus, ts and ms vectors must have the same length");
    
    transmission_time& tt = *(ttr.get());
    const std::size_t n = taus.size();
    writable::doubles r;
    r.reserve(n);
    for(std::size_t i=0; i < n; ++i) {
        r.push_back(taus[i] >= 0
                      ? tt.density(ts[i] + taus[i], ms[i]) /
                        tt.survivalprobability(ts[i], 0.0, ms[i])
                      : 0);
    }
    return r;
}

[[cpp11::register]]
doubles nextnetR_time_hazardrate(const transmission_time_R& ttr, doubles taus) {
    if (!ttr) throw std::runtime_error("time distribution cannot be NULL"); 

    transmission_time& tt = *(ttr.get());
    const std::size_t n = taus.size();
    writable::doubles r;
    r.reserve(n);
    for(std::size_t i=0; i < n; ++i)
        r.push_back(tt.hazardrate(taus[i]));
    return r;
}

[[cpp11::register]]
doubles nextnetR_time_survivalprobability(const transmission_time_R& ttr, doubles taus, doubles ts, doubles ms) {
    if (!ttr)
      throw std::runtime_error("time distribution cannot be NULL"); 
    if ((taus.size() != ts.size()) || (taus.size() != ms.size()))
      throw std::runtime_error("taus, ts and ms vectors must have the same length");

    transmission_time& tt = *(ttr.get());
    const std::size_t n = taus.size();
    writable::doubles r;
    r.reserve(n);
    for(std::size_t i=0; i < n; ++i) {
      const double tau = taus[i];
      const double t = ts[i];
      const double m = ms[i];
      if ((t == 0.0) && (m == 1))
        r.push_back(tt.survivalprobability(tau));
      else
        r.push_back(tt.survivalprobability(tau, t, m));
    }
    return r;
}

[[cpp11::register]]
doubles nextnetR_time_survivalquantile(const transmission_time_R& ttr, doubles ps, doubles ts, doubles ms) {
  if (!ttr)
    throw std::runtime_error("time distribution cannot be NULL"); 
  if ((ps.size() != ts.size()) || (ps.size() != ms.size()))
    throw std::runtime_error("taus, ts and ms vectors must have the same length");
  
  transmission_time& tt = *(ttr.get());
  const std::size_t n = ps.size();
  writable::doubles r;
  r.reserve(n);
  for(std::size_t i=0; i < n; ++i) {
    const double p = ps[i];
    const double t = ts[i];
    const double m = ms[i];
    if ((t == 0.0) && (m == 1))
      r.push_back(tt.survivalquantile(p));
    else
      r.push_back(tt.survivalquantile(p, t, m));
  }
  return r;
}

[[cpp11::register]]
transmission_time_R nextnetR_exponential_time(double lambda, double pinf = 0.0) {
    return new transmission_time_exponential_pinf(lambda, pinf);
}

[[cpp11::register]]
transmission_time_R nextnetR_lognormal_time(double mean, double var, double pinf = 0.0) {
    return new transmission_time_lognormal(mean, var, pinf);
}

[[cpp11::register]]
transmission_time_R nextnetR_gamma_time(double mean, double var, double pinf = 0.0) {
    return new transmission_time_gamma(mean, var, pinf);
}

[[cpp11::register]]
transmission_time_R nextnetR_weibull_time(double shape, double scale, double pinf = 0.0) {
  return new transmission_time_weibull(shape, scale, pinf);
}

[[cpp11::register]]
transmission_time_R nextnetR_polynomial_rate_time(doubles coeffs) {
  return new transmission_time_polynomial_rate(coeffs.begin(), coeffs.end());
}

[[cpp11::register]]
transmission_time_R nextnetR_deterministic_time(double tau, double pinf = 0.0) {
  return new transmission_time_deterministic(tau);
}

struct rfunction_transmission_time : public transmission_time {
    rfunction_transmission_time(double pinf)
        :transmission_time(pinf)
    {}

    virtual interval_t sample(rng_t& e, interval_t t, double m) const {
        if (sample_rf && (sample_is_trinary || ((t == 0.0) && (m == 1.0)))) {
            // Since the sample() R code will likely use the RNG relinquish before calling
            R_rng_relinquish rngrel;
            return as_cpp<double>(sample_is_trinary ? (*sample_rf)(t, m) : (*sample_rf)());
        }
        else {
            return transmission_time::sample(e, t, m);
        }
    }
    
    virtual double density(interval_t tau) const {
        return as_cpp<double>(density_rf(tau));
    }
    
    virtual double survivalprobability(interval_t tau) const {
        if (survival_is_trinary)
            return as_cpp<double>(survival_rf(tau, 0.0, 1.0));
        else
            return as_cpp<double>(survival_rf(tau));
    }
    
    virtual double survivalprobability(interval_t tau, interval_t t, double m) const {
        if (survival_is_trinary)
            return as_cpp<double>(survival_rf(tau, t, m));
        else
            return transmission_time::survivalprobability(tau, t, m);
    }
    
    virtual double survivalquantile(double u) const {
        if (survivalquantile_rf && quantile_is_trinary)
            return as_cpp<double>((*survivalquantile_rf)(u, 0.0, 1.0));
        else if (survivalquantile_rf && !quantile_is_trinary)
            return as_cpp<double>((*survivalquantile_rf)(u));
        else
            return transmission_time::survivalquantile(u);
    }
    
    virtual double survivalquantile(double u, interval_t t, double m) const {
        if (survivalquantile_rf && quantile_is_trinary)
            return as_cpp<double>((*survivalquantile_rf)(u, t, m));
        else
            return transmission_time::survivalquantile(u, t, m);
    }

    std::optional<function> sample_rf;
    bool sample_is_trinary = false;
    function density_rf = R_NilValue;
    function survival_rf = R_NilValue;
    bool survival_is_trinary = false;
    std::optional<function> survivalquantile_rf;
    bool quantile_is_trinary = false;
};

[[cpp11::register]]
transmission_time_R nextnetR_userdefined_time(
    SEXP sample, bool sample_is_trinary,
    SEXP survival, bool survival_is_trinary,
    SEXP density,
    SEXP survivalquantile, bool quantile_is_trinary,
    double pinfinity = 0.0)
{
    std::unique_ptr<rfunction_transmission_time> r(new rfunction_transmission_time(pinfinity));
    
    if (sample != R_NilValue)
      r->sample_rf = sample;
    r->sample_is_trinary = sample_is_trinary;
    r->survival_rf = survival;
    r->survival_is_trinary = survival_is_trinary;
    r->density_rf = density;
    if (survivalquantile != R_NilValue)
        r->survivalquantile_rf = survivalquantile;
    r->quantile_is_trinary = quantile_is_trinary;
        
    return r.release();
}
