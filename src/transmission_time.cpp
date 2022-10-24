#include <optional>

#include <cpp11.hpp>
#include <cpp11/function.hpp>
#include <cpp11/external_pointer.hpp>

#include "episimR_types.h"
#include "rng.h"

#include "epidemics/types.h"
#include "epidemics/random.h"

using namespace episimR;

using namespace cpp11;
namespace writable = cpp11::writable;

[[cpp11::register]]
doubles episimR_time_sample(int n, const transmission_time_R& ttr, interval_t t, int m) {
    RNG_SCOPE;
    transmission_time& tt = *(ttr.get());
    writable::doubles r;
    r.reserve(n);
    for(int i=0; i < n; ++i)
        r.push_back(tt.sample(rng_engine(), t, m));
    return r;
}

[[cpp11::register]]
doubles episimR_time_density(const transmission_time_R& ttr, doubles a1) {
    transmission_time& tt = *(ttr.get());
    const std::size_t n = a1.size();
    writable::doubles r;
    r.reserve(n);
    for(std::size_t i=0; i < n; ++i)
        r.push_back(tt.density(a1[i]));
    return r;
}

[[cpp11::register]]
doubles episimR_time_hazardrate(const transmission_time_R& ttr, doubles a1) {
  transmission_time& tt = *(ttr.get());
  const std::size_t n = a1.size();
  writable::doubles r;
  r.reserve(n);
  for(std::size_t i=0; i < n; ++i)
    r.push_back(tt.hazardrate(a1[i]));
  return r;
}

[[cpp11::register]]
doubles episimR_time_survivalprobability(const transmission_time_R& ttr, doubles a1) {
  transmission_time& tt = *(ttr.get());
  const std::size_t n = a1.size();
  writable::doubles r;
  r.reserve(n);
  for(std::size_t i=0; i < n; ++i)
    r.push_back(tt.survivalprobability(a1[i]));
  return r;
}

[[cpp11::register]]
doubles episimR_time_survivalquantile(const transmission_time_R& ttr, doubles a1) {
  transmission_time& tt = *(ttr.get());
  const std::size_t n = a1.size();
  writable::doubles r;
  r.reserve(n);
  for(std::size_t i=0; i < n; ++i)
    r.push_back(tt.survivalquantile(a1[i]));
  return r;
}

[[cpp11::register]]
transmission_time_R episimR_exponential_time(double lambda) {
    return new transmission_time_exponential(lambda);
}

[[cpp11::register]]
transmission_time_R episimR_lognormal_time(double mean, double var, double pinf = 0.0) {
    return new transmission_time_lognormal(mean, var, pinf);
}

[[cpp11::register]]
transmission_time_R episimR_gamma_time(double mean, double var, double pinf = 0.0) {
    return new transmission_time_gamma(mean, var, pinf);
}

struct rfunction_transmission_time : public transmission_time {
    rfunction_transmission_time(double pinf)
        :transmission_time(pinf)
    {}

    virtual interval_t sample(rng_t& e, interval_t t, int m) {
        if (sample_rf) {
            // Since the sample() R code will likely use the RNG relinquish before calling
            R_rng_relinquish rngrel;
            const interval_t r = as_cpp<double>((*sample_rf)(t, m));
            return r;
        }
        else {
            return transmission_time::sample(e, t, m);
        }
    }
    
    virtual double density(interval_t tau) {
        return as_cpp<double>(density_rf(tau));
    }
    
    virtual double survivalprobability(interval_t tau) {
        if (survivalprobability_is_trinary)
            return as_cpp<double>(survivalprobability_rf(tau, 0, 1));
        else
            return as_cpp<double>(survivalprobability_rf(tau));
    }
    
    virtual double survivalprobability(interval_t tau, interval_t t, int m) {
        if (survivalprobability_is_trinary)
            return as_cpp<double>(survivalprobability_rf(tau, t, m));
        else
            return transmission_time::survivalprobability(tau, t, m);
    }
    
    virtual double survivalquantile(double u) {
        if (survivalquantile_rf && survivalquantile_is_trinary)
            return as_cpp<double>((*survivalquantile_rf)(u, 0, 1));
        else if (survivalquantile_rf && !survivalquantile_is_trinary)
            return as_cpp<double>((*survivalquantile_rf)(u));
        else
            return transmission_time::survivalquantile(u);
    }
    
    virtual double survivalquantile(double u, interval_t t, int m) {
        if (survivalquantile_rf && survivalquantile_is_trinary)
            return as_cpp<double>((*survivalquantile_rf)(u, t, m));
        else
            return transmission_time::survivalquantile(u, t, m);
    }

    std::optional<function> sample_rf;
    function density_rf = R_NilValue;
    function survivalprobability_rf = R_NilValue;
    bool survivalprobability_is_trinary = false;
    std::optional<function> survivalquantile_rf;
    bool survivalquantile_is_trinary = false;
};

[[cpp11::register]]
transmission_time_R episimR_generic_time(
    SEXP density,
    SEXP survivalprobability, bool probability_is_trinary,
    SEXP survivalquantile, bool quantile_is_trinary,
    SEXP sample, double pinfinity = 0.0)
{
    std::unique_ptr<rfunction_transmission_time> r(new rfunction_transmission_time(pinfinity));
    
    r->density_rf = density;
    r->survivalprobability_rf = survivalprobability;
    r->survivalprobability_is_trinary = probability_is_trinary;
    if (survivalquantile != R_NilValue)
        r->survivalquantile_rf = survivalquantile;
    r->survivalquantile_is_trinary = quantile_is_trinary;
    if (sample != R_NilValue)
        r->sample_rf = sample;
        
    return r.release();
}
