#include <optional>
#include <random>

#include <R_ext/Random.h>
#include <cpp11.hpp>

#include "rng.h"

using namespace cpp11;

namespace episimR {

/******************************
 * R_rng_state_tracker
 ******************************/

void R_rng_state_tracker::exit_scope() {
  if (depth == 0)
    throw std::logic_error("RNG exit_scope() without preceeding enter_scope()");
  else if (depth == 1)
    PutRNGstate();
  --depth;
}

void R_rng_state_tracker::enter_scope() {
  if (depth == 0)
    GetRNGstate();
  else if (depth + 1 < depth)
    throw std::logic_error("Too many RNG enter_scope() calls, maximal depth exceeded");
  ++depth;
}

void R_rng_state_tracker::push() {
  if (depth > 0)
    PutRNGstate();
  stack.push(depth);
  depth = 0;
}

void R_rng_state_tracker::pop() {
  if (stack.empty())
    throw std::logic_error("RNG scope pop() without preceeding pop()");
  
  if ((depth > 0) && (stack.top() == 0))
    PutRNGstate();
  else if ((depth == 0) && (stack.top() > 0))
    GetRNGstate();
  
  depth = stack.top();
  stack.pop();
}

void R_rng_state_tracker::ensure_scope_was_entered() {
  if (depth == 0)
    throw std::logic_error("not within an RNG scope");
}

R_rng_state_tracker r_rng_state;

/******************************
 * rng_engine
 ******************************/

#if (RNG != RNG_CUSTOM)
std::optional<rng_t> engine;
#else
rng_t engine;
#endif

rng_t& rng_engine() {
#if (RNG != RNG_CUSTOM)
  /* Use the R RNG to initialize the C++ RNG once */
  if (!engine) {
    R_rng_scope
    std::seed_seq seed = R_rng_adapter().generate_seed<16>();
    engine = rng_t(seed);
  }
  return *engine;
#else
  /* Use R RNG adapter as the C++ RNG */
  return engine;
#endif
}

} /* namespace episimR */

/******************************
 * episimR_R_rng
 ******************************/

using namespace episimR;

/* Implemenets operator() for the rng_t type defined by epidemic's types.h */
rng_t::result_type rng_t::operator()() {
  r_rng_state.ensure_scope_was_entered();
  return R_unif_index((double)max() + 1.0);
}
