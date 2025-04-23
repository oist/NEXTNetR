#include <optional>
#include <random>

#include <R_ext/Random.h>
#include <cpp11.hpp>

#include "rng.h"

using namespace cpp11;

namespace nextnetR {

/******************************
 * factory rng_engine()
 ******************************/

rng_t engine;

rng_t& rng_engine() {
  nextnetR::r_rng_state.ensure_scope_was_entered();
  return engine;
}

/******************************
 * R_rng_state_tracker
 ******************************/

void R_rng_state_tracker::exit_scope() {
  /* Move RNG state to R if leaving the outermost scope */
  if (depth == 0)
    throw std::logic_error("RNG exit_scope() without preceeding enter_scope()");
  else if (depth == 1)
    PutRNGstate();
  --depth;
}

void R_rng_state_tracker::enter_scope() {
  /* Move RNG state to C when entering the outermost scope or when re-entering after relinquish */
  if (depth + 1 < depth)
    throw std::logic_error("Too many RNG enter_scope() calls, maximal depth exceeded");
  if (depth == 0)
    GetRNGstate();
  ++depth;
#if (RNG != RNG_CUSTOM)
  /* Initialize C++ rng when entering outermost scope but not when re-entering afer relinquish */
  if ((depth == 1) && (stack.size() == 0)) {
    std::seed_seq seed = R_rng_adapter().generate_seed<16>();
    engine.seed(seed);
  }
#endif
}

void R_rng_state_tracker::push() {
  /* Relinquish RNG to R, used before calling R functions from C++ code */
  if (depth > 0)
    PutRNGstate();
  stack.push(depth);
  depth = 0;
}

void R_rng_state_tracker::pop() {
  /* Restore rng scope nesting depth and make sure the RNG state matches it */ 
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
  /* Make sure the RNG state was moved to C */
  if (depth == 0)
    throw std::logic_error("not within an RNG scope");
}

R_rng_state_tracker r_rng_state;

/******************************
 * R_rng_adapter
 ******************************/

R_rng_adapter::result_type R_rng_adapter::operator()() {
  r_rng_state.ensure_scope_was_entered();
  return R_unif_index((double)max() + 1.0);
}

} /* namespace nextnetR */

/******************************************
 * rng_t::operator() [ epidemics/types.h ]
 ******************************************/

using namespace nextnetR;

#if RNG == RNG_CUSTOM
/* Implemenets operator() for the rng_t type defined by epidemic's types.h */
rng_t::result_type rng_t::operator()() {
  return R_rng_adapter()();
}
#endif
