#pragma once

#include <cstdint>
#include <random>
#include <deque>
#include <stack>

#include <cpp11.hpp>
#include <Rmath.h>

#include "nextnet/types.h"

#define RNG_SCOPE_IF_NECESSARY R_rng_scope rngscope

namespace nextnetR {

/**
 * @brief State tracker for R's RNG state
 * 
 * R's RNG state can either be stored on the R side of things, or the C side.
 * Before random number generation function can be called from C, the state
 * must be moved to the C side with GetRNGstate(), and must be moved back with
 * PutRNGState() before R code invokes the RNG again.
 * 
 * This state tracker simplifies stat. enter_scope() moves the state to the
 * C side, unless it is already there. Similarly, exit_scope() moves the state
 * back once the last scope is exited (i.e. depth reaches zero).
 * 
 * push() pushes the current depth onto the depth stack, then reset the depth
 * to zero and moves the state to the R side. pop() restores the most recently
 * pushed depth and moves the state accordingly. push() and pop() allow C
 * code to safely call back into R code that might invoke the rng; such C code
 * should relinquish the RNG by calling push() before invoking R code, and
 * should restore the 
 */
struct R_rng_state_tracker {
  void enter_scope();

  void exit_scope();
  
  void push();

  void pop();
  
  void ensure_scope_was_entered();
  
  unsigned int current_depth() { return depth; }

private:
  unsigned int depth = 0;
  std::stack<unsigned int> stack;
};

extern R_rng_state_tracker r_rng_state;

/**
 * @brief RAII helper to enter an RNG scope, i.e. to move the state to C
 */
struct R_rng_scope {
  R_rng_scope() :depth(r_rng_state.current_depth()) { r_rng_state.enter_scope(); }
  
  ~R_rng_scope() noexcept(false) {
    r_rng_state.exit_scope();
    if (depth != r_rng_state.current_depth())
      throw std::logic_error("invalid nesting of R_rng_scope lifetimes");
  }
  
private:
  unsigned int depth;
};

/**
 * @brief RAII helper to relinquish the RNG before calling back into R
 */
struct R_rng_relinquish {
  R_rng_relinquish() :depth(r_rng_state.current_depth()) { r_rng_state.push(); }
  
  ~R_rng_relinquish() noexcept(false) {
    r_rng_state.pop();
    if (depth != r_rng_state.current_depth())
      throw std::logic_error("invalid nesting of R_rng_relinquish lifetimes");
  }
  
private:
  unsigned int depth;
};

/**
 * @brief Returns the RNG engine singleton
 */
rng_t& rng_engine();

/**
 * @brief Adapter to use R's RNG with C++11's random distributions.
 * 
 * Implement an UniformRandomBitGenerator which can be passed as
 * a generator to the distribution's operator()
 * 
 * This implementation is currently quite inefficient because it
 * enters and exits an RNGScope (i.e. calls GetRNGState and
 * PutRNGState) every time operator() is called. However, the
 * confusing state of affairs regarding RNGScope makes it hard to
 * do better without making a lot of assumptions about how this
 * RNG is used.
 *
 * The generate_seed() method can be used to seed another RNG with
 * the values produced by R's RNG; one way to circumvent to problems
 * outlined above is to only use R's RNG to seed a native C++ RNG.
 */
struct R_rng_adapter {
  typedef std::uint32_t result_type;
  
  static constexpr result_type min() { return 0; }
  
  static constexpr result_type max() { return UINT32_MAX; }
  
  result_type operator()();
  
  template<std::size_t L>
  std::seed_seq generate_seed() {
    nextnetR::r_rng_state.ensure_scope_was_entered();
    // Generate specified number of uint32 values
    std::array<std::uint32_t, L> s;
    for(std::size_t i=0; i < L; ++i)
      s[i] = R_unif_index((double)max() + 1.0);
    // Convert to std::seed_seq
    return std::seed_seq(s.begin(), s.end());
  }
};

} /* namespace nextnetR */
