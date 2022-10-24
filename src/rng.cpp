#include <cpp11.hpp>

#include "episimR_types.h"

using namespace cpp11;

rng_t engine;

rng_t& episimR_rng() {
    return engine;
}
