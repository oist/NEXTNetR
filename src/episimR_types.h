#include <cpp11/external_pointer.hpp>

#include "epidemics/types.h"
#include "epidemics/random.h"

typedef cpp11::external_pointer<transmission_time> transmission_time_R;

rng_t& episimR_rng();
