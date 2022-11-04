#pragma once

#include <cpp11/external_pointer.hpp>

#include "safe_external_ptr.hpp"

#include "epidemics/types.h"
#include "epidemics/random.h"
#include "epidemics/graph.h"
#include "epidemics/algorithm.h"

namespace episimR {

typedef safe_external_pointer<transmission_time> transmission_time_R;
typedef safe_external_pointer<graph> graph_R;
typedef safe_external_pointer<simulation_algorithm> simulation_R;

} /* namespace episimR */

using namespace episimR;
