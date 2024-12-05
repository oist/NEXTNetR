#pragma once

#include "safe_external_ptr.hpp"

#include "nextnet/types.h"
#include "nextnet/random.h"
#include "nextnet/network.h"
#include "nextnet/algorithm.h"

namespace nextnetR {

typedef safe_external_pointer<transmission_time> transmission_time_R;
typedef safe_external_pointer<network> network_R;
typedef safe_external_pointer<simulation_algorithm> simulation_R;
typedef safe_external_pointer<simulate_on_temporal_network> sotn_R;

} /* namespace nextnetR */

using namespace nextnetR;
