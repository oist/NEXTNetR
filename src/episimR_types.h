#include <cpp11/external_pointer.hpp>

#include "epidemics/types.h"
#include "epidemics/random.h"
#include "epidemics/graph.h"
#include "epidemics/algorithm.h"

namespace episimR {

typedef cpp11::external_pointer<transmission_time> transmission_time_R;
typedef cpp11::external_pointer<graph> graph_R;

struct simulation_holder {
    transmission_time_R transmission_time = nullptr;
    transmission_time_R reset_time = nullptr;
    graph_R graph = nullptr;
    std::unique_ptr<simulation_algorithm> simulation = nullptr;
};

typedef cpp11::external_pointer<simulation_holder> simulation_R;

} /* namespace episimR */

using namespace episimR;
