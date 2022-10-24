#include <cpp11/external_pointer.hpp>

#include "epidemics/types.h"
#include "epidemics/random.h"
#include "epidemics/graph.h"

namespace episimR {

typedef cpp11::external_pointer<transmission_time> transmission_time_R;
typedef cpp11::external_pointer<graph> graph_R;

} /* namespace episimR */

using namespace episimR;
