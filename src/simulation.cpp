#include <memory>
#include <optional>

#include <cpp11.hpp>
#include <cpp11/function.hpp>
#include <cpp11/external_pointer.hpp>

#include "episimR_types.h"
#include "rng.h"

#include "epidemics/types.h"
#include "epidemics/algorithm.h"
#include "epidemics/NextReaction.h"

using namespace episimR;

using namespace cpp11;
namespace writable = cpp11::writable;

[[cpp11::register]]
simulation_R episimR_nextreaction_simulation(
    graph_R nw, transmission_time_R psi, SEXP rho_
) {
    if (!nw) throw std::runtime_error("graph cannot be NULL"); 
    if (!psi) throw std::runtime_error("transmission time distribution cannot be NULL");

    transmission_time* rho = nullptr;
    if (rho_ != R_NilValue)
        rho = transmission_time_R(rho).get();

    return { new simulate_next_reaction(*nw.get(), *psi.get(), rho),
             writable::list({"nw"_nm = nw, "psi"_nm = psi, "rho"_nm = rho_ }),
             true, true };
}

[[cpp11::register]]
transmission_time_R episimR_simulation_transmissiontime(const simulation_R& sim) {
    if (!sim) throw std::runtime_error("simulation cannot be NULL");
    
    return ((list)sim.protected_data())["psi"];
}

[[cpp11::register]]
SEXP episimR_simulation_resettime(const simulation_R& sim) {
    if (!sim) throw std::runtime_error("simulation cannot be NULL");
    
    return ((list)sim.protected_data())["rho"];
}

[[cpp11::register]]
graph_R episimR_simulation_graph(const simulation_R& sim) {
    if (!sim) throw std::runtime_error("simulation cannot be NULL");
    
    return ((list)sim.protected_data())["nw"];
}

[[cpp11::register]]
bool episimR_simulation_isinfected(const simulation_R& sim, integers nodes) {
    if (!sim) throw std::runtime_error("simulation cannot be NULL");
    
    /* Create output */
    const std::size_t l = nodes.size();
    writable::logicals r;
    r.reserve(l);
    
    /* Fill */
    for(int j = 0; j < l; ++j)
        r.push_back(sim->is_infected(nodes[j] + 1));
    
    return r;
}

[[cpp11::register]]
void episimR_simulation_addinfections(const simulation_R& sim, integers nodes, doubles times) {
    if (!sim) throw std::runtime_error("simulation cannot be NULL");

    if (nodes.size() != times.size())
        throw std::runtime_error("number of nodes and times must agree");
    
    /* Convert to vector of pairs */
    const std::size_t l = nodes.size();
    std::vector<std::pair<node_t, absolutetime_t>> v;
    for(int j = 0; j < l; ++j) {
        const node_t n = nodes[j];
        if ((n < 1) || (n == NA_INTEGER))
            throw std::runtime_error("invalid node");
        
        const double t = times[j];
        if (!std::isfinite(t))
            throw std::runtime_error("infection times must be finite");
        
        v.push_back({n - 1, t});
    }
    
    sim->add_infections(v);
}

[[cpp11::register]]
data_frame episimR_simulation_step(const simulation_R& sim_, int steps) {
    if (!sim_) throw std::runtime_error("simulation cannot be NULL");
    simulation_algorithm& sim = *sim_.get();

    RNG_SCOPE_IF_NECESSARY;
    writable::doubles times;
    times.reserve(steps);
    writable::strings kinds;
    kinds.reserve(steps);
    writable::integers nodes;
    nodes.reserve(steps);

    for(int i = 0; i < steps; ++i) {
        const std::optional<event_t> ev_opt = sim.step(rng_engine());
        if (!ev_opt) break;
        const event_t ev = *ev_opt;
        
        times.push_back(ev.time);
        switch (ev.kind) {
            case event_kind::infection: kinds.push_back("infection"); break;
            case event_kind::reset: kinds.push_back("reset"); break;
            default: kinds.push_back(NA_STRING); break;
        }
        nodes.push_back(ev.node + 1);
    }
    
    return writable::data_frame({
        "time"_nm = times,
        "kind"_nm = kinds,
        "node"_nm = nodes
    });
}
