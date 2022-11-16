#include "compiler.h"

#include <utility>
#include <memory>
#include <optional>

WARNINGS_DISABLE
#include <boost/scope_exit.hpp>
WARNINGS_ENABLE

#include <cpp11.hpp>
#include <cpp11/function.hpp>
#include <cpp11/external_pointer.hpp>

#include "episimR_types.h"
#include "rng.h"
#include "options.h"

#include "epidemics/types.h"
#include "epidemics/algorithm.h"
#include "epidemics/NextReaction.h"
#include "epidemics/NextReactionMeanField.h"
#include "epidemics/nMGA.h"

using namespace std::literals;
using namespace episimR;
using namespace cpp11;
namespace writable = cpp11::writable;

[[cpp11::register]]
simulation_R episimR_nextreaction_simulation(
    graph_R nw, transmission_time_R psi, sexp rho_, list opts
) {
    if (!nw) throw std::runtime_error("graph cannot be NULL"); 
    if (!psi) throw std::runtime_error("transmission time distribution cannot be NULL");

    // Reset time rho is optional, translate R_NilValue to nullptr
    transmission_time* rho = nullptr;
    if (rho_ != R_NilValue)
        rho = transmission_time_R(rho_).get();

    // Decoded options
    bool shuffle_neighbours = true;
    bool edges_concurrent = false;
    const list opts_out = process_options(opts,
        option("shuffle_neighbours", shuffle_neighbours),
        option("edges_concurrent", edges_concurrent));

    return { new simulate_next_reaction(*nw.get(), *psi.get(), rho,
                                        shuffle_neighbours, edges_concurrent),
             writable::list({"nw"_nm = nw, "psi"_nm = psi, "rho"_nm = rho_,
                             "opts"_nm = opts_out,
                             "cinf"_nm = writable::doubles { 0 },
                             "tinf"_nm = writable::doubles { 0 },
                             "trst"_nm = writable::doubles { 0 }}),
             true, true };
}

[[cpp11::register]]
simulation_R episimR_nextreaction_simulation_meanfield(
    int N, double R0, transmission_time_R psi, sexp rho_, list opts
) {
    if (!psi) throw std::runtime_error("transmission time distribution cannot be NULL");

    // Reset time rho is optional, translate R_NilValue to nullptr
    transmission_time* rho = nullptr;
    if (rho_ != R_NilValue)
        rho = transmission_time_R(rho_).get();

    // Decoded options
    const list opts_out = process_options(opts);

    return { new simulate_next_reaction_mean_field(N, R0, *psi.get(), rho),
             writable::list({"nw"_nm = R_NilValue, "psi"_nm = psi, "rho"_nm = rho_,
                             "opts"_nm = opts_out,
                             "cinf"_nm = writable::doubles { 0 },
                             "tinf"_nm = writable::doubles { 0 },
                             "trst"_nm = writable::doubles { 0 }}),
             true, true };
}

[[cpp11::register]]
simulation_R episimR_nmga_simulation(
    graph_R nw, transmission_time_R psi, sexp rho_, list opts
) {
    if (!nw) throw std::runtime_error("graph cannot be NULL"); 
    if (!psi) throw std::runtime_error("transmission time distribution cannot be NULL");

    // Reset time rho is optional, translate R_NilValue to nullptr
    transmission_time* rho = nullptr;
    if (rho_ != R_NilValue)
        rho = transmission_time_R(rho_).get();

    // Decoded options
    int approx_threshold = 100;
    double max_dt = NAN;
    double tauprec = 1e-6;
    const list opts_out = process_options(opts,
        option("approx_threshold", approx_threshold),
        option("max_dt", max_dt),
        option("tauprec", tauprec));
        
    return { new simulate_nmga(*nw.get(), *psi.get(), rho,
                               false, approx_threshold,
                               max_dt, tauprec),
             writable::list({"nw"_nm = nw, "psi"_nm = psi, "rho"_nm = rho_,
                             "opts"_nm = opts_out,
                             "cinf"_nm = writable::doubles { 0 },
                             "tinf"_nm = writable::doubles { 0 },
                             "trst"_nm = writable::doubles { 0 }}),
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
list episimR_simulation_options(const simulation_R& sim) {
    if (!sim) throw std::runtime_error("simulation cannot be NULL");
    
    return ((list)sim.protected_data())["opts"];
}

[[cpp11::register]]
list episimR_simulation_ninfected(const simulation_R& sim) {
  if (!sim) throw std::runtime_error("simulation cannot be NULL");
  
  const list sim_data = sim.protected_data();
  return writable::list({
      "total_infected"_nm = sim_data["tinf"],
      "total_reset"_nm = sim_data["trst"],
      "infected"_nm = sim_data["cinf"]
  });
}

[[cpp11::register]]
logicals episimR_simulation_isinfected(const simulation_R& sim, integers nodes) {
    if (!sim) throw std::runtime_error("simulation cannot be NULL");
    
    /* Create output */
    const std::size_t l = nodes.size();
    writable::logicals r;
    r.reserve(l);
    
    /* Fill */
    for(std::size_t j = 0; j < l; ++j)
        r.push_back(sim->is_infected(nodes[j] - 1));

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
    for(std::size_t j = 0; j < l; ++j) {
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
    
    /* Get current infection & reset counters, arrange for them to be updated at the end */
    writable::list sim_data(std::move(sim_.protected_data())); // move means modify in place
    double infected_ = ((doubles)sim_data["cinf"])[0];
    double total_infected_ = ((doubles)sim_data["tinf"])[0];
    double total_reset_ = ((doubles)sim_data["trst"])[0];
    BOOST_SCOPE_EXIT(&sim_data, &infected_, &total_infected_, &total_reset_) {
        sim_data["cinf"] = writable::doubles { infected_ };
        sim_data["tinf"] = writable::doubles { total_infected_ };
        sim_data["trst"] = writable::doubles { total_reset_ };
    } BOOST_SCOPE_EXIT_END
      
    /* Prepare output columns */
    writable::doubles times;
    times.reserve(steps);

    writable::integers kinds;
    kinds.reserve(steps);
    // Make it a factor vector, not plain integers
    writable::strings kinds_levels;
    for(int i=0; name((event_kind)i) != NULL; ++i)
      kinds_levels.push_back(name((event_kind)i));
    kinds.attr("class") = writable::strings { "factor" };
    kinds.attr("levels") = kinds_levels;

    writable::integers nodes;
    nodes.reserve(steps);

    writable::doubles total_infected;
    total_infected.reserve(steps);

    writable::doubles total_reset;
    total_reset.reserve(steps);
    
    writable::doubles infected;
    infected.reserve(steps);
    
    /* Execute steps */
    for(int i = 0; i < steps; ++i) {
        const std::optional<event_t> ev_opt = sim.step(rng_engine());
        if (!ev_opt) break;
        const event_t ev = *ev_opt;
        
        times.push_back(ev.time);
        switch (ev.kind) {
            case event_kind::outside_infection:
            case event_kind::infection:
                ++total_infected_;
                ++infected_;
                break;
            case event_kind::reset:
                ++total_reset_;
                --infected_;
                break;
            default:
              break;
        }
        kinds.push_back(name(ev.kind) ? (int)(ev.kind) + 1 : NA_INTEGER);
        nodes.push_back(ev.node + 1);
        total_infected.push_back(total_infected_);
        total_reset.push_back(total_reset_);
        infected.push_back(infected_);
    }
    
    /* Return data frame */
    return writable::data_frame({
        "time"_nm = times,
        "kind"_nm = kinds,
        "node"_nm = nodes,
        "total_infected"_nm = total_infected,
        "total_reset"_nm = total_reset,
        "infected"_nm = infected
    });
}
