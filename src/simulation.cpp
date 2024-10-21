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
    simulate_next_reaction::params p;
    bool static_network = false;
    const list opts_out = process_options(opts,
        option("shuffle_neighbours", p.shuffle_neighbours),
        option("edges_concurrent", p.edges_concurrent),
        option("SIR", p.SIR),
        option("static_network", static_network));

    // Create simulation algorithm
    simulation_algorithm* sim = new simulate_next_reaction(*nw.get(), *psi.get(), rho, p);

    // Create dynamic network simulator unless network is static, or override was set
    const bool use_sodn = (dynamic_cast<dynamic_network*>(nw.get()) != nullptr) && !static_network;
    simulate_on_dynamic_network* sodn = use_sodn ? new simulate_on_dynamic_network(*sim) : nullptr;

    // Return simulation algorithm, store dynamic network simulator in meta-data list
    return { sim,
             writable::list({"nw"_nm = nw, "psi"_nm = psi, "rho"_nm = rho_,
                             "sodn"_nm = use_sodn ? (sexp)sodn_R(sodn) : (sexp)R_NilValue,
                             "opts"_nm = opts_out,
                             "step"_nm = writable::integers { 0 },
                             "nwst"_nm = writable::integers { 0 },
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
                             "sodn"_nm = R_NilValue,
                             "opts"_nm = opts_out,
                             "step"_nm = writable::integers { 0 },
                             "nwst"_nm = writable::integers { 0 },
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
    simulate_nmga::params params;
    const list opts_out = process_options(opts,
        option("approx_threshold", params.approximation_threshold),
        option("max_dt", params.maximal_dt),
        option("tauprec", params.tau_precision),
        option("SIR", params.SIR));
        
    return { new simulate_nmga(*nw.get(), *psi.get(), rho, params),
             writable::list({"nw"_nm = nw, "psi"_nm = psi, "rho"_nm = rho_,
                             "sodn"_nm = R_NilValue,
                             "opts"_nm = opts_out,
                             "step"_nm = writable::integers { 0 },
                             "nwst"_nm = writable::integers { 0 },
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
data_frame episimR_simulation_run(const simulation_R& sim_, list stop, list opts) {
    RNG_SCOPE_IF_NECESSARY;

    if (!sim_) throw std::runtime_error("simulation cannot be NULL");

    // Decode stopping condition
    int stop_epidemic_steps = NA_INTEGER;
    int stop_network_steps = NA_INTEGER;
    double stop_time = R_PosInf;
    double stop_infected = R_PosInf;
    double stop_total_infected = R_PosInf;
    double stop_total_reset = R_PosInf;
    process_options(stop,
        option("epidemic_steps", stop_epidemic_steps),
        option("network_steps", stop_network_steps),
        option("time", stop_time),
        option("infected", stop_infected),
        option("total_infected", stop_total_infected),
        option("total_reset", stop_total_reset));

    // Decode options
    bool epidemic_events = true;
    bool network_events = false;
    process_options(opts,
        option("network_events", network_events),
        option("epidemic_events", epidemic_events));

    // Get simulation alorithm
    simulation_algorithm& sim = *sim_.get();
    writable::list sim_data(std::move(sim_.protected_data())); // move means modify in place

    // Get dynamic network simulation algorithm
    const sexp sodn_sexp = (sexp)sim_data["sodn"];
    simulate_on_dynamic_network* const sodn = (
        (sodn_sexp != R_NilValue) ? sodn_R(sodn_sexp).get() : nullptr);

    // Get current infection & reset counters, arrange for them to be updated at the end
    double time = R_NegInf;
    int epidemic_step_ = ((integers)sim_data["step"])[0];
    int network_step_ = ((integers)sim_data["nwst"])[0];
    double infected_ = ((doubles)sim_data["cinf"])[0];
    double total_infected_ = ((doubles)sim_data["tinf"])[0];
    double total_reset_ = ((doubles)sim_data["trst"])[0];
    BOOST_SCOPE_EXIT(&sim_data, &epidemic_step_, &network_step_, &infected_, &total_infected_, &total_reset_) {
        sim_data["step"] = writable::integers { epidemic_step_ };
        sim_data["nwst"] = writable::integers { network_step_ };
        sim_data["cinf"] = writable::doubles { infected_ };
        sim_data["tinf"] = writable::doubles { total_infected_ };
        sim_data["trst"] = writable::doubles { total_reset_ };
    } BOOST_SCOPE_EXIT_END
      
    // Prepare output columns
    writable::doubles times;
    writable::integers epidemic_steps;
    writable::integers network_steps;
    writable::integers nodes;
    writable::integers neighbours;
    writable::doubles total_infected;
    writable::doubles total_reset;
    writable::doubles infected;

    // kinds are a factor vector
    writable::integers kinds;
    kinds.reserve(1); // make sure the SEXP is not null
    writable::strings kinds_levels;
    int network_event_offset = 0;
    for(int i=0; name((event_kind)i) != NULL; ++i, ++network_event_offset)
      kinds_levels.push_back(name((event_kind)i));
    for(int i=0; name((network_event_kind)i) != NULL; ++i)
      kinds_levels.push_back(name((network_event_kind)i));

    // Execute steps
    while(true) {
        // Check stopping conditions
        check_user_interrupt();
        if ((stop_epidemic_steps != NA_INTEGER) && (epidemic_step_ >= stop_epidemic_steps))
            break;
        if ((stop_network_steps != NA_INTEGER) && (network_step_ >= stop_network_steps))
            break;
        if (time > stop_time)
            break;
        if (infected_ >= stop_infected)
            break;
        if (total_infected_ >= stop_total_infected)
            break;
        if (total_reset_ >= stop_total_reset)
            break;

        // Perform step, use dynamic network simulator if available
        std::optional<network_or_epidemic_event_t> any_ev_opt;
        if (sodn != nullptr)
            any_ev_opt = sodn->step(rng_engine());
        else
            any_ev_opt = sim.step(rng_engine());

        // Stop if there are no more events
        if (!any_ev_opt)
            break;
        const network_or_epidemic_event_t any_ev = *any_ev_opt;

        // Current row
        int kind;
        int node;
        int neighbour = NA_INTEGER;

        if (std::holds_alternative<event_t>(any_ev)) {
            // Epidemic event
            const auto ev = std::get<event_t>(any_ev);
            // Update state
            time = ev.time;
            ++epidemic_step_;
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
            // Report event?
            if (!epidemic_events)
                continue;
            // Fill row
            kind = name(ev.kind) ? (int)(ev.kind) + 1 : NA_INTEGER;
            node = ev.node + 1;
        }
        else if (std::holds_alternative<network_event_t>(any_ev)) {
            // Network event
            const auto ev = std::get<network_event_t>(any_ev);
            // Update state
            time = ev.time;
            ++network_step_;
            // Report event?
            if (!network_events)
                continue;
            // Fill row
            kind = name(ev.kind) ? (int)(ev.kind) + network_event_offset + 1 : NA_INTEGER;
            node = ev.source_node + 1;
            neighbour = ev.target_node + 1;
        }
        else throw std::logic_error("unknown event type");

        // Append row
        times.push_back(time);
        epidemic_steps.push_back(epidemic_step_);
        network_steps.push_back(network_step_);
        kinds.push_back(kind);
        nodes.push_back(node);
        neighbours.push_back(neighbour);
        total_infected.push_back(total_infected_);
        total_reset.push_back(total_reset_);
        infected.push_back(infected_);
    }
    
    // Make kinds a factor
    kinds.attr("class") = writable::strings { "factor" };
    kinds.attr("levels") = kinds_levels;

    // Return data frame
    return writable::data_frame({
        "time"_nm = times,
        "epidemic_step"_nm = epidemic_steps,
        "network_step"_nm = network_steps,
        "kind"_nm = kinds,
        "node"_nm = nodes,
        "neighbour"_nm = neighbours,
        "total_infected"_nm = total_infected,
        "total_reset"_nm = total_reset,
        "infected"_nm = infected
    });
}
