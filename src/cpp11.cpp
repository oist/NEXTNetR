// Generated by cpp11: do not edit by hand
// clang-format off

#include "episimR_types.h"
#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// graph.cpp
integers episimR_graph_outdegree(const graph_R& nw, integers nodes);
extern "C" SEXP _episimR_episimR_graph_outdegree(SEXP nw, SEXP nodes) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_graph_outdegree(cpp11::as_cpp<cpp11::decay_t<const graph_R&>>(nw), cpp11::as_cpp<cpp11::decay_t<integers>>(nodes)));
  END_CPP11
}
// graph.cpp
integers episimR_graph_neighbour(const graph_R& nw, integers nodes, integers indices);
extern "C" SEXP _episimR_episimR_graph_neighbour(SEXP nw, SEXP nodes, SEXP indices) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_graph_neighbour(cpp11::as_cpp<cpp11::decay_t<const graph_R&>>(nw), cpp11::as_cpp<cpp11::decay_t<integers>>(nodes), cpp11::as_cpp<cpp11::decay_t<integers>>(indices)));
  END_CPP11
}
// graph.cpp
list episimR_graph_adjacencylist(const graph_R& nw);
extern "C" SEXP _episimR_episimR_graph_adjacencylist(SEXP nw) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_graph_adjacencylist(cpp11::as_cpp<cpp11::decay_t<const graph_R&>>(nw)));
  END_CPP11
}
// graph.cpp
list episimR_graph_bounds(const graph_R& nw);
extern "C" SEXP _episimR_episimR_graph_bounds(SEXP nw) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_graph_bounds(cpp11::as_cpp<cpp11::decay_t<const graph_R&>>(nw)));
  END_CPP11
}
// graph.cpp
doubles_matrix<> episimR_graph_coordinates(const graph_R& nw, integers nodes);
extern "C" SEXP _episimR_episimR_graph_coordinates(SEXP nw, SEXP nodes) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_graph_coordinates(cpp11::as_cpp<cpp11::decay_t<const graph_R&>>(nw), cpp11::as_cpp<cpp11::decay_t<integers>>(nodes)));
  END_CPP11
}
// graph.cpp
graph_R episimR_erdos_reyni_graph(int size, double avg_degree);
extern "C" SEXP _episimR_episimR_erdos_reyni_graph(SEXP size, SEXP avg_degree) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_erdos_reyni_graph(cpp11::as_cpp<cpp11::decay_t<int>>(size), cpp11::as_cpp<cpp11::decay_t<double>>(avg_degree)));
  END_CPP11
}
// graph.cpp
graph_R episimR_fully_connected_graph(int size);
extern "C" SEXP _episimR_episimR_fully_connected_graph(SEXP size) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_fully_connected_graph(cpp11::as_cpp<cpp11::decay_t<int>>(size)));
  END_CPP11
}
// graph.cpp
graph_R episimR_acyclic_graph(int size, double avg_degree, bool reduced_root_degree);
extern "C" SEXP _episimR_episimR_acyclic_graph(SEXP size, SEXP avg_degree, SEXP reduced_root_degree) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_acyclic_graph(cpp11::as_cpp<cpp11::decay_t<int>>(size), cpp11::as_cpp<cpp11::decay_t<double>>(avg_degree), cpp11::as_cpp<cpp11::decay_t<bool>>(reduced_root_degree)));
  END_CPP11
}
// graph.cpp
graph_R episimR_configmodel_graph(integers degrees);
extern "C" SEXP _episimR_episimR_configmodel_graph(SEXP degrees) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_configmodel_graph(cpp11::as_cpp<cpp11::decay_t<integers>>(degrees)));
  END_CPP11
}
// graph.cpp
graph_R episimR_configmodel_clustered_alpha_graph(integers degrees, double alpha, double beta);
extern "C" SEXP _episimR_episimR_configmodel_clustered_alpha_graph(SEXP degrees, SEXP alpha, SEXP beta) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_configmodel_clustered_alpha_graph(cpp11::as_cpp<cpp11::decay_t<integers>>(degrees), cpp11::as_cpp<cpp11::decay_t<double>>(alpha), cpp11::as_cpp<cpp11::decay_t<double>>(beta)));
  END_CPP11
}
// graph.cpp
graph_R episimR_configmodel_clustered_ck_graph(integers degrees, SEXP ck, double beta);
extern "C" SEXP _episimR_episimR_configmodel_clustered_ck_graph(SEXP degrees, SEXP ck, SEXP beta) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_configmodel_clustered_ck_graph(cpp11::as_cpp<cpp11::decay_t<integers>>(degrees), cpp11::as_cpp<cpp11::decay_t<SEXP>>(ck), cpp11::as_cpp<cpp11::decay_t<double>>(beta)));
  END_CPP11
}
// graph.cpp
graph_R episimR_configmodel_clustered_triangles_graph(integers degrees, integers triangles, double beta);
extern "C" SEXP _episimR_episimR_configmodel_clustered_triangles_graph(SEXP degrees, SEXP triangles, SEXP beta) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_configmodel_clustered_triangles_graph(cpp11::as_cpp<cpp11::decay_t<integers>>(degrees), cpp11::as_cpp<cpp11::decay_t<integers>>(triangles), cpp11::as_cpp<cpp11::decay_t<double>>(beta)));
  END_CPP11
}
// graph.cpp
graph_R episimR_barabasialbert_graph(int size, int m);
extern "C" SEXP _episimR_episimR_barabasialbert_graph(SEXP size, SEXP m) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_barabasialbert_graph(cpp11::as_cpp<cpp11::decay_t<int>>(size), cpp11::as_cpp<cpp11::decay_t<int>>(m)));
  END_CPP11
}
// graph.cpp
graph_R episimR_cubiclattice2d_graph(int edge_length);
extern "C" SEXP _episimR_episimR_cubiclattice2d_graph(SEXP edge_length) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_cubiclattice2d_graph(cpp11::as_cpp<cpp11::decay_t<int>>(edge_length)));
  END_CPP11
}
// graph.cpp
graph_R episimR_cubiclattice3d_graph(int edge_length);
extern "C" SEXP _episimR_episimR_cubiclattice3d_graph(SEXP edge_length) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_cubiclattice3d_graph(cpp11::as_cpp<cpp11::decay_t<int>>(edge_length)));
  END_CPP11
}
// graph.cpp
graph_R episimR_cubiclattice4d_graph(int edge_length);
extern "C" SEXP _episimR_episimR_cubiclattice4d_graph(SEXP edge_length) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_cubiclattice4d_graph(cpp11::as_cpp<cpp11::decay_t<int>>(edge_length)));
  END_CPP11
}
// graph.cpp
graph_R episimR_cubiclattice5d_graph(int edge_length);
extern "C" SEXP _episimR_episimR_cubiclattice5d_graph(SEXP edge_length) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_cubiclattice5d_graph(cpp11::as_cpp<cpp11::decay_t<int>>(edge_length)));
  END_CPP11
}
// graph.cpp
graph_R episimR_cubiclattice6d_graph(int edge_length);
extern "C" SEXP _episimR_episimR_cubiclattice6d_graph(SEXP edge_length) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_cubiclattice6d_graph(cpp11::as_cpp<cpp11::decay_t<int>>(edge_length)));
  END_CPP11
}
// graph.cpp
graph_R episimR_cubiclattice7d_graph(int edge_length);
extern "C" SEXP _episimR_episimR_cubiclattice7d_graph(SEXP edge_length) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_cubiclattice7d_graph(cpp11::as_cpp<cpp11::decay_t<int>>(edge_length)));
  END_CPP11
}
// graph.cpp
graph_R episimR_cubiclattice8d_graph(int edge_length);
extern "C" SEXP _episimR_episimR_cubiclattice8d_graph(SEXP edge_length) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_cubiclattice8d_graph(cpp11::as_cpp<cpp11::decay_t<int>>(edge_length)));
  END_CPP11
}
// graph.cpp
graph_R episimR_brownian_proximity_dyngraph(int size, double avg_degree, double radius, double D, SEXP dt);
extern "C" SEXP _episimR_episimR_brownian_proximity_dyngraph(SEXP size, SEXP avg_degree, SEXP radius, SEXP D, SEXP dt) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_brownian_proximity_dyngraph(cpp11::as_cpp<cpp11::decay_t<int>>(size), cpp11::as_cpp<cpp11::decay_t<double>>(avg_degree), cpp11::as_cpp<cpp11::decay_t<double>>(radius), cpp11::as_cpp<cpp11::decay_t<double>>(D), cpp11::as_cpp<cpp11::decay_t<SEXP>>(dt)));
  END_CPP11
}
// graph.cpp
graph_R episimR_empirical_dyngraph(std::string file, bool finite_duration, double dt);
extern "C" SEXP _episimR_episimR_empirical_dyngraph(SEXP file, SEXP finite_duration, SEXP dt) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_empirical_dyngraph(cpp11::as_cpp<cpp11::decay_t<std::string>>(file), cpp11::as_cpp<cpp11::decay_t<bool>>(finite_duration), cpp11::as_cpp<cpp11::decay_t<double>>(dt)));
  END_CPP11
}
// graph.cpp
graph_R episimR_stored_graph(r_string filename);
extern "C" SEXP _episimR_episimR_stored_graph(SEXP filename) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_stored_graph(cpp11::as_cpp<cpp11::decay_t<r_string>>(filename)));
  END_CPP11
}
// graph.cpp
graph_R episimR_userdefined_graph(list input_al);
extern "C" SEXP _episimR_episimR_userdefined_graph(SEXP input_al) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_userdefined_graph(cpp11::as_cpp<cpp11::decay_t<list>>(input_al)));
  END_CPP11
}
// simulation.cpp
simulation_R episimR_nextreaction_simulation(graph_R nw, transmission_time_R psi, sexp rho_, list opts);
extern "C" SEXP _episimR_episimR_nextreaction_simulation(SEXP nw, SEXP psi, SEXP rho_, SEXP opts) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_nextreaction_simulation(cpp11::as_cpp<cpp11::decay_t<graph_R>>(nw), cpp11::as_cpp<cpp11::decay_t<transmission_time_R>>(psi), cpp11::as_cpp<cpp11::decay_t<sexp>>(rho_), cpp11::as_cpp<cpp11::decay_t<list>>(opts)));
  END_CPP11
}
// simulation.cpp
simulation_R episimR_nextreaction_simulation_meanfield(int N, double R0, transmission_time_R psi, sexp rho_, list opts);
extern "C" SEXP _episimR_episimR_nextreaction_simulation_meanfield(SEXP N, SEXP R0, SEXP psi, SEXP rho_, SEXP opts) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_nextreaction_simulation_meanfield(cpp11::as_cpp<cpp11::decay_t<int>>(N), cpp11::as_cpp<cpp11::decay_t<double>>(R0), cpp11::as_cpp<cpp11::decay_t<transmission_time_R>>(psi), cpp11::as_cpp<cpp11::decay_t<sexp>>(rho_), cpp11::as_cpp<cpp11::decay_t<list>>(opts)));
  END_CPP11
}
// simulation.cpp
simulation_R episimR_nmga_simulation(graph_R nw, transmission_time_R psi, sexp rho_, list opts);
extern "C" SEXP _episimR_episimR_nmga_simulation(SEXP nw, SEXP psi, SEXP rho_, SEXP opts) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_nmga_simulation(cpp11::as_cpp<cpp11::decay_t<graph_R>>(nw), cpp11::as_cpp<cpp11::decay_t<transmission_time_R>>(psi), cpp11::as_cpp<cpp11::decay_t<sexp>>(rho_), cpp11::as_cpp<cpp11::decay_t<list>>(opts)));
  END_CPP11
}
// simulation.cpp
transmission_time_R episimR_simulation_transmissiontime(const simulation_R& sim);
extern "C" SEXP _episimR_episimR_simulation_transmissiontime(SEXP sim) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_simulation_transmissiontime(cpp11::as_cpp<cpp11::decay_t<const simulation_R&>>(sim)));
  END_CPP11
}
// simulation.cpp
SEXP episimR_simulation_resettime(const simulation_R& sim);
extern "C" SEXP _episimR_episimR_simulation_resettime(SEXP sim) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_simulation_resettime(cpp11::as_cpp<cpp11::decay_t<const simulation_R&>>(sim)));
  END_CPP11
}
// simulation.cpp
graph_R episimR_simulation_graph(const simulation_R& sim);
extern "C" SEXP _episimR_episimR_simulation_graph(SEXP sim) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_simulation_graph(cpp11::as_cpp<cpp11::decay_t<const simulation_R&>>(sim)));
  END_CPP11
}
// simulation.cpp
list episimR_simulation_options(const simulation_R& sim);
extern "C" SEXP _episimR_episimR_simulation_options(SEXP sim) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_simulation_options(cpp11::as_cpp<cpp11::decay_t<const simulation_R&>>(sim)));
  END_CPP11
}
// simulation.cpp
list episimR_simulation_ninfected(const simulation_R& sim);
extern "C" SEXP _episimR_episimR_simulation_ninfected(SEXP sim) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_simulation_ninfected(cpp11::as_cpp<cpp11::decay_t<const simulation_R&>>(sim)));
  END_CPP11
}
// simulation.cpp
logicals episimR_simulation_isinfected(const simulation_R& sim, integers nodes);
extern "C" SEXP _episimR_episimR_simulation_isinfected(SEXP sim, SEXP nodes) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_simulation_isinfected(cpp11::as_cpp<cpp11::decay_t<const simulation_R&>>(sim), cpp11::as_cpp<cpp11::decay_t<integers>>(nodes)));
  END_CPP11
}
// simulation.cpp
void episimR_simulation_addinfections(const simulation_R& sim, integers nodes, doubles times);
extern "C" SEXP _episimR_episimR_simulation_addinfections(SEXP sim, SEXP nodes, SEXP times) {
  BEGIN_CPP11
    episimR_simulation_addinfections(cpp11::as_cpp<cpp11::decay_t<const simulation_R&>>(sim), cpp11::as_cpp<cpp11::decay_t<integers>>(nodes), cpp11::as_cpp<cpp11::decay_t<doubles>>(times));
    return R_NilValue;
  END_CPP11
}
// simulation.cpp
data_frame episimR_simulation_run(const simulation_R& sim_, list stop, list opts);
extern "C" SEXP _episimR_episimR_simulation_run(SEXP sim_, SEXP stop, SEXP opts) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_simulation_run(cpp11::as_cpp<cpp11::decay_t<const simulation_R&>>(sim_), cpp11::as_cpp<cpp11::decay_t<list>>(stop), cpp11::as_cpp<cpp11::decay_t<list>>(opts)));
  END_CPP11
}
// transmission_time.cpp
doubles episimR_time_sample(int n, const transmission_time_R& ttr, interval_t t, int m);
extern "C" SEXP _episimR_episimR_time_sample(SEXP n, SEXP ttr, SEXP t, SEXP m) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_time_sample(cpp11::as_cpp<cpp11::decay_t<int>>(n), cpp11::as_cpp<cpp11::decay_t<const transmission_time_R&>>(ttr), cpp11::as_cpp<cpp11::decay_t<interval_t>>(t), cpp11::as_cpp<cpp11::decay_t<int>>(m)));
  END_CPP11
}
// transmission_time.cpp
doubles episimR_time_density(const transmission_time_R& ttr, doubles taus);
extern "C" SEXP _episimR_episimR_time_density(SEXP ttr, SEXP taus) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_time_density(cpp11::as_cpp<cpp11::decay_t<const transmission_time_R&>>(ttr), cpp11::as_cpp<cpp11::decay_t<doubles>>(taus)));
  END_CPP11
}
// transmission_time.cpp
doubles episimR_time_hazardrate(const transmission_time_R& ttr, doubles taus);
extern "C" SEXP _episimR_episimR_time_hazardrate(SEXP ttr, SEXP taus) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_time_hazardrate(cpp11::as_cpp<cpp11::decay_t<const transmission_time_R&>>(ttr), cpp11::as_cpp<cpp11::decay_t<doubles>>(taus)));
  END_CPP11
}
// transmission_time.cpp
doubles episimR_time_survivalprobability(const transmission_time_R& ttr, doubles taus, doubles ts, integers ms);
extern "C" SEXP _episimR_episimR_time_survivalprobability(SEXP ttr, SEXP taus, SEXP ts, SEXP ms) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_time_survivalprobability(cpp11::as_cpp<cpp11::decay_t<const transmission_time_R&>>(ttr), cpp11::as_cpp<cpp11::decay_t<doubles>>(taus), cpp11::as_cpp<cpp11::decay_t<doubles>>(ts), cpp11::as_cpp<cpp11::decay_t<integers>>(ms)));
  END_CPP11
}
// transmission_time.cpp
doubles episimR_time_survivalquantile(const transmission_time_R& ttr, doubles ps, doubles ts, integers ms);
extern "C" SEXP _episimR_episimR_time_survivalquantile(SEXP ttr, SEXP ps, SEXP ts, SEXP ms) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_time_survivalquantile(cpp11::as_cpp<cpp11::decay_t<const transmission_time_R&>>(ttr), cpp11::as_cpp<cpp11::decay_t<doubles>>(ps), cpp11::as_cpp<cpp11::decay_t<doubles>>(ts), cpp11::as_cpp<cpp11::decay_t<integers>>(ms)));
  END_CPP11
}
// transmission_time.cpp
transmission_time_R episimR_exponential_time(double lambda);
extern "C" SEXP _episimR_episimR_exponential_time(SEXP lambda) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_exponential_time(cpp11::as_cpp<cpp11::decay_t<double>>(lambda)));
  END_CPP11
}
// transmission_time.cpp
transmission_time_R episimR_lognormal_time(double mean, double var, double pinf);
extern "C" SEXP _episimR_episimR_lognormal_time(SEXP mean, SEXP var, SEXP pinf) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_lognormal_time(cpp11::as_cpp<cpp11::decay_t<double>>(mean), cpp11::as_cpp<cpp11::decay_t<double>>(var), cpp11::as_cpp<cpp11::decay_t<double>>(pinf)));
  END_CPP11
}
// transmission_time.cpp
transmission_time_R episimR_gamma_time(double mean, double var, double pinf);
extern "C" SEXP _episimR_episimR_gamma_time(SEXP mean, SEXP var, SEXP pinf) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_gamma_time(cpp11::as_cpp<cpp11::decay_t<double>>(mean), cpp11::as_cpp<cpp11::decay_t<double>>(var), cpp11::as_cpp<cpp11::decay_t<double>>(pinf)));
  END_CPP11
}
// transmission_time.cpp
transmission_time_R episimR_generic_time(SEXP density, SEXP survivalprobability, bool probability_is_trinary, SEXP survivalquantile, bool quantile_is_trinary, SEXP sample, double pinfinity);
extern "C" SEXP _episimR_episimR_generic_time(SEXP density, SEXP survivalprobability, SEXP probability_is_trinary, SEXP survivalquantile, SEXP quantile_is_trinary, SEXP sample, SEXP pinfinity) {
  BEGIN_CPP11
    return cpp11::as_sexp(episimR_generic_time(cpp11::as_cpp<cpp11::decay_t<SEXP>>(density), cpp11::as_cpp<cpp11::decay_t<SEXP>>(survivalprobability), cpp11::as_cpp<cpp11::decay_t<bool>>(probability_is_trinary), cpp11::as_cpp<cpp11::decay_t<SEXP>>(survivalquantile), cpp11::as_cpp<cpp11::decay_t<bool>>(quantile_is_trinary), cpp11::as_cpp<cpp11::decay_t<SEXP>>(sample), cpp11::as_cpp<cpp11::decay_t<double>>(pinfinity)));
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_episimR_episimR_acyclic_graph",                         (DL_FUNC) &_episimR_episimR_acyclic_graph,                         3},
    {"_episimR_episimR_barabasialbert_graph",                  (DL_FUNC) &_episimR_episimR_barabasialbert_graph,                  2},
    {"_episimR_episimR_brownian_proximity_dyngraph",           (DL_FUNC) &_episimR_episimR_brownian_proximity_dyngraph,           5},
    {"_episimR_episimR_configmodel_clustered_alpha_graph",     (DL_FUNC) &_episimR_episimR_configmodel_clustered_alpha_graph,     3},
    {"_episimR_episimR_configmodel_clustered_ck_graph",        (DL_FUNC) &_episimR_episimR_configmodel_clustered_ck_graph,        3},
    {"_episimR_episimR_configmodel_clustered_triangles_graph", (DL_FUNC) &_episimR_episimR_configmodel_clustered_triangles_graph, 3},
    {"_episimR_episimR_configmodel_graph",                     (DL_FUNC) &_episimR_episimR_configmodel_graph,                     1},
    {"_episimR_episimR_cubiclattice2d_graph",                  (DL_FUNC) &_episimR_episimR_cubiclattice2d_graph,                  1},
    {"_episimR_episimR_cubiclattice3d_graph",                  (DL_FUNC) &_episimR_episimR_cubiclattice3d_graph,                  1},
    {"_episimR_episimR_cubiclattice4d_graph",                  (DL_FUNC) &_episimR_episimR_cubiclattice4d_graph,                  1},
    {"_episimR_episimR_cubiclattice5d_graph",                  (DL_FUNC) &_episimR_episimR_cubiclattice5d_graph,                  1},
    {"_episimR_episimR_cubiclattice6d_graph",                  (DL_FUNC) &_episimR_episimR_cubiclattice6d_graph,                  1},
    {"_episimR_episimR_cubiclattice7d_graph",                  (DL_FUNC) &_episimR_episimR_cubiclattice7d_graph,                  1},
    {"_episimR_episimR_cubiclattice8d_graph",                  (DL_FUNC) &_episimR_episimR_cubiclattice8d_graph,                  1},
    {"_episimR_episimR_empirical_dyngraph",                    (DL_FUNC) &_episimR_episimR_empirical_dyngraph,                    3},
    {"_episimR_episimR_erdos_reyni_graph",                     (DL_FUNC) &_episimR_episimR_erdos_reyni_graph,                     2},
    {"_episimR_episimR_exponential_time",                      (DL_FUNC) &_episimR_episimR_exponential_time,                      1},
    {"_episimR_episimR_fully_connected_graph",                 (DL_FUNC) &_episimR_episimR_fully_connected_graph,                 1},
    {"_episimR_episimR_gamma_time",                            (DL_FUNC) &_episimR_episimR_gamma_time,                            3},
    {"_episimR_episimR_generic_time",                          (DL_FUNC) &_episimR_episimR_generic_time,                          7},
    {"_episimR_episimR_graph_adjacencylist",                   (DL_FUNC) &_episimR_episimR_graph_adjacencylist,                   1},
    {"_episimR_episimR_graph_bounds",                          (DL_FUNC) &_episimR_episimR_graph_bounds,                          1},
    {"_episimR_episimR_graph_coordinates",                     (DL_FUNC) &_episimR_episimR_graph_coordinates,                     2},
    {"_episimR_episimR_graph_neighbour",                       (DL_FUNC) &_episimR_episimR_graph_neighbour,                       3},
    {"_episimR_episimR_graph_outdegree",                       (DL_FUNC) &_episimR_episimR_graph_outdegree,                       2},
    {"_episimR_episimR_lognormal_time",                        (DL_FUNC) &_episimR_episimR_lognormal_time,                        3},
    {"_episimR_episimR_nextreaction_simulation",               (DL_FUNC) &_episimR_episimR_nextreaction_simulation,               4},
    {"_episimR_episimR_nextreaction_simulation_meanfield",     (DL_FUNC) &_episimR_episimR_nextreaction_simulation_meanfield,     5},
    {"_episimR_episimR_nmga_simulation",                       (DL_FUNC) &_episimR_episimR_nmga_simulation,                       4},
    {"_episimR_episimR_simulation_addinfections",              (DL_FUNC) &_episimR_episimR_simulation_addinfections,              3},
    {"_episimR_episimR_simulation_graph",                      (DL_FUNC) &_episimR_episimR_simulation_graph,                      1},
    {"_episimR_episimR_simulation_isinfected",                 (DL_FUNC) &_episimR_episimR_simulation_isinfected,                 2},
    {"_episimR_episimR_simulation_ninfected",                  (DL_FUNC) &_episimR_episimR_simulation_ninfected,                  1},
    {"_episimR_episimR_simulation_options",                    (DL_FUNC) &_episimR_episimR_simulation_options,                    1},
    {"_episimR_episimR_simulation_resettime",                  (DL_FUNC) &_episimR_episimR_simulation_resettime,                  1},
    {"_episimR_episimR_simulation_run",                        (DL_FUNC) &_episimR_episimR_simulation_run,                        3},
    {"_episimR_episimR_simulation_transmissiontime",           (DL_FUNC) &_episimR_episimR_simulation_transmissiontime,           1},
    {"_episimR_episimR_stored_graph",                          (DL_FUNC) &_episimR_episimR_stored_graph,                          1},
    {"_episimR_episimR_time_density",                          (DL_FUNC) &_episimR_episimR_time_density,                          2},
    {"_episimR_episimR_time_hazardrate",                       (DL_FUNC) &_episimR_episimR_time_hazardrate,                       2},
    {"_episimR_episimR_time_sample",                           (DL_FUNC) &_episimR_episimR_time_sample,                           4},
    {"_episimR_episimR_time_survivalprobability",              (DL_FUNC) &_episimR_episimR_time_survivalprobability,              4},
    {"_episimR_episimR_time_survivalquantile",                 (DL_FUNC) &_episimR_episimR_time_survivalquantile,                 4},
    {"_episimR_episimR_userdefined_graph",                     (DL_FUNC) &_episimR_episimR_userdefined_graph,                     1},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_episimR(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
