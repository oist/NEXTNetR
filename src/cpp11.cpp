// Generated by cpp11: do not edit by hand
// clang-format off

#include "NEXTNetR_types.h"
#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// network.cpp
int nextnetR_network_size(const network_R& nw);
extern "C" SEXP _NEXTNetR_nextnetR_network_size(SEXP nw) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_network_size(cpp11::as_cpp<cpp11::decay_t<const network_R&>>(nw)));
  END_CPP11
}
// network.cpp
bool nextnetR_network_is_undirected(const network_R& nw);
extern "C" SEXP _NEXTNetR_nextnetR_network_is_undirected(SEXP nw) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_network_is_undirected(cpp11::as_cpp<cpp11::decay_t<const network_R&>>(nw)));
  END_CPP11
}
// network.cpp
integers nextnetR_network_outdegree(const network_R& nw, integers nodes);
extern "C" SEXP _NEXTNetR_nextnetR_network_outdegree(SEXP nw, SEXP nodes) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_network_outdegree(cpp11::as_cpp<cpp11::decay_t<const network_R&>>(nw), cpp11::as_cpp<cpp11::decay_t<integers>>(nodes)));
  END_CPP11
}
// network.cpp
integers nextnetR_network_neighbour(const network_R& nw, integers nodes, integers indices);
extern "C" SEXP _NEXTNetR_nextnetR_network_neighbour(SEXP nw, SEXP nodes, SEXP indices) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_network_neighbour(cpp11::as_cpp<cpp11::decay_t<const network_R&>>(nw), cpp11::as_cpp<cpp11::decay_t<integers>>(nodes), cpp11::as_cpp<cpp11::decay_t<integers>>(indices)));
  END_CPP11
}
// network.cpp
list nextnetR_network_adjacencylist(const network_R& nw);
extern "C" SEXP _NEXTNetR_nextnetR_network_adjacencylist(SEXP nw) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_network_adjacencylist(cpp11::as_cpp<cpp11::decay_t<const network_R&>>(nw)));
  END_CPP11
}
// network.cpp
list nextnetR_weighted_network_adjacencylist(const network_R& nw);
extern "C" SEXP _NEXTNetR_nextnetR_weighted_network_adjacencylist(SEXP nw) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_weighted_network_adjacencylist(cpp11::as_cpp<cpp11::decay_t<const network_R&>>(nw)));
  END_CPP11
}
// network.cpp
list nextnetR_network_bounds(const network_R& nw);
extern "C" SEXP _NEXTNetR_nextnetR_network_bounds(SEXP nw) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_network_bounds(cpp11::as_cpp<cpp11::decay_t<const network_R&>>(nw)));
  END_CPP11
}
// network.cpp
doubles_matrix<> nextnetR_network_coordinates(const network_R& nw, integers nodes);
extern "C" SEXP _NEXTNetR_nextnetR_network_coordinates(SEXP nw, SEXP nodes) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_network_coordinates(cpp11::as_cpp<cpp11::decay_t<const network_R&>>(nw), cpp11::as_cpp<cpp11::decay_t<integers>>(nodes)));
  END_CPP11
}
// network.cpp
list nextnetR_reproduction_matrix(const network_R& nw);
extern "C" SEXP _NEXTNetR_nextnetR_reproduction_matrix(SEXP nw) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_reproduction_matrix(cpp11::as_cpp<cpp11::decay_t<const network_R&>>(nw)));
  END_CPP11
}
// network.cpp
network_R nextnetR_empirical_network(r_string filename);
extern "C" SEXP _NEXTNetR_nextnetR_empirical_network(SEXP filename) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_empirical_network(cpp11::as_cpp<cpp11::decay_t<r_string>>(filename)));
  END_CPP11
}
// network.cpp
network_R nextnetR_erdos_renyi_network(int size, double avg_degree);
extern "C" SEXP _NEXTNetR_nextnetR_erdos_renyi_network(SEXP size, SEXP avg_degree) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_erdos_renyi_network(cpp11::as_cpp<cpp11::decay_t<int>>(size), cpp11::as_cpp<cpp11::decay_t<double>>(avg_degree)));
  END_CPP11
}
// network.cpp
network_R nextnetR_fully_connected_network(int size);
extern "C" SEXP _NEXTNetR_nextnetR_fully_connected_network(SEXP size) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_fully_connected_network(cpp11::as_cpp<cpp11::decay_t<int>>(size)));
  END_CPP11
}
// network.cpp
network_R nextnetR_acyclic_network(int size, double avg_degree, bool reduced_root_degree);
extern "C" SEXP _NEXTNetR_nextnetR_acyclic_network(SEXP size, SEXP avg_degree, SEXP reduced_root_degree) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_acyclic_network(cpp11::as_cpp<cpp11::decay_t<int>>(size), cpp11::as_cpp<cpp11::decay_t<double>>(avg_degree), cpp11::as_cpp<cpp11::decay_t<bool>>(reduced_root_degree)));
  END_CPP11
}
// network.cpp
network_R nextnetR_configmodel_network(integers degrees);
extern "C" SEXP _NEXTNetR_nextnetR_configmodel_network(SEXP degrees) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_configmodel_network(cpp11::as_cpp<cpp11::decay_t<integers>>(degrees)));
  END_CPP11
}
// network.cpp
network_R nextnetR_configmodel_clustered_alpha_network(integers degrees, double alpha, double beta);
extern "C" SEXP _NEXTNetR_nextnetR_configmodel_clustered_alpha_network(SEXP degrees, SEXP alpha, SEXP beta) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_configmodel_clustered_alpha_network(cpp11::as_cpp<cpp11::decay_t<integers>>(degrees), cpp11::as_cpp<cpp11::decay_t<double>>(alpha), cpp11::as_cpp<cpp11::decay_t<double>>(beta)));
  END_CPP11
}
// network.cpp
network_R nextnetR_configmodel_clustered_ck_network(integers degrees, SEXP ck, double beta);
extern "C" SEXP _NEXTNetR_nextnetR_configmodel_clustered_ck_network(SEXP degrees, SEXP ck, SEXP beta) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_configmodel_clustered_ck_network(cpp11::as_cpp<cpp11::decay_t<integers>>(degrees), cpp11::as_cpp<cpp11::decay_t<SEXP>>(ck), cpp11::as_cpp<cpp11::decay_t<double>>(beta)));
  END_CPP11
}
// network.cpp
network_R nextnetR_configmodel_clustered_triangles_network(integers degrees, integers triangles, double beta);
extern "C" SEXP _NEXTNetR_nextnetR_configmodel_clustered_triangles_network(SEXP degrees, SEXP triangles, SEXP beta) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_configmodel_clustered_triangles_network(cpp11::as_cpp<cpp11::decay_t<integers>>(degrees), cpp11::as_cpp<cpp11::decay_t<integers>>(triangles), cpp11::as_cpp<cpp11::decay_t<double>>(beta)));
  END_CPP11
}
// network.cpp
network_R nextnetR_watts_strogatz_network(int size, int k, double p);
extern "C" SEXP _NEXTNetR_nextnetR_watts_strogatz_network(SEXP size, SEXP k, SEXP p) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_watts_strogatz_network(cpp11::as_cpp<cpp11::decay_t<int>>(size), cpp11::as_cpp<cpp11::decay_t<int>>(k), cpp11::as_cpp<cpp11::decay_t<double>>(p)));
  END_CPP11
}
// network.cpp
network_R nextnetR_barabasialbert_network(int size, int m);
extern "C" SEXP _NEXTNetR_nextnetR_barabasialbert_network(SEXP size, SEXP m) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_barabasialbert_network(cpp11::as_cpp<cpp11::decay_t<int>>(size), cpp11::as_cpp<cpp11::decay_t<int>>(m)));
  END_CPP11
}
// network.cpp
network_R nextnetR_cubiclattice2d_network(int edge_length);
extern "C" SEXP _NEXTNetR_nextnetR_cubiclattice2d_network(SEXP edge_length) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_cubiclattice2d_network(cpp11::as_cpp<cpp11::decay_t<int>>(edge_length)));
  END_CPP11
}
// network.cpp
network_R nextnetR_cubiclattice3d_network(int edge_length);
extern "C" SEXP _NEXTNetR_nextnetR_cubiclattice3d_network(SEXP edge_length) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_cubiclattice3d_network(cpp11::as_cpp<cpp11::decay_t<int>>(edge_length)));
  END_CPP11
}
// network.cpp
network_R nextnetR_cubiclattice4d_network(int edge_length);
extern "C" SEXP _NEXTNetR_nextnetR_cubiclattice4d_network(SEXP edge_length) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_cubiclattice4d_network(cpp11::as_cpp<cpp11::decay_t<int>>(edge_length)));
  END_CPP11
}
// network.cpp
network_R nextnetR_cubiclattice5d_network(int edge_length);
extern "C" SEXP _NEXTNetR_nextnetR_cubiclattice5d_network(SEXP edge_length) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_cubiclattice5d_network(cpp11::as_cpp<cpp11::decay_t<int>>(edge_length)));
  END_CPP11
}
// network.cpp
network_R nextnetR_cubiclattice6d_network(int edge_length);
extern "C" SEXP _NEXTNetR_nextnetR_cubiclattice6d_network(SEXP edge_length) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_cubiclattice6d_network(cpp11::as_cpp<cpp11::decay_t<int>>(edge_length)));
  END_CPP11
}
// network.cpp
network_R nextnetR_cubiclattice7d_network(int edge_length);
extern "C" SEXP _NEXTNetR_nextnetR_cubiclattice7d_network(SEXP edge_length) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_cubiclattice7d_network(cpp11::as_cpp<cpp11::decay_t<int>>(edge_length)));
  END_CPP11
}
// network.cpp
network_R nextnetR_cubiclattice8d_network(int edge_length);
extern "C" SEXP _NEXTNetR_nextnetR_cubiclattice8d_network(SEXP edge_length) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_cubiclattice8d_network(cpp11::as_cpp<cpp11::decay_t<int>>(edge_length)));
  END_CPP11
}
// network.cpp
network_R nextnetR_brownian_proximity_temporalnetwork(int size, double avg_degree, double radius, double D0, double D1, double gamma, SEXP dt);
extern "C" SEXP _NEXTNetR_nextnetR_brownian_proximity_temporalnetwork(SEXP size, SEXP avg_degree, SEXP radius, SEXP D0, SEXP D1, SEXP gamma, SEXP dt) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_brownian_proximity_temporalnetwork(cpp11::as_cpp<cpp11::decay_t<int>>(size), cpp11::as_cpp<cpp11::decay_t<double>>(avg_degree), cpp11::as_cpp<cpp11::decay_t<double>>(radius), cpp11::as_cpp<cpp11::decay_t<double>>(D0), cpp11::as_cpp<cpp11::decay_t<double>>(D1), cpp11::as_cpp<cpp11::decay_t<double>>(gamma), cpp11::as_cpp<cpp11::decay_t<SEXP>>(dt)));
  END_CPP11
}
// network.cpp
network_R nextnetR_empirical_temporalnetwork(std::string file, bool finite_duration, double dt);
extern "C" SEXP _NEXTNetR_nextnetR_empirical_temporalnetwork(SEXP file, SEXP finite_duration, SEXP dt) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_empirical_temporalnetwork(cpp11::as_cpp<cpp11::decay_t<std::string>>(file), cpp11::as_cpp<cpp11::decay_t<bool>>(finite_duration), cpp11::as_cpp<cpp11::decay_t<double>>(dt)));
  END_CPP11
}
// network.cpp
network_R nextnetR_sirx_temporalnetwork(const network_R& nw, double kappa0, double kappa);
extern "C" SEXP _NEXTNetR_nextnetR_sirx_temporalnetwork(SEXP nw, SEXP kappa0, SEXP kappa) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_sirx_temporalnetwork(cpp11::as_cpp<cpp11::decay_t<const network_R&>>(nw), cpp11::as_cpp<cpp11::decay_t<double>>(kappa0), cpp11::as_cpp<cpp11::decay_t<double>>(kappa)));
  END_CPP11
}
// network.cpp
network_R nextnetR_adjacencylist_network(list input_al, bool is_undirected);
extern "C" SEXP _NEXTNetR_nextnetR_adjacencylist_network(SEXP input_al, SEXP is_undirected) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_adjacencylist_network(cpp11::as_cpp<cpp11::decay_t<list>>(input_al), cpp11::as_cpp<cpp11::decay_t<bool>>(is_undirected)));
  END_CPP11
}
// network.cpp
network_R nextnetR_weighted_adjacencylist_network(list input_al, bool is_undirected);
extern "C" SEXP _NEXTNetR_nextnetR_weighted_adjacencylist_network(SEXP input_al, SEXP is_undirected) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_weighted_adjacencylist_network(cpp11::as_cpp<cpp11::decay_t<list>>(input_al), cpp11::as_cpp<cpp11::decay_t<bool>>(is_undirected)));
  END_CPP11
}
// simulation.cpp
simulation_R nextnetR_nextreaction_simulation(network_R nw, transmission_time_R psi, sexp rho_, list opts);
extern "C" SEXP _NEXTNetR_nextnetR_nextreaction_simulation(SEXP nw, SEXP psi, SEXP rho_, SEXP opts) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_nextreaction_simulation(cpp11::as_cpp<cpp11::decay_t<network_R>>(nw), cpp11::as_cpp<cpp11::decay_t<transmission_time_R>>(psi), cpp11::as_cpp<cpp11::decay_t<sexp>>(rho_), cpp11::as_cpp<cpp11::decay_t<list>>(opts)));
  END_CPP11
}
// simulation.cpp
simulation_R nextnetR_nextreaction_simulation_meanfield(int N, double R0, transmission_time_R psi, sexp rho_, list opts);
extern "C" SEXP _NEXTNetR_nextnetR_nextreaction_simulation_meanfield(SEXP N, SEXP R0, SEXP psi, SEXP rho_, SEXP opts) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_nextreaction_simulation_meanfield(cpp11::as_cpp<cpp11::decay_t<int>>(N), cpp11::as_cpp<cpp11::decay_t<double>>(R0), cpp11::as_cpp<cpp11::decay_t<transmission_time_R>>(psi), cpp11::as_cpp<cpp11::decay_t<sexp>>(rho_), cpp11::as_cpp<cpp11::decay_t<list>>(opts)));
  END_CPP11
}
// simulation.cpp
simulation_R nextnetR_nmga_simulation(network_R nw, transmission_time_R psi, sexp rho_, list opts);
extern "C" SEXP _NEXTNetR_nextnetR_nmga_simulation(SEXP nw, SEXP psi, SEXP rho_, SEXP opts) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_nmga_simulation(cpp11::as_cpp<cpp11::decay_t<network_R>>(nw), cpp11::as_cpp<cpp11::decay_t<transmission_time_R>>(psi), cpp11::as_cpp<cpp11::decay_t<sexp>>(rho_), cpp11::as_cpp<cpp11::decay_t<list>>(opts)));
  END_CPP11
}
// simulation.cpp
transmission_time_R nextnetR_simulation_transmissiontime(const simulation_R& sim);
extern "C" SEXP _NEXTNetR_nextnetR_simulation_transmissiontime(SEXP sim) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_simulation_transmissiontime(cpp11::as_cpp<cpp11::decay_t<const simulation_R&>>(sim)));
  END_CPP11
}
// simulation.cpp
SEXP nextnetR_simulation_resettime(const simulation_R& sim);
extern "C" SEXP _NEXTNetR_nextnetR_simulation_resettime(SEXP sim) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_simulation_resettime(cpp11::as_cpp<cpp11::decay_t<const simulation_R&>>(sim)));
  END_CPP11
}
// simulation.cpp
network_R nextnetR_simulation_network(const simulation_R& sim);
extern "C" SEXP _NEXTNetR_nextnetR_simulation_network(SEXP sim) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_simulation_network(cpp11::as_cpp<cpp11::decay_t<const simulation_R&>>(sim)));
  END_CPP11
}
// simulation.cpp
list nextnetR_simulation_options(const simulation_R& sim);
extern "C" SEXP _NEXTNetR_nextnetR_simulation_options(SEXP sim) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_simulation_options(cpp11::as_cpp<cpp11::decay_t<const simulation_R&>>(sim)));
  END_CPP11
}
// simulation.cpp
list nextnetR_simulation_ninfected(const simulation_R& sim);
extern "C" SEXP _NEXTNetR_nextnetR_simulation_ninfected(SEXP sim) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_simulation_ninfected(cpp11::as_cpp<cpp11::decay_t<const simulation_R&>>(sim)));
  END_CPP11
}
// simulation.cpp
logicals nextnetR_simulation_isinfected(const simulation_R& sim, integers nodes);
extern "C" SEXP _NEXTNetR_nextnetR_simulation_isinfected(SEXP sim, SEXP nodes) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_simulation_isinfected(cpp11::as_cpp<cpp11::decay_t<const simulation_R&>>(sim), cpp11::as_cpp<cpp11::decay_t<integers>>(nodes)));
  END_CPP11
}
// simulation.cpp
void nextnetR_simulation_addinfections(const simulation_R& sim, integers nodes, doubles times);
extern "C" SEXP _NEXTNetR_nextnetR_simulation_addinfections(SEXP sim, SEXP nodes, SEXP times) {
  BEGIN_CPP11
    nextnetR_simulation_addinfections(cpp11::as_cpp<cpp11::decay_t<const simulation_R&>>(sim), cpp11::as_cpp<cpp11::decay_t<integers>>(nodes), cpp11::as_cpp<cpp11::decay_t<doubles>>(times));
    return R_NilValue;
  END_CPP11
}
// simulation.cpp
data_frame nextnetR_simulation_run(const simulation_R& sim_, list stop, list opts);
extern "C" SEXP _NEXTNetR_nextnetR_simulation_run(SEXP sim_, SEXP stop, SEXP opts) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_simulation_run(cpp11::as_cpp<cpp11::decay_t<const simulation_R&>>(sim_), cpp11::as_cpp<cpp11::decay_t<list>>(stop), cpp11::as_cpp<cpp11::decay_t<list>>(opts)));
  END_CPP11
}
// transmission_time.cpp
doubles nextnetR_time_sample(int n, const transmission_time_R& ttr, interval_t t, int m);
extern "C" SEXP _NEXTNetR_nextnetR_time_sample(SEXP n, SEXP ttr, SEXP t, SEXP m) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_time_sample(cpp11::as_cpp<cpp11::decay_t<int>>(n), cpp11::as_cpp<cpp11::decay_t<const transmission_time_R&>>(ttr), cpp11::as_cpp<cpp11::decay_t<interval_t>>(t), cpp11::as_cpp<cpp11::decay_t<int>>(m)));
  END_CPP11
}
// transmission_time.cpp
doubles nextnetR_time_density(const transmission_time_R& ttr, doubles taus);
extern "C" SEXP _NEXTNetR_nextnetR_time_density(SEXP ttr, SEXP taus) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_time_density(cpp11::as_cpp<cpp11::decay_t<const transmission_time_R&>>(ttr), cpp11::as_cpp<cpp11::decay_t<doubles>>(taus)));
  END_CPP11
}
// transmission_time.cpp
doubles nextnetR_time_hazardrate(const transmission_time_R& ttr, doubles taus);
extern "C" SEXP _NEXTNetR_nextnetR_time_hazardrate(SEXP ttr, SEXP taus) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_time_hazardrate(cpp11::as_cpp<cpp11::decay_t<const transmission_time_R&>>(ttr), cpp11::as_cpp<cpp11::decay_t<doubles>>(taus)));
  END_CPP11
}
// transmission_time.cpp
doubles nextnetR_time_survivalprobability(const transmission_time_R& ttr, doubles taus, doubles ts, integers ms);
extern "C" SEXP _NEXTNetR_nextnetR_time_survivalprobability(SEXP ttr, SEXP taus, SEXP ts, SEXP ms) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_time_survivalprobability(cpp11::as_cpp<cpp11::decay_t<const transmission_time_R&>>(ttr), cpp11::as_cpp<cpp11::decay_t<doubles>>(taus), cpp11::as_cpp<cpp11::decay_t<doubles>>(ts), cpp11::as_cpp<cpp11::decay_t<integers>>(ms)));
  END_CPP11
}
// transmission_time.cpp
doubles nextnetR_time_survivalquantile(const transmission_time_R& ttr, doubles ps, doubles ts, integers ms);
extern "C" SEXP _NEXTNetR_nextnetR_time_survivalquantile(SEXP ttr, SEXP ps, SEXP ts, SEXP ms) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_time_survivalquantile(cpp11::as_cpp<cpp11::decay_t<const transmission_time_R&>>(ttr), cpp11::as_cpp<cpp11::decay_t<doubles>>(ps), cpp11::as_cpp<cpp11::decay_t<doubles>>(ts), cpp11::as_cpp<cpp11::decay_t<integers>>(ms)));
  END_CPP11
}
// transmission_time.cpp
transmission_time_R nextnetR_exponential_time(double lambda);
extern "C" SEXP _NEXTNetR_nextnetR_exponential_time(SEXP lambda) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_exponential_time(cpp11::as_cpp<cpp11::decay_t<double>>(lambda)));
  END_CPP11
}
// transmission_time.cpp
transmission_time_R nextnetR_lognormal_time(double mean, double var, double pinf);
extern "C" SEXP _NEXTNetR_nextnetR_lognormal_time(SEXP mean, SEXP var, SEXP pinf) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_lognormal_time(cpp11::as_cpp<cpp11::decay_t<double>>(mean), cpp11::as_cpp<cpp11::decay_t<double>>(var), cpp11::as_cpp<cpp11::decay_t<double>>(pinf)));
  END_CPP11
}
// transmission_time.cpp
transmission_time_R nextnetR_gamma_time(double mean, double var, double pinf);
extern "C" SEXP _NEXTNetR_nextnetR_gamma_time(SEXP mean, SEXP var, SEXP pinf) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_gamma_time(cpp11::as_cpp<cpp11::decay_t<double>>(mean), cpp11::as_cpp<cpp11::decay_t<double>>(var), cpp11::as_cpp<cpp11::decay_t<double>>(pinf)));
  END_CPP11
}
// transmission_time.cpp
transmission_time_R nextnetR_generic_time(SEXP density, SEXP survivalprobability, bool probability_is_trinary, SEXP survivalquantile, bool quantile_is_trinary, SEXP sample, double pinfinity);
extern "C" SEXP _NEXTNetR_nextnetR_generic_time(SEXP density, SEXP survivalprobability, SEXP probability_is_trinary, SEXP survivalquantile, SEXP quantile_is_trinary, SEXP sample, SEXP pinfinity) {
  BEGIN_CPP11
    return cpp11::as_sexp(nextnetR_generic_time(cpp11::as_cpp<cpp11::decay_t<SEXP>>(density), cpp11::as_cpp<cpp11::decay_t<SEXP>>(survivalprobability), cpp11::as_cpp<cpp11::decay_t<bool>>(probability_is_trinary), cpp11::as_cpp<cpp11::decay_t<SEXP>>(survivalquantile), cpp11::as_cpp<cpp11::decay_t<bool>>(quantile_is_trinary), cpp11::as_cpp<cpp11::decay_t<SEXP>>(sample), cpp11::as_cpp<cpp11::decay_t<double>>(pinfinity)));
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_NEXTNetR_nextnetR_acyclic_network",                         (DL_FUNC) &_NEXTNetR_nextnetR_acyclic_network,                         3},
    {"_NEXTNetR_nextnetR_adjacencylist_network",                   (DL_FUNC) &_NEXTNetR_nextnetR_adjacencylist_network,                   2},
    {"_NEXTNetR_nextnetR_barabasialbert_network",                  (DL_FUNC) &_NEXTNetR_nextnetR_barabasialbert_network,                  2},
    {"_NEXTNetR_nextnetR_brownian_proximity_temporalnetwork",      (DL_FUNC) &_NEXTNetR_nextnetR_brownian_proximity_temporalnetwork,      7},
    {"_NEXTNetR_nextnetR_configmodel_clustered_alpha_network",     (DL_FUNC) &_NEXTNetR_nextnetR_configmodel_clustered_alpha_network,     3},
    {"_NEXTNetR_nextnetR_configmodel_clustered_ck_network",        (DL_FUNC) &_NEXTNetR_nextnetR_configmodel_clustered_ck_network,        3},
    {"_NEXTNetR_nextnetR_configmodel_clustered_triangles_network", (DL_FUNC) &_NEXTNetR_nextnetR_configmodel_clustered_triangles_network, 3},
    {"_NEXTNetR_nextnetR_configmodel_network",                     (DL_FUNC) &_NEXTNetR_nextnetR_configmodel_network,                     1},
    {"_NEXTNetR_nextnetR_cubiclattice2d_network",                  (DL_FUNC) &_NEXTNetR_nextnetR_cubiclattice2d_network,                  1},
    {"_NEXTNetR_nextnetR_cubiclattice3d_network",                  (DL_FUNC) &_NEXTNetR_nextnetR_cubiclattice3d_network,                  1},
    {"_NEXTNetR_nextnetR_cubiclattice4d_network",                  (DL_FUNC) &_NEXTNetR_nextnetR_cubiclattice4d_network,                  1},
    {"_NEXTNetR_nextnetR_cubiclattice5d_network",                  (DL_FUNC) &_NEXTNetR_nextnetR_cubiclattice5d_network,                  1},
    {"_NEXTNetR_nextnetR_cubiclattice6d_network",                  (DL_FUNC) &_NEXTNetR_nextnetR_cubiclattice6d_network,                  1},
    {"_NEXTNetR_nextnetR_cubiclattice7d_network",                  (DL_FUNC) &_NEXTNetR_nextnetR_cubiclattice7d_network,                  1},
    {"_NEXTNetR_nextnetR_cubiclattice8d_network",                  (DL_FUNC) &_NEXTNetR_nextnetR_cubiclattice8d_network,                  1},
    {"_NEXTNetR_nextnetR_empirical_network",                       (DL_FUNC) &_NEXTNetR_nextnetR_empirical_network,                       1},
    {"_NEXTNetR_nextnetR_empirical_temporalnetwork",               (DL_FUNC) &_NEXTNetR_nextnetR_empirical_temporalnetwork,               3},
    {"_NEXTNetR_nextnetR_erdos_renyi_network",                     (DL_FUNC) &_NEXTNetR_nextnetR_erdos_renyi_network,                     2},
    {"_NEXTNetR_nextnetR_exponential_time",                        (DL_FUNC) &_NEXTNetR_nextnetR_exponential_time,                        1},
    {"_NEXTNetR_nextnetR_fully_connected_network",                 (DL_FUNC) &_NEXTNetR_nextnetR_fully_connected_network,                 1},
    {"_NEXTNetR_nextnetR_gamma_time",                              (DL_FUNC) &_NEXTNetR_nextnetR_gamma_time,                              3},
    {"_NEXTNetR_nextnetR_generic_time",                            (DL_FUNC) &_NEXTNetR_nextnetR_generic_time,                            7},
    {"_NEXTNetR_nextnetR_lognormal_time",                          (DL_FUNC) &_NEXTNetR_nextnetR_lognormal_time,                          3},
    {"_NEXTNetR_nextnetR_network_adjacencylist",                   (DL_FUNC) &_NEXTNetR_nextnetR_network_adjacencylist,                   1},
    {"_NEXTNetR_nextnetR_network_bounds",                          (DL_FUNC) &_NEXTNetR_nextnetR_network_bounds,                          1},
    {"_NEXTNetR_nextnetR_network_coordinates",                     (DL_FUNC) &_NEXTNetR_nextnetR_network_coordinates,                     2},
    {"_NEXTNetR_nextnetR_network_is_undirected",                   (DL_FUNC) &_NEXTNetR_nextnetR_network_is_undirected,                   1},
    {"_NEXTNetR_nextnetR_network_neighbour",                       (DL_FUNC) &_NEXTNetR_nextnetR_network_neighbour,                       3},
    {"_NEXTNetR_nextnetR_network_outdegree",                       (DL_FUNC) &_NEXTNetR_nextnetR_network_outdegree,                       2},
    {"_NEXTNetR_nextnetR_network_size",                            (DL_FUNC) &_NEXTNetR_nextnetR_network_size,                            1},
    {"_NEXTNetR_nextnetR_nextreaction_simulation",                 (DL_FUNC) &_NEXTNetR_nextnetR_nextreaction_simulation,                 4},
    {"_NEXTNetR_nextnetR_nextreaction_simulation_meanfield",       (DL_FUNC) &_NEXTNetR_nextnetR_nextreaction_simulation_meanfield,       5},
    {"_NEXTNetR_nextnetR_nmga_simulation",                         (DL_FUNC) &_NEXTNetR_nextnetR_nmga_simulation,                         4},
    {"_NEXTNetR_nextnetR_reproduction_matrix",                     (DL_FUNC) &_NEXTNetR_nextnetR_reproduction_matrix,                     1},
    {"_NEXTNetR_nextnetR_simulation_addinfections",                (DL_FUNC) &_NEXTNetR_nextnetR_simulation_addinfections,                3},
    {"_NEXTNetR_nextnetR_simulation_isinfected",                   (DL_FUNC) &_NEXTNetR_nextnetR_simulation_isinfected,                   2},
    {"_NEXTNetR_nextnetR_simulation_network",                      (DL_FUNC) &_NEXTNetR_nextnetR_simulation_network,                      1},
    {"_NEXTNetR_nextnetR_simulation_ninfected",                    (DL_FUNC) &_NEXTNetR_nextnetR_simulation_ninfected,                    1},
    {"_NEXTNetR_nextnetR_simulation_options",                      (DL_FUNC) &_NEXTNetR_nextnetR_simulation_options,                      1},
    {"_NEXTNetR_nextnetR_simulation_resettime",                    (DL_FUNC) &_NEXTNetR_nextnetR_simulation_resettime,                    1},
    {"_NEXTNetR_nextnetR_simulation_run",                          (DL_FUNC) &_NEXTNetR_nextnetR_simulation_run,                          3},
    {"_NEXTNetR_nextnetR_simulation_transmissiontime",             (DL_FUNC) &_NEXTNetR_nextnetR_simulation_transmissiontime,             1},
    {"_NEXTNetR_nextnetR_sirx_temporalnetwork",                    (DL_FUNC) &_NEXTNetR_nextnetR_sirx_temporalnetwork,                    3},
    {"_NEXTNetR_nextnetR_time_density",                            (DL_FUNC) &_NEXTNetR_nextnetR_time_density,                            2},
    {"_NEXTNetR_nextnetR_time_hazardrate",                         (DL_FUNC) &_NEXTNetR_nextnetR_time_hazardrate,                         2},
    {"_NEXTNetR_nextnetR_time_sample",                             (DL_FUNC) &_NEXTNetR_nextnetR_time_sample,                             4},
    {"_NEXTNetR_nextnetR_time_survivalprobability",                (DL_FUNC) &_NEXTNetR_nextnetR_time_survivalprobability,                4},
    {"_NEXTNetR_nextnetR_time_survivalquantile",                   (DL_FUNC) &_NEXTNetR_nextnetR_time_survivalquantile,                   4},
    {"_NEXTNetR_nextnetR_watts_strogatz_network",                  (DL_FUNC) &_NEXTNetR_nextnetR_watts_strogatz_network,                  3},
    {"_NEXTNetR_nextnetR_weighted_adjacencylist_network",          (DL_FUNC) &_NEXTNetR_nextnetR_weighted_adjacencylist_network,          2},
    {"_NEXTNetR_nextnetR_weighted_network_adjacencylist",          (DL_FUNC) &_NEXTNetR_nextnetR_weighted_network_adjacencylist,          1},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_NEXTNetR(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
