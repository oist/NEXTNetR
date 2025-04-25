#include <utility>
#include <limits>

#include <cpp11.hpp>
#include <cpp11/function.hpp>
#include <cpp11/external_pointer.hpp>

#include "NEXTNetR_types.h"
#include "rng.h"

#include "nextnet/types.h"
#include "nextnet/network.h"
#include "nextnet/temporal_network.h"
#include "nextnet/weighted_network.h"
#include "nextnet/brownian_proximity_network.h"
#include "nextnet/pstream/pstream.h"

using namespace nextnetR;

using namespace cpp11;
namespace writable = cpp11::writable;

[[cpp11::register]]
int nextnetR_network_size(const network_R& nw) {
    if (!nw) stop("network cannot be NULL"); 
    
    /* Must enter RNG scope since networks may generate their topology on the fly */
    RNG_SCOPE_IF_NECESSARY;
    
    return nw->nodes();
}

[[cpp11::register]]
bool nextnetR_network_is_undirected(const network_R& nw) {
    if (!nw) stop("network cannot be NULL"); 
    
    /* Must enter RNG scope since networks may generate their topology on the fly */
    RNG_SCOPE_IF_NECESSARY;
    
    return nw->is_undirected();
}

[[cpp11::register]]
bool nextnetR_network_is_simple(const network_R& nw) {
  if (!nw) stop("network cannot be NULL"); 
  
  /* Must enter RNG scope since networks may generate their topology on the fly */
  RNG_SCOPE_IF_NECESSARY;
  
  return nw->is_simple();
}

[[cpp11::register]]
bool nextnetR_network_is_weighted(const network_R& nw) {
  if (!nw) stop("network cannot be NULL"); 
  
  return (as_weighted_network(nw.get()) != nullptr);
}

[[cpp11::register]]
bool nextnetR_network_is_temporal(const network_R& nw) {
  if (!nw) stop("network cannot be NULL"); 
  
  return (dynamic_cast<temporal_network*>(nw.get()) != nullptr);
}

[[cpp11::register]]
integers nextnetR_network_outdegree(const network_R& nw, integers nodes) {
    if (!nw) stop("network cannot be NULL"); 
    
    /* Must enter RNG scope since networks may generate their topology on the fly */
    RNG_SCOPE_IF_NECESSARY;
    
    const std::size_t l = nodes.size();
    writable::integers r;
    r.reserve(l);
    
    const node_t N = nw->nodes();
    for(std::size_t j = 0; j < l; ++j) {
        const int n = nodes[j];
        r.push_back(((n >= 1) && (n <= (node_t)N)) ? nw->outdegree(n - 1) : NA_INTEGER);
    }
    return r;
}

[[cpp11::register]]
integers nextnetR_network_neighbour(const network_R& nw, integers nodes, integers indices) {
    if (!nw) stop("network cannot be NULL"); 

    /* Must enter RNG scope since networks may generate their topology on the fly */
    RNG_SCOPE_IF_NECESSARY;
    
    if (nodes.size() != indices.size())
        stop("number of nodes and number of indices must agree");
    
    /* Create output */
    const std::size_t l = nodes.size();
    writable::integers r;
    r.reserve(l);
    
    /* Fill */
    const node_t N = nw->nodes();
    for(std::size_t j = 0; j < l; ++j) {
        const node_t n = nodes[j];
        if ((n < 1) || (n > (node_t)N)) {
            r.push_back(NA_INTEGER);
            continue;
        }
        
        const int i = indices[j];
        const int m = nw->neighbour(n - 1, i - 1);
        r.push_back((m >= 0) ? (m + 1) : NA_INTEGER);
    }
    return r;
}

[[cpp11::register]]
list nextnetR_network_neighbour_weight(const network_R& nw, integers nodes, integers indices) {
    if (!nw) stop("network cannot be NULL"); 

    weighted_network* const wnw = as_weighted_network(nw.get());
    if (wnw == nullptr)
        stop("network is not weighted");

    /* Must enter RNG scope since networks may generate their topology on the fly */
    RNG_SCOPE_IF_NECESSARY;
    
    if (nodes.size() != indices.size())
        stop("number of nodes and number of indices must agree");
    
    /* Create outputs */
    const std::size_t l = nodes.size();
    writable::integers ns;
    writable::doubles ws;
    ns.reserve(l);
    ws.reserve(l);
    
    /* Fill */
    const node_t N = wnw->nodes();
    for(std::size_t j = 0; j < l; ++j) {
        const node_t n = nodes[j];
        if ((n < 1) || (n > (node_t)N)) {
            ns.push_back(NA_INTEGER);
            continue;
        }
        
        const int i = indices[j];
        double w = NA_REAL;
        const int m = wnw->neighbour(n - 1, i - 1, &w);
        ns.push_back((m >= 0) ? (m + 1) : NA_INTEGER);
        ws.push_back((m >= 0) ? w : NA_REAL);
    }

    return writable::list({
        "n"_nm = ns,
        "w"_nm = ws
    });
}

[[cpp11::register]]
list nextnetR_network_adjacencylist(const network_R& nw, bool above_diagonal) {
    if (!nw) stop("network cannot be NULL"); 

    const bool is_undirected = nw->is_undirected();
    if (!is_undirected && above_diagonal)
      stop("above_diagonal=TRUE is invalid for directed networks");
    
    /* Must enter RNG scope since networks may generate their topology on the fly */
    RNG_SCOPE_IF_NECESSARY;

    /* Allocate output, a vector of node labels and a list of neighbour vectors */
    const int l = nw->nodes();
    writable::list al;
    al.reserve(l);

    /* Interate over nodes and generate output */
    for(int n = 0; n < l; ++n) {
        const int d = nw->outdegree(n);
        writable::integers node_al;
        node_al.reserve(d);
        for(int i=0; i < d; ++i) {
            const int m = nw->neighbour(n, i);
            if (is_undirected && above_diagonal && (m <n))
              continue;
            node_al.push_back((m >= 0) ? (m + 1) : NA_INTEGER);
        }
        
        /* Append to output */
        al.push_back(node_al);
    }

    return al;
}

[[cpp11::register]]
list nextnetR_weighted_network_adjacencylist(const network_R& nw, bool above_diagonal) {
    if (!nw) stop("network cannot be NULL"); 

    weighted_network* const nw_weighted = as_weighted_network(nw.get());
    if (nw_weighted == nullptr)
        stop("network is not weighted");

    const bool is_undirected = nw->is_undirected();
    if (!is_undirected && above_diagonal)
      stop("above_diagonal=TRUE is invalid for directed networks");
    
    /* Must enter RNG scope since networks may generate their topology on the fly */
    RNG_SCOPE_IF_NECESSARY;

    /* Allocate output, a vector of node labels and a list of neighbour vectors */
    const int l = nw_weighted->nodes();
    writable::list al;
    al.reserve(l);

    /* Iterate over nodes and generate output */
    for(int n = 0; n < l; ++n) {
        const int d = nw_weighted->outdegree(n);
        writable::integers node_al;
        writable::doubles node_weights;
        node_al.reserve(d);
        for(int i=0; i < d; ++i) {
            double w = NAN;
            const int m = nw_weighted->neighbour(n, i, &w);
            if (is_undirected && above_diagonal && (m <n))
              continue;
            node_al.push_back((m >= 0) ? (m + 1) : NA_INTEGER);
            node_weights.push_back((m >= 0) ? w : NA_INTEGER);
        }
        
        /* Append to output */
        al.push_back(writable::list {
            "n"_nm = node_al,
            "w"_nm = node_weights
        });
    }

    return al;
}


[[cpp11::register]]
list nextnetR_network_bounds(const network_R& nw) {
    if (!nw) stop("network cannot be NULL"); 

    RNG_SCOPE_IF_NECESSARY;
    
    network_embedding* const ebd = dynamic_cast<network_embedding*>(nw.get());
    if (ebd == nullptr)
        stop("network is not embedded into R^d");
    
    std::vector<double> a, b;
    ebd->bounds(a, b);

    return writable::list({
        writable::doubles(a.begin(), a.end()),
        writable::doubles(b.begin(), b.end())
    });
}

[[cpp11::register]]
doubles_matrix<> nextnetR_network_coordinates(const network_R& nw, integers nodes) {
    if (!nw) stop("network cannot be NULL"); 

    RNG_SCOPE_IF_NECESSARY;
    
    network_embedding* const ebd = dynamic_cast<network_embedding*>(nw.get());
    if (ebd == nullptr)
        stop("network is not embedded into R^d");
    
    const std::size_t d = ebd->dimensionality();
    writable::doubles_matrix<> results(nodes.size(), d);
    std::vector<double> v(d, 0.0);
    for(std::size_t i=0, l=nodes.size(); i < l; ++i) {
        if (ebd->coordinates(nodes[i]-1, v)) {
            for(std::size_t j=0; j < d; ++j)
                results(i, j) = v[j];
        }
    } 
    return results;
}

[[cpp11::register]]
list nextnetR_reproduction_matrix(const network_R& nw) {
    if (!nw) stop("network cannot be NULL"); 

    RNG_SCOPE_IF_NECESSARY;


    double r = NAN, c = NAN, k1 = NAN, k2 = NAN, k3 = NAN, m_bar = NAN, R0 = NAN, R_r = NAN, R_pert = NAN;
    std::vector<std::vector<double>> Mkk = reproduction_matrix(*nw.get(), 3, &r, &c, &k1, &k2, &k3,
                                                               &m_bar, &R0, &R_r, &R_pert);

    // copy result into an R matrix
    writable::doubles_matrix<> M(Mkk.size(), Mkk.size());
    for(std::size_t i=0; i < Mkk.size(); ++i) {
        for(std::size_t j=0; j < Mkk[i].size(); ++j)
            M(i, j) = Mkk[i][j];
    }


    return writable::list({
        "M"_nm = M,
        "r"_nm = r,
        "c"_nm = c,
        "k1"_nm = k1,
        "k2"_nm = k2,
        "k3"_nm = k3,
        "m_bar"_nm = m_bar,
        "R0"_nm = R0,
        "R_r"_nm = R_r,
        "R_pert"_nm = R_pert
    });
}

[[cpp11::register]]
network_R nextnetR_empirical_network(
  strings path, bool undirected, bool simplify, node_t idxbase, strings sep, bool gzip
) {
    if (path.size() != 1)
        stop("expected a single path");
    if (sep.size() != 1)
        stop("expected a single separator");
    
    std::string path_ = (std::string)(path[0]);
    std::string sep_ = (std::string)(sep[0]);
    if (sep_.size() != 1)
        stop("expected a single character as separator");

    redi::ipstream gzfile;
    std::ifstream plainfile;
    std::istream& file = gzip ? (std::istream&)gzfile : (std::istream&)plainfile;
    
    if (gzip)
        gzfile.open("gzip",  std::vector<std::string> { "gzip", "-cd", path_ },
                    redi::pstreambuf::pstdout);
    else
        plainfile.open(path_);
    
    network_R nw = new empirical_network(file, undirected, simplify, idxbase, sep_[0]);
    if (!file.eof())
        stop("failed to read %s", path_.c_str());
    return nw;
}

[[cpp11::register]]
network_R nextnetR_empirical_weightednetwork(
  strings path, bool undirected, bool simplify, node_t idxbase, strings csep, strings wsep,
  bool gzip
) {
    if (path.size() != 1)
        stop("expected a single path");
    if (csep.size() != 1)
        stop("expected a single column separator");
    if (wsep.size() != 1)
        stop("expected a single weight separator");
    
    std::string path_ = (std::string)(path[0]);
    std::string csep_ = (std::string)(csep[0]);
    if (csep_.size() != 1)
        stop("expected a single character as column separator");
    std::string wsep_ = (std::string)(wsep[0]);
    if (wsep_.size() != 1)
        stop("expected a single character as weight separator");

    redi::ipstream gzfile;
    std::ifstream plainfile;
    std::istream& file = gzip ? (std::istream&)gzfile : (std::istream&)plainfile;
    
    if (gzip)
        gzfile.open("gzcat",  std::vector<std::string> { "gzcat", path_ },
                    redi::pstreambuf::pstdout);
    else
        plainfile.open(path_);
    
    network_R nw = new weighted_empirical_network(
      file, undirected, simplify, idxbase, csep_[0], wsep_[0]);
    if (!file.eof())
        stop("failed to read %s", path_.c_str());
    return nw;
}

[[cpp11::register]]
network_R nextnetR_erdos_renyi_network(int size, double avg_degree) {
    RNG_SCOPE_IF_NECESSARY;
    return new erdos_renyi(size, avg_degree, rng_engine());
}

[[cpp11::register]]
network_R nextnetR_fully_connected_network(int size) {
    RNG_SCOPE_IF_NECESSARY;
    return new fully_connected(size, rng_engine());
}

[[cpp11::register]]
network_R nextnetR_acyclic_network(int size, double avg_degree, bool reduced_root_degree) {
    RNG_SCOPE_IF_NECESSARY;
    return new acyclic(avg_degree, reduced_root_degree, rng_engine());
}

[[cpp11::register]]
network_R nextnetR_configmodel_network(integers degrees) {
    RNG_SCOPE_IF_NECESSARY;
    return new config_model(std::vector<int>(degrees.begin(), degrees.end()), rng_engine());
}

[[cpp11::register]]
network_R nextnetR_configmodel_clustered_alpha_network(integers degrees, double alpha, double beta) {
    RNG_SCOPE_IF_NECESSARY;
    return new config_model_clustered_serrano(std::vector<int>(degrees.begin(), degrees.end()),
                                              alpha, beta, rng_engine());
}

[[cpp11::register]]
network_R nextnetR_configmodel_clustered_ck_network(integers degrees, SEXP ck, double beta) {
    RNG_SCOPE_IF_NECESSARY;
    function ck_rf = ck;
    auto ck_lambda = [ck_rf] (int k) { return as_cpp<double>(ck_rf(k)); };
    return new config_model_clustered_serrano(std::vector<int>(degrees.begin(), degrees.end()),
                                              ck_lambda, beta, rng_engine());
}

[[cpp11::register]]
network_R nextnetR_configmodel_clustered_triangles_network(integers degrees, integers triangles, double beta) {
    RNG_SCOPE_IF_NECESSARY;
    return new config_model_clustered_serrano(std::vector<int>(degrees.begin(), degrees.end()),
                                              std::vector<int>(triangles.begin(), triangles.end()),
                                              beta, rng_engine());
}


[[cpp11::register]]
network_R nextnetR_watts_strogatz_network(int size, int k, double p) {
    RNG_SCOPE_IF_NECESSARY;
    return new watts_strogatz(size, k, p, rng_engine());
}

[[cpp11::register]]
network_R nextnetR_barabasialbert_network(int size, int m) {
    RNG_SCOPE_IF_NECESSARY;
    return new barabasi_albert(size, rng_engine(), m);
}

[[cpp11::register]]
network_R nextnetR_cubiclattice2d_network(int edge_length) {
    return new cubic_lattice_2d(edge_length);
}

[[cpp11::register]]
network_R nextnetR_cubiclattice3d_network(int edge_length) {
    return new cubic_lattice_3d(edge_length);
}

[[cpp11::register]]
network_R nextnetR_cubiclattice4d_network(int edge_length) {
    return new cubic_lattice_4d(edge_length);
}

[[cpp11::register]]
network_R nextnetR_cubiclattice5d_network(int edge_length) {
    return new cubic_lattice_5d(edge_length);
}

[[cpp11::register]]
network_R nextnetR_cubiclattice6d_network(int edge_length) {
    return new cubic_lattice_6d(edge_length);
}

[[cpp11::register]]
network_R nextnetR_cubiclattice7d_network(int edge_length) {
    return new cubic_lattice_7d(edge_length);
}

[[cpp11::register]]
network_R nextnetR_cubiclattice8d_network(int edge_length) {
    return new cubic_lattice_8d(edge_length);
}

/*
 * nextnetR_adjacencylist_network
 */

[[cpp11::register]]
network_R nextnetR_adjacencylist_network(list input_al, bool is_undirected, bool above_diagonal) {
    if (!is_undirected && above_diagonal)
        stop("above_diagonal=TRUE only supported for undirected networks");

    const std::size_t n = input_al.size();
    if (n > std::numeric_limits<node_t>::max())
        stop("too many nodes");

    /* directed edges whose counterpart hasn't been observed */
    std::unordered_map<edge_t, std::ptrdiff_t, pair_hash> directed_edges = {};

    std::vector<std::vector<node_t>> adjacencylist(n);
    std::unordered_set<node_t> seen;
    bool is_simple = true;
    for(node_t u = 0; u < (node_t)n; ++u) {
        /* Append adjacency list for node u */
        std::vector<node_t> &u_adj = adjacencylist.at(u);

        /* Fill adjacencylist for node u */
        const integers input_u_adj = input_al[u];
        const std::size_t k = input_u_adj.size();
        u_adj.reserve(k);
        seen.clear();
        for(std::size_t i = 0; i < k; ++i) {
            /* Find target node v */
            node_t v = ((integers)input_u_adj)[i];

            /* Check validity and change index base of node indices to zero */
            if ((v < 1) || (v > (node_t)n))
                stop("nodes must be labelled consecutively from 1 to n");
            v -= 1;

            /* Track whether the network is simple */
            const bool is_selfedge = (u == v);
            const bool is_multiedge = !seen.insert(v).second;
            is_simple = is_simple && !is_selfedge && !is_multiedge;
  
            /* If claimed to be undirected, edge (u,v) exists iff (v,u) exits.
            * We track that by keep a balance of forward minus reverse multiplicity
            * of every undirected edge.
            */
            if (is_undirected && !above_diagonal && !is_selfedge) {
              const edge_t e  = (u < v) ? edge_t(u, v) : edge_t(v, u);
              auto i = directed_edges.find(e);
              if (i == directed_edges.end())
                  i = directed_edges.insert({ e, 0 }).first;
              i->second += (u < v) ? 1 : -1;
              if (i->second == 0)
                  directed_edges.erase(i);
            }

            /* If undirected and asked to add reversed edges only edges u <= v should exist */
            if (is_undirected && above_diagonal && (u > v)) 
                stop("unexpected edge (%d -> %d), %d > %d", u+1, v+1, u+1, v+1);

            /* Add edge (and possibly reverse edge) */  
            u_adj.push_back(v);
            if (above_diagonal && !is_selfedge)
              adjacencylist.at(v).push_back(u);
        }
    }

    /* Check that the network was indeed undirected */
    if (is_undirected && !directed_edges.empty())
        stop("%d edges have non-matching forward and reverse edges but is_undirected=TRUE",
             directed_edges.size());

    return new adjacencylist_network(std::move(adjacencylist),
                               is_undirected, is_simple);
}

/*
 * TEMPORAL NETWORKS
 */

[[cpp11::register]]
network_R nextnetR_erdos_renyi_temporalnetwork(int size, double avg_degree, double timescale) {
  RNG_SCOPE_IF_NECESSARY;
  
  return new temporal_erdos_renyi(size, avg_degree, timescale, rng_engine());
}

[[cpp11::register]]
network_R nextnetR_brownian_proximity_temporalnetwork(int size, double avg_degree, double radius,
                                            double D0, double D1, double gamma, SEXP dt) {
    RNG_SCOPE_IF_NECESSARY;
    if (dt == R_NilValue)
        return new brownian_proximity_network(size, avg_degree, radius, D0, D1, gamma,
                                              rng_engine());
    else
        return new brownian_proximity_network(size, avg_degree, radius, D0, D1, gamma,
                                              as_cpp<double>(dt), rng_engine());
}

[[cpp11::register]]
network_R nextnetR_empirical_contact_temporalnetwork(std::string path, bool finite_duration, double dt) {
    RNG_SCOPE_IF_NECESSARY;
    std::ifstream file(path);
    if (finite_duration)
        return new empirical_contact_network(file, empirical_contact_network::finite_duration, dt);
    else
        return new empirical_contact_network(file, empirical_contact_network::infitesimal_duration, dt);
}

[[cpp11::register]]
network_R nextnetR_sirx_temporalnetwork(const network_R& nw, double kappa0, double kappa) {
    if (!nw) stop("underlying network cannot be null"); 

    RNG_SCOPE_IF_NECESSARY;
    
    /* temporal_sirx_network stores the underlying network by reference, so
     * add it to the externalptr metadata to extend its lifetime
    */
    return { new temporal_sirx_network(*nw.get(), kappa0, kappa, rng_engine()),
             writable::list({"nw"_nm = nw}),
             true, true };
}

[[cpp11::register]]
network_R nextnetR_activity_driven_temporalnetwork(
    doubles activities_, int m,
    double eta_sus, double eta_inf, double b_sus, double b_inf
) {
    RNG_SCOPE_IF_NECESSARY;

    std::vector<double> activities(activities_.begin(), activities_.end());    
    return new activity_driven_network(
        std::move(activities), m, eta_sus, eta_inf, b_sus, b_inf, rng_engine());
}


/*
 * WEIGHTED NETWORKS
 */

[[cpp11::register]]
network_R nextnetR_erdos_renyi_weightednetwork(int size, double avg_degree, doubles weights, doubles probabilities) {
  RNG_SCOPE_IF_NECESSARY;

  std::vector<double> weights_stdvec(weights.begin(), weights.end());
  std::vector<double> probabilities_stdvec(probabilities.begin(), probabilities.end());
  return new weighted_erdos_renyi(size, avg_degree, weights_stdvec, probabilities_stdvec, rng_engine());
}

/*
 * nextnetR_adjacencylist_weightednetwork
 */

[[cpp11::register]]
network_R nextnetR_adjacencylist_weightednetwork(list input_al, bool is_undirected, bool above_diagonal) {
    if (!is_undirected && above_diagonal)
        stop("above_diagonal=TRUE only supported for undirected networks");
    
    const std::size_t n = input_al.size();
    if (n > std::numeric_limits<node_t>::max())
        stop("too many nodes");

    /* directed edges whose counterpart hasn't been observed */
    std::unordered_map<edge_t, double, pair_hash> directed_edges = {};
    
    std::vector<std::vector<std::pair<node_t, double>>> adjacencylist(n);
    std::unordered_set<node_t> seen;
    bool is_simple = true;
    for(node_t u = 0; u < (node_t)n; ++u) {
        /* Append adjacency list for node u */
        std::vector<std::pair<node_t, double>> &u_adj = adjacencylist.at(u);
        
        /* Fill adjacency list for node u */
        const integers input_u_adj = ((list)input_al[u])["n"];
        const doubles input_u_weights = ((list)input_al[u])["w"];
        if (input_u_adj.size() != input_u_weights.size())
            stop("sizes of m (neighbours) and w (weights) vectors must agree");
        const std::size_t k = input_u_adj.size();
        u_adj.reserve(k);
        seen.clear();
        for(std::size_t i = 0; i < k; ++i) {
            /* Find target node v */
            node_t v = ((integers)input_u_adj)[i];
            const double w = ((doubles)input_u_weights)[i];
            
            /* Check validity and change index base of node indices to zero */
            if ((v < 1) || (v > (node_t)n))
              stop("nodes must be labelled consecutively from 1 to n");
            v -= 1;
            
            /* Multi-edges are not allowed for weighted networks */
            if (seen.find(v) != seen.end())
                stop("multi-edges are not supported (" +
                                         std::to_string(u+1) + " -> " + std::to_string(v+1) +
                                         "already seen)");
            seen.insert(v);
            
            /* Track whether the network is simple */
            const bool is_selfedge = (u == v);
            is_simple = is_simple && !is_selfedge;
                
            /* If claimed to be undirected, edge (u,v) exists iff (v,u) exits.
             * We track that by keep a balance of forward minus reverse total
             * weight of every undirected edge. Note that we rely on floating
             * point addition/substraction to cancel perfectly here.
             */
            if (is_undirected && !above_diagonal && !is_selfedge) {
                const edge_t e  = (u < v) ? edge_t(u, v) : edge_t(v, u);
                auto i = directed_edges.find(e);
                if (i == directed_edges.end())
                    i = directed_edges.insert({ e, 0 }).first;
                i->second += (u < v) ? w : -w;
                if (i->second == 0.0)
                    directed_edges.erase(i);
            }
            
            /* If undirected and asked to add reversed edges only edges u <= v should exist */
            if (is_undirected && above_diagonal && (u > v)) 
                stop("unexpected edge %d -> %d, %d > %v", u+1, v+1, u+1, v+1);
              
            /* Add edge */
            u_adj.emplace_back(v, w);
            if (above_diagonal && !is_selfedge)
                adjacencylist.at(v).emplace_back(u, w);
        }
    }
    
    /* Check that the network was indeed undirected */
    /* Check that the network was indeed undirected */
    if (is_undirected && !directed_edges.empty())
        stop("%d edges have non-matching forward and reverse edges but is_undirected=TRUE",
             directed_edges.size());
    
    return new weighted_adjacencylist_network(std::move(adjacencylist),
                                              is_undirected, is_simple);
}
