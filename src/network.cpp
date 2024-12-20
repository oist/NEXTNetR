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

using namespace nextnetR;

using namespace cpp11;
namespace writable = cpp11::writable;

[[cpp11::register]]
int nextnetR_network_size(const network_R& nw) {
    if (!nw) throw std::runtime_error("graph cannot be NULL"); 
    
    /* Must enter RNG scope since graphs may generate their topology on the fly */
    RNG_SCOPE_IF_NECESSARY;
    
    return nw->nodes();
}

[[cpp11::register]]
bool nextnetR_network_is_undirected(const network_R& nw) {
    if (!nw) throw std::runtime_error("graph cannot be NULL"); 
    
    /* Must enter RNG scope since graphs may generate their topology on the fly */
    RNG_SCOPE_IF_NECESSARY;
    
    return nw->is_undirected();
}

[[cpp11::register]]
integers nextnetR_network_outdegree(const network_R& nw, integers nodes) {
    if (!nw) throw std::runtime_error("graph cannot be NULL"); 
    
    /* Must enter RNG scope since graphs may generate their topology on the fly */
    RNG_SCOPE_IF_NECESSARY;
    
    const std::size_t l = nodes.size();
    writable::integers r;
    r.reserve(l);
    for(std::size_t j = 0; j < l; ++j) {
        const int n = nodes[j];
        r.push_back(((n >= 1) && (n <= (node_t)l)) ? nw->outdegree(n - 1) : NA_INTEGER);
    }
    return r;
}

[[cpp11::register]]
integers nextnetR_network_neighbour(const network_R& nw, integers nodes, integers indices) {
    if (!nw) throw std::runtime_error("graph cannot be NULL"); 

    /* Must enter RNG scope since graphs may generate their topology on the fly */
    RNG_SCOPE_IF_NECESSARY;
    
    if (nodes.size() != indices.size())
        throw std::runtime_error("number of nodes and number of indices must agree");
    
    /* Create output */
    const std::size_t l = nodes.size();
    writable::integers r;
    r.reserve(l);
    
    /* Fill */
    for(std::size_t j = 0; j < l; ++j) {
        const node_t n = nodes[j];
        if ((n < 1) || (n > (node_t)l)) {
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
list nextnetR_network_adjacencylist(const network_R& nw) {
    if (!nw) throw std::runtime_error("graph cannot be NULL"); 

    /* Must enter RNG scope since graphs may generate their topology on the fly */
    RNG_SCOPE_IF_NECESSARY;

    /* TODO: For instances of network_adjacencylist, this could be done more efficiently */

    /* Allocate output, a vector of node labels and a list of neighbour vectors */
    const int l = nw->nodes();
    writable::integers nodes;
    writable::list neighbours;
    nodes.reserve(l);
    neighbours.reserve(l);

    /* Interate over nodes and generate output */
    for(int n = 0; n < l; ++n) {
        const int d = nw->outdegree(n);
        writable::integers node_neighbours;
        node_neighbours.reserve(d);
        for(int i=0; i < d; ++i) {
            const int m = nw->neighbour(n, i);
            node_neighbours.push_back((m >= 0) ? (m + 1) : NA_INTEGER);
        }
        
        /* Append to output */
        nodes.push_back(n + 1);
        neighbours.push_back(node_neighbours);
    }

    return writable::list({
        "nodes"_nm = nodes,
        "neighbours"_nm = neighbours
    });
}

[[cpp11::register]]
list nextnetR_weighted_network_adjacencylist(const network_R& nw) {
    if (!nw) throw std::runtime_error("network cannot be NULL"); 

    weighted_network* const nw_weighted = dynamic_cast<weighted_network*>(nw.get());
    if (nw_weighted == nullptr)
        throw std::runtime_error("network is not weighted");

    /* Must enter RNG scope since graphs may generate their topology on the fly */
    RNG_SCOPE_IF_NECESSARY;

    /* TODO: For instances of weighted_network_adjacencylist, this could be done more efficiently */

    /* Allocate output, a vector of node labels and a list of neighbour vectors */
    const int l = nw_weighted->nodes();
    writable::integers nodes;
    writable::list neighbours;
    nodes.reserve(l);
    neighbours.reserve(l);

    /* Interate over nodes and generate output */
    for(int n = 0; n < l; ++n) {
        const int d = nw_weighted->outdegree(n);
        writable::integers node_neighbours;
        writable::doubles node_weights;
        node_neighbours.reserve(d);
        for(int i=0; i < d; ++i) {
            double w = NAN;
            const int m = nw_weighted->neighbour(n, i, &w);
            node_neighbours.push_back((m >= 0) ? (m + 1) : NA_INTEGER);
            node_weights.push_back((m >= 0) ? w : NA_INTEGER);
        }
        
        /* Append to output */
        nodes.push_back(n + 1);
        neighbours.push_back(writable::list {
            "n"_nm = node_neighbours,
            "w"_nm = node_weights
        });
    }

    return writable::list({
        "nodes"_nm = nodes,
        "neighbours"_nm = neighbours
    });
}


[[cpp11::register]]
list nextnetR_network_bounds(const network_R& nw) {
    if (!nw) throw std::runtime_error("graph cannot be NULL"); 

    RNG_SCOPE_IF_NECESSARY;
    
    network_embedding* const ebd = dynamic_cast<network_embedding*>(nw.get());
    if (ebd == nullptr)
        throw std::runtime_error("graph is not embedded into R^d");
    
    std::vector<double> a, b;
    ebd->bounds(a, b);

    return writable::list({
        writable::doubles(a.begin(), a.end()),
        writable::doubles(b.begin(), b.end())
    });
}

[[cpp11::register]]
doubles_matrix<> nextnetR_network_coordinates(const network_R& nw, integers nodes) {
    if (!nw) throw std::runtime_error("graph cannot be NULL"); 

    RNG_SCOPE_IF_NECESSARY;
    
    network_embedding* const ebd = dynamic_cast<network_embedding*>(nw.get());
    if (ebd == nullptr)
        throw std::runtime_error("graph is not embedded into R^d");
    
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
    if (!nw) throw std::runtime_error("graph cannot be NULL"); 

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
network_R nextnetR_empirical_network(r_string filename) {
    return new empirical_network((std::string)filename);
}

[[cpp11::register]]
network_R nextnetR_erdos_reyni_network(int size, double avg_degree) {
    RNG_SCOPE_IF_NECESSARY;
    return new erdos_reyni(size, avg_degree, rng_engine());
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
    return new watts_strogatz(size, k, p, rng_engine());
}

[[cpp11::register]]
network_R nextnetR_barabasialbert_network(int size, int m) {
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
network_R nextnetR_empirical_temporalnetwork(std::string file, bool finite_duration, double dt) {
    RNG_SCOPE_IF_NECESSARY;
    if (finite_duration)
        return new empirical_temporal_network(file, empirical_temporal_network::finite_duration, dt);
    else
        return new empirical_temporal_network(file, empirical_temporal_network::infitesimal_duration, dt);
}

[[cpp11::register]]
network_R nextnetR_sirx_temporalnetwork(const network_R& nw, double kappa0, double kappa) {
    if (!nw) throw std::runtime_error("underlying graph cannot be null"); 

    RNG_SCOPE_IF_NECESSARY;
    
    /* temporal_sirx_network stores the underlying network by reference, so
     * add it to the externalptr metadata to extend its lifetime
    */
    return { new temporal_sirx_network(*nw.get(), kappa0, kappa),
             writable::list({"nw"_nm = nw}),
             true, true };
}

namespace {

class nextnetR_adjacencylist_network_impl : public adjacencylist_network {
public:
    nextnetR_adjacencylist_network_impl(std::vector<std::vector<node_t>>&& al, bool undirected_ = false)
        :undirected(undirected_)
    {
        adjacencylist = std::move(al);
    }
    
    virtual bool is_undirected() {
        return undirected;
    }
    
private:
    bool undirected;
};
  
}


[[cpp11::register]]
network_R nextnetR_adjacencylist_network(list input_al, bool is_undirected) {
    const std::size_t n = input_al.size();
    if (n > std::numeric_limits<node_t>::max())
      throw std::runtime_error("too many nodes");

    /* directed edges whose counterpart hasn't been observed */
    std::set<std::size_t> directed_edges;
    
    std::vector<std::vector<node_t>> adjacencylist;
    adjacencylist.reserve(n);
    for(node_t u = 0; u < (node_t)n; ++u) {
        /* Append adjacency list for node u */
        adjacencylist.emplace_back();
        std::vector<node_t> &u_adj = adjacencylist.back();
        
        /* Fill adjacencylist for node u */
        const integers input_u_adj = input_al[u];
        const std::size_t k = input_u_adj.size();
        u_adj.reserve(k);
        std::set<node_t> seen;
        for(std::size_t i = 0; i < k; ++i) {
            /* Find target node v */
            const node_t v = ((integers)input_u_adj)[i];
            
            /* Check validity */
            if ((v < 1) || (v > (node_t)n))
                throw std::runtime_error("nodes must be labelled consecutively from 1 to n");
            if (seen.find(v) != seen.end())
                throw std::runtime_error("multi-edges are not supported (" +
                                         std::to_string(u) + " -> " + std::to_string(v) +
                                         "already seen)");
            seen.insert(v);
                
            /* If claimed to be undirected, edge (u,v) exists iff (v,u) exits.
             * We insert into directed_edges upon seeing {u,v} for the first time,
             * and delete when second a second time.
             */
            if (is_undirected) {
                const std::size_t e = edge_index_undirected(u, v);
                const auto i = directed_edges.find(e);
                if (i == directed_edges.end())
                    directed_edges.insert(e);
                else
                    directed_edges.erase(i);
            }
              
            /* Add edge */  
            u_adj.push_back(v - 1);
        }
    }
    
    /* Check that the network was indeed undirected */
    if (is_undirected && !directed_edges.empty())
        throw std::runtime_error(std::to_string(directed_edges.size()) + " directed edges found " +
                                 "in network claimed to be undirected");
    
    return new nextnetR_adjacencylist_network_impl(std::move(adjacencylist), is_undirected);
}

namespace {

class nextnetR_weighted_adjacencylist_network_impl : public weighted_adjacencylist_network {
public:
    nextnetR_weighted_adjacencylist_network_impl(std::vector<std::vector<std::pair<node_t, double>>>&& al, bool undirected_ = false)
        :undirected(undirected_)
    {
        adjacencylist = std::move(al);
    }
    
    virtual bool is_undirected() {
        return undirected;
    }
    
private:
    bool undirected;
};
  
}

[[cpp11::register]]
network_R nextnetR_weighted_adjacencylist_network(list input_al, bool is_undirected) {
    const std::size_t n = input_al.size();
    if (n > std::numeric_limits<node_t>::max())
      throw std::runtime_error("too many nodes");

    /* directed edges whose counterpart hasn't been observed */
    std::set<std::size_t> directed_edges;
    
    std::vector<std::vector<std::pair<node_t, double>>> adjacencylist;
    adjacencylist.reserve(n);
    for(node_t u = 0; u < (node_t)n; ++u) {
        /* Append adjacency list for node u */
        adjacencylist.emplace_back();
        std::vector<std::pair<node_t, double>> &u_adj = adjacencylist.back();
        
        /* Fill adjacencylist for node u */
        const integers input_u_adj = ((list)input_al[u])["n"];
        const doubles input_u_weights = ((list)input_al[u])["w"];
        if (input_u_adj.size() != input_u_weights.size())
            throw std::runtime_error("sizes of m (neighbours) and w (weights) vectors must agree");
        const std::size_t k = input_u_adj.size();
        u_adj.reserve(k);
        std::set<node_t> seen;
        for(std::size_t i = 0; i < k; ++i) {
            /* Find target node v */
            const node_t v = ((integers)input_u_adj)[i];
            const double w = ((doubles)input_u_weights)[i];
            
            /* Check validity */
            if ((v < 1) || (v > (node_t)n))
                throw std::runtime_error("nodes must be labelled consecutively from 1 to n");
            if (seen.find(v) != seen.end())
                throw std::runtime_error("multi-edges are not supported (" +
                                         std::to_string(u) + " -> " + std::to_string(v) +
                                         "already seen)");
            seen.insert(v);
                
            /* If claimed to be undirected, edge (u,v) exists iff (v,u) exits.
             * We insert into directed_edges upon seeing {u,v} for the first time,
             * and delete when second a second time.
             */
            if (is_undirected) {
                const std::size_t e = edge_index_undirected(u, v);
                const auto i = directed_edges.find(e);
                if (i == directed_edges.end())
                    directed_edges.insert(e);
                else
                    directed_edges.erase(i);
            }
              
            /* Add edge */  
            u_adj.emplace_back(v - 1, w);
        }
    }
    
    /* Check that the network was indeed undirected */
    if (is_undirected && !directed_edges.empty())
        throw std::runtime_error(std::to_string(directed_edges.size()) + " directed edges found " +
                                 "in network claimed to be undirected");
    
    return new nextnetR_weighted_adjacencylist_network_impl(std::move(adjacencylist), is_undirected);
}
