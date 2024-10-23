#include <utility>
#include <limits>

#include <cpp11.hpp>
#include <cpp11/function.hpp>
#include <cpp11/external_pointer.hpp>

#include "episimR_types.h"
#include "rng.h"

#include "epidemics/types.h"
#include "epidemics/graph.h"
#include "epidemics/dynamic_graph.h"
#include "epidemics/brownian_proximity_graph.h"

using namespace episimR;

using namespace cpp11;
namespace writable = cpp11::writable;

[[cpp11::register]]
int episimR_graph_size(const graph_R& nw) {
    if (!nw) throw std::runtime_error("graph cannot be NULL"); 
    
    /* Must enter RNG scope since graphs may generate their topology on the fly */
    RNG_SCOPE_IF_NECESSARY;
    
    return nw->nodes();
}

[[cpp11::register]]
bool episimR_graph_is_undirected(const graph_R& nw) {
    if (!nw) throw std::runtime_error("graph cannot be NULL"); 
    
    /* Must enter RNG scope since graphs may generate their topology on the fly */
    RNG_SCOPE_IF_NECESSARY;
    
    return nw->is_undirected();
}

[[cpp11::register]]
integers episimR_graph_outdegree(const graph_R& nw, integers nodes) {
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
integers episimR_graph_neighbour(const graph_R& nw, integers nodes, integers indices) {
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
list episimR_graph_adjacencylist(const graph_R& nw) {
    if (!nw) throw std::runtime_error("graph cannot be NULL"); 

    /* Must enter RNG scope since graphs may generate their topology on the fly */
    RNG_SCOPE_IF_NECESSARY;

    /* TODO: For instances of graph_adjacencylist, this could be done more efficiently */

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
list episimR_graph_bounds(const graph_R& nw) {
    if (!nw) throw std::runtime_error("graph cannot be NULL"); 

    RNG_SCOPE_IF_NECESSARY;
    
    graph_embedding* const ebd = dynamic_cast<graph_embedding*>(nw.get());
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
doubles_matrix<> episimR_graph_coordinates(const graph_R& nw, integers nodes) {
    if (!nw) throw std::runtime_error("graph cannot be NULL"); 

    RNG_SCOPE_IF_NECESSARY;
    
    graph_embedding* const ebd = dynamic_cast<graph_embedding*>(nw.get());
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
graph_R episimR_stored_graph(r_string filename) {
    return new imported_network((std::string)filename);
}

[[cpp11::register]]
graph_R episimR_erdos_reyni_graph(int size, double avg_degree) {
    RNG_SCOPE_IF_NECESSARY;
    return new erdos_reyni(size, avg_degree, rng_engine());
}

[[cpp11::register]]
graph_R episimR_fully_connected_graph(int size) {
    RNG_SCOPE_IF_NECESSARY;
    return new fully_connected(size, rng_engine());
}

[[cpp11::register]]
graph_R episimR_acyclic_graph(int size, double avg_degree, bool reduced_root_degree) {
    RNG_SCOPE_IF_NECESSARY;
    return new acyclic(avg_degree, reduced_root_degree, rng_engine());
}

[[cpp11::register]]
graph_R episimR_configmodel_graph(integers degrees) {
    RNG_SCOPE_IF_NECESSARY;
    return new config_model(std::vector<int>(degrees.begin(), degrees.end()), rng_engine());
}

[[cpp11::register]]
graph_R episimR_configmodel_clustered_alpha_graph(integers degrees, double alpha, double beta) {
    RNG_SCOPE_IF_NECESSARY;
    return new config_model_clustered_serrano(std::vector<int>(degrees.begin(), degrees.end()),
                                              alpha, beta, rng_engine());
}

[[cpp11::register]]
graph_R episimR_configmodel_clustered_ck_graph(integers degrees, SEXP ck, double beta) {
    RNG_SCOPE_IF_NECESSARY;
    function ck_rf = ck;
    auto ck_lambda = [ck_rf] (int k) { return as_cpp<double>(ck_rf(k)); };
    return new config_model_clustered_serrano(std::vector<int>(degrees.begin(), degrees.end()),
                                              ck_lambda, beta, rng_engine());
}

[[cpp11::register]]
graph_R episimR_configmodel_clustered_triangles_graph(integers degrees, integers triangles, double beta) {
    RNG_SCOPE_IF_NECESSARY;
    return new config_model_clustered_serrano(std::vector<int>(degrees.begin(), degrees.end()),
                                              std::vector<int>(triangles.begin(), triangles.end()),
                                              beta, rng_engine());
}

[[cpp11::register]]
graph_R episimR_barabasialbert_graph(int size, int m) {
    return new barabasi_albert(size, rng_engine(), m);
}

[[cpp11::register]]
graph_R episimR_cubiclattice2d_graph(int edge_length) {
    return new cubic_lattice_2d(edge_length);
}

[[cpp11::register]]
graph_R episimR_cubiclattice3d_graph(int edge_length) {
    return new cubic_lattice_3d(edge_length);
}

[[cpp11::register]]
graph_R episimR_cubiclattice4d_graph(int edge_length) {
    return new cubic_lattice_4d(edge_length);
}

[[cpp11::register]]
graph_R episimR_cubiclattice5d_graph(int edge_length) {
    return new cubic_lattice_5d(edge_length);
}

[[cpp11::register]]
graph_R episimR_cubiclattice6d_graph(int edge_length) {
    return new cubic_lattice_6d(edge_length);
}

[[cpp11::register]]
graph_R episimR_cubiclattice7d_graph(int edge_length) {
    return new cubic_lattice_7d(edge_length);
}

[[cpp11::register]]
graph_R episimR_cubiclattice8d_graph(int edge_length) {
    return new cubic_lattice_8d(edge_length);
}

[[cpp11::register]]
graph_R episimR_brownian_proximity_dyngraph(int size, double avg_degree, double radius, double D, SEXP dt) {
    RNG_SCOPE_IF_NECESSARY;
    if (dt == R_NilValue)
        return new brownian_proximity_graph(size, avg_degree, radius, D, rng_engine());
    else
        return new brownian_proximity_graph(size, avg_degree, radius, D, as_cpp<double>(dt), rng_engine());
}

[[cpp11::register]]
graph_R episimR_empirical_dyngraph(std::string file, bool finite_duration, double dt) {
    RNG_SCOPE_IF_NECESSARY;
    if (finite_duration)
        return new dynamic_empirical_network(file, dynamic_empirical_network::finite_duration, dt);
    else
        return new dynamic_empirical_network(file, dynamic_empirical_network::infitesimal_duration, dt);
}

[[cpp11::register]]
graph_R episimR_sirx_dyngraph(const graph_R& nw, double kappa0, double kappa) {
    if (!nw) throw std::runtime_error("underlying graph cannot be null"); 

    RNG_SCOPE_IF_NECESSARY;
    
    /* dynamic_sirx_networ stores the underlying network by reference, so
     * add it to the externalptr metadata to extend its lifetime
    */
    return { new dynamic_sirx_network(*nw.get(), kappa0, kappa),
             writable::list({"nw"_nm = nw}),
             true, true };
}

namespace {

class graph_userdefined : public graph_adjacencylist {
public:
    graph_userdefined(std::vector<std::vector<node_t>>&& al, bool undirected_ = false)
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
graph_R episimR_userdefined_graph(list input_al, bool is_undirected) {
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
    
    return new graph_userdefined(std::move(adjacencylist));
}
