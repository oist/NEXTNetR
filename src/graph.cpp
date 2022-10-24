#include <optional>

#include <cpp11.hpp>
#include <cpp11/function.hpp>
#include <cpp11/external_pointer.hpp>

#include "episimR_types.h"
#include "rng.h"

#include "epidemics/types.h"
#include "epidemics/graph.h"

using namespace episimR;

using namespace cpp11;
namespace writable = cpp11::writable;

[[cpp11::register]]
integers episimR_graph_outdegree(const graph_R& nw, integers nodes) {
    if (!nw) throw std::runtime_error("graph cannot be NULL"); 
    
    /* Must enter RNG scope since graphs may generate their topology on the fly */
    RNG_SCOPE_IF_NECESSARY;
    
    const std::size_t l = nodes.size();
    writable::integers r;
    r.reserve(l);
    for(int j = 0; j < l; ++j) {
        const int n = nodes[j];
        r.push_back(((n >= 1) && (n <= l)) ? nw->outdegree(n - 1) : NA_INTEGER);
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
    for(int j = 0; j < l; ++j) {
        const node_t n = nodes[j];
        if ((n < 1) || (n > l)) {
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
graph_R episimR_scalefree_graph(int size, bool assortative) {
    RNG_SCOPE_IF_NECESSARY;
    return new scale_free(size, rng_engine());
}

[[cpp11::register]]
graph_R episimR_stored_graph(r_string filename, bool assortative) {
    return new imported_network((std::string)filename);
}
