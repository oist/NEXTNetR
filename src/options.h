#pragma once

#include <memory>
#include <optional>
#include <unordered_map>

#include <cpp11.hpp>
#include <cpp11/function.hpp>
#include <cpp11/external_pointer.hpp>

#include "episimR_types.h"
#include "rng.h"

#include "epidemics/types.h"
#include "epidemics/algorithm.h"
#include "epidemics/NextReaction.h"

using namespace std::literals;
using namespace episimR;
using namespace cpp11;
namespace writable = cpp11::writable;

namespace episimR {

typedef std::unordered_map<std::string, sexp> options_collection_type;

/**
 * @brief Handles a named option with C++ type T
 */
template<typename T>
struct option_handler {
    option_handler(std::string name_, T& dst_) :name(name_), dst(dst_) {}

    /**
     * @brief If the name matches, convert option value to C++, assign
     * to destination, and return assigned value re-wrapped as an R expression
     */
    void operator() (options_collection_type& opts_unprocessed, options_collection_type& opts_out) {
        auto i = opts_unprocessed.find(name);
        if (i != opts_unprocessed.end()) {
            // Found matching option
            sexp v = std::move(i->second);
            opts_unprocessed.erase(i);
            try {
                dst = as_cpp<T>(v);
            }
            catch (const std::exception &e) {
                throw std::runtime_error("invalid value for option "s + name + ": " + e.what());
            }
        }

        // Output option value
        opts_out.insert({ name, as_sexp(dst) });
    }

    std::string name;
    T& dst;
};

/**
 * @brief Declare a named option whose value is stored in a specific destination
 */
template<typename T>
option_handler<T> option(std::string name, T& dst)
{
    return option_handler<T>(name, dst);
}

/**
 * @brief Run the option handlers in order to process unprocessed options
 */
inline void option_handlers_run(options_collection_type& opts_unprocessed, options_collection_type& opts_out)
{
    /* Base case, do nothing */
}
template<typename Handler, typename ...Handlers>
inline void option_handlers_run(options_collection_type& opts_unprocessed, options_collection_type& opts_out,
                                Handler handler, Handlers... tail)
{
    /* Recursive case, handle head, recurse for tail */
    handler(opts_unprocessed, opts_out);
    option_handlers_run(opts_unprocessed, opts_out, tail...);
}

/**
 * @brief Process list of options with the specified list of handlers
 * 
 * Handlers should be values retured by `option(name, destination)`.
 */
template<typename ...Handlers>
list process_options(const list& opts, Handlers... handlers)
{
    if (!opts.empty() && !opts.named())
        throw std::runtime_error("options must be named");

    // Convert list of option values to map
    options_collection_type opts_unprocessed;
    const strings opts_names = opts.names();
    for(R_xlen_t i=0; i < opts.size(); ++i)
        opts_unprocessed.insert({ opts_names[i], opts[i] });

    // Process options
    options_collection_type opts_out;
    option_handlers_run(opts_unprocessed, opts_out, handlers...);

    // Convert output into a named list
    writable::list opts_out_list;
    writable::strings opts_out_names;
    for(const auto& i: opts_out) {
        opts_out_names.push_back(i.first);
        opts_out_list.push_back(i.second);
    }
    Rf_setAttrib(opts_out_list.data(), R_NamesSymbol, opts_out_names.data());

    // Warn about unprocessed options
    for(const auto& i: opts_unprocessed)
        warning("ignoring unknown option "s + i.first);

    return opts_out_list;
}

} /* namespace episimR */
