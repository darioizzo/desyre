// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSYRE_MES4DSYRE_H
#define DSYRE_MES4DSYRE_H

#include <random>
#include <vector>

#include <pagmo/population.hpp>
#include <pagmo/types.hpp>

#include <dsyre/detail/visibility.hpp>

namespace dsyre
{
/// Memetic Evolutionary Strategy for dsyre
/**
 *
 * The term Memetic is widely used, in the context of meta-heuristic search, to indicate a synergy between any
 * population-based approach with local improvement procedures. The resulting algorithms are also referred to, in the
 * literature, as Baldwinian evolutionary algorithms (EAs), Lamarckian EAs, cultural algorithms, or genetic local
 * searches. The very same approach, is seen by many just as an hybridization of a global search technique with a
 * local search technique. Regardless of the terminology and point of view, a memetic approach is applicable to symbolic
 * regression tasks using the dsyre encoding and able to improve considerably on the long standing issue of finding
 * constants in Genetic Programming.
 *
 * @see Izzo, Dario, Francesco Biscani, and Alessio Mereta. "Differentiable genetic programming." In European Conference
 * on Genetic Programming, pp. 35-51. Springer, 2017.
 *
 * In this C++ class we offer an UDA (User Defined Algorithm for the pagmo optimization suite) hybridizing a classic
 * Evolutionary Strategy with a second order Newton
 * search step able to help finding the best values for the ephemeral constants. The resulting algorithm is
 * outlined by the following pseudo-algorithm:
 *
 * @code{.unparsed}
 * > Start from a population (pop) of dimension N
 * > while i < gen
 * > > Mutation: create a new population pop2 mutating the best individual (only the integer part is affected).
 * > > Life long learning: change the constant values using the gradient and hessians of the phenotype.
 * > > Reinsertion: set pop to contain the best N individuals taken from pop and pop2.
 * @endcode
 *
 */
class DSYRE_DLL_PUBLIC mes4dsyre
{
public:
    /// Single entry of the log (gen, fevals, best, constants, formula)
    typedef std::tuple<unsigned, unsigned long long, double, pagmo::vector_double, std::string> log_line_type;
    /// The log
    typedef std::vector<log_line_type> log_type;

    /// Constructor
    mes4dsyre(unsigned gen = 1u, unsigned max_mut = 5u, double ftol = 1e-10, unsigned seed = std::random_device{}());

    /// Algorithm evolve method
    pagmo::population evolve(pagmo::population pop) const;

    /// Sets the seed
    void set_seed(unsigned seed);

    /// Gets the seed
    unsigned get_seed() const;

    /// Sets the algorithm verbosity
    void set_verbosity(unsigned level);

    /// Gets the verbosity level
    unsigned get_verbosity() const
    {
        return m_verbosity;
    }

    /// Algorithm name
    std::string get_name() const;

    /// Extra info
    std::string get_extra_info() const;

    /// Get log
    const log_type &get_log() const;

private:
    // This prints to screen and logs one single line.
    void log_single_line(unsigned gen, unsigned long long fevals, double best_f, const std::string &formula,
                         const std::vector<double> &x, unsigned ncon) const;

public:
    /// Object serialization
    template <typename Archive>
    void serialize(Archive &ar, unsigned)
    {
        ar &m_gen;
        ar &m_max_mut;
        ar &m_ftol;
        ar &m_e;
        ar &m_seed;
        ar &m_verbosity;
        ar &m_log;
    }

private:
    unsigned m_gen;
    unsigned m_max_mut;
    double m_ftol;
    mutable pagmo::detail::random_engine_type m_e;
    unsigned m_seed;
    unsigned m_verbosity;
    mutable log_type m_log;
};
} // namespace dsyre
#endif