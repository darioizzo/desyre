// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSYRE_SYMBOLIC_REGRESSION_H
#define DSYRE_SYMBOLIC_REGRESSION_H

#include <random>
#include <string>
#include <vector>

#include <pagmo/types.hpp>

#include <dsyre/expression.hpp>

namespace dsyre
{
class DSYRE_DLL_PUBLIC sr_problem
{
public:
    // Default constructor. Requirement on pagmo problems.
    sr_problem();

    // The constructor
    sr_problem(const std::vector<std::vector<double>> &points, const std::vector<double> &labels,
                        unsigned length = 20u, std::vector<std::string> kernels = {"sum", "diff", "mul", "div"},
                        unsigned n_con = 0u,
                        bool multi_objective = false // when true the fitness also returns the formula complexity
    );

    /// Number of objectives
    pagmo::vector_double::size_type get_nobj() const;

    /// Fitness computation
    pagmo::vector_double fitness(const pagmo::vector_double &x) const;

    /// Box-bounds
    std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const;

    /// Integer dimension
    pagmo::vector_double::size_type get_nix() const;

    /// Problem name
    std::string get_name() const;

    /// Extra info
    std::string get_extra_info() const;

    /// Human-readable representation of a decision vector. (unsimplified)
    std::string pretty(const pagmo::vector_double &x) const;

    /// Human-readable representation of a decision vector. (simplified using symengine)
    std::string prettier(const pagmo::vector_double &x) const;

    /// Gets the inner expression object
    const expression &get_expression() const;

private:
    // Sanity checks for the dataset
    void sanity_checks(const std::vector<std::vector<double>> &points, const std::vector<double> &labels) const;

public:
    /// Object serialization
    template <typename Archive>
    void serialize(Archive &ar, unsigned)
    {
        ar &m_points;
        ar &m_labels;
        ar &m_length;
        ar &m_kernels;
        ar &m_ncon;
        ar &m_nvar;
        ar &m_multi_objective;
        ar &m_ex;
    }

private:
    // The dataset
    std::vector<std::vector<double>> m_points;
    std::vector<double> m_labels;
    // The length of the us
    unsigned m_length;
    // The kernels requested
    std::vector<std::string> m_kernels;
    // Maximum number of constants allowed
    unsigned m_ncon;
    // Number of variables
    unsigned m_nvar;
    // Multiobjective problem
    bool m_multi_objective;
    // The expression object
    expression m_ex;
};
} // namespace dsyre
#endif