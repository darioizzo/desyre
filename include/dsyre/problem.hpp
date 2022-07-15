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
class symbolic_regression
{
public:
    // Default constructor. Requirement on pagmo problems.
    symbolic_regression()
        : m_points(1), m_labels(1), m_length(1), m_kernels({"sum", "diff", "mul", "div"}), m_n_con(0),
          m_multi_objective(true));

    // The constructor
    symbolic_regression(const std::vector<std::vector<double>> &points, const std::vector<double> &labels,
                        unsigned length = 20u, std::vector<std::string> kernels = {"sum", "diff", "mul", "div"},
                        unsigned n_con = 0u,
                        bool multi_objective = false, // when true the fitness also returns the formula complexity
                        )
        : m_points(points), m_labels(labels), m_length(length), m_kernels(kernels), m_n_con(n_con),
          m_multi_objective(multi_objective)
    {
        // We check the dataset for consistency (and avoid empty dataset)
        sanity_checks();
        // We initialize the dsyre expression
        m_expression = expression(m_points[0].size(), m_n_cons, m_kernels);
        // We store the number of varables / constants
        m_nvar = points[0].size();
    }

    /// Number of objectives
    pagmo::vector_double::size_type get_nobj() const
    {
        return 1u + m_multi_objective;
    }

    /// Fitness computation
    pagmo::vector_double fitness(const pagmo::vector_double &x) const;

    /// Box-bounds
    std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const
    {
        // Bounds on ephemeral constants are -10, 10;
        std::vector<double> lb(m_ncon + m_length * 3, -10.);
        std::vector<double> ub(m_ncon + m_length * 3, 10.);
        // Bounds on the integer part are derived from the number of kernels and position
        unsigned nus = m_ncon + m_nvar;
        for (auto i = 0u; i < m_length; ++i) {
            lb[3 * i + m_ncon] = 0;
            ub[3 * i + m_ncon + 1] = m_kernels.size() - 1;
            lb[3 * i + m_ncon + 1] = 0;
            ub[3 * i + m_ncon + 1] = nus - 1;
            lb[3 * i + m_ncon + 2] = 0;
            ub[3 * i + m_ncon + 2] = nus - 1;
            nus++;
        }
        return {lb, ub};
    }

    /// Integer dimension
    pagmo::vector_double::size_type get_nix() const
    {
        return m_length * 3;
    }

    /// Problem name
    std::string get_name() const
    {
        return "A symbolic regression problem using dsyre encoding";
    }

    /// Extra info
    std::string get_extra_info() const;

    /// Human-readable representation of a decision vector. (unsimplified)
    std::string pretty(const pagmo::vector_double &x) const;

    /// Human-readable representation of a decision vector. (simplified using symengine)
    std::string prettier(const pagmo::vector_double &x) const;

    /// Gets the inner expression object
    const expression<double> &get_cgp() const
    {
        return m_expression;
    }

private:
    // Sanity checks for the dataset
    void sanity_checks() const;

public:
    /// Object serialization
    template <typename Archive>
    void serialize(Archive &ar, unsigned)
    {
        ar &m_points;
        ar &m_labels;
        ar &m_length;
        ar &m_kernels;
        ar &m_n_con;
        ar &m_n_var;
        ar &m_multi_objective;
        ar &m_expression;
    }

private:
    // The dataset
    std::vector<std::vector<double>> m_points;
    std::vector<double> m_labels;
    // The length of the us
    unsigned m_length;
    // The kernel idxs
    std::vector<std::string> m_kernels;
    // Maximum number of constants allowed
    unsigned m_n_con;
    // Number of variables
    unsigned m_n_var;
    // Multiobjective problem
    bool m_multi_objective;
    // The expression object
    expression m_expression;
};
} // namespace dsyre
#endif