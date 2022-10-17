// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSYRE_SYMBOLIC_REGRESSION_H
#define DSYRE_SYMBOLIC_REGRESSION_H

#include <functional>
#include <random>
#include <string>
#include <vector>

#include <pagmo/problem.hpp>
#include <pagmo/s11n.hpp>
#include <pagmo/types.hpp>

#include <dsyre/detail/visibility.hpp>
#include <dsyre/expression.hpp>

namespace dsyre
{
class DSYRE_DLL_PUBLIC sr_problem
{
public:
    // Default constructor. Requirement on pagmo problems.
    sr_problem();

    // The constructor
    sr_problem(const std::vector<std::vector<double>> &points, const std::vector<double> &labels, unsigned length = 20u,
               std::vector<std::string> kernels = {"sum", "diff", "mul", "div"}, unsigned n_con = 0u,
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

    /// Gets the dataset points
    const std::vector<std::vector<double>> &get_points() const;

    /// Gets the dataset labels
    const std::vector<double> &get_labels() const;

    /// Transforms the chromosome into dsyre structures
    void pagmo2dsyre(std::vector<unsigned> &geno, std::vector<double> &cons, const pagmo::vector_double &x) const;

    /// Transforms the dsyre structures into a pagmo chromosome
    void dsyre2pagmo(pagmo::vector_double &x, const std::vector<unsigned> &geno, const std::vector<double> &cons) const;

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

    // Some mutable containers to avoid allocating in fitness
    mutable std::vector<double> m_cons;
    mutable std::vector<unsigned> m_geno;
    mutable std::vector<double> m_mse;
    mutable std::vector<unsigned> m_complexity;
};

namespace details
{
// This function is a global symbol put in the namespace. Its purpose is
// to be overridden in the python bindings so that it can extract from a py::object a
// c++ dsyre::sr_problem. Its use is in the UDAs evolve to access (both in C++ and python)
// the correct UDP.
inline std::function<const dsyre::sr_problem *(const pagmo::problem &)> extract_sr_cpp_py
    = [](const pagmo::problem &p) { return p.extract<dsyre::sr_problem>(); };
} // namespace details
} // namespace dsyre

// Why is this here?
PAGMO_S11N_PROBLEM_EXPORT_KEY(dsyre::sr_problem)
#endif