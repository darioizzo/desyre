// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <random>
#include <string>
#include <vector>

#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <pagmo/types.hpp>
#include <symengine/expression.h>

#include <dsyre/expression.hpp>
#include <dsyre/sr_problem.hpp>

namespace dsyre
{
sr_problem::sr_problem()
    : m_points(1), m_labels(1), m_length(1), m_kernels({"sum", "diff", "mul", "div"}), m_ncon(0),
      m_multi_objective(false)
{
    m_points[0] = {1.};
    m_labels[0] = 1.;
    m_ex = expression(1, 0u, m_kernels);
    m_nvar = 1u;
}

sr_problem::sr_problem(const std::vector<std::vector<double>> &points, const std::vector<double> &labels,
                       unsigned length, std::vector<std::string> kernels, unsigned n_con,
                       bool multi_objective // when true the fitness also returns the formula complexity
                       )
    : m_points(points), m_labels(labels), m_length(length), m_kernels(kernels), m_ncon(n_con),
      m_multi_objective(multi_objective)
{
    // We check the dataset for consistency (and avoid empty dataset)
    sanity_checks(points, labels);
    // We initialize the dsyre expression
    m_ex = expression(m_points[0].size(), m_ncon, m_kernels);
    // We store the number of varables / constants
    m_nvar = points[0].size();
}

pagmo::vector_double::size_type sr_problem::get_nobj() const
{
    return 1u + m_multi_objective;
}

pagmo::vector_double sr_problem::fitness(const pagmo::vector_double &x) const
{
    std::vector<double> retval(1u + m_multi_objective, 0.);
    m_cons.resize(m_ncon);
    m_geno.resize(3 * m_length);
    // We need to copy the pagmo chromosome into dsyre genotype and constants
    std::copy(x.begin(), x.begin() + m_ncon, m_cons.begin());
    std::copy(x.begin() + m_ncon, x.end(), m_geno.begin());
    // Here we compute the mse for all us.
    m_ex.mse(m_mse, m_geno, m_cons, m_points, m_labels);
    // Here we select the smallest.
    auto best_it = std::min_element(m_mse.begin(), m_mse.end());
    retval[0] = *best_it;
    if (m_multi_objective) {
        auto best_idx = std::distance(m_mse.begin(), best_it);
        m_ex.complexity(m_complexity, m_geno);
        retval[1] = m_complexity[best_idx];
    }
    return retval;
}

std::pair<pagmo::vector_double, pagmo::vector_double> sr_problem::get_bounds() const
{
    // Bounds on ephemeral constants are -10, 10;
    std::vector<double> lb(m_ncon + m_length * 3, -10.);
    std::vector<double> ub(m_ncon + m_length * 3, 10.);
    // Bounds on the integer part are derived from the number of kernels and position
    unsigned nus = m_ncon + m_nvar;
    for (auto i = 0u; i < m_length; ++i) {
        lb[3 * i + m_ncon] = 0;
        ub[3 * i + m_ncon] = m_kernels.size() - 1;
        lb[3 * i + m_ncon + 1] = 0;
        ub[3 * i + m_ncon + 1] = nus - 1;
        lb[3 * i + m_ncon + 2] = 0;
        ub[3 * i + m_ncon + 2] = nus - 1;
        nus++;
    }
    return {lb, ub};
}

pagmo::vector_double::size_type sr_problem::get_nix() const
{
    return m_length * 3;
}

std::string sr_problem::get_name() const
{
    return "A symbolic regression problem using dsyre encoding";
}

const expression &sr_problem::get_expression() const
{
    return m_ex;
}

const std::vector<std::vector<double>> &sr_problem::get_points() const
{
    return m_points;
}

const std::vector<double> &sr_problem::get_labels() const
{
    return m_labels;
}

void sr_problem::sanity_checks(const std::vector<std::vector<double>> &points, const std::vector<double> &labels) const
{
    if (points.size() != labels.size()) {
        throw std::invalid_argument(
            "Dataset malformed: the number of labels provided does not match the number of points");
    }
    if (points.size() == 0) {
        throw std::invalid_argument(
            "Dataset malformed: the size of points cannot be zero when constructiong a sr_problem");
    }
    for (auto &item : points) {
        if (item.size() != points[0].size()) {
            throw std::invalid_argument("Dataset malformed: check the dimensions of all points.");
        }
    }
}

void sr_problem::pagmo2dsyre(std::vector<unsigned> &geno, std::vector<double> &cons,
                             const pagmo::vector_double &x) const
{
    geno.resize(x.size() - m_ncon);
    cons.resize(m_ncon);
    std::copy(x.begin(), x.begin() + m_ncon, cons.begin());
    std::copy(x.begin() + m_ncon, x.end(), geno.begin());
}

void sr_problem::dsyre2pagmo(pagmo::vector_double &x, const std::vector<unsigned> &geno,
                             const std::vector<double> &cons) const
{
    x.resize(geno.size() + cons.size());
    std::copy(geno.begin(), geno.end(), x.begin() + m_ncon);
    std::copy(cons.begin(), cons.end(), x.begin());
}

// Pretty (from dsyre) - parenthesis will be unnecessarily used
std::string sr_problem::pretty(const pagmo::vector_double &x) const
{
    std::vector<unsigned> geno;
    std::vector<double> cons;
    std::vector<double> mse;
    std::vector<std::string> retval;

    pagmo2dsyre(geno, cons, x);

    // We compute the idx of the best u in the phenotype
    m_ex.mse(mse, geno, cons, m_points, m_labels);
    auto best_it = std::min_element(mse.begin(), mse.end());
    auto best_idx = std::distance(mse.begin(), best_it);
    // We compute the symbolic representation of all us in the phenotype
    m_ex.sphenotype(retval, geno);

    // We return the best
    return retval[best_idx];
}

// Prettier (from symengine)
std::string sr_problem::prettier(const pagmo::vector_double &x) const
{
    SymEngine::Expression symengine_ex(pretty(x));
    auto retval = fmt::format("{}", symengine_ex);
    return retval;
}

/// Extra info
std::string sr_problem::get_extra_info() const
{
    std::string retval;
    retval += fmt::format("\tNumber of variables: {}\n", m_points[0].size());
    retval += fmt::format("\tNumber of constants: {}\n", m_ncon);
    retval += fmt::format("\tDataset dimension {}\n", m_points.size());
    retval += fmt::format("\tKernels: {}\n", m_kernels);
    return retval;
}
} // namespace dsyre