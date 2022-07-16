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

#include <pagmo/types.hpp>

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
    m_expression = expression(1, 0u, m_kernels);
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
    m_expression = expression(m_points[0].size(), m_ncon, m_kernels);
    // We store the number of varables / constants
    m_nvar = points[0].size();
}

pagmo::vector_double::size_type sr_problem::get_nobj() const
{
    return 1u + m_multi_objective;
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
        ub[3 * i + m_ncon + 1] = m_kernels.size() - 1;
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
    return m_expression;
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
} // namespace dsyre