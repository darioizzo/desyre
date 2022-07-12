// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSYRE_UPDATE_CONSTANTS_HPP
#define DSYRE_UPDATE_CONSTANTS_HPP

#include <Eigen/Dense>
#include <vector>

#include <dsyre/detail/visibility.hpp>


namespace dsyre
{
// Eigen stores indexes and sizes as signed types, while we
// use STL containers thus sizes and indexes are unsigned. To
// make the conversion as painless as possible this template is provided
// allowing, for example, syntax of the type D(_(i),_(j)) to adress an Eigen matrix
// when i and j are unsigned
template <typename I>
static Eigen::DenseIndex _(I n)
{
    return static_cast<Eigen::DenseIndex>(n);
}

DSYRE_DLL_PUBLIC void update_constants(std::vector<double> &constants, const std::vector<double> &function,
                      const std::vector<std::vector<double>> &gradient,
                      const std::vector<std::vector<double>> &hessian);
} // namespace dsyre

#endif
