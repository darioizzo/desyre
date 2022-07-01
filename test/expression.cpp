// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <sstream>

#include <dsyre/expression.hpp>

#include "catch.hpp"

using namespace dsyre;

TEST_CASE("construction")
{
    {
        REQUIRE_NOTHROW(expression(1u, 1u, {0u, 1u, 2u, 3u, 4u, 5u}));
    }
    {
        REQUIRE_THROWS_AS(expression(1u, 1u, {0u, 1u, 2u, 3u, 24u, 5u}), std::invalid_argument);
    }
}

TEST_CASE("random_constants")
{
    {
        auto n_con = 1u;
        double lb = -0.1;
        double ub = 3.;
        expression ex(1u, n_con, {0u, 1u, 2u, 3u, 4u, 5u});
        auto values = ex.random_constants(lb, ub);
        REQUIRE(values.size() == n_con);
        REQUIRE(std::all_of(values.begin(), values.end(), [lb](double x){return x > lb;}));
        REQUIRE(std::all_of(values.begin(), values.end(), [ub](double x){return x < ub;}));
    }
    {
        auto n_con = 7u;
        double lb = -103.2;
        double ub = 1e2;
        expression ex(1u, n_con, {0u, 1u, 2u, 3u, 4u, 5u});
        auto values = ex.random_constants(lb, ub);
        REQUIRE(values.size() == n_con);
        REQUIRE(std::all_of(values.begin(), values.end(), [lb](double x){return x > lb;}));
        REQUIRE(std::all_of(values.begin(), values.end(), [ub](double x){return x < ub;}));
    }
}
