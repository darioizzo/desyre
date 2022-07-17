// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <sstream>
#include <stdexcept>

#include <dsyre/sr_problem.hpp>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <pagmo/problem.hpp>

#include "catch.hpp"

using namespace dsyre;

TEST_CASE("construction")
{
    REQUIRE_NOTHROW(sr_problem({{1., 2}, {3., 4.}}, {1., 2.}, 20u, {"sum", "mul", "inv"}, 1u, false));
    sr_problem udp({{1., 2}, {3., 4.}}, {1., 2.}, 20u, {"sum", "mul", "inv"}, 1u, false);
    REQUIRE_NOTHROW(sr_problem{});
    REQUIRE_NOTHROW(pagmo::problem(udp));

    // Dataset malformed
    REQUIRE_THROWS_AS(sr_problem({{1., 2., 2.}, {3., 4.}}, {1., 2.}, 20u, {"sum", "mul", "inv"}, 1u, false),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(sr_problem({{1., 2.}, {3., 4.}}, {1., 2., 3.}, 20u, {"sum", "mul", "inv"}, 1u, false),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(sr_problem({{1., 2.}, {3., 4.}, {0., -2.}}, {1., 2.}, 20u, {"sum", "mul", "inv"}, 1u, false),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(sr_problem({{}}, {}, 20u, {"sum", "mul", "inv"}, 1u, false), std::invalid_argument);
}

TEST_CASE("fitness")
{
    // Single objective
    {
        sr_problem udp({{1.}, {2.}}, {1., 0.}, 6u, {"sum", "diff", "mul", "div"}, 2u, false);
        std::vector<double> target_fit = {0.5};
        std::vector<double> chromosome = {2.,1., 2, 0, 1, 2, 3, 2};
        REQUIRE(udp.fitness(chromosome) == target_fit);
    }
    // Multi objective
    {
        sr_problem udp({{1.}, {2.}}, {1., 0.}, 6u, {"sum", "diff", "mul", "div"}, 2u, true);
        std::vector<double> target_fit = {0.5, 1};
        std::vector<double> chromosome = {2.,1., 2, 0, 1, 2, 3, 2};
        REQUIRE(udp.fitness(chromosome) == target_fit);
    }
}