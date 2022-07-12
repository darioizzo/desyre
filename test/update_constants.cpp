// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <sstream>
#include <vector>

#include <dsyre/update_constants.hpp>

#include "catch.hpp"

using namespace dsyre;

TEST_CASE("update_constants")
{
    // This tests the first expression as predicted best: depending on both constants
    {
        std::vector<double> cons = {0., 0.};
        std::vector<double> mse = {-100., 0., 0., 0.};
        std::vector<std::vector<double>> grad = {{1., 0., 0., 3.}, {-1., -2., 0., 2.}};
        std::vector<std::vector<double>> hess = {{1., 2., 3., 4.}, {-1., -2, -3., 2.}, {3., 2., 1., -2.}};
        update_constants(cons, mse, grad, hess);
        REQUIRE(cons[0] == Approx(-1.0));
        REQUIRE(cons[1] == Approx(0.0).margin(1e-15));
    }
    // This tests the second expression as predicted best: depending on one constant
    {
        std::vector<double> cons = {0., 0.};
        std::vector<double> mse = {0., -100., 0., 0.};
        std::vector<std::vector<double>> grad = {{1., 0., 0., 3.}, {-1., -2., 0., 2.}};
        std::vector<std::vector<double>> hess = {{1., 2., 3., 4.}, {-1., -2, -3., 2.}, {3., 2., 1., -2.}};
        update_constants(cons, mse, grad, hess);
        REQUIRE(cons[0] == Approx(0.0).margin(1e-15));
        REQUIRE(cons[1] == Approx(1.0));
    }
    // This tests the third expression as predicted best: depending on no constant
    {
        std::vector<double> cons = {0., 0.};
        std::vector<double> mse = {0., 0., -100., 0.};
        std::vector<std::vector<double>> grad = {{1., 0., 0., 3.}, {-1., -2., 0., 2.}};
        std::vector<std::vector<double>> hess = {{1., 2., 3., 4.}, {-1., -2, -3., 2.}, {3., 2., 1., -2.}};
        update_constants(cons, mse, grad, hess);
        REQUIRE(cons[0] == Approx(0.0).margin(1e-15));
        REQUIRE(cons[1] == Approx(0.0).margin(1e-15));
    }
    // This tests the first expression as predicted best: depending again on both constants
    {
        std::vector<double> cons = {0., 0.};
        std::vector<double> mse = {0., 0., 0., -100.};
        std::vector<std::vector<double>> grad = {{1., 0., 0., 3.}, {-1., -2., 0., 2.}};
        std::vector<std::vector<double>> hess = {{1., 2., 3., 4.}, {-1., -2, -3., 2.}, {3., 2., 1., -2.}};
        update_constants(cons, mse, grad, hess);
        REQUIRE(cons[0] == Approx(-0.8333333333333333));
        REQUIRE(cons[1] == Approx(0.16666666666666663));
    }
}