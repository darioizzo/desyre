// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <sstream>

#include <dsyre/expression.hpp>
#include <dsyre/kernels.hpp>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <symengine/expression.h>

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
        CHECK(std::all_of(values.begin(), values.end(), [lb](double x) { return x > lb; }));
        CHECK(std::all_of(values.begin(), values.end(), [ub](double x) { return x < ub; }));
    }
    {
        auto n_con = 7u;
        double lb = -103.2;
        double ub = 1e2;
        expression ex(1u, n_con, {0u, 1u, 2u, 3u, 4u, 5u});
        auto values = ex.random_constants(lb, ub);
        REQUIRE(values.size() == n_con);
        CHECK(std::all_of(values.begin(), values.end(), [lb](double x) { return x > lb; }));
        CHECK(std::all_of(values.begin(), values.end(), [ub](double x) { return x < ub; }));
    }
    // We test for seed control
    {
        auto n_con = 7u;
        double lb = -103.2;
        double ub = 1e2;
        expression ex1(1u, n_con, {0u, 1u, 2u, 3u, 4u, 5u}, 7824323u);
        expression ex2(1u, n_con, {0u, 1u, 2u, 3u, 4u, 5u}, 7824323u);
        auto values1 = ex1.random_constants(lb, ub);
        auto values2 = ex2.random_constants(lb, ub);
        REQUIRE(values1 == values2);
    }
}

TEST_CASE("random_genotype")
{
    {
        auto n_con = 0u;
        auto n_var = 1u;
        unsigned length = 20u;
        std::vector<unsigned> kernels = {0u, 1u, 2u, 3u, 4u, 5u};
        expression ex(n_var, n_con, kernels);
        auto geno = ex.random_genotype(length);
        REQUIRE(geno.size() == length * 3u);
        for (decltype(geno.size()) i = 0u; i < geno.size(); ++i) {
            if (i % 3 == 0u) {
                REQUIRE(std::any_of(kernels.begin(), kernels.end(),
                                    [i, &geno](unsigned kid) { return kid < geno.size(); }));
            } else {
                REQUIRE(geno[i] < i / 3 + n_con + n_var);
            }
        }
    }
    {
        auto n_con = 3u;
        auto n_var = 2u;
        unsigned length = 20u;
        std::vector<unsigned> kernels = {2u, 4u, 1u, 0u, 2u, 5u};
        expression ex(n_var, n_con, kernels);
        auto geno = ex.random_genotype(length);
        REQUIRE(geno.size() == length * 3u);
        for (decltype(geno.size()) i = 0u; i < geno.size(); ++i) {
            if (i % 3 == 0u) {
                REQUIRE(std::any_of(kernels.begin(), kernels.end(),
                                    [i, &geno](unsigned kid) { return kid < geno.size(); }));
            } else {
                REQUIRE(geno[i] < i / 3 + n_con + n_var);
            }
        }
    }
}

TEST_CASE("phenotype")
{
    {
        auto n_con = 1u;
        auto n_var = 1u;
        std::vector<unsigned> kernels = {0u, 1u, 2u, 3u, 5u};
        expression ex(n_var, n_con, kernels);
        // Manually assemble (x + sin(cx)) / x
        std::vector<unsigned> geno = {2, 0, 1, 4, 2, 0, 0, 0, 3, 3, 4, 0};
        std::vector<std::string> objective = {"x", "c", "(x*c)", "sin((x*c))", "(x+sin((x*c)))", "((x+sin((x*c)))/x)"};
        auto sphen = ex.sphenotype(geno, {"x"}, {"c"});

        for (decltype(sphen.size()) i = 0u; i < sphen.size(); ++i) {
            REQUIRE(sphen[i] == objective[i]);
        }
    }
}