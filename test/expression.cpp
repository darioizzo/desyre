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
        REQUIRE_NOTHROW(expression(1u, 1u));
        REQUIRE_NOTHROW(expression(1u, 1u, {"sum", "exp", "sin", "diff"}));
    }
    {
        REQUIRE_THROWS_AS(expression(1u, 1u, {"diff", "cos", "div", "gigietta", "mul"}), std::invalid_argument);
    }
}

TEST_CASE("random_constants")
{
    {
        auto n_con = 1u;
        double lb = -0.1;
        double ub = 3.;
        expression ex(1u, n_con);
        auto values = ex.random_constants(lb, ub);
        REQUIRE(values.size() == n_con);
        CHECK(std::all_of(values.begin(), values.end(), [lb](double x) { return x > lb; }));
        CHECK(std::all_of(values.begin(), values.end(), [ub](double x) { return x < ub; }));
    }
    {
        auto n_con = 7u;
        double lb = -103.2;
        double ub = 1e2;
        expression ex(1u, n_con);
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
        expression ex1(1u, n_con, {"sum", "diff", "mul", "div"}, 7824323u);
        expression ex2(1u, n_con, {"sum", "diff", "mul", "div"}, 7824323u);
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
        std::vector<std::string> kernels = {"sum", "diff", "mul", "div", "sin", "cos"};
        expression ex(n_var, n_con, kernels);
        auto geno = ex.random_genotype(length);
        REQUIRE(geno.size() == length * 3u);
        for (decltype(geno.size()) i = 0u; i < geno.size(); ++i) {
            if (i % 3 == 0u) {
                REQUIRE(std::any_of(ex.get_kernels_idx().begin(), ex.get_kernels_idx().end(),
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
        std::vector<std::string> kernels = {"mul", "sin", "cos", "div", "diff", "sum"};
        expression ex(n_var, n_con, kernels);
        auto geno = ex.random_genotype(length);
        REQUIRE(geno.size() == length * 3u);
        for (decltype(geno.size()) i = 0u; i < geno.size(); ++i) {
            if (i % 3 == 0u) {
                REQUIRE(std::any_of(ex.get_kernels_idx().begin(), ex.get_kernels_idx().end(),
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
        std::vector<std::string> kernels = {"sum", "diff", "mul", "div", "sin"};
        expression ex(n_var, n_con, kernels);
        // Manually assemble (x + sin(cx)) / x
        std::vector<unsigned> geno = {2, 0, 1, 4, 2, 0, 0, 0, 3, 3, 4, 0};
        std::vector<std::string> objective = {"x", "c", "(x*c)", "sin((x*c))", "(x+sin((x*c)))", "((x+sin((x*c)))/x)"};
        std::vector<std::string> sphen;
        ex.sphenotype(sphen, geno, {"x"}, {"c"});

        for (decltype(sphen.size()) i = 0u; i < sphen.size(); ++i) {
            REQUIRE(sphen[i] == objective[i]);
        }
    }
}

TEST_CASE("ddmse")
{
    {
        auto n_con = 1u;
        auto n_var = 1u;
        std::vector<std::string> kernels = {"sum", "diff", "mul", "div"};
        expression ex(n_var, n_con, kernels);
        // Manually assemble [x, c, xc]
        std::vector<unsigned> geno = {2, 0, 1};
        // Manually create a dataset
        std::vector<std::vector<double> > xs = {{1.}};
        std::vector<double> ys = {1.};
        std::vector<double> cons = {2.};
        // Return values
        std::vector<double> mse;
        std::vector<std::vector<double>> dmse, ddmse;
        // Call
        ex.ddmse(geno, cons, xs, ys, mse, dmse, ddmse);
        REQUIRE(mse == std::vector<double>{0, 1, 1});
        REQUIRE(dmse[0] == std::vector<double>{0, 2, 2});
        REQUIRE(ddmse[0] == std::vector<double>{0, 2, 2});
    }
    {
        auto n_con = 1u;
        auto n_var = 1u;
        std::vector<std::string> kernels = {"sum", "diff", "mul", "div"};
        expression ex(n_var, n_con, kernels);
        // Manually assemble [x, c, xc]
        std::vector<unsigned> geno = {2, 0, 1};
        // Manually create a dataset
        std::vector<std::vector<double>> xs = {{1.}, {2.}};
        std::vector<double> ys = {1., 0.};
        std::vector<double> cons = {2.};
        // Return values
        std::vector<double> mse;
        std::vector<std::vector<double>> dmse, ddmse;
        // Call
        ex.ddmse(geno, cons, xs, ys, mse, dmse, ddmse);
        REQUIRE(mse == std::vector<double>{2, 2.5, 8.5});
        REQUIRE(dmse[0] == std::vector<double>{0, 3, 9});
        REQUIRE(ddmse[0] == std::vector<double>{0, 2, 5});
    }
    {
        auto n_con = 2u;
        auto n_var = 1u;
        std::vector<std::string> kernels = {"sum", "diff", "mul", "div"};
        expression ex(n_var, n_con, kernels);
        // Manually assemble [x, c0, c1, xc0, xc0c1]
        std::vector<unsigned> geno = {2, 0, 1, 2, 3, 2};
        // Manually create a dataset
        std::vector<std::vector<double>> xs = {{1.}, {2.}};
        std::vector<double> ys = {1., 0.};
        std::vector<double> cons = {2., 1.};
        // Return values
        std::vector<double> mse;
        std::vector<std::vector<double>> grad, hess;
        // Call
        ex.ddmse(geno, cons, xs, ys, mse, grad, hess);
        // Test
        REQUIRE(mse == std::vector<double>{2, 2.5, 0.5, 8.5, 8.5});

        REQUIRE(grad[0] == std::vector<double>{0, 3, 0, 9, 9});
        REQUIRE(grad[1] == std::vector<double>{0, 0, 1, 0, 18});

        REQUIRE(hess[0] == std::vector<double>{0, 2, 0, 5, 5});
        REQUIRE(hess[1] == std::vector<double>{0, 0, 0, 0, 19});
        REQUIRE(hess[2] == std::vector<double>{0, 0, 2, 0, 20});
    }
}