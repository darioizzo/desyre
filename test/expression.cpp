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
        std::vector<unsigned> geno;
        ex.random_genotype(geno, length);
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
        std::vector<unsigned> geno;
        ex.random_genotype(geno, length);
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

TEST_CASE("remove_nesting")
{
    {
        auto n_con = 1u;
        auto n_var = 1u;
        std::vector<std::string> kernels = {"sum", "mul", "sin", "inv", "cos"};
        expression ex(n_var, n_con, kernels);
        std::vector<unsigned> genotype_not_nested = {1, 0, 1, 4, 1, 1, 3, 1, 2};
        std::vector<unsigned> genotype_nested = {1, 0, 1, 4, 1, 1, 4, 3, 1};
        std::vector<unsigned> tmp, tmp2;
        // 1 - a non nested genotype remains unchanged
        tmp = genotype_not_nested;
        ex.remove_nesting(tmp);
        REQUIRE(tmp == genotype_not_nested);
        // 2 - a nested genotype is changed
        tmp = genotype_nested;
        ex.remove_nesting(tmp);
        REQUIRE(tmp != genotype_nested);
        // 3 - a changed nested genotype is no longer nested
        tmp2 = tmp;
        ex.remove_nesting(tmp);
        REQUIRE(tmp == tmp2);
    }
    {
        auto n_con = 3u;
        auto n_var = 2u;
        std::vector<std::string> kernels = {"sum", "mul", "sin", "inv", "cos"};
        expression ex(n_var, n_con, kernels);
        std::vector<unsigned> genotype_not_nested = {1, 0, 1, 4, 1, 1, 3, 1, 2};
        std::vector<unsigned> genotype_nested = {1, 0, 1, 4, 1, 1, 4, 3, 1, 3, 5, 6, 4, 6, 2};
        std::vector<unsigned> tmp, tmp2;
        // 1 - a non nested genotype remains unchanged
        tmp = genotype_not_nested;
        ex.remove_nesting(tmp);
        REQUIRE(tmp == genotype_not_nested);
        // 2 - a nested genotype is changed
        tmp = genotype_nested;
        ex.remove_nesting(tmp);
        REQUIRE(tmp != genotype_nested);
        // 3 - a changed nested genotype is no longer nested
        tmp2 = tmp;
        ex.remove_nesting(tmp);
        REQUIRE(tmp == tmp2);
    }
}

TEST_CASE("phenotype")
{
    // Test throws
    {
        auto n_con = 1u;
        auto n_var = 1u;
        std::vector<std::string> kernels = {"sum", "diff", "mul", "div", "sin"};
        expression ex(n_var, n_con, kernels);
        // Manually assemble[x, c, cx, sin(cx), x + sin(cx), (x + sin(cx)) / x]
        std::vector<unsigned> geno = {2, 0, 1, 4, 2, 0, 0, 0, 3, 3, 4, 0};
        std::vector<double> phen;
        REQUIRE_THROWS_AS(ex.phenotype(phen, geno, {0.2, 1.2, 3.4}, {0.1}), std::invalid_argument);
        REQUIRE_THROWS_AS(ex.phenotype(phen, geno, {0.0}, {0.2, 1.2, 3.4}), std::invalid_argument);
    }
    {
        auto n_con = 1u;
        auto n_var = 1u;
        std::vector<std::string> kernels = {"sum", "diff", "mul", "div", "sin"};
        expression ex(n_var, n_con, kernels);
        // Manually assemble [x, c, cx, sin(cx), x + sin(cx), (x + sin(cx)) / x]
        std::vector<unsigned> geno = {2, 0, 1, 4, 2, 0, 0, 0, 3, 3, 4, 0};
        double x = 3.1232133;
        double c = -0.21312123;
        std::vector<double> objective = {x, c, x * c, std::sin(x * c), x + std::sin(x * c), (x + std::sin(x * c)) / x};
        std::vector<double> phen;
        ex.phenotype(phen, geno, {x}, {c});
        REQUIRE(phen == objective);
    }
}

TEST_CASE("phenotype_and_complexity")
{
    // Test throws
    {
        auto n_con = 1u;
        auto n_var = 1u;
        std::vector<std::string> kernels = {"sum", "diff", "mul", "div", "sin"};
        expression ex(n_var, n_con, kernels);
        // Manually assemble[x, c, cx, sin(cx), x + sin(cx), (x + sin(cx)) / x]
        std::vector<unsigned> geno = {2, 0, 1, 4, 2, 0, 0, 0, 3, 3, 4, 0}, complexity;
        std::vector<double> phen;
        REQUIRE_THROWS_AS(ex.phenotype_and_complexity(phen, complexity, geno, {0.2, 1.2, 3.4}, {0.1}),
                          std::invalid_argument);
        REQUIRE_THROWS_AS(ex.phenotype_and_complexity(phen, complexity, geno, {0.0}, {0.2, 1.2, 3.4}),
                          std::invalid_argument);
    }
    {
        auto n_con = 1u;
        auto n_var = 1u;
        std::vector<std::string> kernels = {"sum", "diff", "mul", "div", "sin"};
        expression ex(n_var, n_con, kernels);
        // Manually assemble [x, c, cx, sin(cx), x + sin(cx), (x + sin(cx)) / x]
        std::vector<unsigned> geno = {2, 0, 1, 4, 2, 0, 0, 0, 3, 3, 4, 0}, complexity;
        double x = 3.1232133;
        double c = -0.21312123;
        std::vector<double> target_phenotype
            = {x, c, x * c, std::sin(x * c), x + std::sin(x * c), (x + std::sin(x * c)) / x};
        std::vector<unsigned> target_complexity = {1, 1, 3, 4, 6, 8};
        std::vector<double> phen;
        ex.phenotype_and_complexity(phen, complexity, geno, {x}, {c});
        REQUIRE(phen == target_phenotype);
        REQUIRE(complexity == target_complexity);
    }
    {
        auto n_con = 1u;
        auto n_var = 1u;
        std::vector<std::string> kernels = {"sum", "diff", "mul", "div", "sin"};
        expression ex(n_var, n_con, kernels);
        std::vector<unsigned> geno
            = {0, 1, 0, 1, 0, 2, 3, 0, 3, 4, 3, 1, 0, 4, 5, 2, 2, 3, 0, 7, 2, 0, 2, 7, 4, 7, 3, 0, 10, 3},
            complexity;
        double x = 3.1232133;
        double c = -0.21312123;
        std::vector<double> target_phenotype
            = {3.1232133, -0.21312123, 2.91009207, 0.21312123000000005, 14.654632483117705, 0.21151153890618252, 14.866144022023887, 0.6202024013716463, 3.530294471371646, 3.530294471371646, 0.5811998787507274, 0.7943211087507275};
        std::vector<unsigned> target_complexity = { 1, 1, 3, 5, 7, 6, 14, 9, 13, 13, 10, 16 };
        std::vector<double> phen;
        ex.phenotype_and_complexity(phen, complexity, geno, {x}, {c});
        REQUIRE_THAT(phen, Catch::Approx(target_phenotype).margin(1e-8));
        REQUIRE(complexity == target_complexity);
    }
}

TEST_CASE("sphenotype")
{
    {
        auto n_con = 1u;
        auto n_var = 1u;
        std::vector<std::string> kernels = {"sum", "diff", "mul", "div", "sin"};
        expression ex(n_var, n_con, kernels);
        // Manually assemble [x, c, cx, sin(cx), x + sin(cx), (x + sin(cx)) / x]
        std::vector<unsigned> geno = {2, 0, 1, 4, 2, 0, 0, 0, 3, 3, 4, 0};
        std::vector<std::string> phen;
        REQUIRE_THROWS_AS(ex.sphenotype(phen, geno, {"x", "y"}, {"c"}), std::invalid_argument);
        REQUIRE_THROWS_AS(ex.sphenotype(phen, geno, {"x"}, {"c1", "c2"}), std::invalid_argument);
    }
    {
        auto n_con = 1u;
        auto n_var = 1u;
        std::vector<std::string> kernels = {"sum", "diff", "mul", "div", "sin"};
        expression ex(n_var, n_con, kernels);
        // Manually assemble[x, c, cx, sin(cx), x + sin(cx), (x + sin(cx)) / x]
        std::vector<unsigned> geno = {2, 0, 1, 4, 2, 0, 0, 0, 3, 3, 4, 0};
        std::vector<std::string> objective = {"x", "c", "(x*c)", "sin((x*c))", "(x+sin((x*c)))", "((x+sin((x*c)))/x)"};
        std::vector<std::string> sphen;
        ex.sphenotype(sphen, geno, {"x"}, {"c"});
        REQUIRE(sphen == objective);
    }
}

TEST_CASE("dphenotype")
{
    // Test throws
    {
        auto n_con = 1u;
        auto n_var = 1u;
        std::vector<std::string> kernels = {"sum", "diff", "mul", "div", "sin"};
        expression ex(n_var, n_con, kernels);
        // Manually assemble[x, c, cx, sin(cx), x + sin(cx), (x + sin(cx)) / x]
        std::vector<unsigned> geno = {2, 0, 1, 4, 2, 0, 0, 0, 3, 3, 4, 0};
        std::vector<double> phen, dphen;
        ex.phenotype(phen, geno, {0.12}, {0.3});
        // derivative index is too high!
        REQUIRE_THROWS_AS(ex.dphenotype(dphen, geno, phen, 34u), std::invalid_argument);
        // wrong phenotype-genotype sizes
        REQUIRE_THROWS_AS(ex.dphenotype(dphen, geno, {0.3, 0.4}, 1u), std::invalid_argument);
    }
    {
        auto n_con = 1u;
        auto n_var = 1u;
        std::vector<std::string> kernels = {"sum", "diff", "mul", "div", "sin"};
        expression ex(n_var, n_con, kernels);
        // Manually assemble [x, c, cx, sin(cx), x + sin(cx), (x + sin(cx)) / x]
        std::vector<unsigned> geno = {2, 0, 1, 4, 2, 0, 0, 0, 3, 3, 4, 0};
        double x = 3.1232133;
        double c = -0.21312123;
        std::vector<double> phen, dphen;
        ex.phenotype(phen, geno, {x}, {c});
        // Test derivation with respect to x
        {
            std::vector<double> objective = {1.,
                                             0.,
                                             c,
                                             std::cos(x * c) * c,
                                             1 + std::cos(x * c) * c,
                                             ((1 + std::cos(x * c) * c) * x - (x + std::sin(c * x))) / x / x};
            ex.dphenotype(dphen, geno, phen, 0u);
            REQUIRE_THAT(dphen, Catch::Approx(objective).margin(1e-14));
        }
        // Test derivation with respect to c
        {
            std::vector<double> objective = {0., 1., x, std::cos(x * c) * x, std::cos(x * c) * x, std::cos(x * c)};
            ex.dphenotype(dphen, geno, phen, 1u);
            REQUIRE_THAT(dphen, Catch::Approx(objective).margin(1e-14));
        }
    }
}

TEST_CASE("ddphenotype")
{
    // Test throws
    {
        auto n_con = 1u;
        auto n_var = 1u;
        std::vector<std::string> kernels = {"sum", "diff", "mul", "div", "sin"};
        expression ex(n_var, n_con, kernels);
        // Manually assemble[x, c, cx, sin(cx), x + sin(cx), (x + sin(cx)) / x]
        std::vector<unsigned> geno = {2, 0, 1, 4, 2, 0, 0, 0, 3, 3, 4, 0};
        std::vector<double> phen, dphen, ddphen;
        ex.phenotype(phen, geno, {0.12}, {0.3});
        // derivate w.r.t. x
        ex.dphenotype(dphen, geno, phen, 0u);
        // inconsistent dphenotype dimensions
        REQUIRE_THROWS_AS(ex.ddphenotype(ddphen, geno, phen, dphen, {1., 2., 3.}), std::invalid_argument);
        // inconsistent phenotype dphenotype dimensions
        REQUIRE_THROWS_AS(ex.ddphenotype(ddphen, geno, {1.2, 3., 2.}, dphen, dphen), std::invalid_argument);
        // inconsistent phenotype genotype  dimensions
        REQUIRE_THROWS_AS(ex.ddphenotype(ddphen, geno, {1.2, 3., 2.}, {1.2, 3., 2.}, {1.2, 3., 2.}),
                          std::invalid_argument);
    }
    {
        auto n_con = 1u;
        auto n_var = 1u;
        std::vector<std::string> kernels = {"sum", "diff", "mul", "div", "sin"};
        expression ex(n_var, n_con, kernels);
        // Manually assemble [x, c, cx, sin(cx), x + sin(cx), (x + sin(cx)) / x]
        std::vector<unsigned> geno = {2, 0, 1, 4, 2, 0, 0, 0, 3, 3, 4, 0};
        double x = 3.1232133;
        double c = -0.21312123;
        std::vector<double> phen, dphen0, dphen1, ddphen;
        ex.phenotype(phen, geno, {x}, {c});
        ex.dphenotype(dphen0, geno, phen, 0u);
        ex.dphenotype(dphen1, geno, phen, 1u);

        // Test derivation with respect to x2
        {
            std::vector<double> objective
                = {0.,
                   0.,
                   0,
                   -std::sin(x * c) * c * c,
                   -std::sin(x * c) * c * c,
                   ((2 - c * c * x * x) * std::sin(c * x) - 2 * c * x * std::cos(c * x)) / x / x / x};
            ex.ddphenotype(ddphen, geno, phen, dphen0, dphen0);
            REQUIRE_THAT(ddphen, Catch::Approx(objective).margin(1e-14));
        }
        // Test derivation with respect to c2
        {
            std::vector<double> objective
                = {0., 0., 0., -std::sin(x * c) * x * x, -std::sin(x * c) * x * x, -std::sin(x * c) * x};
            ex.ddphenotype(ddphen, geno, phen, dphen1, dphen1);
            REQUIRE_THAT(ddphen, Catch::Approx(objective).margin(1e-14));
        }
        // Test derivation with respect to c x
        {
            std::vector<double> objective = {0.,
                                             0.,
                                             1.,
                                             std::cos(x * c) - std::sin(x * c) * x * c,
                                             std::cos(x * c) - std::sin(x * c) * x * c,
                                             -std::sin(x * c) * c};
            ex.ddphenotype(ddphen, geno, phen, dphen1, dphen0);
            REQUIRE_THAT(ddphen, Catch::Approx(objective).margin(1e-14));
        }
        // Test derivation with respect to x c
        {
            std::vector<double> objective = {0.,
                                             0.,
                                             1.,
                                             std::cos(x * c) - std::sin(x * c) * x * c,
                                             std::cos(x * c) - std::sin(x * c) * x * c,
                                             -std::sin(x * c) * c};
            ex.ddphenotype(ddphen, geno, phen, dphen0, dphen1);
            REQUIRE_THAT(ddphen, Catch::Approx(objective).margin(1e-14));
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
        std::vector<std::vector<double>> xs = {{1.}};
        std::vector<double> ys = {1.};
        std::vector<double> cons = {2.};
        // Return values
        std::vector<double> mse;
        std::vector<std::vector<double>> dmse, ddmse;
        // Call
        ex.ddmse(geno, cons, xs, ys, mse, dmse, ddmse);
        // Test
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
        // Test
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