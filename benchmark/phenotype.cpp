// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <chrono>
#include <iostream>
#include <random>

#include <fmt/ostream.h>
#include <fmt/ranges.h>

#include <dsyre/expression.hpp>

using namespace dsyre;
using namespace std::chrono;

// In this benchmark

std::vector<std::vector<double>> generate_random_data(unsigned N, unsigned dim, unsigned seed)
{
    std::vector<double> phenotype(dim);
    std::vector<std::vector<double>> retval(N, phenotype);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<> dis(-10., 10.);

    for (auto &v_item : retval) {
        for (auto &item : v_item) {
            item = dis(rng);
        }
    }
    return retval;
}

void perform_speed_test(unsigned n_vars, unsigned n_cons, unsigned length, unsigned N, unsigned seed)
{
    // We construct the expression (default arithmetic kernels)
    expression ex(n_vars, n_cons, {"sum", "diff", "mul", "div", "sin", "cos"});

    // We generte the randomly sampled data points
    auto vars = generate_random_data(N, n_vars, seed);
    auto cons = generate_random_data(N, n_cons, seed-10u);

    // Here is where the phenotype will be written
    std::vector<double> phenotype;

    // We log progress
    fmt::print("{} vars, {} cons, {} length, on {} data points: ", n_vars, n_cons, length, N);
    auto genotype = ex.random_genotype(length);
    auto start = high_resolution_clock::now();
    for (decltype(vars.size()) i = 0u; i < vars.size(); ++i) {
        ex.phenotype(phenotype, genotype, vars[i], cons[i]);
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    fmt::print("{}s\n", duration.count() / 1e6);
}
int main()
{
    unsigned seed = 7898935u;
    fmt::print("We change the length of the genotype (numbers of decompositions u):\n");
    perform_speed_test(1, 1, 100, 1000, seed);
    perform_speed_test(1, 1, 1000, 1000, seed);
    perform_speed_test(1, 1, 10000, 1000, seed);
    fmt::print("We change the number of points:\n");
    perform_speed_test(1, 1, 100, 1000, seed);
    perform_speed_test(1, 1, 100, 10000, seed);
    perform_speed_test(1, 1, 100, 100000, seed);

}