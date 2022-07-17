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

#include <dsyre/sr_problem.hpp>

using namespace dsyre;
using namespace std::chrono;

// In this benchmark

std::vector<std::vector<double>> generate_random_points(unsigned N, unsigned dim, std::mt19937 &rng)
{
    std::vector<double> phenotype(dim);
    std::vector<std::vector<double>> retval(N, phenotype);
    std::uniform_real_distribution<> dis(-10., 10.);

    for (auto &v_item : retval) {
        for (auto &item : v_item) {
            item = dis(rng);
        }
    }
    return retval;
}
std::vector<double> generate_random_labels(unsigned N, std::mt19937 &rng)
{
    std::vector<double> retval(N);
    std::uniform_real_distribution<> dis(-10., 10.);

    for (auto &item : retval) {
        item = dis(rng);
    }
    return retval;
}

void perform_speed_test(unsigned n_vars, unsigned n_cons, unsigned length, unsigned N, unsigned seed)
{
    // The rng
    std::mt19937 rng(12201220u);

    // We generte a radnndom dataset
    auto points = generate_random_points(N, n_vars, rng);
    auto labels = generate_random_labels(N, rng);
    auto cons = generate_random_points(1, n_cons, rng);

    // We instantiate the udp
    sr_problem udp(points, labels, length, {"sum", "diff", "mul", "div", "sin", "cos", "inv"}, n_cons, true);

    // We log progress
    fmt::print("{} vars, {} cons, {} length, on {} data points: ", n_vars, n_cons, length, N);
    std::vector<unsigned> genotype;
    std::vector<double> chromosome(length * 3 + n_cons);
    // We construct an expression just to be able to create random genotypes
    expression ex(n_vars, n_cons, {"sum", "diff", "mul", "div", "sin", "cos", "inv"});
    ex.random_genotype(genotype, length, rng);
    // The genotype is always the same
    std::copy(genotype.begin(), genotype.end(), chromosome.begin() + n_cons);
    std::copy(cons[0].begin(), cons[0].end(), chromosome.begin());
    // One only call to fitness will compute the whole phenotype on all dataset
    // Hence this is comparable to the performances of the benchmark phenotype
    auto start = high_resolution_clock::now();
    udp.fitness(chromosome);
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