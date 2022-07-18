// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <algorithm>
#include <assert.h> /* assert */
#include <iostream>
#include <random>
#include <vector>

#include <boost/math/constants/constants.hpp>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <pagmo/algorithm.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <symengine/expression.h>

#include <dsyre/mes4dsyre.hpp>
#include <dsyre/gym/gym.hpp>
#include <dsyre/sr_problem.hpp>

using namespace fmt;
using namespace dsyre;

// Usage ./main n_trials restart verbosity
int main(int argc, char *argv[])
{
    auto n_trials = std::atoi(argv[1]);
    auto restart = std::atoi(argv[2]);
    auto verbosity = std::atoi(argv[3]);

    std::random_device rd;  // only used once to initialise (seed) engine
    std::mt19937 rng(rd()); // random-number engine used (Mersenne-Twister in this case)

    // Generate data
    std::vector<std::vector<double>> xs;
    std::vector<double> ys;
    // generate_1d_data(xs, ys, 10u, -1., 1.);
    dsyre::gym::generate_P1(xs, ys);
    // Other hyperparameters
    auto length = 20u;
    auto max_mut = 15u;
    auto n_con = 1;
    auto popsize = 4u;
    auto gen = restart / popsize;
    // We instantiate the problem
    sr_problem udp(xs, ys, length, {"sum", "mul", "diff", "div"}, n_con, false);
    // We instantiate the algorithm
    mes4dsyre uda(gen, max_mut, 1e-10, 12345u);
    // Pagmo problem
    pagmo::problem prob(udp);
    // Pagmo algorithm
    pagmo::algorithm algo(uda);
    algo.set_verbosity(verbosity);
    // We run n_trials experiments
    auto ERT = 0u;
    auto n_success = 0u;
    for (auto j = 0u; j < n_trials; ++j) {
        // Pagmo population
        pagmo::population pop(prob, popsize);
        // Evolution!
        pop = algo.evolve(pop);
        ERT += pop.get_problem().get_fevals();
        if (pop.champion_f()[0] < 1e-10) {
            n_success++;
            fmt::print(".");
            fflush(stdout);
        } else {
            fmt::print("x");
            fflush(stdout);
        }
    }
    if (n_success > 0u) {
        print("\nERT is {}\n", ERT / n_success);
        print("Successful runs {}\n", n_success);
    } else {
        print("\nNo success, restart less frequently?\n");
    }
    return 0;
}
