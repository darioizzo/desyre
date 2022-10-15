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

#include <dsyre/expression.hpp>
#include <dsyre/formatters.hpp>
#include <dsyre/gym/gym.hpp>
#include <dsyre/symengine.hpp>
#include <dsyre/update_constants.hpp>

using namespace fmt;
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
    dsyre::gym::generate_P1(xs, ys);
    // Allocate some stuff
    auto length = 20u;
    auto n_var = xs[0].size();
    auto n_con = 1u;
    std::vector<double> mse;
    std::vector<double> predicted_mse(n_con + n_var + length, 0.);
    std::vector<std::vector<double>> grad, hess;

    // The expression system
    // dsyre::expression ex(n_var, n_con, {"sum", "mul", "sin", "cos", "exp", "inv"});
    // dsyre::expression ex(n_var, n_con, {"sum", "mul", "diff", "div", "sin", "cos"});
    dsyre::expression ex(n_var, n_con, {"sum", "mul", "diff", "div"});

    // Run the evolution
    // We run n_trials experiments
    auto ERT = 0u;
    auto n_success = 0u;
    std::vector<unsigned> best_x;
    std::vector<double> best_c;
    double best_f, new_f;

    std::vector<unsigned> active;

    for (auto j = 0u; j < n_trials; ++j) {
        // We let each run to convergence
        ex.random_genotype(best_x, length, rng);
        best_c = ex.random_constants(-1., 1., rng);
        // We compute the initial fitness
        ex.mse(mse, best_x, best_c, xs, ys);
        best_f = *std::min_element(mse.begin(), mse.end());
        // Iterations
        auto count = 0u;
        ERT++;
        while (count < restart) {
            for (auto i = 0u; i < 4u; ++i) {
                auto new_x = ex.mutate(best_x, 3 * i + 3, rng);
                ex.remove_nesting(new_x, rng);
                // We now have a new candidate genotype new_x and see what a Newton step could produce.
                // 1 - We compute the mse its gradient and hessian
                ex.ddmse(mse, grad, hess, new_x, best_c, xs, ys);
                // 2 - We compute the new constants based on the best predicted u
                auto new_c = best_c;
                dsyre::update_constants(new_c, mse, grad, hess);
                // 3 - We compute the fitness of the new genotype with the updated constants
                ex.mse(mse, new_x, new_c, xs, ys);
                new_f = *std::min_element(mse.begin(), mse.end());
                count++;
                ERT++;
                // 4 - Reinsertion if better or equal
                if (new_f <= best_f) {
                    best_x = new_x;
                    best_f = new_f;
                    best_c = new_c;
                    // Only if verbosity is > 0
                    if (verbosity > 0) {
                        print("New Best is {} at {} fevals: c value {})\n", best_f, count, best_c);
                    }
                }
            }
            if (best_f < 1e-10) {
                n_success++;
                fmt::print(".");
                fflush(stdout);
                break;
            }
        }
        if (best_f > 1e-10) {
            fmt::print("x");
        }
        fflush(stdout);
    }
    if (n_success > 0u) {
        print("\nERT is {}\n", ERT / n_success);
        print("Successful runs {}\n", n_success);
    } else {
        print("\nNo success, restart less frequently?\n");
    }
    std::vector<std::string> final_best;
    ex.sphenotype(final_best, best_x);
    print("Best phenotype: {}\n", final_best);
    std::vector<SymEngine::Expression> exs;
    for (auto const &raw : final_best) {
        exs.emplace_back(raw);
    }
    print("Best prettied phenotype: {}\n", exs);

    std::vector<double> phen;
    ex.mse(mse, best_x, best_c, xs, ys);
    auto tmp = std::min_element(mse.begin(), mse.end());
    auto idx = std::distance(mse.begin(), tmp);
    print("Best prettied phenotype: {:e}\n", exs[idx]);
    print("Best prettied phenotype: {:v}\n", exs[idx]);

    return 0;
}
