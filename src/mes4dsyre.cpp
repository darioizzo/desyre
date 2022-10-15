// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>

#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>

#include <dsyre/expression.hpp>
#include <dsyre/mes4dsyre.hpp>
#include <dsyre/sr_problem.hpp>
#include <dsyre/update_constants.hpp>

namespace dsyre
{
mes4dsyre::mes4dsyre(unsigned gen, unsigned max_mut, double ftol, unsigned seed)
    : m_gen(gen), m_max_mut(max_mut), m_ftol(ftol), m_e(seed), m_seed(seed), m_verbosity(0u)
{
    if (max_mut == 0u) {
        throw std::invalid_argument("The number of active mutations is zero, it must be at least 1.");
    }
    if (ftol < 0.) {
        throw std::invalid_argument("The ftol is negative, it must be positive or zero.");
    }
}

pagmo::population mes4dsyre::evolve(pagmo::population pop) const
{
    const auto &prob = pop.get_problem();
    auto n_obj = prob.get_nobj();
    auto NP = pop.size();
    auto fevals0 = prob.get_fevals(); // fevals already made
    auto count = 1u;                  // regulates the screen output

    // TODO: We should not use directly the pagmo::problem::extract as otherwise we could not override it in the
    // python bindings. Using a global function, instead, may allow its implementation to be overridden in the
    // bindings.
    auto udp_ptr = prob.extract<sr_problem>();
    // PREAMBLE-------------------------------------------------------------------------------------------------
    // Check whether the problem is suitable for mes4cgp
    // If the UDP in pop is not a symbolic_regression UDP, udp_ptr will be NULL
    if (!udp_ptr) {
        throw std::invalid_argument(prob.get_name() + " does not seem to be a symbolic regression problem. "
                                    + get_name() + " can only be used on problems of the type dsyre::sr_problem ");
    }
    if (n_obj > 1) {
        throw std::invalid_argument(prob.get_name() + " has multiple objectives. " + get_name()
                                    + " can only be used on problems that are single objective.");
    }
    if (NP < 2u) {
        throw std::invalid_argument(get_name() + " needs at least 2 individuals in the population, "
                                    + std::to_string(NP) + " detected");
    }
    // Get out if there is nothing to do.
    if (m_gen == 0u) {
        return pop;
    }
    // ---------------------------------------------------------------------------------------------------------

    // No throws, all valid: we clear the logs
    m_log.clear();
    // We access the inner expression valuable objects
    auto ex = udp_ptr->get_expression();
    auto points = udp_ptr->get_points();
    auto labels = udp_ptr->get_labels();
    // Alias for the number of constants in the expressions.
    auto ncon = prob.get_ncx();
    // We init the best/worst chromosome in the population.
    auto best_idx = pop.best_idx();
    auto worst_idx = pop.worst_idx();
    auto best_x = pop.get_x()[best_idx];
    auto new_x = best_x;
    auto new_f = pop.get_f()[best_idx];
    // Some further allocations (dsyre stuff)
    std::vector<unsigned> best_geno(best_x.size() - ncon), new_geno(best_x.size() - ncon);
    std::vector<double> best_c(ncon), new_c(ncon);
    std::vector<double> mse;
    std::vector<std::vector<double>> grad, hess;

    // Number of mutations will be randomly distributed in 1...m_max)mut
    std::uniform_int_distribution<unsigned> dis(1, m_max_mut);

    // We split the chromosome into integer and float part as to match the API for dsyre expression
    udp_ptr->pagmo2dsyre(best_geno, best_c, best_x);

    // Main loop
    for (decltype(m_gen) gen = 1u; gen <= m_gen; ++gen) {
        // Logs and prints (verbosity modes > 1: a line is added every m_verbosity generations)
        if (m_verbosity > 0u) {
            // Every m_verbosity generations print a log line
            if (gen % m_verbosity == 1u || m_verbosity == 1u) {
                // Every 50 lines print the column names
                if (count % 50u == 1u) {
                    fmt::print("{:<7} {:<7} {:<11} {:<50} \n", "Gen:", "Fevals:", "Best:", "Constants and Formula:");
                }
                auto formula = udp_ptr->prettier(best_x);
                log_single_line(gen - 1, prob.get_fevals() - fevals0, pop.champion_f()[0], formula, best_x, ncon);
                ++count;
            }
        }

        for (decltype(NP) i = 0u; i < NP; ++i) {
            // 1 - We create a mutated genotype
            new_geno = ex.mutate(best_geno, dis(m_e), m_e); // TODO
            ex.remove_nesting(new_geno, m_e);
            // 2 - We compute its mse its gradient and hessian
            ex.ddmse(mse, grad, hess, new_geno, best_c, points, labels);
            // 3 - We update the constant values based on the mse, grad and hess (life long learning)
            new_c = best_c;
            update_constants(new_c, mse, grad, hess);
            // 4 - We compute the new fitness (need to put back into the pagmo format)
            udp_ptr->dsyre2pagmo(new_x, new_geno, new_c);
            new_f = prob.fitness(new_x);
            // 5 - We reinsert (in steady state)
            if (new_f[0] <= pop.get_f()[worst_idx][0]) {
                pop.set_xf(worst_idx, new_x, new_f);
                worst_idx = pop.worst_idx();
                best_geno = new_geno;
                best_c = new_c;
                best_x = new_x;
            }
        }

        // 4 - Exit if ftol is reached
        if (pop.champion_f()[0] <= m_ftol) {
            if (m_verbosity > 0u) {
                auto formula = udp_ptr->prettier(best_x);
                log_single_line(gen, prob.get_fevals() - fevals0, pop.champion_f()[0], formula, best_x, ncon);
                fmt::print("Exit condition -- ftol < {}\n", m_ftol);
            }
            return pop;
        }
    }

    if (m_verbosity > 0u) {
        auto formula = udp_ptr->prettier(best_x);
        log_single_line(m_gen, prob.get_fevals() - fevals0, pop.champion_f()[0], formula, best_x, ncon);
        fmt::print("Exit condition -- max generations = {}", m_gen);
    }
    return pop;
}

void mes4dsyre::set_seed(unsigned seed)
{
    m_e.seed(seed);
    m_seed = seed;
}

unsigned mes4dsyre::get_seed() const
{
    return m_seed;
}

void mes4dsyre::set_verbosity(unsigned level)
{
    m_verbosity = level;
}

unsigned mes4dsyre::get_verbosity() const
{
    return m_verbosity;
}

std::string mes4dsyre::get_name() const
{
    return "M-ES for DSYRE: A memetic Evolutionary Strategy for the dsyre encoding";
}

std::string mes4dsyre::get_extra_info() const
{
    std::string retval;
    retval += fmt::format("\tMaximum number of generations: {}", m_gen);
    retval += fmt::format("\n\tNumber of active mutations: {}", m_max_mut);
    retval += fmt::format("\n\tExit condition of the final loss (ftol): {}", m_ftol);
    retval += fmt::format("\n\tVerbosity: {}", m_verbosity);
    retval += fmt::format("\n\tSeed: {}", m_seed);
    return retval;
}

const mes4dsyre::log_type &mes4dsyre::get_log() const
{
    return m_log;
}

void mes4dsyre::log_single_line(unsigned gen, unsigned long long fevals, double best_f, const std::string &formula,
                                const std::vector<double> &x, unsigned ncon) const
{
    std::vector<double> cons(ncon);
    std::copy(x.begin(), x.begin() + ncon, cons.begin());
    if (formula.length() > 40) {
        fmt::print("{:<7} {:<7} {:.5e} [{:.5e}] {} ...\n", gen, fevals, best_f, fmt::join(cons, ","),
                   formula.substr(0, 40));
    } else {
        fmt::print("{:<7} {:<7} {:.5e} [{:.5e}] {} \n", gen, fevals, best_f, fmt::join(cons, ","), formula);
    }
    m_log.emplace_back(gen, fevals, best_f, cons, formula);
}
} // namespace dsyre