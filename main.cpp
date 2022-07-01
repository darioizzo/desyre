#include <algorithm>
#include <assert.h> /* assert */
#include <iostream>
#include <random>
#include <vector>

#include <fmt/ranges.h>
#include <fmt/ostream.h>
#include <symengine/expression.h>

#include <dsyre/expression.hpp>

using namespace fmt;

inline double P1(const std::vector<double> &x)
{
    return std::pow(x[0], 5) - 3.14151617 * std::pow(x[0], 3) + x[0];
}

inline double P2(const std::vector<double> &x)
{
    return std::pow(x[0], 5) - 3.14151617 * std::pow(x[0], 3) + 2 * 3.14151617 / x[0];
}

inline double P3(const std::vector<double> &x)
{
    return (2.7182 * std::pow(x[0], 5) + std::pow(x[0], 3)) / (x[0] + 1);
}

inline double P4(const std::vector<double> &x)
{
    return std::sin(3.14151617 * x[0]) + 1. / x[0];
}

void generate_data(std::vector<std::vector<double>> &xs, std::vector<double> &ys, unsigned N, double lb, double ub)
{
    xs.resize(N);
    ys.resize(N);
    // i must be double
    for (double i = 0.; i < N; ++i) {
        xs[i] = {lb + i / (N - 1) * (ub - lb)};
        ys[i] = P4(xs[i]);
    }
}

using namespace fmt;
// Usage ./main n_trials restart verbosity
int main(int argc, char *argv[])
{
    auto n_trials = std::atoi(argv[1]);
    auto restart = std::atoi(argv[2]);
    auto verbosity = std::atoi(argv[3]);

    std::random_device rd;  // only used once to initialise (seed) engine
    std::mt19937 rng(rd()); // random-number engine used (Mersenne-Twister in this case)

    // Generate P1 data
    std::vector<std::vector<double>> xs;
    std::vector<double> ys;
    generate_data(xs, ys, 10, -1., 1.);

    // Allocate some stuff
    auto length = 20u;
    std::vector<double> mse(length + 2, 0.), dmse(length + 2, 0.), ddmse(length + 2, 0.), predicted_mse(length + 2, 0.);

    // The expression system 1 var 1 constant +,-,*,/, sin, cos
    dsyre::expression ex(1, 1, {0, 1, 2, 3, 4, 5});

    // Run the evolution
    // We run n_trials experiments
    auto ERT = 0u;
    auto n_success = 0u;
    auto best_x = ex.random_genotype(length);

    for (auto j = 0u; j < n_trials; ++j) {
        // We let each run to convergence
        best_x = ex.random_genotype(length);
        auto best_c = ex.random_constants(-10., 10.);
        auto best_f = ex.fitness(best_x, best_c, xs, ys);
        auto count = 0u;
        count++;
        while (count < restart) {
            for (auto i = 0u; i < 4u; ++i) {
                auto new_x = ex.mutation(best_x, 5u);
                auto new_c = best_c;
                // We now have a new candidate genotype new_x and see what a Newton step could produce.
                // 1 - We compute the mse and its derivatives w.r.t. c
                ex.ddmse(new_x, best_c, xs, ys, mse, dmse, ddmse);
                // 2 - We assume a second order model mse = mse + dmse dc + 1/2 ddmse dc^2 and find the best
                // genotype.
                for (decltype(predicted_mse.size()) i = 0u; i < predicted_mse.size(); ++i) {
                    double dc = 0.;
                    if (std::isfinite(dmse[i]) && std::isfinite(ddmse[i]) && ddmse[i] != 0) {
                        dc = -dmse[i] / ddmse[i];
                    }
                    predicted_mse[i] = mse[i] + dmse[i] * dc + 0.5 * ddmse[i] * dc * dc;
                }
                auto tmp = std::min_element(predicted_mse.begin(), predicted_mse.end());
                auto idx = std::distance(predicted_mse.begin(), tmp);
                // 3 - The new constants are those of the best genotype if not nans
                if (std::isfinite(dmse[idx]) && std::isfinite(ddmse[idx]) && ddmse[idx] != 0) {
                    new_c[0] -= dmse[idx] / ddmse[idx];
                }
                // 4 - We compute the fitness of the new genotype with those constants
                auto new_f = ex.fitness(new_x, new_c, xs, ys);
                count++;
                ERT++;
                if (new_f[0] <= best_f[0]) {
                    best_x = new_x;
                    best_f = new_f;
                    best_c = new_c;
                    // Only if verbosity is > 0
                    if (verbosity > 0) {
                        print("New Best is {} at {} fevals: c value {})\n", best_f, count, best_c[0]);
                    }
                }
            }
            if (best_f[0] < 1e-10) {
                n_success++;
                break;
            }
        }
    }
    if (n_success > 0u) {
        print("ERT is {}\n", ERT / n_success);
        print("Successful runs {}\n", n_success);
    } else {
        print("No success, restart less frequently?\n");
    }
    auto final_best = ex.sphenotype(best_x, {"x"}, {"c"});
    print("Best phenotype: {}\n", final_best);
    std::vector<SymEngine::Expression> exs;
    for (auto const &raw : final_best) {
        exs.emplace_back(raw);
    }
   print("Best prettied phenotype: {}\n", exs);
    return 0;
}