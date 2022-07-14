#include <algorithm>
#include <assert.h> /* assert */
#include <iostream>
#include <random>
#include <vector>

#include <boost/math/constants/constants.hpp>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <symengine/expression.h>

#include <dsyre/expression.hpp>
#include <dsyre/update_constants.hpp>

constexpr double pi = boost::math::constants::pi<double>();
constexpr double napier = boost::math::constants::e<double>();

using namespace fmt;

inline double P1(const std::vector<double> &x)
{
    return std::pow(x[0], 5) - pi * std::pow(x[0], 3) + x[0];
}

inline double P2(const std::vector<double> &x)
{
    return std::pow(x[0], 5) - pi * std::pow(x[0], 3) + 2 * pi / x[0];
}

inline double P3(const std::vector<double> &x)
{
    return (2.7182 * std::pow(x[0], 5) + std::pow(x[0], 3)) / (x[0] + 1);
}

inline double P4(const std::vector<double> &x)
{
    return std::sin(pi * x[0]) + 1. / x[0];
}

inline double P5(const std::vector<double> &x)
{
    return napier * std::pow(x[0], 5) - pi * std::pow(x[0], 3) + x[0];
}

inline double P6(const std::vector<double> &x)
{
    return (napier * std::pow(x[0], 2) - 1.) / (pi * (x[0] + 2));
}

inline double P7(const std::vector<double> &x)
{
    return std::cos(pi * x[0]) + std::sin(napier * x[0]);
}

inline double P8(const std::vector<double> &x)
{
    return 2.5382 * std::cos(x[3]) + x[0] * x[0] - 0.5;
}

void generate_1d_data(std::vector<std::vector<double>> &xs, std::vector<double> &ys, unsigned N, double lb, double ub)
{
    xs.resize(N);
    ys.resize(N);
    // i must be double
    for (double i = 0.; i < N; ++i) {
        xs[i] = {lb + i / (N - 1) * (ub - lb)};
        ys[i] = P6(xs[i]);
    }
}

void generate_md_data(std::vector<std::vector<double>> &xs, std::vector<double> &ys, unsigned N, unsigned m)
{
    std::random_device rd;  // only used once to initialise (seed) engine
    std::mt19937 rng(rd()); // random-number engine used (Mersenne-Twister in this case)
    std::normal_distribution<double> normal(0., 1.);
    xs.resize(N);
    ys.resize(N);
    for (auto &x : xs) {
        x.resize(m);
        for (auto &item : x) {
            item = 2 * normal(rng);
        }
    }
    for (auto i = 0u; i < N; ++i) {
        ys[i] = P8(xs[i]);
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

    // Generate data
    std::vector<std::vector<double>> xs;
    std::vector<double> ys;
    // generate_md_data(xs, ys, 100u, 5u);
    generate_1d_data(xs, ys, 10u, -2.1, 1.);
    //    Allocate some stuff
    auto length = 20u;
    auto n_var = xs[0].size();
    auto n_con = 1u;
    std::vector<double> mse;
    std::vector<double> predicted_mse(n_con + n_var + length, 0.);
    std::vector<std::vector<double>> grad, hess;

    // The expression system
    // dsyre::expression ex(n_var, n_con, {"sum", "mul", "sin", "cos", "exp", "inv"});
    //dsyre::expression ex(n_var, n_con, {"sum", "mul", "diff", "div", "sin", "cos"});
    dsyre::expression ex(n_var, n_con, {"sum", "mul", "diff", "div"});

    // Run the evolution
    // We run n_trials experiments
    auto ERT = 0u;
    auto n_success = 0u;
    std::vector<unsigned> best_x;
    std::vector<double> best_c;

    std::vector<unsigned> active;

    for (auto j = 0u; j < n_trials; ++j) {
        // We let each run to convergence
        ex.random_genotype(best_x, length);
        best_c = ex.random_constants(-1., 1.);
        auto best_f = ex.fitness(best_x, best_c, xs, ys, mse);
        auto count = 0u;
        ERT++;
        while (count < restart) {
            for (auto i = 0u; i < 4u; ++i) {
                auto new_x = ex.mutation3(best_x, mse, 3 * i + 3);
                ex.remove_nesting(new_x);
                // We now have a new candidate genotype new_x and see what a Newton step could produce.
                // 1 - We compute the mse its gradient and hessian
                ex.ddmse(new_x, best_c, xs, ys, mse, grad, hess);
                // 2 - We compute the new constants based on the best predicted u
                auto new_c = best_c;
                dsyre::update_constants(new_c, mse, grad, hess);
                // 3 - We compute the fitness of the new genotype with the updated constants
                auto new_f = ex.fitness(new_x, new_c, xs, ys, mse);
                count++;
                ERT++;
                // 4 - Reinsertion if better or equal
                if (new_f[0] <= best_f[0]) {
                    best_x = new_x;
                    best_f = new_f;
                    best_c = new_c;
                    // Only if verbosity is > 0
                    if (verbosity > 0) {
                        print("New Best is {} at {} fevals: c value {})\n", best_f, count, best_c);
                    }
                }
            }
            if (best_f[0] < 1e-10) {
                n_success++;
                fmt::print(".");
                fflush(stdout);
                break;
            }
        }
        if (best_f[0] > 1e-10) {
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
    mse = ex.mse(best_x, best_c, xs, ys);
    auto tmp = std::min_element(mse.begin(), mse.end());
    auto idx = std::distance(mse.begin(), tmp);
    print("Best prettied phenotype: {}\n", exs[idx]);

    return 0;
}
