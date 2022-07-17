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
#include <dsyre/sr_problem.hpp>

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
        ys[i] = P1(xs[i]);
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
    generate_1d_data(xs, ys, 10u, 1., 3.);
    // Other hyperparameters
    auto length = 20u;
    auto max_mut = 15u;
    auto n_con = 1.;
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
