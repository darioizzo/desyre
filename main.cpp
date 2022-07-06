#include <algorithm>
#include <assert.h> /* assert */
#include <iostream>
#include <random>
#include <vector>

#include <fmt/ostream.h>
#include <fmt/ranges.h>
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

inline double P5(const std::vector<double> &x)
{
    return 2 * std::cos(x[3]) + x[0] * x[0] - 3.123;
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
        ys[i] = P5(xs[i]);
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
    generate_md_data(xs, ys, 100u, 5u);
    //generate_1d_data(xs, ys, 10u, 1., 3.);
    //   Allocate some stuff
    auto length = 20u;
    auto n_var = xs[0].size();
    auto n_con = 2u;
    std::vector<double> mse(length + n_var + n_con, 0.), dmse(length + n_var + n_con, 0.),
        ddmse(length + n_var + n_con, 0.), predicted_mse(length + n_var + n_con, 0.);

    // The expression system 1 var 1 constant +,-,*,/, sin, cos
    dsyre::expression ex(n_var, n_con, {"sum", "mul", "sin", "cos", "exp", "inv"});
    //dsyre::expression ex(n_var, n_con, {"sum", "mul", "diff", "div"});

    // Run the evolution
    // We run n_trials experiments
    auto ERT = 0u;
    auto n_success = 0u;
    std::vector<unsigned> best_x;
    std::vector<double> best_c;

    for (auto j = 0u; j < n_trials; ++j) {
        // We let each run to convergence
        best_x = ex.random_genotype(length);
        // best_x =
        // {2,3,0,1,0,0,1,5,7,0,8,9,0,10,6,2,3,0,1,0,0,1,5,7,0,8,9,0,10,6,2,3,0,1,0,0,1,5,7,0,8,9,0,10,6,2,3,0,1,0,0,1,5,7,0,8,9,0,10,6};
        best_c = ex.random_constants(-1., 1.);
        auto best_f = ex.fitness(best_x, best_c, xs, ys);
        auto count = 0u;
        ERT++;
        while (count < restart) {
            for (auto i = 0u; i < 4u; ++i) {
                auto new_x = ex.mutation2(best_x, 3 * i + 3);
                // auto new_x = best_x;
                ex.remove_nesting(new_x);
                auto new_c = best_c;
                // We now have a new candidate genotype new_x and see what a Newton step could produce.
                // 1 - We compute the mse and its derivatives w.r.t. ca random constant
                std::uniform_int_distribution<unsigned> dis(n_var, n_con + n_var - 1);
                auto c_idx = dis(rng);
                ex.ddmse(new_x, best_c, c_idx, xs, ys, mse, dmse, ddmse);
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
                    new_c[c_idx - n_var] -= dmse[idx] / ddmse[idx];
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
                        print("New Best is {} at {} fevals: c value {})\n", best_f, count, best_c);
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