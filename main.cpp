#include <algorithm>
#include <assert.h> /* assert */
#include <iostream>
#include <random>
#include <vector>

#include <Eigen/Dense>
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
    return 2.12 * std::cos(x[3]) + x[0] * x[0] - 3.123;
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

// Eigen stores indexes and sizes as signed types, while we
// use STL containers thus sizes and indexes are unsigned. To
// make the conversion as painless as possible this template is provided
// allowing, for example, syntax of the type D(_(i),_(j)) to adress an Eigen matrix
// when i and j are unsigned
template <typename I>
static Eigen::DenseIndex _(I n)
{
    return static_cast<Eigen::DenseIndex>(n);
}

void update_constants(std::vector<double> &cons, const std::vector<double> &mse,
                      const std::vector<std::vector<double>> &grad, const std::vector<std::vector<double>> &hess)
{
    auto n_con = cons.size();
    std::vector<double> predicted_mse(mse.size());
    std::vector<std::vector<unsigned>> active(mse.size());

    // We compute for each idx_u the active constants
    for (decltype(predicted_mse.size()) idx_u = 0u; idx_u < predicted_mse.size(); ++idx_u) {
        for (auto j = 0u; j < n_con; ++j) {
            if (grad[j][idx_u] != 0.) {
                active[idx_u].push_back(j);
            }
        }
    }
    // Eigen stuff
    // Hessian
    Eigen::MatrixXd fullH;
    // (active) Hessian
    Eigen::MatrixXd H;
    // (active) Gradient
    Eigen::MatrixXd G;
    // (active) Constants
    Eigen::MatrixXd C;
    // The increments computed for the various idcx_u
    std::vector<Eigen::MatrixXd> dC(mse.size());
    // Full pivoting LU decomposition to check that H is invertible, find the inverse
    // and see if H is positive definite
    Eigen::FullPivLU<Eigen::MatrixXd> fullpivlu;

    // 1 - We assume a second order model mse = mse + dmse dc + 1/2 ddmse dc^2 and find the best
    // genotype.

    // We loop over all expressions (us) represented in the phenotype
    for (decltype(predicted_mse.size()) idx_u = 0u; idx_u < predicted_mse.size(); ++idx_u) {
        // If the current model idx_u contains a single ephemeral constant we avoid to call the Eigen machinery. This
        // makes the code more readable and results in a few lines rather than tens of.
        if (n_con == 1u) {
            double dc = 0.;
            if (std::isfinite(grad[0][idx_u]) && std::isfinite(hess[0][idx_u]) && hess[0][idx_u] != 0.) {
                dc = -grad[0][idx_u] / hess[0][idx_u];
            }
            predicted_mse[idx_u] = mse[idx_u] + grad[0][idx_u] * dc + 0.5 * hess[0][idx_u] * dc * dc;
        } else {
            // We have at least two potential constants and thus may need to deal with matrices.
            // We find out how many constants are actually in the current expression idx_u
            // collecting the indices of non-zero gradients
            auto n_active = active[idx_u].size();
            if (n_active == 0u) {
                predicted_mse[idx_u] = mse[idx_u];
                continue;
            }
            if (n_active == 1u) {
                // Only one ephemeral constant is in the expression, as above, we avoid to call the
                // Eigen machinery. This makes the code more readable.
                double dc = 0.;
                // For the gradient the index is simply active[idx_u][0]
                auto grad_idx = active[idx_u][0];
                // For the Hessian we must pick the diagonal entry which is 0,2,5,9,14, ....
                auto hess_idx = (grad_idx + 1u);
                hess_idx = (hess_idx + 1) * hess_idx / 2 - 1;
                // Then business as usual!
                if (std::isfinite(grad[grad_idx][idx_u]) && std::isfinite(hess[hess_idx][idx_u])
                    && hess[hess_idx][idx_u] != 0.) {
                    dc = -grad[grad_idx][idx_u] / hess[hess_idx][idx_u];
                }
                predicted_mse[idx_u] = mse[idx_u] + grad[grad_idx][idx_u] * dc + 0.5 * hess[hess_idx][idx_u] * dc * dc;
            } else {
                // The full Hessian
                fullH = Eigen::MatrixXd::Zero(_(n_con), _(n_con));
                // (active) Hessian
                H.setZero(n_active, n_active);
                // (active) Gradient stored in a column vector
                G.setZero(n_active, 1);
                // (active) Constants stored in a column vector
                C.setZero(n_active, 1);
                // (variations for the active constants)
                dC[idx_u].setZero(n_active, 1);
                // Here we assemble the full hessian
                auto jk = 0;
                for (decltype(n_con) j = 0u; j < n_con; ++j) {
                    for (decltype(j) k = 0u; k <= j; ++k) {
                        fullH(_(j), _(k)) = hess[jk][idx_u];
                        if (j != k) {
                            fullH(_(k), _(j)) = hess[jk][idx_u];
                        }
                        jk++;
                    }
                }
                // Here we assemble the reduced Hessian
                // The following nested fors loops copy into the active hessian all cols and rows tht are non zero, i.e.
                // corresponding to eph constants that actually are in the expression
                // Probably this can be coded more efficiently using operations on rows and columns instead.
                auto new_row_idx = 0u;
                for (decltype(n_active) row_idx = 0u; row_idx < n_active; ++row_idx) {
                    decltype(n_active) new_col_idx = 0u;
                    for (auto col_idx = 0u; col_idx < n_active; ++col_idx) {
                        H(_(new_row_idx), _(new_col_idx)) = fullH(_(active[idx_u][row_idx]), _(active[idx_u][col_idx]));
                        new_col_idx++;
                    }
                    new_row_idx++;
                }

                // We now construct the (active) gradient and the (active) constant column vectors.
                for (decltype(n_active) i = 0u; i < n_active; ++i) {
                    G(_(i), 0) = grad[active[idx_u][i]][idx_u];
                    C(_(i), 0) = cons[active[idx_u][i]];
                }

                // And we perform a Newton step
                if (G.array().isFinite().all()) {
                    fullpivlu.compute(H);
                    // NOTE: (see e.g.
                    // https://math.stackexchange.com/questions/621040/compute-eigenvalues-of-a-square-matrix-given-lu-decomposition)
                    if (fullpivlu.isInvertible()) {
                        // H is invertible the inverse is defined
                        // it can however contain infinities in some elements
                        // so we check that all elements of the inverse are finite and only then do a Newton
                        // step
                        auto H_inv = fullpivlu.inverse();
                        if (H_inv.array().isFinite().all()) {
                            dC[idx_u] = -H_inv * G;
                            // compute the quadratic model
                            predicted_mse[idx_u] = mse[idx_u] + (G.transpose() * dC[idx_u])(0, 0)
                                                   + 0.5 * (dC[idx_u].transpose() * H * dC[idx_u])(0, 0);
                        } else {
                            predicted_mse[idx_u] = mse[idx_u];
                        }
                    } else {
                        // should never go here though
                        predicted_mse[idx_u] = mse[idx_u];
                    }
                }
            }
        }
    }
    auto tmp = std::min_element(predicted_mse.begin(), predicted_mse.end());
    auto idx_u_best = std::distance(predicted_mse.begin(), tmp);
    // 3 - The new constants are those of the best genotype if not nans
    if (n_con == 1u) {
        if (std::isfinite(grad[0][idx_u_best]) && std::isfinite(hess[0][idx_u_best]) && hess[0][idx_u_best] != 0) {
            cons[0] -= grad[0][idx_u_best] / hess[0][idx_u_best];
        }
    } else {
        if (active[idx_u_best].size() == 0u) {
            return;
        }
        if (active[idx_u_best].size() == 1u) {
            auto grad_idx = active[idx_u_best][0];
            // For the Hessian we must pick the diagonal entry which is 0,2,5,9,14, ....
            auto hess_idx = (grad_idx + 1u);
            hess_idx = (hess_idx + 1) * hess_idx / 2 - 1;
            cons[grad_idx] -= grad[grad_idx][idx_u_best] / hess[hess_idx][idx_u_best];
        } else {
            for (auto i = 0u; i < active[idx_u_best].size(); ++i) {
                cons[active[idx_u_best][i]] += dC[idx_u_best](_(i));
            }
        }
    }
}

void update_constantsOLD(std::vector<double> &c, const std::vector<double> &mse,
                         const std::vector<std::vector<double>> &grad, const std::vector<std::vector<double>> &hess,
                         std::mt19937 &rng)
{
    auto n_con = c.size();
    std::vector<double> predicted_mse(mse.size());
    std::vector<std::vector<unsigned>> active(mse.size());

    for (decltype(predicted_mse.size()) idx_u = 0u; idx_u < predicted_mse.size(); ++idx_u) {
        for (auto j = 0u; j < n_con; ++j) {
            if (grad[j][idx_u] != 0.) {
                active[idx_u].push_back(j);
            }
        }
    }

    // 1 - We assume a second order model mse = mse + dmse dc + 1/2 ddmse dc^2 and find the best
    // genotype.

    // We loop over all expressions (us) represented in the phenotype
    for (decltype(predicted_mse.size()) idx_u = 0u; idx_u < predicted_mse.size(); ++idx_u) {
        // If the current model idx_u contains a single ephemeral constant we avoid to call the Eigen machinery. This
        // makes the code more readable and results in a few lines rather than tens of.
        if (n_con == 1u) {
            double dc = 0.;
            if (std::isfinite(grad[0][idx_u]) && std::isfinite(hess[0][idx_u]) && hess[0][idx_u] != 0.) {
                dc = -grad[0][idx_u] / hess[0][idx_u];
            }
            predicted_mse[idx_u] = mse[idx_u] + grad[0][idx_u] * dc + 0.5 * hess[0][idx_u] * dc * dc;
        } else {
            // We have at least two potential constants and thus may need to deal with matrices.
            // We find out how many constants are actually in the current expression idx_u
            // collecting the indices of non-zero gradients
            if (active[idx_u].size() == 1u) {
                // Only one ephemeral constant is in the expression, as above, we avoid to call the
                // Eigen machinery. This makes the code more readable.
                double dc = 0.;
                // For the gradient the index is simply active[idx_u][0]
                auto grad_idx = active[idx_u][0];
                // For the Hessian we must pick the diagonal entry which is 0,2,5,9,14, ....
                auto hess_idx = (grad_idx + 1u);
                hess_idx = (hess_idx + 1) * hess_idx / 2 - 1;
                // Then business as usual!
                if (std::isfinite(grad[grad_idx][idx_u]) && std::isfinite(hess[hess_idx][idx_u])
                    && hess[hess_idx][idx_u] != 0.) {
                    dc = -grad[grad_idx][idx_u] / hess[hess_idx][idx_u];
                }
                predicted_mse[idx_u] = mse[idx_u] + grad[grad_idx][idx_u] * dc + 0.5 * hess[hess_idx][idx_u] * dc * dc;
            } else {
                // Only one ephemeral constant is in the expression, as above, we avoid to call the
                // Eigen machinery. This makes the code more readable.
                double dc = 0.;
                // For the gradient the index is simply active[idx_u][0]
                std::uniform_int_distribution<unsigned> dis(0, n_con - 1);
                auto grad_idx = dis(rng);
                // For the Hessian we must pick the diagonal entry which is 0,2,5,9,14, ....
                auto hess_idx = (grad_idx + 1u);
                hess_idx = (hess_idx + 1) * hess_idx / 2 - 1;
                // Then business as usual!
                if (std::isfinite(grad[grad_idx][idx_u]) && std::isfinite(hess[hess_idx][idx_u])
                    && hess[hess_idx][idx_u] != 0.) {
                    dc = -grad[grad_idx][idx_u] / hess[hess_idx][idx_u];
                }
                predicted_mse[idx_u] = mse[idx_u] + grad[grad_idx][idx_u] * dc + 0.5 * hess[hess_idx][idx_u] * dc * dc;
            }
        }
    }
    auto tmp = std::min_element(predicted_mse.begin(), predicted_mse.end());
    auto idx_best = std::distance(predicted_mse.begin(), tmp);
    // 3 - The new constants are those of the best genotype if not nans
    if (n_con == 1u) {
        if (std::isfinite(grad[0][idx_best]) && std::isfinite(hess[0][idx_best]) && hess[0][idx_best] != 0) {
            c[0] -= grad[0][idx_best] / hess[0][idx_best];
        }
    } else {
        if (active[idx_best].size() == 1u) {
            auto grad_idx = active[idx_best][0];
            // For the Hessian we must pick the diagonal entry which is 0,2,5,9,14, ....
            auto hess_idx = (grad_idx + 1u);
            hess_idx = (hess_idx + 1) * hess_idx / 2 - 1;
            c[grad_idx] -= grad[grad_idx][idx_best] / hess[hess_idx][idx_best];
        } else {
        }
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
    // generate_1d_data(xs, ys, 10u, 1., 3.);
    //    Allocate some stuff
    auto length = 20u;
    auto n_var = xs[0].size();
    auto n_con = 2u;
    std::vector<double> mse;
    std::vector<double> predicted_mse(n_con + n_var + length, 0.);
    std::vector<std::vector<double>> grad, hess;

    // The expression system
    dsyre::expression ex(n_var, n_con, {"sum", "mul", "sin", "cos", "exp", "inv"});
    // dsyre::expression ex(n_var, n_con, {"sum", "mul", "diff", "div"});

    // Run the evolution
    // We run n_trials experiments
    auto ERT = 0u;
    auto n_success = 0u;
    std::vector<unsigned> best_x;
    std::vector<double> best_c;

    std::vector<unsigned> active;

    for (auto j = 0u; j < n_trials; ++j) {
        // We let each run to convergence
        best_x = ex.random_genotype(length);
        // best_x =
        // {2,3,0,1,0,0,1,5,7,0,8,9,0,10,6,2,3,0,1,0,0,1,5,7,0,8,9,0,10,6,2,3,0,1,0,0,1,5,7,0,8,9,0,10,6,2,3,0,1,0,0,1,5,7,0,8,9,0,10,6};
        best_c = ex.random_constants(-1., 1.);
        auto best_f = ex.fitness(best_x, best_c, xs, ys, mse);
        auto count = 0u;
        ERT++;
        while (count < restart) {
            for (auto i = 0u; i < 4u; ++i) {
                auto new_x = ex.mutation3(best_x, mse, 3 * i + 3);
                // auto new_x = best_x;
                ex.remove_nesting(new_x);
                // We now have a new candidate genotype new_x and see what a Newton step could produce.
                // 1 - We compute the mse its gradient and hessian
                ex.ddmse(new_x, best_c, xs, ys, mse, grad, hess);
                // 2 - We compute the new constants based on the best predicted u
                auto new_c = best_c;
                // fmt::print("before: {}\n", new_c);
                update_constants(new_c, mse, grad, hess);
                // fmt::print("after: {}\n", new_c);

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

// int main()
//{
//     std::vector<double> cons = {0., 0.};
//     std::vector<double> mse = {0., -10., -2., -100.};
//     std::vector<std::vector<double>> grad = {{1., 0., 0., 3.}, {-1., -2., 0., 2.}};
//     std::vector<std::vector<double>> hess = {{1., 2., 3., 4.}, {-1., -2, -3., 2.}, {3., 2., 1., -2.}};
//     update_constants(cons, mse, grad, hess);
//     fmt::print("new: {}\n", cons);
//
//     return 0;
// }