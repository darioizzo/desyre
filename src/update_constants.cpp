// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <vector>
#include <Eigen/Dense>

#include <dsyre/update_constants.hpp>

namespace dsyre
{
// takes cons and updates it.
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
} // namespace dsyre