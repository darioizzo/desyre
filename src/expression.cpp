// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <algorithm>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include <dsyre/detail/visibility.hpp>
#include <dsyre/expression.hpp>
#include <dsyre/kernels.hpp>

namespace dsyre
{

using ukernel_f_ptr = double (*)(double);
using pkernel_f_ptr = std::string (*)(std::string);

// Global array of function pointers for the unary kernels
ukernel_f_ptr ukernel_list[] = {inv, cos, sin, exp};
ukernel_f_ptr dukernel_list[] = {dinv, dcos, dsin, dexp};
ukernel_f_ptr ddukernel_list[] = {ddinv, ddcos, ddsin, ddexp};
pkernel_f_ptr pkernel_list[] = {pinv, pcos, psin, pexp};

// Number of binary operations (+,-,*,/)
unsigned n_binary = 4u;
unsigned n_unary = std::size(ukernel_list);

// Global map between kernel names and an unsigned (must correspond to the order in the global arrays as its used
// in the switch cases.
std::unordered_map<std::string, unsigned> kernel_map{{"sum", 0},
                                                     {"diff", 1},
                                                     {"mul", 2},
                                                     {"div", 3},
                                                     {"inv", n_binary},
                                                     {"cos", n_binary + 1},
                                                     {"sin", n_binary + 2},
                                                     {"exp", n_binary + 3}};

// Constructor
expression::expression(unsigned nvar, unsigned ncon, std::vector<std::string> kernels,
                       decltype(std::random_device{}()) seed)
    : m_nvar(nvar), m_ncon(ncon), m_rng(seed)
{
    m_kernels.resize(kernels.size());
    // We initialize the internal kernel id (unsigned) with the user requests (strings)
    for (decltype(kernels.size()) i = 0u; i < kernels.size(); ++i) {
        if (kernel_map.find(kernels[i]) != kernel_map.end()) {
            m_kernels[i] = kernel_map[kernels[i]];
        } else {
            throw std::invalid_argument("The requested kernel (" + kernels[i] + ") has not been implemeted");
        }
    }
    m_nker = m_kernels.size();
};

// Methods
std::vector<double> expression::random_constants(double lb, double ub)
{
    std::uniform_real_distribution<> dis(lb, ub);
    std::vector<double> retval(m_ncon, 0.);
    for (auto &el : retval) {
        el = dis(m_rng);
    }
    return retval;
}

// Assuming 0:+, 1:-, 2:* 3:/
std::vector<unsigned> expression::random_genotype(unsigned length)
{
    std::vector<unsigned> retval(3 * length);
    unsigned nus = m_nvar + m_ncon;
    std::uniform_int_distribution<int> uni_ker(0, m_nker - 1);
    for (auto i = 0u; i < length; ++i) {
        // Lets pick a random kernel
        auto funidx = uni_ker(m_rng);
        // Lets pick a random u1
        std::uniform_int_distribution<int> uni_us(0, nus - 1);
        unsigned u1idx, u2idx;
        while (true) { // TODO -> CHANGE THIS LOGIC!!!! infinite loops and unclear
            u1idx = uni_us(m_rng);
            if (m_kernels[funidx] < n_binary) { // binary operator
                break;
            } else { // unary operator->do not select a constant
                if (u1idx < m_nvar || u1idx >= m_nvar + m_ncon) {
                    break;
                }
            }
        }
        // ... and a random u1 not equal to u2 if the kernel is diff
        while (true) { // TODO -> CHANGE THIS LOGIC!!!! infinite loops and unclear
            u2idx = uni_us(m_rng);
            if (m_kernels[funidx] != 1u) { // if diff force u1 and u2 to be different
                break;
            } else {
                if (u1idx != u2idx || (nus==1)) {
                    break;
                }
            }
        }
        retval[3 * i] = funidx;
        retval[3 * i + 1] = u1idx;
        retval[3 * i + 2] = u2idx;
        nus++;
    }
    return retval;
}

void expression::phenotype_impl(std::vector<double> &retval, const std::vector<unsigned> &genotype,
                                const std::vector<double> &vars, const std::vector<double> &cons, bool check)
{
    auto n_terminals = m_nvar + m_ncon;
    if (check) {
        check_genotype(genotype);
        if (vars.size() != m_nvar) {
            throw std::invalid_argument(
                "When calling phenotype_impl the number of variables to compute this phenotype is wrong.");
        }
        if (cons.size() != m_ncon) {
            throw std::invalid_argument(
                "When calling phenotype_impl the number of constants to compute this phenotype is wrong.");
        }
    }

    // Size will be the vars+constant values (x) and then the number of triplets F u0 u1
    auto n_triplets = genotype.size() / 3;
    retval.resize(n_terminals + n_triplets);
    // The u0, u1, ... are the values of variables and constants
    std::copy(vars.begin(), vars.end(), retval.begin());
    std::copy(cons.begin(), cons.end(), retval.begin() + m_nvar);

    // We loop and for each triplet compute the corresponding function
    for (decltype(n_triplets) i = 0u; i < n_triplets; ++i) {
        auto u0 = retval[genotype[3 * i + 1]];
        auto u1 = retval[genotype[3 * i + 2]];
        auto fidx = genotype[3 * i];
        switch (m_kernels[fidx]) {
            case 0:
                retval[i + n_terminals] = u0 + u1;
                break;
            case 1:
                retval[i + n_terminals] = u0 - u1;
                break;
            case 2:
                retval[i + n_terminals] = u0 * u1;
                break;
            case 3:
                retval[i + n_terminals] = u0 / u1;
                break;
            // non arithmetic kernels (unary only all assumed before binary ones)
            default:
                retval[i + n_terminals] = ukernel_list[m_kernels[fidx] - n_binary](u0);
        }
    }
}

void expression::phenotype(std::vector<double> &retval, const std::vector<unsigned> &genotype,
                           const std::vector<double> &vars, const std::vector<double> &cons)
{
    phenotype_impl(retval, genotype, vars, cons, true);
}

void expression::sphenotype_impl(std::vector<std::string> &retval, const std::vector<unsigned> &genotype,
                                 const std::vector<std::string> &vars, const std::vector<std::string> &cons, bool check)
{
    auto n_terminals = m_nvar + m_ncon;
    auto vnames = vars;
    auto cnames = cons;
    // If not passed we create boring variable names
    if (vnames.size() == 0) {
        for (auto i = 0u; i < m_nvar; ++i) {
            vnames.push_back("x" + std::to_string(i));
        }
    }
    if (cnames.size() == 0) {
        for (auto i = 0u; i < m_ncon; ++i) {
            cnames.push_back("c" + std::to_string(i));
        }
    }

    if (check) {
        check_genotype(genotype);
        if (vnames.size() != m_nvar) {
            throw std::invalid_argument(
                "When calling sphenotype_impl the number of variables to compute this phenotype is wrong.");
        }
        if (cnames.size() != m_ncon) {
            throw std::invalid_argument(
                "When calling sphenotype_impl the number of constants to compute this phenotype is wrong.");
        }
    }

    // Size will be the vars+constant values (x) and then the number of triplets F u0 u1
    auto n_triplets = genotype.size() / 3;
    retval.resize(n_terminals + n_triplets);
    // The u0, u1, ... are the values of variables and constants
    std::copy(vnames.begin(), vnames.end(), retval.begin());
    std::copy(cnames.begin(), cnames.end(), retval.begin() + m_nvar);

    // We loop and for each triplet compute the corresponding function
    for (decltype(n_triplets) i = 0u; i < n_triplets; ++i) {
        auto u0 = retval[genotype[3 * i + 1]];
        auto u1 = retval[genotype[3 * i + 2]];
        auto fidx = genotype[3 * i];
        switch (m_kernels[fidx]) {
            case 0:
                retval[i + n_terminals] = "(" + u0 + "+" + u1 + ")";
                break;
            case 1:
                retval[i + n_terminals] = "(" + u0 + "-" + u1 + ")";
                break;
            case 2:
                retval[i + n_terminals] = "(" + u0 + "*" + u1 + ")";
                break;
            case 3:
                retval[i + n_terminals] = "(" + u0 + "/" + u1 + ")";
                break;
            // non arithmetic kernels (unary only all assumed before binary ones)
            default:
                retval[i + n_terminals] = pkernel_list[m_kernels[fidx] - n_binary](u0);
        }
    }
}

void expression::sphenotype(std::vector<std::string> &retval, const std::vector<unsigned> &genotype,
                            const std::vector<std::string> &vars, const std::vector<std::string> &cons)
{
    sphenotype_impl(retval, genotype, vars, cons, true);
}

// First order derivatives
void expression::dphenotype_impl(std::vector<double> &retval, const std::vector<unsigned> &genotype,
                                 const std::vector<double> &phenotype, unsigned idx, bool check)
{
    if (check) {
        check_genotype(genotype);
        if (idx >= m_nvar + m_ncon) {
            throw std::invalid_argument("When calling dphenotype_impl the idx of the derivation variable is larger "
                                        "than the sum n_vars + n_cons");
        }
    }
    // Number of terminals (vars and cons)
    unsigned n_terminals = m_nvar + m_ncon;
    // Number of triplets (F idx0, idx1 in the chromosome)
    auto n_triplets = genotype.size() / 3;
    // Size of the return value will be the same as phenotype
    retval.resize(phenotype.size());
    std::fill(retval.data(), retval.data() + m_nvar + m_ncon, 0.);
    // The du0, du1, ... for terminals are all zeros except the idx
    retval[idx] = 1.;
    // We loop and for each triplet compute the corresponding function
    for (decltype(n_triplets) i = 0u; i < n_triplets; ++i) {
        // Retrieve the values
        auto u0 = phenotype[genotype[3 * i + 1]];
        auto u1 = phenotype[genotype[3 * i + 2]];
        // Retrieve the derivatives
        auto d_u0 = retval[genotype[3 * i + 1]];
        auto d_u1 = retval[genotype[3 * i + 2]];
        // Retrieve the function
        auto fidx = genotype[3 * i];
        switch (m_kernels[fidx]) {
            case 0:
                retval[i + n_terminals] = d_u0 + d_u1;
                break;
            case 1:
                retval[i + n_terminals] = d_u0 - d_u1;
                break;
            case 2:
                retval[i + n_terminals] = u0 * d_u1 + d_u0 * u1;
                break;
            case 3:
                retval[i + n_terminals] = (d_u0 * u1 - d_u1 * u0) / (u1 * u1);
                break;
            // non arithmetic kernels (unary only all assumed before binary ones)
            default:
                retval[i + n_terminals] = dukernel_list[m_kernels[fidx] - n_binary](u0) * d_u0;
        }
    }
}

void expression::dphenotype(std::vector<double> &retval, const std::vector<unsigned> &genotype,
                            const std::vector<double> &phenotype, unsigned idx)
{
    dphenotype_impl(retval, genotype, phenotype, idx, true);
}
// Second order derivative
void expression::ddphenotype_impl(std::vector<double> &retval, const std::vector<unsigned> &genotype,
                                  const std::vector<double> &phenotype, const std::vector<double> &d0phenotype,
                                  const std::vector<double> &d1phenotype, bool check)
{
    assert(d0phenotype.size() == d1phenotype.size());
    assert(phenotype.size() == d1phenotype.size());
    if (check) {
        check_genotype(genotype);
        if (d0phenotype.size() != d1phenotype.size()) {
            throw std::invalid_argument(
                "When calling ddphenotype_impl the sizes of the gradients (d2phen, d1phen) is to be equal");
        }
        if (phenotype.size() != d1phenotype.size()) {
            throw std::invalid_argument("When calling ddphenotype_impl the sizes of the gradients (d1phen) is to be "
                                        "equal to the size of the phenotype");
        }
    }

    // Number of terminals (vars and cons)
    unsigned n_terminals = m_nvar + m_ncon;
    // Number of triplets (F idx0, idx1 in the chromosome)
    auto n_triplets = genotype.size() / 3;
    // Size of the return value will be the same as phenotype
    retval.resize(phenotype.size());
    std::fill(retval.data(), retval.data() + m_nvar + m_ncon, 0.);
    // We loop and for each triplet compute the corresponding function
    for (decltype(n_triplets) i = 0u; i < n_triplets; ++i) {
        // Retrieve the values
        auto u0 = phenotype[genotype[3 * i + 1]];
        auto u1 = phenotype[genotype[3 * i + 2]];
        // Retrieve the derivatives
        auto d0_u0 = d0phenotype[genotype[3 * i + 1]];
        auto d0_u1 = d0phenotype[genotype[3 * i + 2]];
        auto d1_u0 = d1phenotype[genotype[3 * i + 1]];
        auto d1_u1 = d1phenotype[genotype[3 * i + 2]];
        // Retrieve the second derivatives
        auto dd_u0 = retval[genotype[3 * i + 1]];
        auto dd_u1 = retval[genotype[3 * i + 2]];
        // Retrieve the function
        auto fidx = genotype[3 * i];
        switch (m_kernels[fidx]) {
            // +
            case 0:
                retval[i + n_terminals] = dd_u0 + dd_u1;
                break;
            // -
            case 1:
                retval[i + n_terminals] = dd_u0 - dd_u1;
                break;
            // *
            case 2:
                retval[i + n_terminals] = dd_u0 * u1 + dd_u1 * u0 + d0_u0 * d1_u1 + d1_u0 * d0_u1;
                break;
            // /
            case 3:
                retval[i + n_terminals] = ((dd_u0 * u1 + d0_u0 * d1_u1 - d0_u1 * d1_u0 - dd_u1 * u0) * u1 * u1
                                           - 2 * u1 * d1_u1 * (d0_u0 * u1 - d0_u1 * u0))
                                          / u1 / u1 / u1 / u1;
                break;
            // non arithmetic kernels (unary only all assumed before binary ones)
            default:
                retval[i + n_terminals] = ddukernel_list[m_kernels[fidx] - n_binary](u0) * d0_u0 * d1_u0
                                          + dukernel_list[m_kernels[fidx] - n_binary](u0) * dd_u0;
        }
    }
}

void expression::ddphenotype(std::vector<double> &retval, const std::vector<unsigned> &genotype,
                             const std::vector<double> &phenotype, const std::vector<double> &d0phenotype,
                             const std::vector<double> &d1phenotype)
{
    ddphenotype_impl(retval, genotype, phenotype, d0phenotype, d1phenotype, true);
}

std::vector<double> expression::mse(const std::vector<unsigned> &genotype, const std::vector<double> &cons,
                                    const std::vector<std::vector<double>> &xs, const std::vector<double> &ys)
{
    auto N = xs.size();
    assert(m_nvar == xs[0].size());
    std::vector<double> retval(m_nvar + m_ncon + genotype.size() / 3, 0u);
    std::vector<double> squared_err;
    for (decltype(xs.size()) i = 0u; i < N; ++i) {
        // compute all values in the phenotype (u0, u1, u2, .... un) at xs[i], cons
        phenotype(squared_err, genotype, xs[i], cons);
        // subtract ys[i] and square
        for (auto &element : squared_err) {
            element -= ys[i];
            element *= element;
        }
        // Add to retval
        std::transform(retval.begin(), retval.end(), squared_err.begin(), retval.begin(), std::plus<double>());
    }
    std::transform(retval.begin(), retval.end(), retval.begin(), [N](double &c) { return c / N; });

    return retval;
}

std::vector<double> expression::fitness(const std::vector<unsigned> &genotype, const std::vector<double> &cons,
                                        const std::vector<std::vector<double>> &xs, const std::vector<double> &ys)
{
    std::vector<double> retval(3);
    auto errors = mse(genotype, cons, xs, ys);
    retval[2] = std::reduce(errors.begin(), errors.end(), 0.0) / errors.size();
    retval[1] = *std::max_element(errors.begin(), errors.end());
    retval[0] = *std::min_element(errors.begin(), errors.end());
    return retval;
}

//// computes mse, dmse and ddmse in one go
// void expression::ddmse(const std::vector<unsigned> &genotype, const std::vector<double> &cons, unsigned c_idx,
//                        const std::vector<std::vector<double>> &xs, const std::vector<double> &ys,
//                        std::vector<double> &mse, std::vector<double> &dmse, std::vector<double> &ddmse)
//{
//     auto N = xs.size();
//
//     assert(mse.size() == genotype.size() / 3 + m_nvar + m_ncon);
//
//     std::vector<double> ph;
//     std::vector<double> dph;
//     std::vector<double> ddph;
//
//     std::fill(mse.begin(), mse.end(), 0.);
//     std::fill(dmse.begin(), dmse.end(), 0.);
//     std::fill(ddmse.begin(), ddmse.end(), 0.);
//
//     // For each point
//     for (decltype(N) i = 0u; i < N; ++i) {
//         // The value of each expression in the single point
//         phenotype_impl(ph, genotype, xs[i], cons, false);
//         // The derivative value w.r.t. a c_idx constant
//         dphenotype_impl(dph, genotype, ph, c_idx, false);
//         // The second derivative value w.r.t. c c
//         ddphenotype_impl(ddph, genotype, ph, dph, dph, false);
//         // We will now store in ph. dph, ddph respectively, the mse, dmse and ddmse
//         // NOTE: The order of the next loops counts and should not be touched.
//         // yi-\hat y_i
//         for (auto j = 0u; j < ph.size(); ++j) {
//             ph[j] -= ys[i];
//         }
//         // 2 ((yi-\hat y_i)d2ydc2+(dy/dc)^2)
//         for (auto j = 0u; j < ddph.size(); ++j) {
//             ddph[j] = 2 * (ph[j] / 2. * ddph[j] + dph[j] * dph[j]);
//         }
//         // 2 (yi-\hat y_i)dydc
//         for (auto j = 0u; j < dph.size(); ++j) {
//             dph[j] = 2 * ph[j] * dph[j];
//         }
//         // finally (yi-\hat y_i)^2
//         for (auto j = 0u; j < ph.size(); ++j) {
//             ph[j] *= ph[j];
//         }
//         // And we accumulate
//         std::transform(mse.begin(), mse.end(), ph.begin(), mse.begin(), std::plus<double>());
//         std::transform(dmse.begin(), dmse.end(), dph.begin(), dmse.begin(), std::plus<double>());
//         std::transform(ddmse.begin(), ddmse.end(), ddph.begin(), ddmse.begin(), std::plus<double>());
//     }
//     std::transform(mse.begin(), mse.end(), mse.begin(), [N](double &c) { return c / N; });
//     std::transform(dmse.begin(), dmse.end(), dmse.begin(), [N](double &c) { return c / N; });
//     std::transform(ddmse.begin(), ddmse.end(), ddmse.begin(), [N](double &c) { return c / N; });
// }

void expression::ddmse(const std::vector<unsigned> &genotype, const std::vector<double> &cons,
                       const std::vector<std::vector<double>> &xs, const std::vector<double> &ys,
                       std::vector<double> &mse, std::vector<std::vector<double>> &gradient,
                       std::vector<std::vector<double>> &hess)
{
    auto N_points = xs.size();
    auto N_us = genotype.size() / 3 + m_nvar + m_ncon;

    // These will store values and derivatives of the expressions (not the loss)
    std::vector<double> ph;
    std::vector<std::vector<double>> dph(m_ncon);
    std::vector<std::vector<double>> ddph(m_ncon * (1 + m_ncon) / 2.);

    // We init all to 0. and resize if needed.
    mse.resize(N_us);
    std::fill(mse.begin(), mse.end(), 0.);

    gradient.resize(m_ncon);
    for (auto &it : gradient) {
        it.resize(N_us);
        std::fill(it.begin(), it.end(), 0.);
    }
    hess.resize(m_ncon * (1 + m_ncon) / 2.);
    for (auto &it : hess) {
        it.resize(N_us);
        std::fill(it.begin(), it.end(), 0.);
    }

    // Loop ove  points
    for (decltype(N_points) i = 0u; i < N_points; ++i) {
        // The value of expression i will be stored in ph[i]
        phenotype_impl(ph, genotype, xs[i], cons, false);
        // The gradient of expression i with respect to constant j will be stored in dph[j][i]
        for (decltype(m_ncon) j = 0u; j < m_ncon; ++j) {
            // The derivative value w.r.t. the jth constant
            dphenotype_impl(dph[j], genotype, ph, j, false);
        }
        // The hessian of expression idx_u with respect to constants jk will be stored in dph[jk][idx_u] - jk =
        // [00,01,02,10,11,20]
        auto jk = 0u;
        for (decltype(m_ncon) j = 0u; j < m_ncon; ++j) {
            for (decltype(j) k = 0u; k <= j; ++k) {
                ddphenotype_impl(ddph[jk], genotype, ph, dph[j], dph[k], false);
                jk++;
            }
        }
        // We now must construct the value, gradient and hessian for the loss (mse) ....
        // 1 - we compute (yi-\hat y_i)
        for (decltype(ph.size()) idx_u = 0u; idx_u < ph.size(); ++idx_u) {
            ph[idx_u] -= ys[i];
        }

        // 2 - we compute the hessian, jk component ((yi-\hat y_i)d2ydcjdck+(dy/dcj)(dy/dck))
        jk = 0u;
        for (decltype(m_ncon) j = 0u; j < m_ncon; ++j) {
            for (decltype(j) k = 0u; k <= j; ++k) {
                for (decltype(ddph[jk].size()) idx_u = 0u; idx_u < ddph[jk].size(); ++idx_u) {
                    ddph[jk][idx_u] = 2 * (ph[idx_u] / 2. * ddph[jk][idx_u] + dph[j][idx_u] * dph[k][idx_u]);
                }
                std::transform(hess[jk].begin(), hess[jk].end(), ddph[jk].begin(), hess[jk].begin(),
                               std::plus<double>());
                jk++;
            }
        }

        // 3 - we compute the gradient 2 (yi-\hat y_i)dydcj
        for (decltype(dph.size()) j = 0u; j < dph.size(); ++j) {
            for (decltype(dph[j].size()) idx_u = 0u; idx_u < dph[j].size(); ++idx_u) {
                dph[j][idx_u] = 2 * ph[idx_u] * dph[j][idx_u];
            }
            std::transform(gradient[j].begin(), gradient[j].end(), dph[j].begin(), gradient[j].begin(),
                           std::plus<double>());
        }
        // 4 - finally we compute (yi-\hat y_i)^2
        for (auto idx_u = 0u; idx_u < ph.size(); ++idx_u) {
            ph[idx_u] *= ph[idx_u];
        }
        std::transform(mse.begin(), mse.end(), ph.begin(), mse.begin(), std::plus<double>());
    }
    // Finlly we divide all by the number of points (can be done before)
    std::transform(mse.begin(), mse.end(), mse.begin(), [N_points](double &c) { return c / N_points; });
    for (auto &it : gradient) {
        std::transform(it.begin(), it.end(), it.begin(), [N_points](double &c) { return c / N_points; });
    }
    for (auto &it : gradient) {
        std::transform(it.begin(), it.end(), it.begin(), [N_points](double &c) { return c / N_points; });
    }
}

std::vector<unsigned> expression::mutation(std::vector<unsigned> genotype, unsigned N)
{
    auto retval = genotype;
    // We generate N randomly selected indexes of the genotype triplets
    auto n_triplets = genotype.size() / 3;
    std::vector<unsigned> choice(n_triplets);
    std::iota(choice.begin(), choice.end(), 0u);
    std::shuffle(choice.begin(), choice.end(), m_rng);
    // We generate a new feasible random genotype
    auto muts = random_genotype(genotype.size());
    // For each selected triplet we use the randomly generated one
    for (auto i = 0u; i < N; ++i) {
        retval[3 * choice[i]] = muts[3 * choice[i]];
        retval[3 * choice[i] + 1] = muts[3 * choice[i] + 1];
        retval[3 * choice[i] + 2] = muts[3 * choice[i] + 2];
    }
    return retval;
}

std::vector<unsigned> expression::mutation2(std::vector<unsigned> genotype, unsigned N)
{
    auto retval = genotype;
    // We generate N randomly selected indexes of the genotype triplets
    std::vector<unsigned> choice(retval.size());
    std::iota(choice.begin(), choice.end(), 0u);
    std::shuffle(choice.begin(), choice.end(), m_rng);
    // For each selected gene, we randomly regenerate a feasible one.
    for (auto i = 0u; i < N; ++i) {
        if (choice[i] % 3 == 0u) {
            std::uniform_int_distribution<unsigned> dis(0, m_kernels.size() - 1);
            retval[choice[i]] = dis(m_rng);
        } else {
            std::uniform_int_distribution<unsigned> dis(0, m_ncon + m_nvar + choice[i] / 3 - 1);
            retval[choice[i]] = dis(m_rng);
        }
    }
    return retval;
}

void expression::remove_nesting(std::vector<unsigned> &g) const
{
    auto n_triplets = g.size() / 3;
    for (auto i = 1u; i < n_triplets; ++i) {   // skip first triplet as connection genes are surely var or con
        if (m_kernels[g[3 * i]] >= n_binary) { // unary operator
            auto u1 = g[3 * i + 1];
            if (u1 > m_ncon + m_nvar) { // not taking in a var or con
                std::uniform_int_distribution<int> random_kernel(0u, m_kernels.size() - 1);
                // this loop is infinite if no binary operators are in the kernels.
                while (m_kernels[g[3 * (u1 - m_ncon - m_nvar)]] >= n_binary) { // and nesting one more unary
                    // we change the operator randomly
                    g[3 * (u1 - m_ncon - m_nvar)] = random_kernel(m_rng);
                }
            }
        }
    }
}

const std::vector<unsigned> &expression::get_kernels_idx() const
{
    return m_kernels;
}

void expression::check_genotype(const std::vector<unsigned> &g) const
{
    // Check that the genotype is made out of triplets.
    if (g.size() % 3 != 0) {
        throw std::invalid_argument("The genotype length must be a multiple of 3.");
    }
    // Check function genes contain valid idxs in the kernels and connections that are acyclic.
    auto n_triplets = g.size() / 3;
    for (auto i = 0u; i < n_triplets; ++i) {
        if (m_kernels[g[3 * i]] >= n_unary + n_binary) {
            throw std::invalid_argument(
                "The genotype contains a function gene id out of the available implemented ones");
        }
        if ((g[3 * i + 1] >= i + m_ncon + m_nvar) || (g[3 * i + 2] >= i + m_ncon + m_nvar)) {
            fmt::print("genotype: {}\n", g);
            fmt::print("position: {}\n", i);
            throw std::invalid_argument("The genotype contains an incompatibe connection gene");
        }
    }
}
} // namespace dsyre