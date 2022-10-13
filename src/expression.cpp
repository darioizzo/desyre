// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#include <dsyre/detail/visibility.hpp>
#include <dsyre/expression.hpp>
#include <dsyre/kernels.hpp>

namespace dsyre
{

using ukernel_f_ptr = double (*)(double);
using pkernel_f_ptr = std::string (*)(std::string);

/** INSTRUCTIONS TO ADD ONE KERNEL (assumes the kernel, its derivative, its second order derivative and its string
representation has been put already in kernels.hpp).

1 - add it to ukernel_list, dukernel_list, ddukernel_list, pkernel_list
2 - update if needed n_binary
3 - add it to the map kernel_map */

// Global array of function pointers for the unary kernels
ukernel_f_ptr ukernel_list[] = {inv, cos, sin, exp};
ukernel_f_ptr dukernel_list[] = {dinv, dcos, dsin, dexp};
ukernel_f_ptr ddukernel_list[] = {ddinv, ddcos, ddsin, ddexp};
pkernel_f_ptr pkernel_list[] = {pinv, pcos, psin, pexp};

// Number of binary operations (+,-,*,/)
unsigned n_binary = 4u;
unsigned n_unary = std::size(ukernel_list);

// Global maps between kernel names and an unsigned (must correspond to the order in the global arrays as its used
// in the switch cases.
std::unordered_map<unsigned, std::string> build_reverse_map(const std::unordered_map<std::string, unsigned> &kernel_map)
{
    std::unordered_map<unsigned, std::string> retval;
    for (auto i = kernel_map.begin(); i != kernel_map.end(); ++i) {
        retval[i->second] = i->first;
    }
    return retval;
}
std::unordered_map<std::string, unsigned> kernel_map{{"sum", 0},
                                                     {"diff", 1},
                                                     {"mul", 2},
                                                     {"div", 3},
                                                     {"inv", n_binary},
                                                     {"cos", n_binary + 1},
                                                     {"sin", n_binary + 2},
                                                     {"exp", n_binary + 3}};

std::unordered_map<unsigned, std::string> reverse_kernel_map = build_reverse_map(kernel_map);

// Constructor
expression::expression(unsigned nvar, unsigned ncon, std::vector<std::string> kernels) : m_nvar(nvar), m_ncon(ncon)
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
std::vector<double> expression::random_constants(double lb, double ub, std::mt19937 &rng) const
{
    std::uniform_real_distribution<> dis(lb, ub);
    std::vector<double> retval(m_ncon, 0.);
    for (auto &el : retval) {
        el = dis(rng);
    }
    return retval;
}

// Assuming 0:+, 1:-, 2:* 3:/
void expression::random_genotype(std::vector<unsigned> &retval, unsigned length, std::mt19937 &rng) const
{
    retval.resize(3 * length);
    unsigned nus = m_nvar + m_ncon;
    std::uniform_int_distribution<int> uni_ker(0, m_nker - 1);
    for (auto i = 0u; i < length; ++i) {
        // Lets pick a random kernel
        auto funidx = uni_ker(rng);
        // Lets pick a random u1
        std::uniform_int_distribution<int> uni_us(0, nus - 1);
        unsigned u1idx, u2idx;
        while (true) { // TODO -> CHANGE THIS LOGIC!!!! infinite loops and unclear
            u1idx = uni_us(rng);
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
            u2idx = uni_us(rng);
            if (m_kernels[funidx] != 1u) { // if diff force u1 and u2 to be different
                break;
            } else {
                if (u1idx != u2idx || (nus == 1)) {
                    break;
                }
            }
        }
        retval[3 * i] = funidx;
        retval[3 * i + 1] = u1idx;
        retval[3 * i + 2] = u2idx;
        nus++;
    }
}

void expression::phenotype_impl(std::vector<double> &retval, const std::vector<unsigned> &genotype,
                                const std::vector<double> &vars, const std::vector<double> &cons) const
{
    auto n_terminals = m_nvar + m_ncon;

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
                           const std::vector<double> &vars, const std::vector<double> &cons) const
{
    check_genotype(genotype);
    if (vars.size() != m_nvar) {
        throw std::invalid_argument("The dimensions of vars to compute this phenotype seems wrong.");
    }
    if (cons.size() != m_ncon) {
        throw std::invalid_argument("The dimensions of cons to compute this phenotype seems wrong.");
    }
    phenotype_impl(retval, genotype, vars, cons);
}

void expression::complexity(std::vector<unsigned> &retval, const std::vector<unsigned> &genotype) const
{
    check_genotype(genotype);
    complexity_impl(retval, genotype);
}

void expression::complexity_impl(std::vector<unsigned> &retval, const std::vector<unsigned> &genotype) const
{
    {
        auto n_terminals = m_nvar + m_ncon;
        // Size will be the vars+constant values (x) and then the number of triplets F u0 u1
        auto n_triplets = genotype.size() / 3;
        retval.resize(n_terminals + n_triplets);
        // The complexity of the model made by only a variable or a constant is one
        std::fill(retval.begin(), retval.begin() + n_terminals, 1u);
        // We loop and for each triplet compute the corresponding function and add to the expression complexity
        for (decltype(n_triplets) i = 0u; i < n_triplets; ++i) {
            auto fidx = genotype[3 * i];
            if (m_kernels[fidx] > n_binary) { // unary
                retval[i + n_terminals] = 1u + retval[genotype[3 * i + 1]];
            } else { // binary
                retval[i + n_terminals] = 1u + retval[genotype[3 * i + 1]] + retval[genotype[3 * i + 2]];
            }
        }
    }
}

void expression::sphenotype_impl(std::vector<std::string> &retval, const std::vector<unsigned> &genotype,
                                 const std::vector<std::string> &vars, const std::vector<std::string> &cons) const
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
                            const std::vector<std::string> &vars, const std::vector<std::string> &cons) const
{
    check_genotype(genotype);
    if (vars.size() != m_nvar && vars.size() > 0) {
        throw std::invalid_argument(
            "When calling sphenotype the number of variables to compute this phenotype is wrong.");
    }
    if (cons.size() != m_ncon && cons.size() > 0) {
        throw std::invalid_argument(
            "When calling sphenotype the number of constants to compute this phenotype is wrong.");
    }

    sphenotype_impl(retval, genotype, vars, cons);
}

// First order derivatives
void expression::dphenotype_impl(std::vector<double> &retval, const std::vector<unsigned> &genotype,
                                 const std::vector<double> &phenotype, unsigned idx) const
{
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
                            const std::vector<double> &phenotype, unsigned idx) const
{
    check_genotype(genotype);
    if (idx >= m_nvar + m_ncon) {
        throw std::invalid_argument("When calling dphenotype the idx of the derivation variable is larger "
                                    "than the sum n_vars + n_cons");
    }
    if ((phenotype.size() - m_nvar - m_ncon) * 3 != genotype.size()) {
        throw std::invalid_argument("When calling dphenotype the phenotype and genotype lengths are not compatible");
    }
    dphenotype_impl(retval, genotype, phenotype, idx);
}
// Second order derivative
void expression::ddphenotype_impl(std::vector<double> &retval, const std::vector<unsigned> &genotype,
                                  const std::vector<double> &phenotype, const std::vector<double> &d0phenotype,
                                  const std::vector<double> &d1phenotype) const
{
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
                             const std::vector<double> &d1phenotype) const
{
    check_genotype(genotype);
    if (d0phenotype.size() != d1phenotype.size()) {
        throw std::invalid_argument(
            "When calling ddphenotype the sizes of the gradients (d2phen, d1phen) is to be equal");
    }
    if (phenotype.size() != d1phenotype.size()) {
        throw std::invalid_argument("When calling ddphenotype the sizes of the gradients (d1phen) is to be "
                                    "equal to the size of the phenotype");
    }
    if ((phenotype.size() - m_nvar - m_ncon) * 3 != genotype.size()) {
        throw std::invalid_argument("When calling ddphenotype the phenotype and genotype lengths are not compatible");
    }

    ddphenotype_impl(retval, genotype, phenotype, d0phenotype, d1phenotype);
}

void expression::mse(std::vector<double> &retval, const std::vector<unsigned> &genotype,
                     const std::vector<double> &cons, const std::vector<std::vector<double>> &xs,
                     const std::vector<double> &ys) const
{
    auto N = xs.size();
    if (xs.size() != ys.size()) {
        throw std::invalid_argument("When computing the mse the data xs, ys seems to be malformed.");
    }

    retval.resize(m_nvar + m_ncon + genotype.size() / 3, 0u);
    std::fill(retval.begin(), retval.end(), 0);
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
}

// This is the main computational routine whose efficiency determines the overall efficiency of dsyre.
// Conventions followed for the grad and hess --------------------------------------------------------
// If we focus (say) on the i constant and some phenotype [u0, u1, ..., un]
// grad[i] -> The derivative of all the phenotype w.r.t. the constant. [du0, du1, ..., dun]
// For the Hessian, the relation is more complex as we flattened a symmetric matrix (i and j) -> ij
// following the convention [00 10 11 20 21 22 ....].
// If we seek the derivative w.r.t the constants i and j, then
// hess[ij] - > The second order derivative w.r.t i and j. [ddu00, ddu10, ddu11, ddu20 ...]
void expression::ddmse(std::vector<double> &mse, std::vector<std::vector<double>> &grad,
                       std::vector<std::vector<double>> &hess, const std::vector<unsigned> &genotype,
                       const std::vector<double> &cons, const std::vector<std::vector<double>> &xs,
                       const std::vector<double> &ys) const
{
    auto N_points = xs.size();
    auto N_us = genotype.size() / 3 + m_nvar + m_ncon;
    if (xs.size() != ys.size()) {
        throw std::invalid_argument("When computing the mse the data xs, ys seems to be malformed.");
    }
    // These will store values and derivatives of the expressions (not the loss)
    std::vector<double> ph;
    std::vector<std::vector<double>> dph(m_ncon);
    std::vector<std::vector<double>> ddph(m_ncon * (1 + m_ncon) / 2.);

    // We init all to 0. and resize if needed.
    mse.resize(N_us);
    std::fill(mse.begin(), mse.end(), 0.);

    grad.resize(m_ncon);
    for (auto &it : grad) {
        it.resize(N_us);
        std::fill(it.begin(), it.end(), 0.);
    }
    hess.resize(m_ncon * (1 + m_ncon) / 2.);
    for (auto &it : hess) {
        it.resize(N_us);
        std::fill(it.begin(), it.end(), 0.);
    }

    // Loop over data points
    for (decltype(N_points) i = 0u; i < N_points; ++i) {
        // The value of expression i will be stored in ph[i]
        phenotype_impl(ph, genotype, xs[i], cons);
        // The gradient of expression i with respect to constant j will be stored in dph[j][i]
        for (decltype(m_ncon) j = 0u; j < m_ncon; ++j) {
            // The derivative value w.r.t. the jth constant
            dphenotype_impl(dph[j], genotype, ph, j + m_nvar);
        }
        // The hessian of expression idx_u with respect to constants jk will be stored in dph[jk][idx_u]
        // jk = [00,10,11,20,21,22,....]
        auto jk = 0u;
        for (decltype(m_ncon) j = 0u; j < m_ncon; ++j) {
            for (decltype(j) k = 0u; k <= j; ++k) {
                ddphenotype_impl(ddph[jk], genotype, ph, dph[j], dph[k]);
                jk++;
            }
        }

        // We now must construct the value, gradient and hessian for the loss (mse) .... (so far we computed it for the
        // phenotype) 1 - we compute (yi-\hat y_i)
        for (decltype(ph.size()) idx_u = 0u; idx_u < ph.size(); ++idx_u) {
            ph[idx_u] -= ys[i];
        }

        // 2 - we compute the hessian, jk component ((yi-\hat y_i)d2ydcjdck+(dy/dcj)(dy/dck))
        jk = 0u;
        for (decltype(m_ncon) j = 0u; j < m_ncon; ++j) {
            for (decltype(j) k = 0u; k <= j; ++k) {
                for (decltype(ddph[jk].size()) idx_u = 0u; idx_u < ddph[jk].size(); ++idx_u) {
                    ddph[jk][idx_u] = 2 * (ph[idx_u] * ddph[jk][idx_u] + dph[j][idx_u] * dph[k][idx_u]);
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
            std::transform(grad[j].begin(), grad[j].end(), dph[j].begin(), grad[j].begin(), std::plus<double>());
        }

        // 4 - and we compute (yi-\hat y_i)^2
        for (auto idx_u = 0u; idx_u < ph.size(); ++idx_u) {
            ph[idx_u] *= ph[idx_u];
        }
        std::transform(mse.begin(), mse.end(), ph.begin(), mse.begin(), std::plus<double>());
    }
    // Finally we divide all by the number of points (can be done before)
    std::transform(mse.begin(), mse.end(), mse.begin(), [N_points](double &c) { return c / N_points; });
    for (auto &it : grad) {
        std::transform(it.begin(), it.end(), it.begin(), [N_points](double &c) { return c / N_points; });
    }
    for (auto &it : hess) {
        std::transform(it.begin(), it.end(), it.begin(), [N_points](double &c) { return c / N_points; });
    }
}

std::vector<unsigned> expression::mutation(const std::vector<unsigned> &genotype, unsigned N, std::mt19937 &rng)
{
    auto retval = genotype;
    // We generate N randomly selected indexes of the genotype triplets
    auto n_triplets = genotype.size() / 3;
    std::vector<unsigned> choice(n_triplets);
    std::iota(choice.begin(), choice.end(), 0u);
    std::shuffle(choice.begin(), choice.end(), rng);
    // We generate a new feasible random genotype
    random_genotype(retval, n_triplets, rng);
    // For each selected triplet we use the randomly generated one
    for (auto i = 0u; i < n_triplets - N; ++i) {
        retval[3 * choice[i]] = genotype[3 * choice[i]];
        retval[3 * choice[i] + 1] = genotype[3 * choice[i] + 1];
        retval[3 * choice[i] + 2] = genotype[3 * choice[i] + 2];
    }
    return retval;
}

std::vector<unsigned> expression::mutation2(const std::vector<unsigned> &genotype, unsigned N, std::mt19937 &rng)
{
    auto retval = genotype;
    // We generate N randomly selected indexes of the genotype triplets
    std::vector<unsigned> choice(retval.size());
    std::iota(choice.begin(), choice.end(), 0u);
    std::shuffle(choice.begin(), choice.end(), rng);
    // For each selected gene, we randomly regenerate a feasible one.
    for (auto i = 0u; i < N; ++i) {
        if (choice[i] % 3 == 0u) {
            std::uniform_int_distribution<unsigned> dis(0, m_kernels.size() - 1);
            retval[choice[i]] = dis(rng);
        } else {
            std::uniform_int_distribution<unsigned> dis(0, m_ncon + m_nvar + choice[i] / 3 - 1);
            retval[choice[i]] = dis(rng);
        }
    }
    return retval;
}

std::vector<unsigned> expression::mutation3(const std::vector<unsigned> &genotype, const std::vector<double> &phenotype,
                                            unsigned N, std::mt19937 &rng)
{
    auto retval = genotype;
    for (auto idx_u = m_nvar + m_ncon; idx_u < phenotype.size(); ++idx_u) {
        if (!std::isfinite(phenotype[idx_u])) {
            std::uniform_int_distribution<unsigned> dis1(0, m_kernels.size() - 1);
            std::uniform_int_distribution<unsigned> dis2(0, idx_u - 1);
            retval[(idx_u - m_nvar - m_ncon) * 3] = dis1(rng);
            retval[(idx_u - m_nvar - m_ncon) * 3 + 1] = dis2(rng);
            retval[(idx_u - m_nvar - m_ncon) * 3 + 2] = dis2(rng);
        }
    }

    // We generate N randomly selected indexes of the genotype triplets
    std::vector<unsigned> choice(retval.size());
    std::iota(choice.begin(), choice.end(), 0u);
    std::shuffle(choice.begin(), choice.end(), rng);
    // For each selected gene, we randomly regenerate a feasible one.
    for (auto i = 0u; i < N; ++i) {
        if (choice[i] % 3 == 0u) {
            std::uniform_int_distribution<unsigned> dis(0, m_kernels.size() - 1);
            retval[choice[i]] = dis(rng);
        } else {
            std::uniform_int_distribution<unsigned> dis(0, m_ncon + m_nvar + choice[i] / 3 - 1);
            retval[choice[i]] = dis(rng);
        }
    }
    return retval;
}

// NOTE: the inv function although unary will not count for nesting as to allow 1/sin.
// this allows for inv(inv) TODO: make a special case
void expression::remove_nesting(std::vector<unsigned> &g, std::mt19937 &rng) const
{
    auto n_triplets = g.size() / 3;
    for (auto i = 1u; i < n_triplets; ++i) {       // skip first triplet as connection genes are surely var or con
        if (m_kernels[g[3 * i]] >= n_binary + 1) { // unary operator (and not inv)
            auto u1 = g[3 * i + 1];
            if (u1 > m_ncon + m_nvar) { // not taking in a var or con
                std::uniform_int_distribution<unsigned> random_kernel(0u, m_kernels.size() - 1);
                // this loop is infinite if no binary operators(or inv) are in the kernels.
                while (m_kernels[g[3 * (u1 - m_ncon - m_nvar)]] >= n_binary + 1) { // and nesting one more unary
                    // we change the operator randomly
                    g[3 * (u1 - m_ncon - m_nvar)] = random_kernel(rng);
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
            std::string err_msg;
            err_msg += fmt::format("genotype: {}\n", g);
            err_msg += fmt::format("position: {}\n", i * 3);
            err_msg += fmt::format("this triplet contains an incompatibe connection gene");
            throw std::invalid_argument(err_msg);
        }
    }
}

/// Streaming operator for dsyre::expression.
/**
 * @param os target stream.
 * @param p the dsyre::expression
 *
 * @return a reference to \p os.
 */
std::ostream &operator<<(std::ostream &os, const expression &p)
{
    os << fmt::format("A differentiable expression using the dsyre encoding.\n");
    os << fmt::format("Number of variables: {}\n", p.m_nvar);
    os << fmt::format("Number of constants: {}\n", p.m_ncon);
    os << fmt::format("Number of kernels: {}\n", p.m_nker);
    std::vector<std::string> kernels_names;
    for (auto item : p.m_kernels) {
        kernels_names.push_back(reverse_kernel_map[item]);
    }
    os << fmt::format("Kernels: {}\n", kernels_names);

    return os;
}

std::unordered_map<std::string, unsigned> get_kernel_map()
{
    return kernel_map;
}

std::unordered_map<unsigned, std::string> get_reverse_kernel_map()
{
    return reverse_kernel_map;
}

} // namespace dsyre