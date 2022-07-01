// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <algorithm>
#include <stdexcept>

#include <dsyre/detail/visibility.hpp>
#include <dsyre/expression.hpp>
#include <dsyre/kernels.hpp>

namespace dsyre
{

using ukernel_f_ptr = double (*)(double);
using pkernel_f_ptr = std::string (*)(std::string);

// Global array of function pointers for the unary kernels
ukernel_f_ptr ukernel_list[] = {cos, sin, exp};
ukernel_f_ptr dukernel_list[] = {dcos, dsin, dexp};
ukernel_f_ptr ddukernel_list[] = {ddcos, ddsin, ddexp};
pkernel_f_ptr pkernel_list[] = {pcos, psin, pexp};

// Constructor
expression::expression(unsigned nvar, unsigned ncon, std::vector<unsigned> kernels,
                       decltype(std::random_device{}()) seed)
    : m_nvar(nvar), m_ncon(ncon), m_kernels(kernels), m_rng(seed)
{
    if (std::any_of(m_kernels.begin(), m_kernels.end(), [](unsigned i) { return i >= n_tot_kernels; })) {
        throw std::invalid_argument("Invalid kernel requested, all ids must be < " + std::to_string(n_ukernel));
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
    unsigned nus = 0u;
    std::uniform_int_distribution<int> uni_ker(0, m_nker - 1);
    for (auto i = 0u; i < length; ++i) {
        // Lets pick a random kernel (but not subtraction as a first pick)
        unsigned funidx = 1;
        while (funidx == 1) {
            funidx = uni_ker(m_rng);
            if (i > 0)
                break;
        }
        // Lets pick a random u1
        std::uniform_int_distribution<int> uni_us(0, m_nvar + m_ncon + nus - 1);
        unsigned u1idx = uni_us(m_rng);
        unsigned u2idx = u1idx;
        // ... and a random u1 not equal to u2 if the kernel is sub
        while (u2idx == u1idx) {
            u2idx = uni_us(m_rng);
            if (funidx != 1)
                break;
        }
        retval[3 * i] = funidx;
        retval[3 * i + 1] = u1idx;
        retval[3 * i + 2] = u2idx;
        nus++;
    }
    return retval;
}

std::vector<std::string> expression::sphenotype(const std::vector<unsigned> &genotype,
                                                const std::vector<std::string> &vars,
                                                const std::vector<std::string> &cons)
{
    assert(m_nvar == vars.size());
    assert(m_ncon == vars.size());
    auto n_terminals = m_nvar + m_ncon;

    // Size will be the vars+constant values (x) and then the number of triplets F u0 u1
    auto n_triplets = genotype.size() / 3;
    std::vector<std::string> sphenotype(n_terminals + n_triplets);
    // The u0, u1, ... are the values of variables and constants
    std::copy(vars.begin(), vars.end(), sphenotype.begin());
    std::copy(cons.begin(), cons.end(), sphenotype.begin() + m_nvar);

    // We loop and for each triplet compute the corresponding function
    for (decltype(n_triplets) i = 0u; i < n_triplets; ++i) {
        auto u0 = sphenotype[genotype[3 * i + 1]];
        auto u1 = sphenotype[genotype[3 * i + 2]];
        auto fidx = genotype[3 * i];
        switch (m_kernels[fidx]) {
            case 0:
                sphenotype[i + n_terminals] = "(" + u0 + "+" + u1 + ")";
                break;
            case 1:
                sphenotype[i + n_terminals] = "(" + u0 + "-" + u1 + ")";
                break;
            case 2:
                sphenotype[i + n_terminals] = "(" + u0 + "*" + u1 + ")";
                break;
            case 3:
                sphenotype[i + n_terminals] = "(" + u0 + "/" + u1 + ")";
                break;
            // non arithmetic kernels (unary only)
            default:
                sphenotype[i + n_terminals] = pkernel_list[m_kernels[fidx] - 4](u0);
        }
    }
    return sphenotype;
}

std::vector<double> expression::phenotype(const std::vector<unsigned> &genotype, const std::vector<double> &vars,
                                          const std::vector<double> &cons)
{
    assert(m_nvar == vars.size());
    assert(m_ncon == vars.size());
    auto n_terminals = m_nvar + m_ncon;

    // Size will be the vars+constant values (x) and then the number of triplets F u0 u1
    auto n_triplets = genotype.size() / 3;
    std::vector<double> phenotype(n_terminals + n_triplets);
    // The u0, u1, ... are the values of variables and constants
    std::copy(vars.begin(), vars.end(), phenotype.begin());
    std::copy(cons.begin(), cons.end(), phenotype.begin() + m_nvar);

    // We loop and for each triplet compute the corresponding function
    for (decltype(n_triplets) i = 0u; i < n_triplets; ++i) {
        auto u0 = phenotype[genotype[3 * i + 1]];
        auto u1 = phenotype[genotype[3 * i + 2]];
        auto fidx = genotype[3 * i];
        switch (m_kernels[fidx]) {
            case 0:
                phenotype[i + n_terminals] = u0 + u1;
                break;
            case 1:
                phenotype[i + n_terminals] = u0 - u1;
                break;
            case 2:
                phenotype[i + n_terminals] = u0 * u1;
                break;
            case 3:
                phenotype[i + n_terminals] = u0 / u1;
                break;
            // non arithmetic kernels (unary only)
            default:
                phenotype[i + n_terminals] = ukernel_list[m_kernels[fidx] - 4](u0);
        }
    }
    return phenotype;
}

// First order derivatives
std::vector<double> expression::dphenotype(const std::vector<unsigned> &genotype, const std::vector<double> &phenotype,
                                           unsigned idx)
{
    assert(idx < m_nvar + m_ncon);
    // Number of terminals (vars and cons)
    unsigned n_terminals = m_nvar + m_ncon;
    // Number of triplets (F idx0, idx1 in the chromosome)
    auto n_triplets = genotype.size() / 3;
    // Size of the return value will be the same as phenotype
    std::vector<double> dphenotype(phenotype.size(), 0.);
    // The du0, du1, ... for terminals are all zeros except the idx
    dphenotype[idx] = 1.;
    // We loop and for each triplet compute the corresponding function
    for (decltype(n_triplets) i = 0u; i < n_triplets; ++i) {
        // Retrieve the values
        auto u0 = phenotype[genotype[3 * i + 1]];
        auto u1 = phenotype[genotype[3 * i + 2]];
        // Retrieve the derivatives
        auto d_u0 = dphenotype[genotype[3 * i + 1]];
        auto d_u1 = dphenotype[genotype[3 * i + 2]];
        // Retrieve the function
        auto fidx = genotype[3 * i];
        switch (m_kernels[fidx]) {
            case 0:
                dphenotype[i + n_terminals] = d_u0 + d_u1;
                break;
            case 1:
                dphenotype[i + n_terminals] = d_u0 - d_u1;
                break;
            case 2:
                dphenotype[i + n_terminals] = u0 * d_u1 + d_u0 * u1;
                break;
            case 3:
                dphenotype[i + n_terminals] = (d_u0 * u1 - d_u1 * u0) / (u1 * u1);
                break;
            // non arithmetic kernels (unary only)
            default:
                dphenotype[i + n_terminals] = dukernel_list[m_kernels[fidx] - 4](u0) * d_u0;
        }
    }
    return dphenotype;
}

// Second order derivative
std::vector<double> expression::ddphenotype(const std::vector<unsigned> &genotype, const std::vector<double> &phenotype,
                                            const std::vector<double> &d0phenotype,
                                            const std::vector<double> &d1phenotype)
{
    assert(d0phenotype.size() == d1phenotype.size());
    assert(phenotype.size() == d1phenotype.size());
    // Number of terminals (vars and cons)
    unsigned n_terminals = m_nvar + m_ncon;
    // Number of triplets (F idx0, idx1 in the chromosome)
    auto n_triplets = genotype.size() / 3;
    // Size of the return value will be the same as phenotype
    std::vector<double> ddphenotype(phenotype.size(), 0.);
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
        auto dd_u0 = ddphenotype[genotype[3 * i + 1]];
        auto dd_u1 = ddphenotype[genotype[3 * i + 2]];
        // Retrieve the function
        auto fidx = genotype[3 * i];
        switch (m_kernels[fidx]) {
            // +
            case 0:
                ddphenotype[i + n_terminals] = dd_u0 + dd_u1;
                break;
            // -
            case 1:
                ddphenotype[i + n_terminals] = dd_u0 - dd_u1;
                break;
            // *
            case 2:
                ddphenotype[i + n_terminals] = dd_u0 * u1 + dd_u1 * u0 + d0_u0 * d1_u1 + d1_u0 * d0_u1;
                break;
            // /
            case 3:
                ddphenotype[i + n_terminals] = ((dd_u0 * u1 + d0_u0 * d1_u1 - d0_u1 * d1_u0 - dd_u1 * u0) * u1 * u1
                                                - 2 * u1 * d1_u1 * (d0_u0 * u1 - d0_u1 * u0))
                                               / u1 / u1 / u1 / u1;
                break;
            // non arithmetic kernels (unary only)
            default:
                ddphenotype[i + n_terminals] = ddukernel_list[m_kernels[fidx] - 4](u0) * d0_u0 * d1_u0
                                               + dukernel_list[m_kernels[fidx] - 4](u0) * dd_u0;
        }
    }
    return ddphenotype;
}

std::vector<double> expression::mse(const std::vector<unsigned> &genotype, const std::vector<double> &cons,
                                    const std::vector<std::vector<double>> &xs, const std::vector<double> &ys)
{
    auto N = xs.size();
    assert(m_nvar == xs[0].size());
    std::vector<double> retval(m_nvar + m_ncon + genotype.size() / 3, 0u);
    for (decltype(xs.size()) i = 0u; i < N; ++i) {
        // compute all values in the phenotype (u0, u1, u2, .... un) at xs[i], cons
        auto squared_err = phenotype(genotype, xs[i], cons);
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

// computes mse, dmse and ddmse in one go
void expression::ddmse(const std::vector<unsigned> &genotype, const std::vector<double> &cons,
                       const std::vector<std::vector<double>> &xs, const std::vector<double> &ys,
                       std::vector<double> &mse, std::vector<double> &dmse, std::vector<double> &ddmse)
{
    auto N = xs.size();

    assert(mse.size() == genotype.size() / 3 + m_nvar + m_ncon);

    std::vector<double> ph;
    std::vector<double> dph;
    std::vector<double> ddph;

    std::fill(mse.begin(), mse.end(), 0.);
    std::fill(dmse.begin(), dmse.end(), 0.);
    std::fill(ddmse.begin(), ddmse.end(), 0.);

    // For each point
    for (decltype(xs.size()) i = 0u; i < xs.size(); ++i) {
        // The value of each expression in the single point
        ph = phenotype(genotype, xs[i], cons);
        // The derivative value w.r.t. c (idx 1)
        dph = dphenotype(genotype, ph, 1);
        // The second derivative value w.r.t. c c
        ddph = ddphenotype(genotype, ph, dph, dph);
        // We will now store in ph. dph, ddph respectively, the mse, dmse and ddmse
        // NOTE: The order of the next loops counts and should not be touched.
        // yi-\hat y_i
        for (auto j = 0u; j < ph.size(); ++j) {
            ph[j] -= ys[i];
        }
        // 2 ((yi-\hat y_i)d2ydc2+(dy/dc)^2)
        for (auto j = 0u; j < ddph.size(); ++j) {
            ddph[j] = 2 * (ph[j] / 2. * ddph[j] + dph[j] * dph[j]);
        }
        // 2 (yi-\hat y_i)dydc
        for (auto j = 0u; j < dph.size(); ++j) {
            dph[j] = 2 * ph[j] * dph[j];
        }
        // finally (yi-\hat y_i)^2
        for (auto j = 0u; j < ph.size(); ++j) {
            ph[j] *= ph[j];
        }
        // And we accumulate
        std::transform(mse.begin(), mse.end(), ph.begin(), mse.begin(), std::plus<double>());
        std::transform(dmse.begin(), dmse.end(), dph.begin(), dmse.begin(), std::plus<double>());
        std::transform(ddmse.begin(), ddmse.end(), ddph.begin(), ddmse.begin(), std::plus<double>());
    }
    std::transform(mse.begin(), mse.end(), mse.begin(), [N](double &c) { return c / N; });
    std::transform(dmse.begin(), dmse.end(), dmse.begin(), [N](double &c) { return c / N; });
    std::transform(ddmse.begin(), ddmse.end(), ddmse.begin(), [N](double &c) { return c / N; });
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
} // namespace dsyre