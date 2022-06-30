// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the heyoka library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DESYRE_EXPRESSION_HPP
#define DESYRE_EXPRESSION_HPP

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <random>
#include <vector>

#include <desyre/detail/visibility.hpp>

namespace desyre
{
class DESYRE_DLL_PUBLIC expression
{
public:
    // Constructor
    expression(unsigned nvar, unsigned ncon, std::vector<unsigned> kernels,
               decltype(std::random_device{}()) seed = std::random_device{}());

    // Generates the constants at random within bounds.
    std::vector<double> random_constants(double lb, double ub);

    // Generates the genotype at random.
    std::vector<unsigned> random_genotype(unsigned length);

    // Computes the phenotype (i.e. the numerical values of all the graph odes)
    std::vector<double> phenotype(const std::vector<unsigned> &genotype, const std::vector<double> &vars,
                                  const std::vector<double> &cons);

    // First order derivatives
    std::vector<double> dphenotype(const std::vector<unsigned> &genotype, const std::vector<double> &phenotype,
                                   unsigned idx);

    // Second order derivatives
    std::vector<double> ddphenotype(const std::vector<unsigned> &genotype, const std::vector<double> &phenotype,
                                    const std::vector<double> &d0phenotype, const std::vector<double> &d1phenotype);

    // Mean Squared Error
    std::vector<double> mse(const std::vector<unsigned> &genotype, const std::vector<double> &cons,
                            const std::vector<std::vector<double>> &xs, const std::vector<double> &ys);

    // Computes mse, dmse and ddmse in one go
    void ddmse(const std::vector<unsigned> &genotype, const std::vector<double> &cons,
               const std::vector<std::vector<double>> &xs, const std::vector<double> &ys, std::vector<double> &mse,
               std::vector<double> &dmse, std::vector<double> &ddmse);

    // Computes the fitness as the best,avg, worst mse over all the nodes
    std::vector<double> fitness(const std::vector<unsigned> &genotype, const std::vector<double> &cons,
                                const std::vector<std::vector<double>> &xs, const std::vector<double> &ys);

    // Mutates the graph
    std::vector<unsigned> mutation(std::vector<unsigned> genotype, unsigned N);

private:
    // Serialization.
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive &ar, unsigned)
    {
        ar &m_nvar;
        ar &m_ncon;
        ar &m_nker;
        ar &m_kernels;
    }

    // data members
    unsigned m_nvar;
    unsigned m_ncon;
    unsigned m_nker;
    std::vector<unsigned> m_kernels;
    std::mt19937 m_rng;
};

} // namespace desyre

#endif