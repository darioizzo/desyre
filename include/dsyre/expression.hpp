// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSYRE_EXPRESSION_HPP
#define DSYRE_EXPRESSION_HPP

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <random>
#include <vector>

#include <dsyre/detail/visibility.hpp>

namespace dsyre
{
class DSYRE_DLL_PUBLIC expression
{
public:
    // Constructor
    expression(unsigned nvar, unsigned ncon, std::vector<std::string> kernels = {"sum", "diff", "mul", "div"},
               decltype(std::random_device{}()) seed = std::random_device{}());

    // Generates the constants at random within bounds.
    std::vector<double> random_constants(double lb, double ub);

    // Generates the genotype at random.
    void random_genotype(std::vector<unsigned> &genotype, unsigned length);

    // Removes nested unary functions (but not inv)
    void remove_nesting(std::vector<unsigned> &g) const;

    // Computes the phenotype (i.e. the numerical values of all the graph odes)
    void phenotype(std::vector<double> &retval, const std::vector<unsigned> &genotype, const std::vector<double> &vars,
                   const std::vector<double> &cons);

    // Computes the symbolic phenotype (i.e. the symbolic expression for all the nodes)
    void sphenotype(std::vector<std::string> &retval, const std::vector<unsigned> &genotype,
                    const std::vector<std::string> &vars = {}, const std::vector<std::string> &cons = {});

    // First order derivatives
    void dphenotype(std::vector<double> &retval, const std::vector<unsigned> &genotype,
                    const std::vector<double> &phenotype, unsigned idx);

    // Second order derivatives
    void ddphenotype(std::vector<double> &retval, const std::vector<unsigned> &genotype,
                     const std::vector<double> &phenotype, const std::vector<double> &d0phenotype,
                     const std::vector<double> &d1phenotype);

    // Mean Squared Error
    std::vector<double> mse(const std::vector<unsigned> &genotype, const std::vector<double> &cons,
                            const std::vector<std::vector<double>> &xs, const std::vector<double> &ys);

    // Computes mse, dmse and ddmse in one go
    void ddmse(const std::vector<unsigned> &genotype, const std::vector<double> &cons,
               const std::vector<std::vector<double>> &xs, const std::vector<double> &ys, std::vector<double> &mse,
               std::vector<std::vector<double>> &dmse, std::vector<std::vector<double>> &ddmse);

    // Computes the fitness as the best,avg, worst mse over all the nodes
    std::vector<double> fitness(const std::vector<unsigned> &genotype, const std::vector<double> &cons,
                                const std::vector<std::vector<double>> &xs, const std::vector<double> &ys,
                                std::vector<double> &errors);

    // Mutates the graph
    std::vector<unsigned> mutation(const std::vector<unsigned> &genotype, unsigned N);
    std::vector<unsigned> mutation2(const std::vector<unsigned> &genotype, unsigned N);
    std::vector<unsigned> mutation3(const std::vector<unsigned> &genotype, const std::vector<double> &phenotype,
                                    unsigned N);

    void check_genotype(const std::vector<unsigned> &genotype) const;

    const std::vector<unsigned> &get_kernels_idx() const;

private:
    // Computes the phenotype (i.e. the numerical values of all the graph odes) - no checks
    void phenotype_impl(std::vector<double> &retval, const std::vector<unsigned> &genotype,
                        const std::vector<double> &vars, const std::vector<double> &cons);
    // Computes the symbolic phenotype (i.e. the symbolic expression for all the nodes) - no checks
    void sphenotype_impl(std::vector<std::string> &retval, const std::vector<unsigned> &genotype,
                         const std::vector<std::string> &vars, const std::vector<std::string> &cons);
    // First order derivatives - no checks
    void dphenotype_impl(std::vector<double> &retval, const std::vector<unsigned> &genotype,
                         const std::vector<double> &phenotype, unsigned idx);
    // Second order derivatives - no checks
    void ddphenotype_impl(std::vector<double> &retval, const std::vector<unsigned> &genotype,
                          const std::vector<double> &phenotype, const std::vector<double> &d0phenotype,
                          const std::vector<double> &d1phenotype);

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
    mutable std::mt19937 m_rng;
};

} // namespace dsyre

#endif