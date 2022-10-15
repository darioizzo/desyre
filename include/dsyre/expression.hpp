// Copyright 2020, 2021, 2022 Francesco Biscani (bluescarni@gmail.com), Dario Izzo (dario.izzo@gmail.com)
//
// This file is part of the dsyre library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSYRE_EXPRESSION_HPP
#define DSYRE_EXPRESSION_HPP

#include <random>
#include <unordered_map>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <dsyre/detail/visibility.hpp>

namespace dsyre
{

class DSYRE_DLL_PUBLIC expression
{
public:
    // Constructor
    expression(unsigned nvar = 1, unsigned ncon = 0, std::vector<std::string> kernels = {"sum", "diff", "mul", "div"});

    // Generates the constants at random within bounds.
    std::vector<double> random_constants(double lb, double ub, std::mt19937 &rng) const;

    // Generates the genotype at random.
    void random_genotype(std::vector<unsigned> &genotype, unsigned length, std::mt19937 &rng) const;

    // Removes nested unary functions (but not inv)
    void remove_nesting(std::vector<unsigned> &g, std::mt19937 &rng) const;

    // Computes the phenotype (i.e. the numerical values of all the graph nodes)
    void phenotype(std::vector<double> &retval, const std::vector<unsigned> &genotype, const std::vector<double> &vars,
                   const std::vector<double> &cons) const;

    // Computes the us complexity
    void complexity(std::vector<unsigned> &retval, const std::vector<unsigned> &geno) const;

    // Computes the symbolic phenotype (i.e. the symbolic expression for all the nodes)
    void sphenotype(std::vector<std::string> &retval, const std::vector<unsigned> &genotype,
                    const std::vector<std::string> &vars = {}, const std::vector<std::string> &cons = {}) const;

    // First order derivatives
    void dphenotype(std::vector<double> &retval, const std::vector<unsigned> &genotype,
                    const std::vector<double> &phenotype, unsigned idx) const;

    // Second order derivatives
    void ddphenotype(std::vector<double> &retval, const std::vector<unsigned> &genotype,
                     const std::vector<double> &phenotype, const std::vector<double> &d0phenotype,
                     const std::vector<double> &d1phenotype) const;

    // Mean Squared Error
    void mse(std::vector<double> &retval, const std::vector<unsigned> &genotype, const std::vector<double> &cons,
             const std::vector<std::vector<double>> &xs, const std::vector<double> &ys) const;

    // Computes mse, dmse and ddmse in one go
    void ddmse(std::vector<double> &mse, std::vector<std::vector<double>> &gradient,
               std::vector<std::vector<double>> &hessian, const std::vector<unsigned> &genotype,
               const std::vector<double> &cons, const std::vector<std::vector<double>> &xs,
               const std::vector<double> &ys) const;

    // Mutates the graph
    std::vector<unsigned> mutate_triplets(const std::vector<unsigned> &genotype, unsigned N, std::mt19937 &rng) const;
    std::vector<unsigned> mutate(const std::vector<unsigned> &genotype, unsigned N, std::mt19937 &rng) const;
    std::vector<unsigned> mutation3(const std::vector<unsigned> &genotype, const std::vector<double> &phenotype,
                                    unsigned N, std::mt19937 &rng);
    void check_genotype(const std::vector<unsigned> &genotype) const;

    const std::vector<unsigned> &get_kernels_idx() const;

    // Streaming operator friendship declaration
    friend std::ostream &operator<<(std::ostream &os, const expression &d);

private:
    // Computes the phenotype (i.e. the numerical values of all the graph nodes) - no checks
    void phenotype_impl(std::vector<double> &retval, const std::vector<unsigned> &genotype,
                        const std::vector<double> &vars, const std::vector<double> &cons) const;
    // Computes the phenotype (i.e. the numerical values of all the graph nodes) - no checks
    void complexity_impl(std::vector<unsigned> &complexity, const std::vector<unsigned> &geno) const;
    // Computes the symbolic phenotype (i.e. the symbolic expression for all the nodes) - no checks
    void sphenotype_impl(std::vector<std::string> &retval, const std::vector<unsigned> &genotype,
                         const std::vector<std::string> &vars, const std::vector<std::string> &cons) const;
    // First order derivatives - no checks
    void dphenotype_impl(std::vector<double> &retval, const std::vector<unsigned> &genotype,
                         const std::vector<double> &phenotype, unsigned idx) const;
    // Second order derivatives - no checks
    void ddphenotype_impl(std::vector<double> &retval, const std::vector<unsigned> &genotype,
                          const std::vector<double> &phenotype, const std::vector<double> &d0phenotype,
                          const std::vector<double> &d1phenotype) const;

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
};

// Streaming operator
DSYRE_DLL_PUBLIC std::ostream &operator<<(std::ostream &os, const expression &d);

// Available kernels. Since the map is built manually in the cpp file, this helper allow access to it.
DSYRE_DLL_PUBLIC std::unordered_map<std::string, unsigned> get_kernel_map();
DSYRE_DLL_PUBLIC std::unordered_map<unsigned, std::string> get_reverse_kernel_map();

} // namespace dsyre

#endif