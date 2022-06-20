#include <algorithm>
#include <assert.h> /* assert */
#include <iostream>
#include <random>
#include <vector>

#include <fmt/ranges.h>

struct expression {
    expression(unsigned nvar, unsigned ncon, std::vector<unsigned> kernels,
               decltype(std::random_device{}()) seed = std::random_device{}())
        : m_nvar(nvar), m_ncon(ncon), m_kernels(kernels), m_rng(seed)
    {
        m_nker = kernels.size();
    };

    std::vector<unsigned> random_genotype(unsigned length)
    {
        std::vector<unsigned> retval(3 * length);
        unsigned nus = 0u;
        std::uniform_int_distribution<int> uni_ker(0, m_nker - 1);
        for (auto i = 0u; i < length; ++i) {
            // Lets pick a random kernel (but not subtraction as a first pick)
            unsigned funidx = 1;
            while (funidx == 1) {
                funidx = uni_ker(m_rng);
                if (i > 0) break;
            }
            // Lets pick a random u0
            std::uniform_int_distribution<int> uni_us(0, m_nvar + nus - 1);
            unsigned u0idx = uni_us(m_rng);
            unsigned u1idx = u0idx;
            // ... and a random u1 not equal to u2 if the kernel is sub
            while (u1idx == u0idx) {
                u1idx = uni_us(m_rng);
                if (funidx != 1) break;
            }
            retval[3 * i] = funidx;
            retval[3 * i + 1] = u0idx;
            retval[3 * i + 2] = u1idx;
            nus++;
        }
        return retval;
    }

    std::vector<double> phenotype(const std::vector<unsigned> &genotype, const std::vector<double> &x)
    {
        unsigned n_terminals = m_nvar + m_ncon;
        assert(n_terminals == x.size());
        // Size will be the vars+constant values (x) and then the number of triplets F u0 u1
        auto n_triplets = genotype.size() / 3;
        std::vector<double> phenotype(n_terminals + n_triplets);
        // The u0, u1, ... are the values of variables and constants
        std::copy(x.begin(), x.end(), phenotype.begin());
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
                case 4:
                    phenotype[i + n_terminals] = std::sin(u0);
                    break;
                default:
                    throw;
            }
        }
        return phenotype;
    }

    // First order derivatives
    std::vector<double> dphenotype(const std::vector<unsigned> &genotype, const std::vector<double> &phenotype,
                                   unsigned idx)
    {
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
                case 4:
                    dphenotype[i + n_terminals] = std::cos(u0) * d_u0;
                    break;
                default:
                    throw;
            }
        }
        return dphenotype;
    }

    // Second order derivative
    std::vector<double> ddphenotype(const std::vector<unsigned> &genotype, const std::vector<double> &phenotype,
                                    const std::vector<double> &d0phenotype, const std::vector<double> &d1phenotype)
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
                // sin
                case 4:
                    ddphenotype[i + n_terminals] = -std::sin(u0) * d0_u0 * d1_u0 + std::cos(u0) * dd_u0;
                    break;
                default:
                    throw;
            }
        }
        return ddphenotype;
    }

    std::vector<double> mse(const std::vector<unsigned> &genotype, const std::vector<std::vector<double>> &xs,
                            const std::vector<double> &ys)
    {
        assert(m_nvar + m_ncon == xs[0].size());
        std::vector<double> retval(m_nvar + m_ncon + genotype.size() / 3, 0u);
        for (decltype(xs.size()) i = 0u; i < xs.size(); ++i) {
            // compute all values in the phenotype (u0, u1, u2, .... un) at xs[i]
            auto squared_err = phenotype(genotype, xs[i]);
            // subtract ys[i] and square
            for (auto &element : squared_err) {
                element -= ys[i];
                element *= element;
            }
            // Add to retval
            std::transform(retval.begin(), retval.end(), squared_err.begin(), retval.begin(), std::plus<double>());
        }
        return retval;
    }

    std::vector<double> fitness(const std::vector<unsigned> &genotype, const std::vector<std::vector<double>> &xs,
                                const std::vector<double> &ys)
    {
        std::vector<double> retval(3);
        auto errors = mse(genotype, xs, ys);
        retval[2] = std::reduce(errors.begin(), errors.end(), 0.0) / errors.size();
        retval[1] = *std::max_element(errors.begin(), errors.end());
        retval[0] = *std::min_element(errors.begin(), errors.end());
        return retval;
    }

    std::vector<unsigned> mutation(std::vector<unsigned> genotype, unsigned N)
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

    unsigned m_nvar;
    unsigned m_ncon;
    unsigned m_nker;
    std::vector<unsigned> m_kernels;
    std::mt19937 m_rng;
};

inline double koza_quintic(const std::vector<double> &x)
{
    return std::pow(x[0], 5) - 2 * std::pow(x[0], 3) + x[0];
}

void generate_data(std::vector<std::vector<double>> &xs, std::vector<double> &ys, unsigned N, double lb, double ub)
{
    xs.resize(N);
    ys.resize(N);
    // i must be double
    for (double i = 0.; i < N; ++i) {
        xs[i] = {lb + i / (N - 1) * (ub - lb)};
        ys[i] = koza_quintic(xs[i]);
    }
}

using namespace fmt;
// Usage ./main n_trials verbosity
int main(int argc, char *argv[])
{
    // auto n_trials = std::atoi(argv[1]);
    // auto verbosity = std::atoi(argv[2]);

    std::random_device rd;  // only used once to initialise (seed) engine
    std::mt19937 rng(rd()); // random-number engine used (Mersenne-Twister in this case)
    fmt::print("Working with f=[x, c, xc, sin(x), x+sin(x), (x+sin(cx))/x] \n");
    // One variable, one constants +,-,*,/, sin
    expression ex(1, 1, {0, 1, 2, 3, 4});
    // Manually assemble (x + sin(cx)) / x
    std::vector<unsigned> genotype = {2, 0, 1, 4, 2, 0, 0, 0, 3, 3, 4, 0};
    // Compute the phenotype at x=1, c=0.5
    auto f = ex.phenotype(genotype, {1., 1. / 2.});
    print("f: {}\n", f);
    // Compute the dphenotype at x=1, c=0.5
    // d dx
    auto f_dx = ex.dphenotype(genotype, f, 0);
    print("dfdx: {}\n", f_dx);
    // d dc
    auto f_dc = ex.dphenotype(genotype, f, 1);
    print("dfdc: {}\n", f_dc);
    // Compute the ddphenotype at x=1, c=0.5
    // d dxdx
    auto f_dx_dx = ex.ddphenotype(genotype, f, f_dx, f_dx);
    print("dfdxdx: {}\n", f_dx_dx);
    // d dcdc
    auto f_dc_dc = ex.ddphenotype(genotype, f, f_dc, f_dc);
    print("dfdcdc: {}\n", f_dc_dc);
    // d dxdc
    auto f_dx_dc = ex.ddphenotype(genotype, f, f_dx, f_dc);
    print("dfdxdc: {}\n", f_dx_dc);
    // d dcdx
    auto f_dc_dx = ex.ddphenotype(genotype, f, f_dc, f_dx);
    print("dfdcdx: {}\n", f_dc_dx);

    return 0;
}
