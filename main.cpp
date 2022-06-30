#include <algorithm>
#include <assert.h> /* assert */
#include <iostream>
#include <random>
#include <vector>

#include <fmt/ranges.h>

using namespace fmt;
using kernel_f_ptr = double (*)(double, double);

// Functions
double my_sin(double a, double b)
{
    return std::sin(a);
}
double my_dsin(double a, double b)
{
    return std::cos(a);
}
double my_ddsin(double a, double b)
{
    return -std::sin(a);
}
double my_cos(double a, double b)
{
    return std::cos(a);
}
double my_dcos(double a, double b)
{
    return -std::sin(a);
}
double my_ddcos(double a, double b)
{
    return -std::cos(a);
}

// Global array of function pointers
kernel_f_ptr kernel_list[] = {my_cos, my_sin};
kernel_f_ptr dkernel_list[] = {my_dcos, my_dsin};
kernel_f_ptr ddkernel_list[] = {my_ddcos, my_ddsin};


struct expression {
    expression(unsigned nvar, unsigned ncon, std::vector<unsigned> kernels,
               decltype(std::random_device{}()) seed = std::random_device{}())
        : m_nvar(nvar), m_ncon(ncon), m_kernels(kernels), m_rng(seed)
    {
        m_nker = kernels.size();
    };

    std::vector<double> random_constants(double lb, double ub)
    {
        std::uniform_real_distribution<> dis(lb, ub);
        std::vector<double> retval(m_ncon, 0.);
        for (auto &el : retval) {
            el = dis(m_rng);
        }
        return retval;
    }

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

    std::vector<double> phenotype(const std::vector<unsigned> &genotype, const std::vector<double> &vars,
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
                default:
                    phenotype[i + n_terminals] = kernel_list[m_kernels[fidx]-4](u0, u1);
            }
        }
        return phenotype;
    }

    // First order derivatives
    std::vector<double> dphenotype(const std::vector<unsigned> &genotype, const std::vector<double> &phenotype,
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
                default:
                    dphenotype[i + n_terminals] = dkernel_list[m_kernels[fidx]-4](u0, u1) * d_u0;
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
                // non arithmetic kernels
                default:
                    ddphenotype[i + n_terminals] = ddkernel_list[m_kernels[fidx]-4](u0, u1) * d0_u0 * d1_u0 + dkernel_list[m_kernels[fidx]-4](u0, u1) * dd_u0;
            }
        }
        return ddphenotype;
    }

    std::vector<double> mse(const std::vector<unsigned> &genotype, const std::vector<double> &cons,
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

    std::vector<double> fitness(const std::vector<unsigned> &genotype, const std::vector<double> &cons,
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
    void dfitness(const std::vector<unsigned> &genotype, const std::vector<double> &cons,
                  const std::vector<std::vector<double>> &xs, const std::vector<double> &ys, std::vector<double> &mse,
                  std::vector<double> &dmse, std::vector<double> &ddmse)
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

void generate_data(std::vector<std::vector<double>> &xs, std::vector<double> &ys, unsigned N, double lb, double ub)
{
    xs.resize(N);
    ys.resize(N);
    // i must be double
    for (double i = 0.; i < N; ++i) {
        xs[i] = {lb + i / (N - 1) * (ub - lb)};
        ys[i] = P4(xs[i]);
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

    // Generate P1 data
    std::vector<std::vector<double>> xs;
    std::vector<double> ys;
    generate_data(xs, ys, 10, -1., 1.);

    // Allocate some stuff
    auto length = 20u;
    std::vector<double> mse(length + 2, 0.), dmse(length + 2, 0.), ddmse(length + 2, 0.), predicted_mse(length + 2, 0.);

    // The expression system 1 var 1 constant , +,-,*,/, sin, cos
    expression ex(1, 1, {0, 1, 2, 3, 4, 5});

    // Run the evolution
    // We run n_trials experiments
    auto ERT = 0u;
    auto n_success = 0u;
    for (auto j = 0u; j < n_trials; ++j) {
        // We let each run to convergence
        auto best_x = ex.random_genotype(length);
        auto best_c = ex.random_constants(-10., 10.);
        auto best_f = ex.fitness(best_x, best_c, xs, ys);
        auto count = 0u;
        count++;
        while (count < restart) {
            for (auto i = 0u; i < 4u; ++i) {
                auto new_x = ex.mutation(best_x, 5);
                auto new_c = best_c;
                // We now have a new candidate genotype new_x and see what a Newton step could produce.
                // 1 - We compute the mse and its derivateives w.r.t. c
                ex.dfitness(new_x, best_c, xs, ys, mse, dmse, ddmse);
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
                    new_c[0] -= dmse[idx] / ddmse[idx];
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
                        print("New Best is {} at {} fevals: c value {})\n", best_f, count, best_c[0]);
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
    return 0;
}