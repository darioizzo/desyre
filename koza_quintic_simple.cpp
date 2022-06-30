#include <algorithm>
#include <assert.h> /* assert */
#include <iostream>
#include <random>
#include <vector>

#include <fmt/ranges.h>

// Binary kernels
double add(double a, double b)
{
    return a + b;
}

double sub(double a, double b)
{
    return a - b;
}

double mul(double a, double b)
{
    return a * b;
}

double div(double a, double b)
{
    return a / b;
}

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
        std::vector<double> phenotype(m_nvar + m_ncon + n_triplets);
        // The u0, u1, ... are the values of variables and constants
        std::copy(x.begin(), x.end(), phenotype.begin());
        // We loop and for each triplet compute the corresponding function
        for (decltype(n_triplets) i = 0u; i < n_triplets; ++i) {
            auto u0 = phenotype[genotype[3 * i + 1]];
            auto u1 = phenotype[genotype[3 * i + 2]];
            auto fidx = genotype[3 * i];
            switch (m_kernels[fidx]) {
                case 0:
                    phenotype[i + n_terminals] = add(u0, u1);
                    break;
                case 1:
                    phenotype[i + n_terminals] = sub(u0, u1);
                    break;
                case 2:
                    phenotype[i + n_terminals] = mul(u0, u1);
                    break;
                case 3:
                    phenotype[i + n_terminals] = div(u0, u1);
                    break;
                default:
                    throw;
            }
        }
        return phenotype;
    }

    // First order derivatives
    std::vector<double> dphenotype(const std::vector<unsigned> &genotype, const std::vector<unsigned> &phenotype,
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
                    dphenotype[i + n_terminals] = add(d_u0, d_u1);
                    break;
                case 1:
                    dphenotype[i + n_terminals] = sub(d_u0, d_u1);
                    break;
                case 2:
                    dphenotype[i + n_terminals] = add(mul(u0, d_u1), mul(d_u0, u1));
                    break;
                case 3:
                    dphenotype[i + n_terminals] = div(sub(mul(d_u0, u1), mul(u0, d_u1)), mul(u1, u1));
                    break;
                default:
                    throw;
            }
        }
        return dphenotype;
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

// Usage ./koza n_trials verbosity
int main(int argc, char *argv[])
{
    auto n_trials = std::atoi(argv[1]);
    auto restart = std::atoi(argv[2]);
    auto verbosity = std::atoi(argv[3]);

    std::random_device rd;  // only used once to initialise (seed) engine
    std::mt19937 rng(rd()); // random-number engine used (Mersenne-Twister in this case)
    // One variable, no constants +,-,*,/
    expression ex(1, 0, {0, 1, 2, 3});
    // Generate Koza data
    std::vector<std::vector<double>> xs;
    std::vector<double> ys;
    generate_data(xs, ys, 10, -3, 3);

    // Run the evolution
    // We run n_trials experiments
    auto ERT = 0u;
    auto n_success = 0u;
    for (auto j = 0u; j < n_trials; ++j) {
        fmt::print(" {}", j);
        fflush(stdout);
        // We let each run to convergence
        auto best_x = ex.random_genotype(50);
        auto best_f = ex.fitness(best_x, xs, ys);
        auto count = 0u;
        count++;
        while (count < restart) {
            for (auto i = 0u; i < 4u; ++i) {
                auto new_x = ex.mutation(best_x, 7);
                auto new_f = ex.fitness(new_x, xs, ys);
                count++;
                if (new_f[0] <= best_f[0]) {
                    best_x = new_x;
                    best_f = new_f;
                    // Only if verbosity is > 0
                    if (verbosity > 0) {
                        fmt::print("New Best is {} at {} fevals)\n", best_f, count);
                    }
                }
            }
            if (best_f[0] < 1e-10) {
                fmt::print(".");
                fflush(stdout);
                n_success++;
                break;
            }
        }
        ERT += count;
    }
    fmt::print("\n\nERT is {}\n", ERT / n_success);
    fmt::print("Number of successful trials is {}\n", n_success);

    return 0;
}
