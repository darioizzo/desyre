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

    std::vector<double> phenotype(const std::vector<double> &x, const std::vector<unsigned> &genotype)
    {
        // Size will be the vars+constant values (x) and then the number of triplets F u0 u1
        auto n_triplets = genotype.size() / 3;
        std::vector<double> retval(m_nvar + m_ncon + n_triplets);
        // The u0, u1, ... are the values of variables and constants
        std::copy(x.begin(), x.end(), retval.begin());
        // We loop and for each triplet compute the corresponding function
        for (decltype(n_triplets) i = 0u; i < n_triplets; ++i) {
            auto u0val = retval[genotype[3 * i + 1]];
            auto u1val = retval[genotype[3 * i + 2]];
            auto fidx = genotype[3 * i];
            switch (fidx) {
                case 0:
                    retval[i + x.size()] = add(u0val, u1val);
                    break;
                case 1:
                    retval[i + x.size()] = sub(u0val, u1val);
                    break;
                case 2:
                    retval[i + x.size()] = mul(u0val, u1val);
                    break;
                case 3:
                    retval[i + x.size()] = div(u0val, u1val);
                    break;
                default:
                    throw;
            }
        }
        return retval;
    }

    std::vector<double> mse(const std::vector<unsigned> &genotype, const std::vector<std::vector<double>> &xs,
                            const std::vector<double> &ys)
    {
        assert(m_nvar + m_ncon == xs[0].size());
        std::vector<double> retval(m_nvar + m_ncon + genotype.size() / 3, 0u);
        for (decltype(xs.size()) i = 0u; i < xs.size(); ++i) {
            // compute all values in the phenotype (u0, u1, u2, .... un) at xs[i]
            auto squared_err = phenotype(xs[i], genotype);
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


// Usage ./main n_trials verbosity
int main(int argc, char *argv[])
{
    auto n_trials = std::atoi(argv[1]);
    auto verbosity = std::atoi(argv[2]);

    std::random_device rd;  // only used once to initialise (seed) engine
    std::mt19937 rng(rd()); // random-number engine used (Mersenne-Twister in this case)
    // One variable, no constants +,-,*,/
    expression ex(1, 0, {0, 1, 2, 3});
    // A random genotype
    auto genotype = ex.random_genotype(5);
    fmt::print("[{}]\n", fmt::join(genotype, ", "));
    // And the phenotype evaluated at x = 1.
    auto phenotype = ex.phenotype({1.}, genotype);
    fmt::print("{}\n", fmt::join(phenotype, ", "));
    // Generate Koza data
    std::vector<std::vector<double>> xs;
    std::vector<double> ys;
    generate_data(xs, ys, 10, -3, 3);
    fmt::print("{}\n", fmt::join(xs, ", "));
    fmt::print("{}\n", fmt::join(ys, ", "));
    // Print the error
    auto err = ex.mse(genotype, xs, ys);
    fmt::print("{}\n", fmt::join(err, ", "));
    // Print the fitness
    auto fit = ex.fitness(genotype, xs, ys);
    fmt::print("{}\n", fmt::join(fit, ", "));
    // Mutate (1,2,3 and check)
    fmt::print("[{}]\n", fmt::join(ex.mutation(genotype, 1), ", "));
    fmt::print("[{}]\n", fmt::join(ex.mutation(genotype, 2), ", "));
    fmt::print("[{}]\n", fmt::join(ex.mutation(genotype, 3), ", "));
    // Run the evolution
    // We run n_trials experiments
    auto ERT= 0u;
    for (auto j = 0u; j < n_trials; ++j) {
        // We let each run to convergence
        auto best_x = ex.random_genotype(20);
        auto best_f = ex.fitness(best_x, xs, ys);
        auto count = 0u;
        count++;
        while (true) {
            for (auto i = 0u; i < 4u; ++i) {
                auto new_x = ex.mutation(best_x, 5);
                auto new_f = ex.fitness(new_x, xs, ys);
                count++;
                if (new_f[0] <= best_f[0]) {
                    best_x = new_x;
                    best_f = new_f;
                    // Only if verbosity is > 0
                    if (verbosity>0) {
                        fmt::print("New Best is {} at {} fevals)\n", best_f, count);
                }
                }
            }
            if (best_f[0] < 1e-10) {
                break;
            }
        }
        ERT+=count;
    }
    fmt::print("ERT is {}\n", ERT/n_trials);
    return 0;
}
