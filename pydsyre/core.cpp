#include <cmath>
#include <random>
#include <string>
#include <vector>

#include <iostream>

#include <boost/optional.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <dsyre/expression.hpp>

#include "docstrings.hpp"

namespace py = pybind11;
using namespace pydsyre;

// We instantiate and init the global random number generator for pydsyre.
std::random_device rd;
std::mt19937 rng(rd());

PYBIND11_MODULE(core, m)
{
    py::options options;
    options.disable_function_signatures();
    m.doc() = module_doc();

    // We expose the global random number generator seeding
    m.def(
        "set_global_seed", [](unsigned seed) { rng.seed(seed); }, set_global_seed_doc().c_str(), py::arg("seed"));
    m.def(
        "rand", []() { return rng(); }, rand_doc().c_str());

    auto exp_
        = py::class_<dsyre::expression>(m, "expression", expression_doc().c_str())
              // Default constructor
              .def(py::init<>())
              // Constructor from number of variables, number of constants, and kernels
              .def(py::init<unsigned, unsigned, std::vector<std::string>>(), py::arg("nvars"), py::arg("ncons"),
                   py::arg("kernels"))
              .def("__repr__",
                   [](const dsyre::expression &ex) -> std::string {
                       std::ostringstream oss;
                       std::cout << ex;
                       return oss.str();
                   })
              .def(
                  "random_constants",
                  [](const dsyre::expression &ex, double lb, double ub) { return ex.random_constants(lb, ub, rng); },
                  py::arg("lb") = -1., py::arg("ub") = 1., expression_random_constants_doc().c_str())
              .def(
                  "random_genotype",
                  [](const dsyre::expression &ex, unsigned length) {
                      std::vector<unsigned> retval;
                      ex.random_genotype(retval, length, rng);
                      return retval;
                  },
                  py::arg("length"), expression_random_genotype_doc().c_str())
              .def(
                  "remove_nesting",
                  [](const dsyre::expression &ex, std::vector<unsigned> &g) {
                      ex.remove_nesting(g, rng);
                      return g;
                  },
                  py::arg("geno"), expression_remove_nesting_doc().c_str())
              .def(
                  "phenotype",
                  [](const dsyre::expression &ex, const std::vector<unsigned> &g, const std::vector<double> &v,
                     const std::vector<double> &c) {
                      std::vector<double> retval;
                      ex.phenotype(retval, g, v, c);
                      return retval;
                  },
                  py::arg("geno"), py::arg("vars"), py::arg("cons"), expression_phenotype_doc().c_str())
              .def(
                  "complexity",
                  [](const dsyre::expression &ex, const std::vector<unsigned> &g) {
                      std::vector<unsigned> retval;
                      ex.complexity(retval, g);
                      return retval;
                  },
                  py::arg("geno"), expression_complexity_doc().c_str())
              .def(
                  "sphenotype",
                  [](const dsyre::expression &ex, const std::vector<unsigned> &g, const std::vector<std::string> &v,
                     const std::vector<std::string> &c) {
                      std::vector<std::string> retval;
                      ex.sphenotype(retval, g, v, c);
                      return retval;
                  },
                  py::arg("geno"), py::arg("vars") = std::vector<std::string>{},
                  py::arg("cons") = std::vector<std::string>{}, expression_sphenotype_doc().c_str())
              .def(
                  "mse",
                  [](const dsyre::expression &ex, const std::vector<unsigned> &geno, const std::vector<double> &cons,
                     const std::vector<std::vector<double>> &xs, const std::vector<double> &ys) {
                      std::vector<double> retval;
                      ex.mse(retval, geno, cons, xs, ys);
                      return retval;
                  },
                  expression_mse_doc().c_str())
              .def(
                  "ddmse",
                  [](const dsyre::expression &ex, const std::vector<unsigned> &geno, const std::vector<double> &cons,
                     const std::vector<std::vector<double>> &xs, const std::vector<double> &ys) {
                      std::vector<double> mse;
                      std::vector<std::vector<double>> grad;
                      std::vector<std::vector<double>> hess;
                      ex.ddmse(mse, grad, hess, geno, cons, xs, ys);
                      py::tuple retval = py::make_tuple(mse, grad, hess);
                      return retval;
                  },
                  expression_ddmse_doc().c_str());
} // namespace details