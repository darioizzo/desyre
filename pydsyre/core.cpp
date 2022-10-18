#include <cmath>
#include <random>
#include <string>
#include <vector>

#include <iostream>

#include <boost/optional.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <dsyre/expression.hpp>
#include <dsyre/gym/gym.hpp>
#include <dsyre/mes4dsyre.hpp>
#include <dsyre/sr_problem.hpp>

#include "common_utils.hpp"
#include "docstrings.hpp"

namespace py = pybind11;
using namespace pydsyre;

// We instantiate and init the global random number generator for pydsyre. The keyword thread_local should make sure
// each thread will maintain a separate copy.
thread_local std::mt19937 rng(std::random_device{}());

using gym_ptr = void (*)(std::vector<std::vector<double>> &, std::vector<double> &);
template <gym_ptr F>
inline void expose_data_from_the_gym(py::module &m, const std::string &name, const std::string &docstring)
{
    m.def(
        name.c_str(),
        []() {
            std::vector<std::vector<double>> xs;
            std::vector<double> ys;
            F(xs, ys);
            return py::make_tuple(xs, ys);
        },
        docstring.c_str());
}

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
    // We expose the kernel maps
    m.def("get_kernel_map", &dsyre::get_kernel_map, get_kernel_map_doc().c_str());
    m.def("get_reverse_kernel_map", &dsyre::get_reverse_kernel_map, get_reverse_kernel_map_doc().c_str());

    // We override the global symbol used to extract a pointer to dsyre::sr_problem from a pagmo::problem
    // This will be called by the UDAs (mes4dsyre) to check for problem compatibility.
    dsyre::details::extract_sr_cpp_py = [](const pagmo::problem &p) -> const dsyre::sr_problem * {
        // We extract a generic python object
        auto py_ptr = p.extract<py::object>();
        // .. if failed, we fail
        if (!py_ptr) {
            return nullptr;
        }
        // .. else we attempt to cast
        const dsyre::sr_problem *retval;
        try {
            retval = py_ptr->cast<const dsyre::sr_problem *>();
        } catch (py::cast_error const &) {
            retval = nullptr;
        }
        return retval;
    };

    // Exposing the expression class
    py::class_<dsyre::expression>(m, "expression", expression_doc().c_str())
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
            "phenotype",
            [](const dsyre::expression &ex, const std::vector<unsigned> &g, const std::vector<std::vector<double>> &xs,
               const std::vector<double> &c) {
                std::vector<std::vector<double>> retval(xs.size());
                for (decltype(xs.size()) i = 0u; i < xs.size(); ++i) {
                    ex.phenotype(retval[i], g, xs[i], c);
                }
                return retval;
            },
            py::arg("geno"), py::arg("xs"), py::arg("cons"))
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
            py::arg("geno"), py::arg("vars") = std::vector<std::string>{}, py::arg("cons") = std::vector<std::string>{},
            expression_sphenotype_doc().c_str())
        .def(
            "mse",
            [](const dsyre::expression &ex, const std::vector<unsigned> &geno, const std::vector<double> &cons,
               const std::vector<std::vector<double>> &xs, const std::vector<double> &ys) {
                std::vector<double> retval;
                ex.mse(retval, geno, cons, xs, ys);
                return retval;
            },
            py::arg("geno"), py::arg("cons"), py::arg("xs"), py::arg("ys"), expression_mse_doc().c_str())
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
            py::arg("geno"), py::arg("cons"), py::arg("xs"), py::arg("ys"), expression_ddmse_doc().c_str())
        .def(
            "mutate",
            [](const dsyre::expression &ex, const std::vector<unsigned> &geno, unsigned N) {
                return ex.mutate(geno, N, rng);
            },
            py::arg("geno"), py::arg("N"), expression_mutate_doc().c_str())
        .def(
            "mutate_triplets",
            [](const dsyre::expression &ex, const std::vector<unsigned> &geno, unsigned N) {
                return ex.mutate(geno, N, rng);
            },
            py::arg("geno"), py::arg("N"), expression_mutate_triplets_doc().c_str())
        .def("get_kernels_idx", &dsyre::expression::get_kernels_idx, expression_get_kernels_idx_doc().c_str())
        .def("check_genotype", &dsyre::expression::check_genotype, py::arg("geno"),
             expression_check_genotype_doc().c_str())
        .def(py::pickle(&pickle_getstate<dsyre::expression>, &pickle_setstate<dsyre::expression>));

    // Exposing the sr_problem class (UDP)
    py::class_<dsyre::sr_problem>(m, "sr_problem", sr_problem_doc().c_str())
        // Default constructor
        .def(py::init<>())
        .def(py::init<const std::vector<std::vector<double>> &, const std::vector<double> &, unsigned,
                      std::vector<std::string>, unsigned, bool>(),
             py::arg("xs"), py::arg("ys"), py::arg("length") = 20,
             py::arg("kernels") = std::vector<std::string>{"sum", "mul", "diff", "inv"}, py::arg("ncons") = 0,
             py::arg("multi_objective") = false)
        .def("get_nobj", &dsyre::sr_problem::get_nobj)
        .def("fitness", &dsyre::sr_problem::fitness)
        .def("get_bounds", &dsyre::sr_problem::get_bounds)
        .def("get_nix", &dsyre::sr_problem::get_nix)
        .def("get_name", &dsyre::sr_problem::get_name)
        .def("pretty", &dsyre::sr_problem::pretty)
        .def("prettier", &dsyre::sr_problem::prettier)
        .def(py::pickle(&pickle_getstate<dsyre::sr_problem>, &pickle_setstate<dsyre::sr_problem>));

    // Exposing the mes4dsyre algorithm class (UDA)
    py::class_<dsyre::mes4dsyre>(m, "mes4dsyre", mes4dsyre_doc().c_str())
        // Default constructor
        .def(py::init<>())
        .def(py::init<unsigned, unsigned, double, unsigned>(), py::arg("gen"), py::arg("max_mut"), py::arg("ftol"),
             py::arg("seed"))
        .def("__repr__", &dsyre::mes4dsyre::get_extra_info)
        .def("evolve", &dsyre::mes4dsyre::evolve)
        .def("set_seed", &dsyre::mes4dsyre::set_seed)
        .def("get_seed", &dsyre::mes4dsyre::get_seed)
        .def("set_verbosity", &dsyre::mes4dsyre::set_verbosity)
        .def("get_verbosity", &dsyre::mes4dsyre::get_verbosity)
        .def("get_name", &dsyre::mes4dsyre::get_name)
        .def("get_log", &generic_log_getter<dsyre::mes4dsyre>, mes4dsyre_get_log_doc().c_str())
        .def(py::pickle(&pickle_getstate<dsyre::mes4dsyre>, &pickle_setstate<dsyre::mes4dsyre>));

    // Making data from the gym available in python
    expose_data_from_the_gym<&dsyre::gym::generate_P0>(m, "generate_P0", generate_koza_quintic_doc());
    // From Our paper
    expose_data_from_the_gym<&dsyre::gym::generate_P1>(m, "generate_P1", generate_P1_doc());
    expose_data_from_the_gym<&dsyre::gym::generate_P2>(m, "generate_P2", generate_P2_doc());
    expose_data_from_the_gym<&dsyre::gym::generate_P3>(m, "generate_P3", generate_P3_doc());
    expose_data_from_the_gym<&dsyre::gym::generate_P4>(m, "generate_P4", generate_P4_doc());
    expose_data_from_the_gym<&dsyre::gym::generate_P5>(m, "generate_P5", generate_P5_doc());
    expose_data_from_the_gym<&dsyre::gym::generate_P6>(m, "generate_P6", generate_P6_doc());
    expose_data_from_the_gym<&dsyre::gym::generate_P7>(m, "generate_P7", generate_P7_doc());
    // From PySR
    expose_data_from_the_gym<&dsyre::gym::generate_P8>(m, "generate_P8", generate_P8_doc());
    // From Vladi paper
    expose_data_from_the_gym<&dsyre::gym::generate_P9>(m, "generate_P9", generate_kotanchek_doc());
    expose_data_from_the_gym<&dsyre::gym::generate_P10>(m, "generate_P10", generate_salutowicz_doc());
    expose_data_from_the_gym<&dsyre::gym::generate_P11>(m, "generate_P11", generate_salutowicz2d_doc());
    expose_data_from_the_gym<&dsyre::gym::generate_P12>(m, "generate_P12", generate_uball5d_doc());
    expose_data_from_the_gym<&dsyre::gym::generate_P13>(m, "generate_P13", generate_ratpol3d_doc());
    expose_data_from_the_gym<&dsyre::gym::generate_P14>(m, "generate_P14", generate_sinecosine_doc());
    expose_data_from_the_gym<&dsyre::gym::generate_P15>(m, "generate_P15", generate_ripple_doc());
    expose_data_from_the_gym<&dsyre::gym::generate_P16>(m, "generate_P16", generate_ratpol2d_doc());
    // NIST data
    expose_data_from_the_gym<&dsyre::gym::generate_P17>(m, "generate_P17", generate_chwirut1_doc());
    expose_data_from_the_gym<&dsyre::gym::generate_P18>(m, "generate_P18", generate_chwirut2_doc());
    expose_data_from_the_gym<&dsyre::gym::generate_P19>(m, "generate_P19", generate_daniel_wood_doc());
    expose_data_from_the_gym<&dsyre::gym::generate_P20>(m, "generate_P20", generate_gauss1_doc());
    expose_data_from_the_gym<&dsyre::gym::generate_P21>(m, "generate_P21", generate_kirby2_doc());
    expose_data_from_the_gym<&dsyre::gym::generate_P22>(m, "generate_P22", generate_lanczos2_doc());
    expose_data_from_the_gym<&dsyre::gym::generate_P23>(m, "generate_P23", generate_misra1b_doc());
    // MISC data
    expose_data_from_the_gym<&dsyre::gym::generate_P24>(m, "generate_P24", generate_luca1_doc());

} // namespace details
