#include <cmath>

#include <boost/optional.hpp>
#include <pybind11/pybind11.h>
#include <tbb/global_control.h>

#include <dsyre/expression.hpp>

#include "docstrings.hpp"

namespace py = pybind11;
using namespace pydsyre;

PYBIND11_MODULE(core, m)
{
    py::options options;
    options.disable_function_signatures();
    m.doc() = module_doc();

    auto exp_ = py::class_<dsyre::expression>(
                    m, "pydsyre_expression",
                    "An object to generate and manipulate symbolic regression expressions using the dsyre encoding.")
                    // Default constructor exposed
                    .def(py::init<>())
                    // From number of variamles, number of constants, kernels
                    .def(py::init<unsigned, unsigned, std::vector<std::string>>(), py::arg("nvar"), py::arg("ncon"),
                         py::arg("kernels"), pydsyre_expression_init_doc().c_str());

} // namespace details
