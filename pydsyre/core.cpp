#include <cmath>
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

PYBIND11_MODULE(core, m)
{
    py::options options;
    options.disable_function_signatures();
    m.doc() = module_doc();

    auto exp_ = py::class_<dsyre::expression>(
                    m, "pydsyre_expression",
                    "An object to generate and manipulate symbolic regression expressions using the dsyre encoding.")
                    // Default constructor
                    .def(py::init<>())
                    // Constructor from number of variables, number of constants, and kernels
                    .def(py::init<unsigned, unsigned, std::vector<std::string>>(), py::arg("nvar"), py::arg("ncon"),
                         py::arg("kernels"), pydsyre_expression_init_doc().c_str())
                    .def("__repr__", [](const dsyre::expression &instance) -> std::string {
                        std::ostringstream oss;
                        std::cout << instance;
                        return oss.str();
                    });

} // namespace details
