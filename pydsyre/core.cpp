#include <cmath>

#include <boost/optional.hpp>
#include <pybind11/pybind11.h>
#include <tbb/global_control.h>

#include <dsyre/expression.hpp>

#include "docstrings.hpp"

namespace py = pybind11;

PYBIND11_MODULE(core, m)
{
    py::options options;
    options.disable_function_signatures();
    m.doc() = pydsyre::module_doc();
} // namespace details
