#ifndef PYDSYRE_COMMON_UTILS_HPP
#define PYDSYRE_COMMON_UTILS_HPP

#include <algorithm>
#include <vector>

#include <boost/numeric/conversion/cast.hpp>
#include <pagmo/types.hpp>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "python_includes.hpp"

namespace py = pybind11;

namespace pydsyre
{
// Throw a Python exception of type "type" with associated
// error message "msg".
[[noreturn]] inline void py_throw(PyObject *type, const char *msg)
{
    PyErr_SetString(type, msg);
    throw py::error_already_set();
}

// Utils to expose algo log.
template <typename UDA>
inline py::list generic_log_getter(const UDA &a)
{
    py::list retval;
    for (const auto &t : a.get_log()) {
        retval.append(t);
    }
    return retval;
}

template <typename UDA>
inline void expose_algo_log(py::class_<UDA> &algo_class, const char *doc)
{
    algo_class.def("get_log", &generic_log_getter<UDA>, doc);
}

// Serialization support
template <typename UDX>
inline py::tuple udx_pickle_getstate(const UDX &udx)
{
    // The idea here is that first we extract a char array
    // into which a has been serialized, then we turn
    // this object into a Python bytes object and return that.
    std::ostringstream oss;
    {
        boost::archive::binary_oarchive oarchive(oss);
        oarchive << udx;
    }
    auto s = oss.str();
    return py::make_tuple(py::bytes(s.data(), boost::numeric_cast<py::size_t>(s.size())));
}

template <typename UDX>
inline UDX udx_pickle_setstate(py::tuple state)
{
    // Similarly, first we extract a bytes object from the Python state,
    // and then we build a C++ string from it. The string is then used
    // to deserialized the object.
    if (py::len(state) != 1) {
        py_throw(PyExc_ValueError, ("the state tuple passed for udp/uda deserialization "
                                    "must have 1 element, but instead it has "
                                    + std::to_string(py::len(state)) + " element(s)")
                                       .c_str());
    }

    auto ptr = PyBytes_AsString(state[0].ptr());
    if (!ptr) {
        py_throw(PyExc_TypeError, "a bytes object is needed in the dcgp serialization setstate");
    }

    std::istringstream iss;
    iss.str(std::string(ptr, ptr + py::len(state[0])));
    UDX udx;
    {
        boost::archive::binary_iarchive iarchive(iss);
        iarchive >> udx;
    }

    return udx;
}

} // namespace pydsyre

#endif
