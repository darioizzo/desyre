#include <string>

#include "docstrings.hpp"

namespace pydsyre
{

std::string module_doc()
{
    return R"(The core functionalities implemented in cpp and exposed to python
)";
}

std::string pydsyre_expression_init_doc()
{
    return R"(Construct a pydsyre expression object from the number of variables, constants and non-linearities (i.e. kernels)
)";
}

}