#include <string>

#include "docstrings.hpp"

namespace pydsyre
{

std::string module_doc()
{
    return R"(The core functionalities implemented in cpp and exposed to python
)";
}

std::string set_global_seed_doc()
{
    return R"(set_global_seed(seed)
        
Sets the seed for the dsyre global random number generator used behind the scenes by the dsyre (C++) code.

Examples:
    >>> import pydsyre as dsy
    >>> dsy.set_global_seed(seed = 45)
)";
}

std::string rand_doc()
{
    return R"(rand()
        
Returns the next random number from the dsyre global random number generator used behind the scenes by the dsyre (C++) code.

Examples:
    >>> import pydsyre as dsy
    >>> random_unsigned = dsy.rand()
)";
}

std::string expression_doc()
{
    return R"(__init__(nvars, ncons, kernels)
        
Main pydsyre expression class. It represents the encoding of generic symbolic expressions as defined by:

u0 = x0
u1 = x1
u2 = c1
u3 = f3(u_i, u_j) .. i,j < 3
u4 = f4(u_i, u_j) .. i,j < 4
...

Args:
    nvars: number of variables in the expression
    ncons: number of constants (parameters) in the expression
    kernels: list of functions to be used as nonlinearities (e.g. ["sum", "diff"])

Raises:
    TypeError: if *nvars*, *ncons* are negative.
    ValueError: if one of the requested kernel is not implemented. 

Examples:
    >>> import pydsyre as dsy
    >>> ex = dsy.pydsyre_expression(nvars=1, ncons=1, kernels=["sum","mul","diff"])
    >>> print(ex)
    A differentiable expression using the dsyre encoding.
    Number of variables: 1
    Number of constants: 1
    Number of kernels: 3
    Kernels: ["sum", "mul", "diff"]
)";
}

std::string expression_random_constants_doc()
{
    return R"(random_constants(lb = -1., ub = 1.)
        
Returns a list of randomly initialized constants in the selected bounds and compatible with the dsyre exression (i.e. of the correct size)

Args:
    nvars: number of variables in the expression.

Returns:
    a list of randomly initialized constants.

Examples:
    >>> import pydsyre as dsy
    >>> ex = dsy.pydsyre_expression(nvars=1, ncons=3, kernels=["sum","mul","diff"])
    >>> cons = ex.random_constants(lb=-1., ub=1.)
    >>> len(cons)
    3
)";
}

std::string expression_random_genotype_doc()
{
    return R"(random_genotype(lenght)
        
Returns a random genotype compatible with the expression.

Args:
    length: number of expressions generated on top on nvars and ncons (u0, u1, u2,....,u_{length+ncons+nvars}).
            Each will consist of a triplet (f,u_i, u_j), hence the genotype dimension will be 3 * length

Raises:
    TypeError: if *length* is negative.

Returns:
    a random genotype (list of unsigned).

Examples:
    >>> import pydsyre as dsy
    >>> ex = dsy.pydsyre_expression(nvars=1, ncons=3, kernels=["sum","mul","diff"])
    >>> geno = ex.random_genotype(length = 10)
    >>> len(geno)
    30
)";
}

}