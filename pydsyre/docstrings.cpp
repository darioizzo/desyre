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

std::string get_kernel_map_doc()
{
    return R"(get_kernel_map()
        
Returns a dictionary mapping kernel names to their global id

Examples:
    >>> import pydsyre as dsy
    >>> kmap = dsy.get_kernel_map()
)";
}

std::string get_reverse_kernel_map_doc()
{
    return R"(get_reverse_kernel_map()
        
Returns a dictionary mapping global ids to the kernel names

Examples:
    >>> import pydsyre as dsy
    >>> kmap = dsy.get_reverse_kernel_map()
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
    >>> ex = dsy.expression(nvars=1, ncons=1, kernels=["sum","mul","diff"])
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

std::string expression_remove_nesting_doc()
{
    return R"(remove_nesting(geno)
        
Modifies geno so that nested non linear expressions such as sin(sin(sin(x))) are no longer in the phenotype.

Args:
    geno: genotype that needs nesting removed from.

Raises:
    TypeError: if length is negative.
    ValueError: if geno is incompatible with the expression.


Returns:
    the genotype with nesting removed. Note that this is not a copy as the input geno will be changed and returned.

Examples:
    >>> import pydsyre as dsy
    >>> ex = dsy.pydsyre_expression(nvars=1, ncons=3, kernels=["sum","mul","diff"])
    >>> geno = ex.random_genotype(length = 10)
    >>> geno = ex.remove_nesting(geno)
)";
}

std::string expression_phenotype_doc()
{
    return R"(phenotype(geno, vars, cons)
        
Computes the numerical value of the phenotype expressed by geno.

Args:
    geno: genotype.
    vars: variables.
    cons: constants.

Raises:
    ValueError: if geno is incompatible with the expression.
    ValueError: if the dimensions of vars and cons are not conistent.

Returns:
    the numeric phenotype.

Examples:
    >>> import pydsyre as dsy
    >>> ex = dsy.pydsyre_expression(nvars=1, ncons=3, kernels=["sum","mul","diff"])
    >>> geno = ex.random_genotype(length = 10)
    >>> cons = ex.random_constants(-1,1)
    >>> phen = ex.phenotype(geno, [1.23], cons)
)";
}

std::string expression_complexity_doc()
{
    return R"(complexity(geno)
        
Computes the complexity of the phenotype expressed by geno.

Args:
    geno: genotype.

Raises:
    ValueError: if geno is incompatible with the expression.

Returns:
    the formulas complexity.

Examples:
    >>> import pydsyre as dsy
    >>> ex = dsy.pydsyre_expression(nvars=1, ncons=3, kernels=["sum","mul","diff"])
    >>> geno = ex.random_genotype(length = 10)
    >>> ex.complexity(geno)
)";
}

std::string expression_sphenotype_doc()
{
    return R"(sphenotype(geno, vars = [], cons = [])
        
Computes the symbolic value of the phenotype expressed by geno.

Args:
    geno: genotype.
    vars: variables. If empty names will be given automatically as "x0" ...
    cons: constants. If empty names will be given automatically as "c0" ...

Raises:
    ValueError: if geno is incompatible with the expression.
    ValueError: if the dimensions of vars and cons are not conistent.

Returns:
    the symbolic phenotype.

Examples:
    >>> import pydsyre as dsy
    >>> ex = dsy.pydsyre_expression(nvars=1, ncons=3, kernels=["sum","mul","diff"])
    >>> geno = ex.random_genotype(length = 10)
    >>> cons = ex.random_constants(-1,1)
    >>> sphen = ex.sphenotype(geno)
)";
}

std::string expression_mse_doc()
{
    return R"(mse(geno, cons, xs, ys)
        
Computes the mean squared error of the phenotype expressed by geno as evaluated on xs.

Args:
    geno: genotype.
    cons: constants value.
    xs: dataset.
    ys: labels.

Raises:
    ValueError: if geno is incompatible with the expression.
    ValueError: if xs or ys are malformed.

Returns:
    the mean squared error on the dataset, labels for each of the phenotype expressions.

Examples:
    >>> import pydsyre as dsy
    >>> import numpy as np
    >>> ex = dsy.expression(nvars=1, ncons=3, kernels=["sum","mul","diff"])
    >>> geno = ex.random_genotype(length = 10)
    >>> cons = ex.random_constants(-1,1)
    >>> xs = np.random.randn(12,1)
    >>> ys = np.random.randn(12)
    >>> ex.mse(geno, cons, xs, ys)
)";
}

std::string expression_ddmse_doc()
{
    return R"(ddmse(geno, cons, xs, ys)
        
Computes the mean squared error of the phenotype expressed by geno as evaluated on xs, its gradient and its hessian with repect to the constants.

Args:
    geno: genotype.
    cons: constants value.
    xs: dataset.
    ys: labels.

Raises:
    ValueError: if geno is incompatible with the expression.
    ValueError: if xs or ys are malformed.

Returns:
    the mean squared error, the gradient and the hessian. If we focus (say) on the i constant
    and indicate the phenotype with [u0, u1, ..., un], then grad[i] will contain the derivative of
    each expression in the the phenotype w.r.t. the constant, i.e. [du0, du1, ..., dun]
    For the Hessian, the relation is more complex as we flatten a symmetric matrix (i and j) -> ij
    following the convention ij = [00 10 11 20 21 22 ....].
    If we seek the derivative w.r.t the constants i and j, then hess[ij] contains the second order
    derivative w.r.t i and j, i.e. [didju0, didju1, didju2, ...]

Examples:
    >>> import pydsyre as dsy
    >>> import numpy as np
    >>> ex = dsy.expression(nvars=1, ncons=3, kernels=["sum","mul","diff"])
    >>> geno = ex.random_genotype(length = 10)
    >>> cons = ex.random_constants(-1,1)
    >>> xs = np.random.randn(12,1)
    >>> ys = np.random.randn(12)
    >>> mse, grad, hess = ex.ddmse(geno, cons, xs, ys)
)";
}

std::string expression_get_kernels_idx_doc()
{
    return R"(get_kernel_idx()
        
Returns:
    the global idx of the kernels used in the expression.

Examples:
    >>> import pydsyre as dsy
    >>> ex = dsy.expression(nvars=1, ncons=3, kernels=["sum","mul","diff"])
    >>> ex.get_kernel_idx()
)";
}

std::string expression_check_genotype_doc()
{
    return R"(check_genotype(geno)
        
Raises:
    ValueError: if geno is incompatible with the expression.

Examples:
    >>> import pydsyre as dsy
    >>> ex = dsy.expression(nvars=1, ncons=3, kernels=["sum","mul","diff"])
    >>> ex.check_genotype([0,1,0,0,1,1,0,4,3])
)";
}

} // namespace pydsyre