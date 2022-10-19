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

.. math::
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
phenotype(geno, xs, cons)
        
Computes the numerical value of the phenotype expressed by geno.

Args:
    geno: genotype.
    vars: variables (single instance).
    xs: variables (multiple instances).
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
        
Computes the complexity of the phenotype expressed by geno. The complexity c(u_k) of each 
expression in the phenotype (i.e. u_k=f(u_i,u_j)) is computed as 1+c(u_i)+c(u_j) for binary operators
and as 1+c(u_i) for unary operators. 

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

std::string expression_mutate_doc()
{
    return R"(mutate(geno, N)
        
Mutates N reandomly selected genes in geno.

Args:
    geno: genotype.
    N: number of genes to mutate.

Returns:
    the mutated genotype.

Examples:
    >>> import pydsyre as dsy
    >>> ex = dsy.pydsyre_expression(nvars=1, ncons=3, kernels=["sum","mul","diff"])
    >>> geno = ex.random_genotype(length = 10)
    >>> ex.mutate(geno, 4)
)";
}

std::string expression_mutate_triplets_doc()
{
    return R"(mutate_triplets(geno, N)
        
Mutates N reandomly selected triplets in geno. Note this will produce three times 
the amount of mutations w.r.t. a normal mutation (using the same N)

Args:
    geno: genotype.
    N: number of triplets to mutate.

Returns:
    the mutated genotype.

Examples:
    >>> import pydsyre as dsy
    >>> ex = dsy.pydsyre_expression(nvars=1, ncons=3, kernels=["sum","mul","diff"])
    >>> geno = ex.random_genotype(length = 10)
    >>> ex.mutate_triplets(geno, 4)
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

std::string sr_problem_doc()
{
    return R"(__init__(xs, ys, ncons, multi_objective)

A symbolic regression problem. This can be used as a UDP in pagmo to search for a
mathematical expression represented in the dsyre encoding and the value of its constants,
so that the input data are reproduced successfully.

The chromosome will start with the constant values and will be followed by the representation of the formula
as triplets of integers.

Args:
    xs: input data. This determines the number of variables in the expression.
    ys: labels.
    ncons: number of constants in the expression.
    multi_objective: when True the problem considers also the formula complexity.

Raises:
    ValueError: if xs and ys are malformed.

)";
}

std::string mes4dsyre_doc()
{
    return R"(mes4dsyre(gen, max_mut, ftol, seed = ranom)

The term Memetic is widely used, in the context of meta-heuristic search, to indicate a synergy between any
population-based approach with local improvement procedures. The resulting algorithms are also referred to, in the
literature, as Baldwinian evolutionary algorithms (EAs), Lamarckian EAs, cultural algorithms, or genetic local
searches. The very same approach, is seen by many just as an hybridization of a global search technique with a
local search technique. Regardless of the terminology and point of view, a memetic approach is applicable to symbolic
regression tasks and able to improve considerably on the long standing issue of finding constants in
Genetic Programming.

In this class we offer an UDA (User Defined Algorithm for the pagmo optimization suite) hybridizing a classic
Evolutionary Strategy with a second order Newton search step able to help finding the best values for the ephemeral constants.

The resulting algorithm is outlined by the following pseudo-algorithm:

* Start from a population (pop) of dimension N
*  while i < gen
*  > > Mutation: create a new population pop2 mutating N times the best individual
*  > > Life long learning: apply a one step of a second order Newton method to each individual (only the continuous part is affected) 
*  > > Reinsertion: set pop to contain the best N individuals taken from pop and pop2

Args:
    gen: number of generations.
    max_mut: maximum number of mutations at each new candidate.
    ftol: the algorithm will exit when the loss is below this tolerance.
    seed: seed used by the internal random number generator (default is random)

Raises:
    ValueError: if  *max_mut* or *gen* are negative, or if *ftol* is negative.

.. note::
    mes4dsyre is tailored to solve :class:`dsyre.symbolic_regression` problems and will not work on different types.
    )";
}

std::string mes4dsyre_get_log_doc()
{
    return R"(    
)";
}

std::string generate_koza_quintic_doc()
{
    return R"(
Generates the data for the classic Koza quintic regression problem.

.. math::
   y = x^5 - 2 x^3 + x

x is sampled in ten equally spaced points in [-3,3].

Returns:
    A tuple (xs,ys) containing problem data.
)";
}

std::string generate_P1_doc()
{
    return R"(Generates the problem P1 from the paper:

Izzo, Dario, Francesco Biscani, and Alessio Mereta. "Differentiable genetic programming." 
European Conference on Genetic Programming. Springer, 2017.

The functional form of such a problem is:

.. math::
     y = x^5 - \pi x^3 + x

x is sampled in ten equally spaced points in [1,3].

Returns:
    A tuple (xs,ys) containing problem data.
)";
}

std::string generate_P2_doc()
{
    return R"(
Generates the problem P2 from the paper:

Izzo, Dario, Francesco Biscani, and Alessio Mereta. "Differentiable genetic programming." 
European Conference on Genetic Programming. Springer, 2017.

The functional form of such a problem is:

.. math::
     y = x^5 - \pi x^3 + \frac{\pi}{x}

x is sampled in ten equally spaced points in [0.1,5].

Returns:
    A tuple (xs,ys) containing problem data.
)";
}

std::string generate_P3_doc()
{
    return R"(
Generates the problem P3 from the paper:

Izzo, Dario, Francesco Biscani, and Alessio Mereta. "Differentiable genetic programming." 
European Conference on Genetic Programming. Springer, 2017.

The functional form of such a problem is:

.. math::
     y = \frac{e x^5 + x^3}{x+1}

x is sampled in ten equally spaced points in [-0.9,1].

Returns:
    A tuple (xs,ys) containing problem data.
)";
}
std::string generate_P4_doc()
{
    return R"(
Generates the problem P4 from the paper:

Izzo, Dario, Francesco Biscani, and Alessio Mereta. "Differentiable genetic programming." 
European Conference on Genetic Programming. Springer, 2017.

The functional form of such a problem is:

.. math::
     y = \sin(\pi x) + \frac 1x

x is sampled in ten equally spaced points in [-1,1].

Returns:
    A tuple (xs,ys) containing problem data.
)";
}
std::string generate_P5_doc()
{
    return R"(
Generates the problem P5 from the paper:

Izzo, Dario, Francesco Biscani, and Alessio Mereta. "Differentiable genetic programming." 
European Conference on Genetic Programming. Springer, 2017.

The functional form of such a problem is:

.. math::
     y = e x^5 - \pi x^3 + x

x is sampled in ten equally spaced points in [1,3].

Returns:
    A tuple (xs,ys) containing problem data.
)";
}
std::string generate_P6_doc()
{
    return R"(
Generates the problem P6 from the paper:

Izzo, Dario, Francesco Biscani, and Alessio Mereta. "Differentiable genetic programming." 
European Conference on Genetic Programming. Springer, 2017.

The functional form of such a problem is:

.. math::
     y = \frac{e x^2 - 1}{\pi (x + 2)}

x is sampled in ten equally spaced points in [-2.1,1].

Returns:
    A tuple (xs,ys) containing problem data.
)";
}
std::string generate_P7_doc()
{
    return R"(
Generates the problem P7 from the paper:

Izzo, Dario, Francesco Biscani, and Alessio Mereta. "Differentiable genetic programming." 
European Conference on Genetic Programming. Springer, 2017.

The functional form of such a problem is:

.. math::
     y = \cos(\pi x) + \sin(e x)

x is sampled in ten equally spaced points in [-1,1].

Returns:
    A tuple (xs,ys) containing problem data.
    )";
}

std::string generate_P8_doc()
{
    return R"(
Generates the problem used as example in the code PySR

The functional form of such a problem is:

.. math::
     y = 2.5382 * \cos(x[3]) + x[0] * x[0] - 0.5;

x is five dimensional and normally sampled.

Returns:
    A tuple (xs,ys) containing problem data.
)";
}

std::string generate_kotanchek_doc()
{
    return R"(
Generates the problem Kotanchek from the paper:
Vladislavleva, Ekaterina J., Guido F. Smits, and Dick Den Hertog.
"Order of nonlinearity as a complexity measure for models generated by symbolic regression via pareto genetic
programming." IEEE Transactions on Evolutionary Computation 13.2 (2008): 333-349. 

The functional form of such a problem is:
.. math::
   y = \frac{e^{-(x_1-1)^2}}{1.2+(x_2-2.5)^2}

:math:`x_1` and :math:`x_2` are sampled in one hundred randomly selected points in [0.3,4]x[0.3,4].

Returns:
    A tuple containing problem data.
)";
}

std::string generate_salutowicz_doc()
{
    return R"(
Generates the problem Salutowicz from the paper:
Vladislavleva, Ekaterina J., Guido F. Smits, and Dick Den Hertog.
"Order of nonlinearity as a complexity measure for models generated by symbolic regression via pareto genetic
programming." IEEE Transactions on Evolutionary Computation 13.2 (2008): 333-349. 

The functional form of such a problem is:
.. math::
   y = e^{-x} x^3 \cos x\sin x (\cos x \sin^2 x - 1)

x is sampled in one hundred points uniformly sampled in [0.5,10].

Returns:
    A tuple containing problem data.
)";
}

std::string generate_salutowicz2d_doc()
{
    return R"(
Generates the problem Salutowicz2D from the paper:
Vladislavleva, Ekaterina J., Guido F. Smits, and Dick Den Hertog.
"Order of nonlinearity as a complexity measure for models generated by symbolic regression via pareto genetic
programming." IEEE Transactions on Evolutionary Computation 13.2 (2008): 333-349. 

The functional form of such a problem is:
.. math::
   y = e^{-x} x^3 \cos x_1\sin x_1 (\cos x_1 \sin^2 x_1 - 1) * (x_2 - 5)

:math:`x_1` and :math:`x_2` are sampled in 601 randomly selected points in [0.05,10]x[0.05,10].

Returns:
    A tuple containing problem data.
)";
}

std::string generate_uball5d_doc()
{
    return R"(
Generates the problem UBall5D from the paper:
Vladislavleva, Ekaterina J., Guido F. Smits, and Dick Den Hertog.
"Order of nonlinearity as a complexity measure for models generated by symbolic regression via pareto genetic
programming." IEEE Transactions on Evolutionary Computation 13.2 (2008): 333-349. 

The functional form of such a problem is:
.. math::
   y = \frac{10}{5 + \sum_{i=1}^5 (x_i-3)^2}

:math:`x_i` are sampled in 1024 randomly selected points in :math:`[0.05,6.05]^5`.

Returns:
    A tuple containing problem data.
)";
}

std::string generate_ratpol3d_doc()
{
    return R"(
Generates the problem RatPol3D from the paper:
Vladislavleva, Ekaterina J., Guido F. Smits, and Dick Den Hertog.
"Order of nonlinearity as a complexity measure for models generated by symbolic regression via pareto genetic
programming." IEEE Transactions on Evolutionary Computation 13.2 (2008): 333-349. 

The functional form of such a problem is:
.. math::
   y = 30 \frac{(x_1 - 3)(x_3 - 1)}{x_2^2(x_1-10)}

:math:`x_1`, :math:`x_2`, :math:`x_3` are sampled in 300 randomly selected points in [0.05,2] x [1,2].

Returns:
    A tuple containing problem data.
)";
}

std::string generate_sinecosine_doc()
{
    return R"(
Generates the problem SineCosine from the paper:
Vladislavleva, Ekaterina J., Guido F. Smits, and Dick Den Hertog.
"Order of nonlinearity as a complexity measure for models generated by symbolic regression via pareto genetic
programming." IEEE Transactions on Evolutionary Computation 13.2 (2008): 333-349. 

The functional form of such a problem is:
.. math::
   y = 6 \sin(x_1)\cos(x_2)

:math:`x_1`, :math:`x_2` are sampled in 30 randomly selected points in [0.1,5.9] x [0.1,5.9].

Returns:
    A tuple containing problem data.
)";
}

std::string generate_ripple_doc()
{
    return R"(
Generates the problem Ripple from the paper:
Vladislavleva, Ekaterina J., Guido F. Smits, and Dick Den Hertog.
"Order of nonlinearity as a complexity measure for models generated by symbolic regression via pareto genetic
programming." IEEE Transactions on Evolutionary Computation 13.2 (2008): 333-349. 

The functional form of such a problem is:
.. math::
   y = (x_1-3)(x_2-3) + 2\sin((x_1-4)(x_2-4))

:math:`x_1`, :math:`x_2` are sampled in 300 randomly selected points in [0.05,6.05] x [0.05,6.05].

Returns:
    A tuple containing problem data.
)";
}

std::string generate_ratpol2d_doc()
{
    return R"(
Generates the problem RatPol2D from the paper:
Vladislavleva, Ekaterina J., Guido F. Smits, and Dick Den Hertog.
"Order of nonlinearity as a complexity measure for models generated by symbolic regression via pareto genetic
programming." IEEE Transactions on Evolutionary Computation 13.2 (2008): 333-349. 

The functional form of such a problem is:
.. math::
   y = \frac{(x_1-3)^4+(x_2-3)^3-(x_2-3)}{(x_2-2)^4+10}

:math:`x_1`, :math:`x_2` are sampled in 50 randomly selected points in [0.05,6.05] x [0.05,6.05].

Returns:
    A tuple containing problem data.
)";
}

std::string generate_chwirut1_doc()
{
    return R"(
These data are the result of a NIST study involving ultrasonic calibration. The response variable is ultrasonic response, 
and the predictor variable is metal distance. (see https://www.itl.nist.gov/div898/strd/nls/data/chwirut1.shtml)

A proposed good model for such a problem is:
.. math::
   y = \frac{e^{-\beta_1 x}}{\beta_2 + \beta_3 x} + \epsilon

Returns:
    A tuple containing problem data.
)";
}

std::string generate_luca1_doc()
{
    return R"(
These data are the result of an industrial study. Since the original data were subject to intellectual property
restrictions the data has been altered adding an unkown amount of noise and artifacts. The abscissa is some unknown
material property while in the y axis we can read a temperature value in Celsius.

Returns:
    A tuple containing problem data.
)";
}

std::string generate_chwirut2_doc()
{
    return R"(
These data are the result of a NIST study involving ultrasonic calibration. The response variable is ultrasonic response, 
and the predictor variable is metal distance. (see https://www.itl.nist.gov/div898/strd/nls/data/chwirut2.shtml)

A proposed good model for such a problem is:
.. math::
   y = \frac{e^{-\beta_1 x}}{\beta_2 + \beta_3 x} + \epsilon
with respect to the problem chwirut1, less points are included here.

Returns:
    A tuple containing problem data.
)";
}

std::string generate_daniel_wood_doc()
{
    return R"(
These data and model are described in Daniel and Wood (1980), and originally published in E.S.Keeping, 
"Introduction to Statistical Inference," Van Nostrand Company, Princeton, NJ, 1962, p. 354. The response variable is energy
radieted from a carbon filament lamp per cm**2 per second, and the predictor variable is the absolute temperature 
of the filament in 1000 degrees Kelvin. (see https://www.itl.nist.gov/div898/strd/nls/data/daniel_wood.shtml)

A proposed good model for such a problem is:
.. math::
   y = \beta_1 x^{\beta_2} + \epsilon

Returns:
    A tuple containing problem data.
)";
}

std::string generate_gauss1_doc()
{
    return R"(
The data are two well-separated Gaussians on a decaying exponential baseline plus normally 
distributed zero-mean noise with variance = 6.25. (see https://www.itl.nist.gov/div898/strd/nls/data/gauss1.shtml)

A proposed good model for such a problem is:
.. math::
   y = \beta_1 e^{-\beta_2 x} + \beta_3 e^{-\frac{(x-\beta_4)^2}{\beta_5^2}} + \beta_6 e^{-\frac{(x-\beta_7)^2}{\beta_8^2}} + \epsilon

Returns:
    A tuple containing problem data.
)";
}

std::string generate_kirby2_doc()
{
    return R"(
These data are the result of a NIST study involving scanning electron microscope line with standards. 151 
observations are included. (see https://www.itl.nist.gov/div898/strd/nls/data/kirby2.shtml)

A proposed good model for such a problem is:
.. math::
   y = \frac{\beta_1 + \beta_2 x + \beta_3 x^2}{1 + \beta_4 x + \beta_5 x^2} + \epsilon

Returns:
    A tuple containing problem data.
)";
}

std::string generate_lanczos2_doc()
{
    return R"(
These data are taken from an example discussed in Lanczos (1956). The data were generated to 6-digits
of accuracy using the formula below. (see https://www.itl.nist.gov/div898/strd/nls/data/lanczos2.shtml)

A good model for such a problem is, trivially:
.. math::
   y = \beta_1 e^{-\beta_2 x} + \beta_3 e^{-\beta_4 x} + \beta_5 e^{-\beta_6 x} + \epsilon

Returns:
    A tuple containing problem data.
)";
}

std::string generate_misra1b_doc()
{
    return R"(
These data are the result of a NIST study involving dental research in monomolecular adsorption. 
The response variable is volume, and the predictor variable is pressure. 14 observations are 
available. (see https://www.itl.nist.gov/div898/strd/nls/data/misra1b.shtml)

A good model for such a problem is:
.. math::
   y = \beta_1 \left( 1 - \frac 1{\left(1 + \beta_1 \frac x2\right)^2}\right) + \epsilon

Returns:
    A tuple containing problem data.
)";
}
} // namespace pydsyre