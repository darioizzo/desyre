#ifndef PYDSYRE_DOCSTRINGS_HPP
#define PYDSYRE_DOCSTRINGS_HPP

#include <string>

namespace pydsyre
{
// Modules
std::string module_doc();

// Utilities
std::string set_global_seed_doc();
std::string rand_doc();
std::string get_kernel_map_doc();
std::string get_reverse_kernel_map_doc();

// Expression
std::string expression_doc();
std::string expression_random_constants_doc();
std::string expression_random_genotype_doc();
std::string expression_remove_nesting_doc();
std::string expression_phenotype_doc();
std::string expression_complexity_doc();
std::string expression_sphenotype_doc();
std::string expression_mse_doc();
std::string expression_ddmse_doc();
std::string expression_mutate_doc();
std::string expression_mutate_triplets_doc();
std::string expression_get_kernels_idx_doc();
std::string expression_check_genotype_doc();

// sr_problem
std::string sr_problem_doc();

// mes4dsyre
std::string mes4dsyre_doc();
std::string mes4dsyre_get_log_doc();
} // namespace pydsyre

#endif