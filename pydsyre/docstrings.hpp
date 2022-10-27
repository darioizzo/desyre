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

// gym
// Classic
std::string generate_koza_quintic_doc();
// From our paper
std::string generate_P1_doc();
std::string generate_P2_doc();
std::string generate_P3_doc();
std::string generate_P4_doc();
std::string generate_P5_doc();
std::string generate_P6_doc();
std::string generate_P7_doc();
// From PySR
std::string generate_P8_doc();
// From Vladi paper
std::string generate_kotanchek_doc();
std::string generate_salutowicz_doc();
std::string generate_salutowicz2d_doc();
std::string generate_uball5d_doc();
std::string generate_ratpol3d_doc();
std::string generate_sinecosine_doc();
std::string generate_ripple_doc();
std::string generate_ratpol2d_doc();
// NIST data
std::string generate_chwirut1_doc();
std::string generate_chwirut2_doc();
std::string generate_daniel_wood_doc();
std::string generate_gauss1_doc();
std::string generate_kirby2_doc();
std::string generate_lanczos2_doc();
std::string generate_misra1b_doc();
// MISC data
std::string generate_luca1_doc();
} // namespace pydsyre

#endif