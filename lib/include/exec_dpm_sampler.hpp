#ifndef EXEC_DPM_SAMPLER_HPP
#define EXEC_DPM_SAMPLER_HPP

#include <string>

int exec_dpm_sampler(const std::string i_filein, const int j_iterations, const double a_alpha, const int t_history, const int R_seed, const double k_cluster_avg_reads, const int K_cluster_start);

#endif