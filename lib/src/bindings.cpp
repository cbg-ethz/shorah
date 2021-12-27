
#include <pybind11/pybind11.h>

#include "fil.hpp"
#include "exec_dpm_sampler.hpp"

namespace py = pybind11;
using namespace py::literals;

PYBIND11_MODULE(libshorah, m) {
    //m.doc() = "pybind11 example plugin"; // optional module docstring TODO

    m.def("fil", &fil, 
        "in_bam"_a, 
        py::arg("sigma") = 0.01, 
        py::arg("max_depth") = 10000, 
        py::arg("amplicon_mode") = false, 
        py::arg("drop_snv") = false, 
    R"pbdoc(
        Args:
            in_bam:
            sigma:
            max_depth:
            amplicon_mode: Toggle modes. Shotgun is the default.
            drop_snv: Drops SNVs that are adjacent to insertions/deletions (alternate behaviour)
        
        Returns:
            0 if successful. 
    )pbdoc");

    m.def("exec_dpm_sampler", &exec_dpm_sampler, 
        "i_filein"_a, 
        "j_iterations"_a,
        "a_alpha"_a,
        "t_history"_a,
        py::arg("R_seed") = 0, 
        py::arg("k_cluster_avg_reads") = 0, 
        py::arg("K_cluster_start") = 0, 
    R"pbdoc(
        Args:
            i_filein:
            j_iterations: sampling iterations
            a_alpha:
            t_history: history time
            R_seed:
            K_cluster_start: Start value for number of clusters, not compatible `k_cluster_avg_reads`
            k_cluster_avg_reads: Average number of reads in each startcluster, not compatible with `k_cluster_avg_reads`
        
        Returns:
            0 if successful. 
    )pbdoc");
}