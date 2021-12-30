#ifndef B2W_HPP
#define B2W_HPP

#include <string>

int b2w(
    const std::string bam_file, const std::string fasta_name,  
    int win, int inc, int min_overlap, int max, int cov_thrd, 
    bool skip_indel, const std::string region_name);

#endif