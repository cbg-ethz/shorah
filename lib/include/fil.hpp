#ifndef FIL_HPP
#define FIL_HPP

#include <string>

int fil(const std::string in_bam, const std::string file_to_append, 
    const std::string out_file_prefix, const float sigma, const int max_depth, 
    const bool amplicon_mode, const bool drop_snv);

#endif