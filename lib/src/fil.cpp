/*
 # ETH Zurich

 # This file is part of ShoRAH.
 # ShoRAH is free software: you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.

 # ShoRAH is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.

 # You should have received a copy of the GNU General Public License
 # along with ShoRAH.  If not, see <http://www.gnu.org/licenses/>.
*/

/******************************
fil.cpp
Kerensa McElroy
ETH Zurich
adapted from samtools/calDep.c
 ******************************/
#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <set>
#include <algorithm>
#include <exception>
#include <cassert>
#include <cmath>
#include <cfenv>
#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "fil.hpp"

#define UNUSED(expr) (void)(expr)

typedef struct {
    int maxdepth;
    double sig;
    bool skip_indel;
} paramstruct_t;

struct StrandCounts {
    StrandCounts() : forward(0), reverse(0) {}

    int forward;
    int reverse;
};

bam_plp_t _fetch_reads(htsFile* inFile, bam_hdr_t* header, hts_idx_t* idx,
                       const std::string& chr, size_t pos_begin, size_t length,
                       const paramstruct_t& params) {
    // based on samtools/bam_aux.c:bam_parse_region
    int targetid; // targetid : target (as header's numeric form)
    if ((targetid = bam_name2id(header, chr.c_str())) < 0) {
        std::cerr << "Invalid region " << chr << std::endl;
        throw std::exception();
    }

    bam_plp_t plp_iter = bam_plp_init(0, 0);
    bam_plp_set_maxcnt(plp_iter, params.maxdepth);

    /*
    * Get all the reads that cover the interesting region (on chromosomes/target "targetid", position "pos"
    * And store them on a pileup
    */
    // based on samtools/bam.c:bam_fetch
    bam1_t *b = bam_init1();
    // pos_begin is 1-based
    hts_itr_t* iter = sam_itr_queryi(idx, targetid, pos_begin-1, pos_begin-1+length);
    while (sam_itr_next(inFile, iter, b) >= 0) {
        // callback for bam_fetch()
        // based on samtools/bam_plbuf.c
        if (bam_plp_push(plp_iter, b) < 0) {
            // TODO trap errors
            std::cerr << "!";
        }
    }
    bam_itr_destroy(iter);
    bam_destroy1(b);

    // based on samtools/bam_plbuf.c
    bam_plp_push(plp_iter, 0);

    return plp_iter;
}

void _update_intersection(std::set<std::string>& covering_reads,
                          std::set<std::string>& covering_reads_pos,
                          bool first_iteration) {
    if (first_iteration) {
        std::swap(covering_reads, covering_reads_pos);
    } else {
        // https://stackoverflow.com/a/13448094
        std::set<std::string> intersection;
        std::set_intersection(
            covering_reads.begin(), covering_reads.end(),
            covering_reads_pos.begin(), covering_reads_pos.end(),
            std::inserter(intersection, intersection.begin()));
        covering_reads = intersection;
    }
}

StrandCounts _coverage(htsFile* inFile, bam_hdr_t* header, hts_idx_t* idx,
                       const std::string& chr, size_t pos_begin, size_t length,
                       const paramstruct_t& params) {
    bam_plp_t plp_iter = _fetch_reads(inFile, header, idx, chr, pos_begin,
                                      length, params);

    // NOTE: using read IDs and separate sets for forward and reverse reads,
    //       since bam1_t* in pileup point to ephemeral memory locations
    std::set<std::string> covering_reads_forward, covering_reads_reverse;
    int n_plp, tid, p_pos;
    const bam_pileup1_t *plp;

    bool first_iteration = true;
    while ((plp = bam_plp_next(plp_iter, &tid, &p_pos, &n_plp)) != 0) {
        ++p_pos;  // convert 0-based to 1-based indexing
        if ((p_pos < pos_begin) || (pos_begin + length - 1 < p_pos))
            continue;

        std::set<std::string> covering_reads_pos_forward, covering_reads_pos_reverse;
        for (int i = 0; i < n_plp; i++) { // take each read in turn
            const bam_pileup1_t* p = plp + i;
            if (bam_is_rev(p->b))
                covering_reads_pos_reverse.insert(
                    std::string(bam_get_qname(p->b)));
            else
                covering_reads_pos_forward.insert(
                    std::string(bam_get_qname(p->b)));
        }

        _update_intersection(covering_reads_forward,
                             covering_reads_pos_forward, first_iteration);
        _update_intersection(covering_reads_reverse,
                             covering_reads_pos_reverse, first_iteration);
        first_iteration = false;
    }
    // bam_plp_reset(plp_iter);
    bam_plp_destroy(plp_iter);

    StrandCounts counts;
    counts.forward = covering_reads_forward.size();
    counts.reverse = covering_reads_reverse.size();
    return counts;
}

StrandCounts _count_matching_deletions(htsFile* inFile, bam_hdr_t* header,
                                       hts_idx_t* idx,
                                       const std::string& chr, size_t pos_begin,
                                       size_t deletion_length,
                                       const paramstruct_t& params) {
    bam_plp_t plp_iter = _fetch_reads(inFile, header, idx, chr, pos_begin, 1,
                                      params);

    StrandCounts counts;
    int n_plp, tid, p_pos;
    const bam_pileup1_t *plp;
    while ((plp = bam_plp_next(plp_iter, &tid, &p_pos, &n_plp)) != 0) {
        if (p_pos + 1 == pos_begin) {  // p_pos 0-based
            for (int i = 0; i < n_plp; i++) { // take each read in turn
                const bam_pileup1_t* p = plp + i;
                // for an explanation of "is_del" and "indel", see 
                // https://www.biostars.org/p/104301/#104311
                if (-p->indel == deletion_length) {
                    if (bam_is_rev(p->b))
                        ++counts.reverse;
                    else
                        ++counts.forward;
                }
            }
        }
    }
    // bam_plp_reset(plp_iter);
    bam_plp_destroy(plp_iter);

    return counts;
}

StrandCounts _count_matching_substitutions(htsFile* inFile, bam_hdr_t* header,
                                           hts_idx_t* idx, const std::string& chr,
                                           size_t pos, char var,
                                           const paramstruct_t& params) {
    bam_plp_t plp_iter = _fetch_reads(inFile, header, idx, chr, pos, 1,
                                      params);

    StrandCounts counts;
    int n_plp, tid, p_pos;
    const bam_pileup1_t *plp;
    while ((plp = bam_plp_next(plp_iter, &tid, &p_pos, &n_plp)) != 0) {
        if (p_pos + 1 == pos) {  // p_pos 0-based
            for (int i = 0; i < n_plp; i++) { // take each read in turn
                const bam_pileup1_t* p = plp + i;
                // for an explanation of "is_del" and "indel", see 
                // https://www.biostars.org/p/104301/#104311
                if (!p->is_del) {
                    char c = seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)];
                    if (c == var) {
                        if (bam_is_rev(p->b))
                            ++counts.reverse;
                        else
                            ++counts.forward;
                    }
                }
            }
        }
    }
    // bam_plp_reset(plp_iter);
    bam_plp_destroy(plp_iter);

    return counts;
}


int fil(const std::string in_bam, const std::string file_to_append, 
    const std::string out_file_prefix, const float sigma, const int max_depth, 
    const bool amplicon_mode, const bool drop_snv
) {
    paramstruct_t params = { 
        max_depth,
        sigma,
        drop_snv
    };
    int amplicon = amplicon_mode == true ? 1 : 0;

    hts_idx_t* idx = NULL;
    htsFile* inFile = NULL;

    inFile = hts_open(in_bam.c_str(), "r");
    idx = sam_index_load(inFile, in_bam.c_str());

    if (inFile == NULL) {
        std::cerr << "Failed to open BAM file  " << in_bam.c_str() << std::endl;
        return 1;
    }
    bam_hdr_t* header = sam_hdr_read(inFile);
    if (idx == NULL) {
        std::cerr << "BAM indexing file is not available." << std::endl;
        return 2;
    }

    std::ifstream fl (file_to_append, std::ios_base::in); 
    if (fl.fail()) {
        std::cerr << "Failed to open SNV file." << std::endl;
        return 3;
    }
    std::string str_buf;

    std::ostringstream oss_fn;
    oss_fn << out_file_prefix << std::fixed << std::setprecision(6) << params.sig << ".tsv";
    std::string filename = oss_fn.str();
    std::ofstream snpsOut;
    snpsOut.open(filename, std::ios::out);

    // process the input line by line (one SNV after the other)
    while (std::getline(fl, str_buf)) {
        /*
         * General workflow
         *
         * for each SNVs, we fetch all the reads that cover this position,
         * push them on a pileup, and finally read the pileup.
         */

        std::string chr;   // chr:   target chromosone in text form
        int pos; // pos: 1-based position of SNV
        std::string ref;     // reference base (allow for deletions with variable length)
        std::string var;   // variant nucleotide (for future insertions)
        std::istringstream iss(str_buf);
        // amplicon only has one window, parses a single frequency and posterior and uses fr1 and p1
        UNUSED(amplicon);
        //  - except that, because we're just spitting the exact same output, and only *adding* new columns,
        //    not modifying the existing ones, well we actually don't need to care.
        // so we just parse the first couple of interesting columns and keep the line as-is for output
        if (iss >> chr >> pos >> ref >> var) { // Get the SNV location
            /// for the frequencies and posteriors :
            //stdstrings fr1, fr2, fr3, p1, p2, p3;
            // for amplicon:
            //  >> fr1 >> p1
            // for everything else:
            //  >> fr1 >> fr2 >> fr3 >> p1 >> p2 >> p2

            StrandCounts st_freq_all;
            StrandCounts st_freq_Mut;
            assert(ref.size() > 0);
            // A deletion is reported at the preceding position w.r.t the
            // reference. The reference column 'ref' contains the reference
            // base followed by the deleted bases.
            if (ref.size() > 1) {
                int deletion_length = ref.size() - 1;
                assert(pos > 0);

                // construct pileup for coverage (pos + 1, pos + 1 + del_length)
                st_freq_all = _coverage(
                    inFile, header, idx, chr, pos + 1, deletion_length, params);
                st_freq_Mut = _count_matching_deletions(
                    inFile, header, idx, chr, pos, deletion_length, params);
            } else {
                // construct pileup at pos
                st_freq_all = _coverage(
                    inFile, header, idx, chr, pos, 1, params);
                st_freq_Mut = _count_matching_substitutions(
                    inFile, header, idx, chr, pos, var[0], params);
            }

            /*
             * Do stats on the counts - is there a strand bias ?
             */
            
            // init sane defaults
            double pval = 1.;

            if (st_freq_all.forward == 0 and st_freq_all.reverse == 0) {
                // historically this could happen due to B2W and FIL counting indels in different way.
                std::cerr << "Bug: No reads found at position " << pos <<" ?!?" << std::endl;
            } else if (st_freq_Mut.forward == 0 and st_freq_Mut.reverse == 0) {
                std::cerr << "Critical: No mutations found at position " << pos <<" ?!?" << std::endl;
            } else {
                double mean = (double)st_freq_all.forward /
                            (st_freq_all.forward +
                            st_freq_all.reverse);  // forward read ratio, all reads at this position
                double fr =
                    (double)st_freq_Mut.forward /
                    (st_freq_Mut.forward +
                    st_freq_Mut.reverse);  // forward read ratio, only reads with variant at this position
                double alpha;
                double beta;
                double sigma = params.sig;
                double prMut;
                int m = st_freq_Mut.forward + st_freq_Mut.reverse;  // total reads with variant
                int k;
                double tail = 0.0;
                alpha = mean / sigma;  // alpha and beta for beta binomial
                beta = (1 - mean) / sigma;
                
                if (alpha < 1E-6 or beta < 1E-6) {
                    // ...at that point it's not even bias, it's a monopoly ! ;-)
                    std::cerr << "Warning: All reads at position " << pos <<" in the same " << ((alpha < 1E-6) ? "reverse" : "forward") << " orientation ?" << std::endl;
                } else {


                    /*!
                        TODO

                        @abstract P(X[ib] == x) with X ~ BetaBinom(n[ib], mean, sigma) -- formula 4 from doi:10.1186/1471-2164-14-501

                        NOTE as of Boost 1.73, there is still no beta-binomial distribution yet
                        see: https://www.boost.org/doc/libs/1_73_0/libs/math/doc/html/math_toolkit/issues.html#math_toolkit.issues.feature_requests

                        NOTE boost::math::lgamma selected as it has a well documented precision and boost is already used in ShoRAH for dpm_sampler anyway
                        see: https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/sf_gamma/lgamma.html
                      */
                    auto pdfbetabinom = [&](const double x, const double n, const double mean, const double sigma) -> double {
                        return std::exp(
                            std::lgamma(n + 1.)
                            - std::lgamma(x + 1.) 
                            - std::lgamma(n - x + 1.)
                            + std::lgamma(1. / sigma)
                            + std::lgamma(x + (mean * (1. / sigma)))
                            + std::lgamma(n + ((1. - mean) / sigma) - x)
                            - std::lgamma(n + (1. / sigma))
                            - std::lgamma(mean * (1. / sigma))
                            - std::lgamma((1. - mean) / sigma)
                        );
                    };
                    std::feclearexcept(FE_ALL_EXCEPT);
                    if (fr < mean) {
                        for (k = 0; k <= st_freq_Mut.forward; k++)
                            tail += pdfbetabinom (k, m, mean, sigma);
                            // calculate cumulative distribution
                    } else {
                        for (k = st_freq_Mut.forward; k <= m; k++)
                            tail += pdfbetabinom (k, m, mean, sigma);
                            // the other tail, if required
                    }
                    if(std::fetestexcept(FE_UNDERFLOW)) {
                        // catch exception
                        std::cerr << "Warning: at position " << pos << " -- underflow while computing CDF BetaBinom(n=" << m << ", mu=" << mean << ", sigma=" << sigma << ")" << std::endl;
                    }
                    prMut = tail * 2;  // two sided test
                    if (prMut > 1) prMut = 1;
                    pval = prMut;  // p value
                }
            }

            // write the results add
            snpsOut << str_buf << '\t' // we will simply append the new column instead of re-writing them

            //  for amplicons :
            //      << chr << '\t' << pos << '\t' << ref << '\t' << var << '\t'
            //      << fr1 << '\t' << p1 << '\t'
            //  for the normal mode :
            //      << chr << '\t' << pos << '\t' << ref << '\t' << var << '\t'
            //      << fr1 << '\t' << fr2 << '\t' << fr3 << '\t'
            //      << p1 << '\t' << p2 << '\t' << p3 << '\t'

                    << st_freq_Mut.forward << '\t' << st_freq_Mut.reverse << '\t'
                    << st_freq_all.forward << '\t' << st_freq_all.reverse << '\t' << pval << '\n';
        }
    }
    snpsOut.close();
    fl.close();
    bam_hdr_destroy(header);
    hts_idx_destroy(idx);
    hts_close(inFile);
    return 0;
}
