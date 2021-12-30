/*
 # Copyright 2007-2012
 # Niko Beerenwinkel,
 # Nicholas Eriksson,
 # Moritz Gerstung,
 # Lukas Geyrhofer,
 # Kerensa McElroy
 # Osvaldo Zagordi,
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
b2w.cpp
Kerensa McElroy
UNSW, Sydney Australia
adapted from samtools/calDep.c
 ******************************/

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <vector>
#include <getopt.h>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "b2w.hpp"

#define UNUSED(expr) (void)(expr)

// data for fetch_func and pileup_func
typedef struct
{
    int win, inc, min_overlap, max, cov_thrd;
    bool skip_indel;
} paramstruct_t;

// to consider, based on samtools/sam.c
//     if (mask < 0) mask = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
//     else mask |= BAM_FUNMAP;
//     while ((ret = samread(fp, b)) >= 0) {
//         // bam_plp_push() itself now filters out unmapped reads only
//         if (b->core.flag & mask) b->core.flag |= BAM_FUNMAP;


// callback for bam_plbuf_init() (for making windows)
static int pileup_func(const uint32_t win_b, const uint32_t win_e, const uint32_t pos, const bam_pileup1_t* pl, char* rd, const bool skip_indel = false)
{
    const int i = pos - win_b;     // base position within window
    char* mm = rd + i;             // location of base
    if (pos >= win_b && pos <= win_e) {  // make sure read position within window
        // for the exact signification of "is_del" and "indel", see https://www.biostars.org/p/104301/#104311
        if (pl->is_del == 1) {                 // remove deletions, replacing with reference
            // char *rb;
            // int td;
            // td = pl->b->core.tid;
            // rb = fai_fetch(tmp->fai, tmp->in->header->target_name[td], &tmp->len);
            *mm = '-';
        } else {  // replace 'N' in rd with current read base
            if ( (! skip_indel) || (pl->indel == 0))
                // based on samtools/bam.h
                *mm = seq_nt16_str[bam_seqi(bam_get_seq(pl->b), pl->qpos)];
        }
    }
    return 0;
}


int b2w(
    const std::string bam_file, const std::string fasta_name,  
    int win, int inc, int min_overlap, int max, int cov_thrd, 
    bool skip_indel, const std::string region_name)
{
    int c = 0;           // for parsing command line arguments

    htsFile* inFile;
    faidx_t* fai;
    std::ofstream reads;
    std::ofstream covOut;  // stores window coverages

    paramstruct_t param = {  // data for callback functions
        // based on default from cli.py and shotgun.py
        win,   // window size
        inc, // size / shifts
        min_overlap, // min_overlap = size * win_min_ext
        max, // max_c = max_coverage / size
        cov_thrd,  // coverage threshold
        skip_indel   // (historically default behaviour of shorah is to never skip insertions)
    };

    if (NULL == (inFile = hts_open(bam_file.c_str(), "r"))) {  // open bam file
        std::fprintf(stderr, "Failed to open BAM file %s\n", bam_file.c_str());
        return 2;
    }

    if (fai_build(fasta_name.c_str())) {  // generate reference index
        std::fprintf(stderr, "Failed to index FASTA file %s\n", fasta_name.c_str());
        return 3;
    }
    if (NULL == (fai = fai_load(fasta_name.c_str()))) {
        std::fprintf(stderr, "Failed to load FASTA file %s\n", fasta_name.c_str());
        return 2;
    }

    const int idx_min_shift = 0; // NOTE BAI sufficient for up to 2^29-1. If we move from virus to organism with longer chromosomes (plants ?), we should switch to HTS_FMT_CSI [default idx_min_shift = 14]
    int iBuild = bam_index_build(bam_file.c_str(), idx_min_shift);  // generate bam index
    UNUSED(iBuild);
    hts_idx_t* idx;
    if (NULL == (idx = sam_index_load(inFile, bam_file.c_str()))) {  // load bam index
        std::fprintf(stderr, "BAM indexing file is not available.\n");
        return 3;
    }

    bam_hdr_t* header = sam_hdr_read(inFile);
    const int32_t n = header->n_targets;          // number of chromosomes
    const uint32_t* ln = header->target_len;    // chromosome lengths
    bam_plp_t plp_iter = bam_plp_init(0, 0);
    covOut.open("coverage.txt", std::ios::out);   // open file to store window coverage
    reads.open("reads.fas", std::ios::out);

    /*!
      @abstract Common code to process a region of interest

      @param  tid   chromosome ID as is defined in the header
      @param  roi_b   regio of interest's start coordinate, 0-based
      @param  roi_e   regio of interest's end coordinate, 0-based

      @param  idx      pointer to the alignment index
      @param  ln       chromosome lengths from header
      @param  header   bam header
      @param  inFile   input bam file
      @param  faidx_t  input indexed fasta file;
      @param  reads    all read in ROI written here
      @param  covOut   file storing window coverage
      @param  plp_iter pileup_iterator
      @param  param    command-line provided parameters
     */
    auto process_ROI = [idx, ln, &header, &inFile, &reads, &covOut, &plp_iter, param] (const int tid, const int roi_b, const uint32_t roi_e) {
        const auto lnth = ln[tid]; // this chromosome ID's length  as is defined in the header

        {
            std::vector<int> rLen(lnth, 0);  // list for read start coverage
            std::vector<char> rd;
            rd.reserve(301);  // reasonable initial buffer size to avoid frequent re-allocation. Most HTS tend to use protocols that yield between 100 and 250 bp
            // based on samtools/bam.c:bam_fetch
            bam1_t *b = bam_init1();
            hts_itr_t* iter = sam_itr_queryi(idx, tid, roi_b, roi_e);
            if (iter == nullptr) {
                std::cerr << "error while processing ROI: " << tid << ":" << roi_b << "-" << roi_e << "\n";
                return;
            }
            while (sam_itr_next(inFile, iter, b) >= 0) {
                // callback two (for printing reads)
                uint32_t Rstart = b->core.pos;
                // based on samtools/bam.h:bam_calend
                int readLen = b->core.n_cigar ? bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b)) : 1;
                uint32_t Rend = Rstart + readLen;

                int* cv = rLen.data() + Rstart;
                *cv += 1;
                if ((*cv < param.max) && (readLen >= param.min_overlap)) {
                    rd.assign(readLen,'N');
                    rd.push_back('\0');
                    // based on samtools/bam_plbuf.c
                    {
                        int n_plp, tid, pos;
                        const bam_pileup1_t *plp;
                        if (0 <= bam_plp_push(plp_iter, b)) {
                            bam_plp_push(plp_iter, 0);
                            while ((plp = bam_plp_next(plp_iter, &tid, &pos, &n_plp)) != 0)
                                pileup_func(Rstart, Rend, pos, plp, rd.data(), param.skip_indel);
                        } // TODO trap errors
                    }
                    bam_plp_reset(plp_iter);

                    reads << bam_get_qname(b) << "\t" << roi_b << "\t" << roi_e << "\t" << Rstart + 1
                            << "\t" << Rend << "\t" << rd.data() << "\n";
                }
            }

            bam_itr_destroy(iter);
            bam_destroy1(b);
        }

        unsigned wSize = param.win;         // window size
        char *rd = new char[wSize + 1];     // store blank read segment of 'N's

        for(uint32_t win_b = roi_b,             // window's start coordinate, 0-based
                     win_e = win_b + wSize - 1; // window's end coordinate, 0-based

            win_e <= roi_e;  // make windows over length of region (or chromosome)

            win_b += param.inc,  // increment windows
            win_e += param.inc
        ) {
            unsigned cov = 0;
            std::vector<int> rLen(lnth, 0);  // list for read start coverage
            char* filename = NULL;  // filename for window files
            asprintf(&filename, "w-%s-%u-%u.reads.fas",  // read window filename
                    header->target_name[tid], win_b + 1, win_e + 1);
            // TODO clean ref_name of special caracters - that would be an alternative to processing everything with regex down the line
            // BUG the solution currently used by ShoRAH can still fail when path '/' (or on windows '\\' and ':') characters are present in the bam's header->target_name

            std::ofstream outFile;  // output file for current window
            outFile.open(filename, std::ios::out);
            {   // based on samtools/bam.c:bam_fetch
                bam1_t *b = bam_init1();
                hts_itr_t* iter = sam_itr_queryi(idx, tid, win_b, win_e);
                if (iter == nullptr) {
                    std::cerr << "error while processing WIN: " << tid << ":" << win_b << "-" << win_e << "\n";
                    return;
                }
                while (sam_itr_next(inFile, iter, b) >= 0) {
                    // callback one for bam_fetch() (for making windows)
                    int overlap = 0;        // read overlap with window
                    uint32_t Rstart = b->core.pos, Rend = Rstart;
                    // based on samtools/bam.h:bam_calend
                    Rend += b->core.n_cigar ? bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b)) : 1;

                    int* cv = rLen.data() + Rstart;  // tracks nr of reads starting at ref pos
                    *cv += 1;

                    if (*cv < param.max) {                              // limit coverage
                        overlap = std::min(Rend, win_e) - std::max(Rstart, win_b) + 1; // calculate overlap
                        if (overlap >= param.min_overlap) {
                            for (unsigned i = 0; i < wSize; i++) {
                                rd[i] = 'N';
                            }
                            rd[wSize] = 0;
                            {
                                int n_plp, tid, pos;
                                const bam_pileup1_t *plp;
                                if (0 <= bam_plp_push(plp_iter, b)) { // push a read to pileup buffer
                                    cov += 1;                 // increment window coverage
                                    bam_plp_push(plp_iter, 0); // but make sure only 1 read at a time
                                    while ((plp = bam_plp_next(plp_iter, &tid, &pos, &n_plp)) != 0)
                                        pileup_func(win_b, win_e, pos, plp, rd, param.skip_indel);
                                } // TODO trap errors
                            }
                            bam_plp_reset(plp_iter);

                            outFile  // write read header
                                << '>' << bam_get_qname(b) << ' ' << Rstart << '\n';
                            outFile  // write read
                                << rd << "\n";
                        }
                    }
                }

                bam_itr_destroy(iter);
                bam_destroy1(b);
            }
            outFile.close();

            if (cov > param.cov_thrd) {     // ignore windows with low-coverage
                covOut << filename << "\t"  // write out coverage information
                        << header->target_name[tid] << "\t" << win_b + 1 << "\t"
                        << win_e + 1 << "\t" << cov << "\n";
            }
            free(filename);

            if (param.inc == 0) {
                break;
            }  // OZ: in amplicon mode, one window only
        }
        delete [] rd;
    };

    if (region_name.empty()) {          // region not specified
        for (int i = 0; i < n; i++)
            process_ROI(i,  // take each chromosome in turn
                    0,      // region of interest's start coordinate, 0-based
                    ln[i]); // chromosome length == region of interest's end
    } else {  // parse specific region
        // based on samtools/bam_aux.c:bam_parse_region
        int ref_id = -1, roi_b = 0, roi_e = std::numeric_limits<decltype(roi_e)>::max();
        {
            const char* name_lim = hts_parse_reg(region_name.c_str(), &roi_b, &roi_e);
            if (name_lim) { // valid name ?
                if (roi_b <= roi_e) { // valid range ?
                    const std::ptrdiff_t nl = name_lim - region_name.c_str(); // get only the 'name' part (everything before the colon ':')
                    assert(nl > 0);
                    char *name = strndup(region_name.c_str(), nl);
                    ref_id = bam_name2id(header, name);
                    free(name);
                }
            } else
                // not parsable as a region, but possibly a sequence named "foo:a"
                ref_id = bam_name2id(header, region_name.c_str());
        }

        if (ref_id < 0) {
            std::fprintf(stderr, "Invalid region %s\n", region_name.c_str());
            return 4;
        }
        roi_b -= 3 * param.inc;  // make sure start and end of region are covered by 3 windows
        roi_e += 3 * param.inc;
        if (roi_b < 0) {
            roi_b = 0;
        }
        if ((uint32_t)roi_e >= ln[ref_id]) {
            roi_e = ln[ref_id] - 1;
        }

        process_ROI(ref_id, roi_b, roi_e);
    }
    bam_plp_destroy(plp_iter);
    bam_hdr_destroy(header);
    hts_idx_destroy(idx); 
    fai_destroy(fai);
    hts_close(inFile);
    covOut.close();
    reads.close();
    return 0;
}
