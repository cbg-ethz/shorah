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

#include <getopt.h>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string.h>
#include <string>
#include <limits>

#include <assert.h>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#define UNUSED(expr) (void)(expr)

// data for fetch_func and pileup_func
typedef struct
{
    int /*b, e,*/ win, inc, min_overlap, /*len, cov,*/ max;
    htsFile* in;
//    uint32_t beg, end;  // b, e store begin and end of entire region of interest
//    bam_plbuf_t* buf;   // beg, end store begin and end of current window
    faidx_t* fai;
//    void* rd;               // stores working read
//    void* rLen;             // stores coverage of read starts on reference
    std::ofstream outFile;  // output file for current window
    std::ofstream reads;
} tmpstruct_t;

// to consider, based on samtools/sam.c
//     if (mask < 0) mask = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
//     else mask |= BAM_FUNMAP;
//     while ((ret = samread(fp, b)) >= 0) {
//         // bam_plp_push() itself now filters out unmapped reads only
//         if (b->core.flag & mask) b->core.flag |= BAM_FUNMAP;


// callback for bam_plbuf_init() (for making windows)
//static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t* pl, void* data)
static int pileup_func(const uint32_t win_b, const uint32_t win_e, const uint32_t pos, const bam_pileup1_t* pl, char* rd)
{
    const int i = pos - win_b;     // base position within window
    char* mm = rd + i;             // location of base
    if (pos >= win_b && pos <= win_e) {  // make sure read position within window
        if (pl->is_del == 1) {                 // remove deletions, replacing with reference
            // char *rb;
            // int td;
            // td = pl->b->core.tid;
            // rb = fai_fetch(tmp->fai, tmp->in->header->target_name[td], &tmp->len);
            *mm = '-';
        } else {  // replace 'N' in tmp.read with current read base
            // based on samtools/bam.h
            *mm = seq_nt16_str[bam_seqi(bam_get_seq(pl->b), pl->qpos)];
        }
    }
    return 0;
}


// main
int main(int argc, char* argv[])
{
    int c = 0;           // for parsing command line arguments
    int ref = 0;         // index for retrieving reference chromosome

    std::ofstream covOut;  // stores window coverages

    tmpstruct_t tmp;  // data for callback functions

    char help_string[] =
        "\nUsage: b2w [options] <in.bam> <in.fasta> region\n\nOptions:\n\t-w: window length "
        "(INT)\n\t-i: increment (INT)\n\t-m: minimum overlap (INT)\n\t-x: max reads starting at a "
        "position (INT)\n\t-h: show this help\n\n";

    while ((c = getopt(argc, argv, "w:i:m:x:h")) != EOF) {
        switch (c) {
            case 'w':
                tmp.win = atof(optarg);
                break;
            case 'i':
                tmp.inc = atof(optarg);
                break;
            case 'm':
                tmp.min_overlap = atof(optarg);
                break;
            case 'x':
                tmp.max = atof(optarg);
                break;
            case 'h':
                fprintf(stdout, "%s", help_string);
                exit(EXIT_SUCCESS);
            default:
                fprintf(stderr, "%s", help_string);
                exit(EXIT_FAILURE);
        }
    }
    if (argc < 11 || argc > 12) {
        fprintf(stderr, "%s", help_string);
        return 1;
    }

    if (NULL == (tmp.in = hts_open(argv[optind], "r"))) {  // open bam file
        fprintf(stderr, "Failed to open BAM file %s\n", argv[optind]); 
        return 2;
    }

    if (fai_build(argv[optind + 1])) {  // generate reference index
        fprintf(stderr, "Failed to index FASTA file %s\n", argv[optind + 1]); 
        return 3;
    }
    if (NULL == (tmp.fai = fai_load(argv[optind + 1]))) {
        fprintf(stderr, "Failed to load FASTA file %s\n", argv[optind + 1]); 
        return 2;
    }

    const int idx_min_shift = 0; // NOTE BAI sufficient for up to 2^29-1. If we move from virus to organism with longer chromosomes (plants ?), we should switch to HTS_FMT_CSI [default idx_min_shift = 14]
    int iBuild = bam_index_build(argv[optind], idx_min_shift);  // generate bam index
    UNUSED(iBuild);
    hts_idx_t* idx;
    if (NULL == (idx = hts_idx_load(argv[optind], idx_min_shift ? HTS_FMT_CSI : HTS_FMT_BAI))) {  // load bam index
        fprintf(stderr, "BAM indexing file is not available.\n");
        return 3;
    }

    
    bam_hdr_t* header = sam_hdr_read(tmp.in);
    const int32_t n = header->n_targets;          // number of chromosomes
    const uint32_t* ln = header->target_len;    // chromosome lengths
//    tmp.buf = bam_plbuf_init(pileup_func, &tmp);  // initiate pileup buffer
    bam_plp_t plp_iter = bam_plp_init(0, 0);
    covOut.open("coverage.txt", std::ios::out);   // open file to store window coverage
    tmp.reads.open("reads.fas", std::ios::out);

#if 0 // from samtools
/*!
      @abstract Retrieve the alignments that are overlapped with the
      specified region.  (For BAM files only; see also samfetch() in sam.h.)

      @discussion A user defined function will be called for each
      retrieved alignment ordered by its start position.

      @param  fp    BAM file handler
      @param  idx   pointer to the alignment index
      @param  tid   chromosome ID as is defined in the header
      @param  beg   start coordinate, 0-based
      @param  end   end coordinate, 0-based
      @param  data  user provided data (will be transferred to func)
      @param  func  user defined function
     */
int bam_fetch(bamFile fp, const bam_index_t *idx, int tid, int beg, int end, void *data, bam_fetch_f func)
{
    int ret;
    bam_iter_t iter;
    bam1_t *b;
    b = bam_init1();
    iter = bam_iter_query(idx, tid, beg, end);
    while ((ret = bam_iter_read(fp, iter, b)) >= 0) func(b, data);
    bam_iter_destroy(iter);
    bam_destroy1(b);
    return (ret == -1)? 0 : ret;
}
#endif



    /*!
      @param  tid   chromosome ID as is defined in the header
      @param  roi_b   regio of interest's start coordinate, 0-based
      @param  roi_e   regio of interest's end coordinate, 0-based
      
      captures:
      @param  idx   pointer to the alignment index
      @param  ln    chromosome lengths from header
      @param  header bam header
      @param  tmp   ref to tmpstruct_t data
     */ 
    auto process_ROI = [idx, ln, &header, &covOut, &plp_iter, &tmp] (const int tid, const int roi_b, const uint32_t roi_e) {
        const auto lnth = ln[tid]; // this chromosome ID's length  as is defined in the header

        {
            int rLen[lnth];  // list for read start coverage
            for (int j = 0; j < lnth; j++) {
                rLen[j] = 0;
            }
//            bam_fetch(tmp.in->x.bam, idx, i, tmp.b, tmp.e, &tmp,
//                     fetch_func2);  // fetch and write reads
            // based on samtools/bam.c
            // These BAM iterator functions work only on BAM files.  To work with either
            // BAM or CRAM files use the sam_index_load() & sam_itr_*() functions.
            bam1_t *b = bam_init1();
            hts_itr_t* iter = bam_itr_queryi(idx, tid, roi_b, roi_e);
            while (hts_itr_next(tmp.in->fp.bgzf, iter, b, 0) >= 0) {
                // callback two (for printing reads)
                uint32_t Rstart = b->core.pos;
                // based on samtools/bam.h
                // Rend = bam_calend(&b->core, bam_get_cigar(b));
                int readLen = b->core.n_cigar ? bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b)) : 1;
                uint32_t Rend = Rstart + readLen;
                
                int* cv = rLen + Rstart;
                *cv += 1;
                if ((*cv < tmp.max) && (readLen >= tmp.min_overlap)) {
                    char rd[readLen + 1];
                    for (int i = 0; i < readLen; i++) {
                        rd[i] = 'N';
                    }
                    rd[readLen] = 0;
#if 0
int bam_plbuf_push(const bam1_t *b, bam_plbuf_t *buf)
{
    int ret, n_plp, tid, pos;
    const bam_pileup1_t *plp;
    ret = bam_plp_push(buf->iter, b);
    if (ret < 0) return ret;
    while ((plp = bam_plp_next(buf->iter, &tid, &pos, &n_plp)) != 0)
        buf->func(tid, pos, n_plp, plp, buf->data);
    return 0;
}
#endif
//                    bam_plbuf_push(b, tmp->buf);
                    {
                        int n_plp, tid, pos;
                        const bam_pileup1_t *plp;
                        if (0 <= bam_plp_push(plp_iter, b)) {
                            while ((plp = bam_plp_next(plp_iter, &tid, &pos, &n_plp)) != 0)
                                pileup_func(Rstart, Rend, pos, plp, rd);
                                //pileup_func(tid, pos, n_plp, plp, &tmp);
//                    bam_plbuf_push(0, tmp->buf);
                            bam_plp_push(plp_iter, 0);
                            while ((plp = bam_plp_next(plp_iter, &tid, &pos, &n_plp)) != 0)
                                pileup_func(Rstart, Rend, pos, plp, rd);
                                //pileup_func(tid, pos, n_plp, plp, &tmp);
                        } // TODO trap errors
                    }
//                    bam_plbuf_reset(tmp->buf);
                    bam_plp_reset(plp_iter);

                    // based on samtools/bam_plbuf.c

                    tmp.reads << bam_get_qname(b) << "\t" << roi_b << "\t" << roi_e << "\t" << Rstart + 1
                            << "\t" << Rend << "\t" << rd << "\n";
                }
            }

            bam_itr_destroy(iter);
            bam_destroy1(b);
        }

        uint32_t win_b = roi_b;                // window's start coordinate, 0-based
        uint32_t win_e = win_b + tmp.win - 1;  // window's end coordinate, 0-based (tmp.win : command line parameter for window length)
        while (win_e <= roi_e) {  // make windows over length of region (or chromosome)
            unsigned cov = 0;
            int rLen[lnth];  // list for read start coverage
            for (int j = 0; j < lnth;
                 j++) {  // restart read start coverage count for each window
                rLen[j] = 0;
            }
            char* filename = NULL;  // filename for window files
            asprintf(&filename, "w-%s-%u-%u.reads.fas",  // read window filename
                    header->target_name[tid], win_b + 1, win_e + 1);
            // TODO clean ref_name of special caracters - that would be an alternative to processing everything with regex down the line
            // BUG the solution currently used by ShoRAH can still fail when path '/' (or on windows '\\' and ':') characters are present in the bam's header->target_name

            tmp.outFile.open(filename, std::ios::out);

//                bam_fetch(tmp.in->x.bam, idx, i, tmp.beg, tmp.end, &tmp,
//                          fetch_func1);         // fetch and write reads
            {   // based on samtools/bam.c
                // These BAM iterator functions work only on BAM files.  To work with either
                // BAM or CRAM files use the sam_index_load() & sam_itr_*() functions.
                bam1_t *b = bam_init1();
                hts_itr_t* iter = bam_itr_queryi(idx, tid, win_b, win_e);
                while (hts_itr_next(tmp.in->fp.bgzf, iter, b, 0) >= 0) {
                    // callback one for bam_fetch() (for making windows)
                    int overlap = 0;        // read overlap with window
                    unsigned wSize;              // window size
                    uint32_t Rstart = b->core.pos, Rend = Rstart;
                    // based on samtools/bam.h
                    // Rend = bam_calend(&b->core, bam_get_cigar(b));
                    Rend += b->core.n_cigar ? bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b)) : 1;
                    wSize = tmp.win;

                    int* cv = rLen + Rstart;  // tracks nr of reads starting at ref pos
                    *cv += 1;

                    if (*cv < tmp.max) {                              // limit coverage
                        // TODO: overlap = std::min(Rend, win_e) - std::max(Rstart, win_b) + 1
                        if (Rstart <= win_b && Rend >= win_e) {  // calculate overlap
                            overlap = win_e - win_b + 1;
                        } else if (Rstart >= win_b && Rend <= win_e) {
                            overlap = Rend - Rstart + 1;
                        } else if (Rstart <= win_b && Rend <= win_e) {
                            overlap = Rend - win_b + 1;
                        } else if (Rstart >= win_b && Rend >= win_e) {
                            overlap = win_e - Rstart + 1;
                        }
                        if (overlap >= tmp.min_overlap) {
                            //char rd[wSize + 1]; store blank read segment of 'N's
                            char *rd = new char[wSize + 1];
                            for (unsigned i = 0; i < wSize; i++) {
                                rd[i] = 'N';
                            }
                            rd[wSize] = 0;
//                        bam_plbuf_push(b, tmp->buf);  // push a read to pileup buffer
                            {
                                int n_plp, tid, pos;
                                const bam_pileup1_t *plp;
                                if (0 <= bam_plp_push(plp_iter, b)) { // push a read to pileup buffer
                                    while ((plp = bam_plp_next(plp_iter, &tid, &pos, &n_plp)) != 0)
                                        pileup_func(win_b, win_e, pos, plp, rd);
                                        //pileup_func(tid, pos, n_plp, plp, &tmp);
                                    cov += 1;                 // increment window coverage
//                      bam_plbuf_push(0, tmp->buf);
                                    bam_plp_push(plp_iter, 0); // but make sure only 1 read at a time
                                    while ((plp = bam_plp_next(plp_iter, &tid, &pos, &n_plp)) != 0)
                                        pileup_func(win_b, win_e, pos, plp, rd);
                                    //pileup_func(tid, pos, n_plp, plp, &tmp);
                                } // TODO trap errors
                            }
//                    bam_plbuf_reset(tmp->buf);
                            bam_plp_reset(plp_iter);

                            tmp.outFile  // write read header
                                << '>' << bam_get_qname(b) << ' ' << Rstart << '\n';
                            tmp.outFile  // write read
                                << rd << "\n";
                            delete [] rd;
                        }
                    }
                }

                bam_itr_destroy(iter);
                bam_destroy1(b);
            }

            if (cov > 0) {             // ignore empty windows
                covOut << filename << "\t"  // write out coverage information
                        << header->target_name[tid] << "\t" << win_b + 1 << "\t"
                        << win_e + 1 << "\t" << cov << "\n";
            }
            tmp.outFile.close();
            free(filename);

            win_b += tmp.inc;  // increment windows
            win_e += tmp.inc;

            if (tmp.inc == 0) {
                break;
            }  // OZ: in amplicon mode, one window only
        }
    };
    
    if (argc == optind + 2) {          // region not specified
        for (int i = 0; i < n; i++) {  // take each chromosome in turn
            process_ROI(i,  //
                    0,      // region of interest's start coordinate, 0-based 
                    ln[i]); // chromosome length == region of interest's end
        }
    } else {  // parse specific region
#if 0
int bam_parse_region(bam_header_t *header, const char *str, int *ref_id, int *beg, int *end)
{
    const char *name_lim = hts_parse_reg(str, beg, end);
    if (name_lim) {
        char *name = malloc(name_lim - str + 1);
        memcpy(name, str, name_lim - str);
        name[name_lim - str] = '\0';
        *ref_id = bam_name2id(header, name);
        free(name);
    }
    else {
        // not parsable as a region, but possibly a sequence named "foo:a"
        *ref_id = bam_name2id(header, str);
        *beg = 0; *end = INT_MAX;
    }
    if (*ref_id == -1) return -1;
    return *beg <= *end? 0 : -1;
}
       
#endif
//        bam_parse_region(tmp.in->header, argv[optind + 2], &ref, &tmp.b, &tmp.e);
        // based on samtools/bam_aux.c
        int ref_id = -1, roi_b = 0, roi_e = std::numeric_limits<typeof(roi_e)>::max();
        const char* roi_str = argv[optind + 2];
        {
            const char* name_lim = hts_parse_reg(roi_str, &roi_b, &roi_e);
            if (name_lim) { // valid name ?
                if (roi_b <= roi_e) { // valid range ?
                    const std::ptrdiff_t nl = name_lim - roi_str; // get only the 'name' part (everything before the colon ':')
                    assert(nl > 0);
                    char *name = strndup(roi_str, nl);
                    ref_id = bam_name2id(header, name);
                    free(name);
                }
            } else
                // not parsable as a region, but possibly a sequence named "foo:a"                
                ref_id = bam_name2id(header, roi_str);
        }

        if (ref_id < 0) {
            fprintf(stderr, "Invalid region %s\n", roi_str);
            return 4;
        }
        roi_b -= 3 * tmp.inc;  // make sure start and end of region are covered by 3 windows
        roi_e += 3 * tmp.inc;
        if (roi_b < 0) {
            roi_b = 0;
        }
        if ((uint32_t)roi_e >= ln[ref_id]) {
            roi_e = ln[ref_id] - 1;
        }

        process_ROI(ref_id, roi_b, roi_e);
#if 0
        int rLen[ln[ref]];
        tmp.len = ln[ref];

        for (uint32_t j = 0; j < ln[ref]; j++) {
            rLen[j] = 0;
        }
        tmp.rLen = rLen;
        bam_fetch(tmp.in->x.bam, idx, ref, tmp.b, tmp.e, &tmp,
                  fetch_func2);  // fetch and write reads
        tmp.beg = tmp.b;
        tmp.end = tmp.beg + tmp.win - 1;

        while (tmp.end < (uint32_t)tmp.e) {  // make windows over region
            tmp.cov = 0;
            for (uint32_t j = 0; j < ln[ref];
                 j++) {  // restart read start coverage count for each window
                rLen[j] = 0;
            }
            tmp.rLen = rLen;
            sprintf(filename, "w-%s-%d-%d.reads.fas",  // read window filename
                    tmp.in->header->target_name[ref], tmp.beg + 1, tmp.end + 1);
            tmp.outFile.open(filename, std::ios::out);
            bam_fetch(tmp.in->x.bam, idx, ref, tmp.beg, tmp.end, &tmp,
                      fetch_func1);  // fetch and write reads
            if (tmp.cov != 0) {      // ignore empty windows
                covOut << filename << "\t" << tmp.in->header->target_name[ref] << "\t"
                       << tmp.beg + 1 << "\t" << tmp.end + 1 << "\t" << tmp.cov << "\n";
            }
            tmp.outFile.close();
            tmp.beg += tmp.inc;  // increment windows
            tmp.end += tmp.inc;
            if (tmp.inc == 0) {
                break;
            }  // OZ: in amplicon mode, one window only
        }
#endif
    }
//    bam_plbuf_destroy(tmp.buf);
    hts_idx_destroy(idx); 
    hts_close(tmp.in);
    covOut.close();
    tmp.reads.close();
    exit(EXIT_SUCCESS);
}
