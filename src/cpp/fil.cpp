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
#include <gsl/gsl_sf.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <functional>
#include <limits>

#include <string.h>
#include <assert.h>

#include <htslib/sam.h>
#include <htslib/faidx.h>

typedef struct
{
    int pos, pos1, SNP, forwM, revM, forwT, revT, maxdepth;
    double sig, pval;
    char nuc;
    htsFile* in;
} tmpstruct_t;

// reinit members of tmpstruct_t
static int re_init_tmp(tmpstruct_t* tmp)
{
    // re-init sane defaults for struct members that should be filled in by pileup callback fonction
    tmp->forwM = tmp->revM = tmp->forwT = tmp->revT = 0;
    tmp->pval = 1.;
    return 0;
}

// callback for bam_plbuf_init()
static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t* pl, void* data)
{
    tmpstruct_t* tmp = (tmpstruct_t*)data;
    std::map<char, int> st_freq_all;
    std::map<char, int> st_freq_Mut;
    st_freq_all['f'] = 0;  // store frequencies of forward and reverse reads
    st_freq_all['r'] = 0;
    st_freq_Mut['f'] = 0;
    st_freq_Mut['r'] = 0;
    if ((int)pos >= tmp->pos && (int)pos < tmp->pos + 1) {
        for (int i = 0; i < n; i++)  // take each read in turn
        {
            char c;
            const bam_pileup1_t* p = pl + i;
            if (p->indel == 0 and p->is_del != 1)
//              c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b),
//                                               p->qpos)];  // get base for this read
                // based on samtools/bam.h
                c = seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)];
            else if (p->is_del == 1)
                c = '-';
            else
                c = '*';
            if (bam_is_rev(p->b) == 0)  // forward read
                st_freq_all['f']++;
            else
                st_freq_all['r']++;  // reverse read
            if (c == tmp->nuc) {     // base is variant
                if (bam_is_rev(p->b) == 0) {
                    st_freq_Mut['f']++;  // forward read
                } else
                    st_freq_Mut['r']++;  // reverse read
            }
        }
        double mean = (double)st_freq_all['f'] /
                      (st_freq_all['f'] +
                       st_freq_all['r']);  // forward read ratio, all reads at this position
        double fr =
            (double)st_freq_Mut['f'] /
            (st_freq_Mut['f'] +
             st_freq_Mut['r']);  // forward read ratio, only reads with variant at this position
        double alpha;
        double beta;
        double sigma = tmp->sig;
        double prMut;
        int m = st_freq_Mut['f'] + st_freq_Mut['r'];  // total reads with variant
        int k;
        double tail = 0.0;
        alpha = mean / sigma;  // alpha and beta for beta binomial
        beta = (1 - mean) / sigma;
        if (alpha < 1E-6 or beta < 1E-6) {
            fprintf(stderr, "All reads in the same orientation?\n");
            return 11;
        }
        if (fr < mean) {
            for (k = 0; k <= st_freq_Mut['f']; k++) {
                tail += gsl_sf_exp((
                    gsl_sf_lngamma((double)m + 1) - gsl_sf_lngamma((double)k + 1) -
                    gsl_sf_lngamma((double)m - (double)k + 1) + gsl_sf_lngamma(1 / sigma) +
                    gsl_sf_lngamma((double)k + (mean * (1 / sigma))) +
                    gsl_sf_lngamma((double)m + ((1 - mean) / sigma) - (double)k) -
                    gsl_sf_lngamma(mean * (1 / sigma)) - gsl_sf_lngamma((1 - mean) / sigma) -
                    gsl_sf_lngamma((double)m + (1 / sigma))));  // calculate cumulative distribution
            }
        } else {
            for (k = st_freq_Mut['f']; k <= m; k++) {
                tail += gsl_sf_exp(
                    (gsl_sf_lngamma((double)m + 1) - gsl_sf_lngamma((double)k + 1) -
                     gsl_sf_lngamma((double)m - (double)k + 1) + gsl_sf_lngamma(1 / sigma) +
                     gsl_sf_lngamma((double)k + (mean * (1 / sigma))) +
                     gsl_sf_lngamma((double)m + ((1 - mean) / sigma) - (double)k) -
                     gsl_sf_lngamma(mean * (1 / sigma)) - gsl_sf_lngamma((1 - mean) / sigma) -
                     gsl_sf_lngamma((double)m + (1 / sigma))));  // the other tail, if required
            }
        }
        tmp->forwM = st_freq_Mut['f'];
        tmp->revM = st_freq_Mut['r'];
        tmp->forwT = st_freq_all['f'];
        tmp->revT = st_freq_all['r'];
        prMut = tail * 2;  // two sided test
        if (prMut > 1) prMut = 1;
        tmp->pval = prMut;  // p value
    }
    return 0;
}

// main
int main(int argc, char* argv[])
{
    tmpstruct_t tmp;
    int c = 0;
    hts_idx_t* idx;
    int amplicon = 0;
    while ((c = getopt(argc, argv, "b:v:x:a")) != EOF) {
        switch (c) {
            case 'b':
                tmp.in = hts_open(optarg, "r");
                idx = hts_idx_load(optarg, HTS_FMT_BAI);
                 // NOTE BAI sufficient for up to 2^29-1. If we move from virus to organism with longer chromosomes (plants ?), we should switch to HTS_FMT_CSI
                break;
            case 'v':
                tmp.sig = atof(optarg);
                break;
            case 'x':
                tmp.maxdepth = atoi(optarg);
                break;
            case 'a':
                amplicon = 1;
        }
    }
    if (tmp.in == 0) {
        fprintf(stderr, "Failed to open BAM file %s\n", argv[1]);
        return 1;
    }
    bam_hdr_t* header = sam_hdr_read(tmp.in);
    if (idx == 0) {
        fprintf(stderr, "BAM indexing file is not available.\n");
        return 2;
    }
//  bam_plbuf_t* buf;
//  buf = bam_plbuf_init(pileup_func, &tmp);
    bam_plp_t plp_iter = bam_plp_init(0, 0);
    bam_plp_set_maxcnt(plp_iter, tmp.maxdepth);
    FILE* fl = fopen("SNV.txt", "rt");
    if (fl == NULL) {
        fprintf(stderr, "Failed to open SNV file.\n");
        return 3;
    }
    // BUG way too many fixed size for my comfort
    char chr[100], ref, var, fr1[7], fr2[7], fr3[7], p1[7], p2[7], p3[7];
    int pos, name;
    char str_buf[200];
    char reg[200];
    char* filename = NULL;
    std::ofstream snpsOut;
    asprintf(&filename, "SNVs_%f.txt", tmp.sig);
    snpsOut.open(filename, std::ios::out);

    // setup reader and writer depending on amplicon mode
    std::function<int(const char*)> input_scanf;
    std::function<void(std::ofstream &)> output_SNPs;
    if (amplicon) {
        // amplicon only has one window, parses a single frequency and posterior and uses fr1 and p1
        // only
        input_scanf = [&] (const char* line) -> bool {
                // BUG scanf formart string completely misses any size-limitation, and is at high-risk of over-flowing (cue in the valgrind-furby)
                return sscanf(line, "%s %d %c %c %s %s", chr, &pos, &ref, &var, fr1, p1) == 6;
            };
        output_SNPs = [&] (std::ofstream &out) {
                out << chr << '\t' << pos << '\t' << ref << '\t' << var << '\t' << fr1 << '\t'
                    << p1 << '\t' << tmp.forwM << '\t' << tmp.revM << '\t' << tmp.forwT << '\t'
                    << tmp.revT << '\t' << tmp.pval << '\n';
            };
    } else {
        input_scanf = [&] (const char* line) -> bool {
                // BUG scanf formart string completely misses any size-limitation, and is at high-risk of over-flowing (cue in the valgrind-furby)
                return sscanf(line, "%s %d %c %c %s %s %s %s %s %s", chr, &pos, &ref, &var, fr1, fr2,
                       fr3, p1, p2, p3) == 10;
            };
        output_SNPs = [&] (std::ofstream &out) {
                out << chr << '\t' << pos << '\t' << ref << '\t' << var << '\t' << fr1 << '\t'
                    << fr2 << '\t' << fr3 << '\t' << p1 << '\t' << p2 << '\t' << p3 << '\t'
                    << tmp.forwM << '\t' << tmp.revM << '\t' << tmp.forwT << '\t' << tmp.revT
                    << '\t' << tmp.pval << '\n';
            };
    }

    // process the input line by line
    while (fgets(str_buf, 200, fl) != NULL) {
        if (input_scanf(str_buf)) {
            sprintf(reg, "%s:%d-%d", chr, pos, pos);
            tmp.nuc = var;

            re_init_tmp(&tmp);
//          bam_parse_region(tmp.in->header, reg, &name, &tmp.pos, &tmp.pos1);
            // based on samtools/bam_aux.c:bam_parse_region
            name = -1; tmp.pos = 0; tmp.pos1 = std::numeric_limits<typeof(tmp.pos1)>::max();
            {
                const char* name_lim = hts_parse_reg(reg, &tmp.pos, &tmp.pos1);
                if (name_lim) { // valid name ?
                    if (tmp.pos <= tmp.pos1) { // valid range ?
                        const std::ptrdiff_t nl = name_lim - reg; // get only the 'name' part (everything before the colon ':')
                        assert(nl > 0);
                        char *tnam = strndup(reg, nl);
                        name = bam_name2id(header, tnam);
                        free(tnam);
                    }
                } else {
                    // not parsable as a region, but possibly a sequence named "foo:a"
                    name = bam_name2id(header, reg);
                    tmp.pos = 0;
                    tmp.pos1 = header->target_len[name] - 1;
                }
            }

            if (name < 0) {
                fprintf(stderr, "Invalid region %s\n", reg);
                return 4;
            }

//          bam_fetch(tmp.in->x.bam, idx, name, tmp.pos, tmp.pos1, buf, fetch_func);
            // based on samtools/bam.c:bam_fetch
            // TODO These BAM iterator functions work only on BAM files.  To work with either BAM or CRAM files use the sam_index_load() & sam_itr_*() functions.
            bam1_t *b = bam_init1();
            hts_itr_t* iter = bam_itr_queryi(idx, name, tmp.pos, tmp.pos1);
            while (hts_itr_next(tmp.in->fp.bgzf, iter, b, 0) >= 0) {
                // callback for bam_fetch()
//              bam_plbuf_push(b, buf);  // push all reads covering variant to pileup buffer
                // based on samtools/bam_plbuf.c
                const bam_pileup1_t *plp;
                if (bam_plp_push(plp_iter, b) < 0) {
                    // TODO trap errors
                    fprintf(stderr, "!");
                } else {
                    // TODO call-back here too ?
                    int n_plp, tid, p_pos;
                    const bam_pileup1_t *plp;
                    while ((plp = bam_plp_next(plp_iter, &tid, &p_pos, &n_plp)) != 0)
                            pileup_func(tid, p_pos, n_plp, plp, &tmp);
                }
            }
            bam_itr_destroy(iter);
            bam_destroy1(b);

//          bam_plbuf_push(0, buf);
            // based on samtools/bam_plbuf.c
            bam_plp_push(plp_iter, 0);
            {
                int n_plp, tid, p_pos;
                const bam_pileup1_t *plp;
                while ((plp = bam_plp_next(plp_iter, &tid, &p_pos, &n_plp)) != 0)
                    pileup_func(tid, p_pos, n_plp, plp, &tmp);
            }
//          bam_plbuf_reset(buf);
            bam_plp_reset(plp_iter);
            output_SNPs(snpsOut);
        }
    }
    snpsOut.close();
    free(filename);
    fclose(fl);
//  bam_plbuf_destroy(buf);
    bam_plp_destroy(plp_iter);
    bam_hdr_destroy(header);
    hts_idx_destroy(idx);
    hts_close(tmp.in);
}
