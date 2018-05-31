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
#include <cstdio>
#include <fstream>
#include <sstream>
#include <map>
#include <string>

#include <gsl/gsl_sf.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>

#define UNUSED(expr) (void)(expr)

typedef struct
{
    int pos, pos1, SNP, forwM, revM, forwT, revT, maxdepth;
    double sig, pval;
    char nuc; // nucleotide that we're looking for in this SNP
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
    tmp.maxdepth = 10000; tmp.sig = 0.01; // cli.py's defaults
    int c = 0;
    hts_idx_t* idx = NULL;
    int amplicon = 0;
    htsFile* inFile = NULL;
    while ((c = getopt(argc, argv, "b:v:x:a")) != EOF) {
        switch (c) {
            case 'b':
                inFile = hts_open(optarg, "r");
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
    if (inFile == NULL) {
        std::fprintf(stderr, "Failed to open BAM file %s\n", argv[1]);
        return 1;
    }
    bam_hdr_t* header = sam_hdr_read(inFile);
    if (idx == NULL) {
        std::fprintf(stderr, "BAM indexing file is not available.\n");
        return 2;
    }
    bam_plp_t plp_iter = bam_plp_init(0, 0);
    bam_plp_set_maxcnt(plp_iter, tmp.maxdepth);
    std::ifstream fl ("SNV.txt", std::ios_base::in);
    if (fl.fail()) {
        std::fprintf(stderr, "Failed to open SNV file.\n");
        return 3;
    }
    std::string str_buf ("");

    char* filename = NULL;
    asprintf(&filename, "SNVs_%f.txt", tmp.sig);
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

        std::string chr; // chr:   target chromosone in text form
        int name, pos;   // name : target (as header's numeric form)  pos: 1-based position of SNV
        char ref, var;   // reference and variant for this SNV
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

            /*
             * Prepare the data
             */
            re_init_tmp(&tmp);
            tmp.nuc = var;
#if 1
            // as seen in htslib/hts.c:hts_parse_reg
            tmp.pos = pos - 1;  // begin is 0-based (so position-1)
            tmp.pos1 = pos;     // end is 1-based (so straight position)
#else       // alternative : based on samtools/bam_aux.c:bam_parse_region
            {
                char* reg = NULL;
                tmp.pos = 0; tmp.pos1 = std::numeric_limits<typeof(tmp.pos1)>::max();
                asprintf(&reg, ":%d-%d", pos, pos);
                hts_parse_reg(reg, &tmp.pos, &tmp.pos1);
                free(reg);
            }
#endif
            // based on samtools/bam_aux.c:bam_parse_region
            if ((name = bam_name2id(header, chr.c_str())) < 0) {
                std::fprintf(stderr, "Invalid region %s\n", chr.c_str());
                return 4;
            }

            /*
             * Phase 1:
             *
             * Get all the reads that cover the interesting region (on chromosomes/target "name", position "pos"
             * And store them on a pileup
             */
            // based on samtools/bam.c:bam_fetch
            // TODO These BAM iterator functions work only on BAM files.  To work with either BAM or CRAM files use the sam_index_load() & sam_itr_*() functions.
            bam1_t *b = bam_init1();
            hts_itr_t* iter = bam_itr_queryi(idx, name, tmp.pos, tmp.pos1);
            while (hts_itr_next(inFile->fp.bgzf, iter, b, 0) >= 0) {
                // callback for bam_fetch()
                // based on samtools/bam_plbuf.c
                if (bam_plp_push(plp_iter, b) < 0) {
                    // TODO trap errors
                    std::fprintf(stderr, "!");
                }
            }
            bam_itr_destroy(iter);
            bam_destroy1(b);

            // based on samtools/bam_plbuf.c
            bam_plp_push(plp_iter, 0);
            {

            /*
             * Phase 2:
             *
             * Get the whole pileup and look at what's happening at the SNVs site.
             */
                int n_plp, tid, p_pos;
                const bam_pileup1_t *plp;
                while ((plp = bam_plp_next(plp_iter, &tid, &p_pos, &n_plp)) != 0) {
                    if (p_pos == tmp.pos) { // only pay attention to pileups in out target position
                        pileup_func(tid, p_pos, n_plp, plp, &tmp);
                        break;
                    }
                }
            }

            // write the results add
            bam_plp_reset(plp_iter);
            snpsOut << str_buf << '\t' // we will simply append the new column instead of re-writing them

            //  for amplicons :
            //      << chr << '\t' << pos << '\t' << ref << '\t' << var << '\t'
            //      << fr1 << '\t' << p1 << '\t'
            //  for the normal mode :
            //      << chr << '\t' << pos << '\t' << ref << '\t' << var << '\t'
            //      << fr1 << '\t' << fr2 << '\t' << fr3 << '\t'
            //      << p1 << '\t' << p2 << '\t' << p3 << '\t'

                    << tmp.forwM << '\t' << tmp.revM << '\t'
                    << tmp.forwT << '\t' << tmp.revT << '\t' << tmp.pval << '\n';
        }
    }
    snpsOut.close();
    std::free(filename);
    fl.close();
    bam_plp_destroy(plp_iter);
    bam_hdr_destroy(header);
    hts_idx_destroy(idx);
    hts_close(inFile);
}
