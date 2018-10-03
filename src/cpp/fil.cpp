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

#include <gsl/gsl_sf.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>

#define UNUSED(expr) (void)(expr)

typedef struct {
    int maxdepth;
    double sig;
    bool skip_indel;
} paramstruct_t;



// main
int main(int argc, char* argv[])
{
    paramstruct_t params = { 10000, 0.01, // cli.py's defaults
            false }; // mimics b2w's behaviour (historically fil used to do the other way around, leading to contradictions in results).
    int c = 0;
    hts_idx_t* idx = NULL;
    int amplicon = 0;
    htsFile* inFile = NULL;
    while ((c = getopt(argc, argv, "b:v:x:ad")) != EOF) {
        switch (c) {
            case 'b':
                inFile = hts_open(optarg, "r");
                idx = hts_idx_load(optarg, HTS_FMT_BAI);
                 // NOTE BAI sufficient for up to 2^29-1. If we move from virus to organism with longer chromosomes (plants ?), we should switch to HTS_FMT_CSI
                break;
            case 'v':
                params.sig = atof(optarg);
                break;
            case 'x':
                params.maxdepth = atoi(optarg);
                break;
            case 'a':
                amplicon = 1;
                break;
            case 'd':
                params.skip_indel = true; // this will bring the old historical behaviour of fil
                break;
        }
    }
    if (inFile == NULL) {
        std::cerr << "Failed to open BAM file  " << argv[1] << std::endl;
        return 1;
    }
    bam_hdr_t* header = sam_hdr_read(inFile);
    if (idx == NULL) {
        std::cerr << "BAM indexing file is not available." << std::endl;
        return 2;
    }
    bam_plp_t plp_iter = bam_plp_init(0, 0);
    bam_plp_set_maxcnt(plp_iter, params.maxdepth);
    std::ifstream fl ("SNV.txt", std::ios_base::in);
    if (fl.fail()) {
        std::cerr << "Failed to open SNV file." << std::endl;
        return 3;
    }
    std::string str_buf;

    std::ostringstream oss_fn;
    oss_fn << "SNVs_" << std::fixed << std::setprecision(6) << params.sig << ".txt";
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
        int targetid, pos; // targetid : target (as header's numeric form)  pos: 1-based position of SNV
        char ref, var;     // reference and variant nucleotide for this SNV
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
            enum class Strand { f, r };
            std::map<Strand, int> st_freq_all;
            std::map<Strand, int> st_freq_Mut;
            // init sane defaults
            st_freq_all[Strand::f] = 0;  // store frequencies of forward and reverse reads
            st_freq_all[Strand::r] = 0;
            st_freq_Mut[Strand::f] = 0;
            st_freq_Mut[Strand::r] = 0;


            // as seen in htslib/hts.c:hts_parse_reg
            const int pos_begin = pos - 1,  // begin is 0-based (so position-1)
                      pos_end = pos;        // end is 1-based (so straight position)
            // based on samtools/bam_aux.c:bam_parse_region
            if ((targetid = bam_name2id(header, chr.c_str())) < 0) {
                std::cerr << "Invalid region " << chr << std::endl;
                return 4;
            }

            /*
             * Phase 1:
             *
             * Get all the reads that cover the interesting region (on chromosomes/target "targetid", position "pos"
             * And store them on a pileup
             */
            // based on samtools/bam.c:bam_fetch
            // TODO These BAM iterator functions work only on BAM files.  To work with either BAM or CRAM files use the sam_index_load() & sam_itr_*() functions.
            bam1_t *b = bam_init1();
            hts_itr_t* iter = bam_itr_queryi(idx, targetid, pos_begin, pos_end);
            while (hts_itr_next(inFile->fp.bgzf, iter, b, 0) >= 0) {
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
            {

            /*
             * Phase 2:
             *
             * Get the whole pileup and look at what's happening at the SNVs site.
             */
                int n_plp, tid, p_pos;
                const bam_pileup1_t *plp;
                while ((plp = bam_plp_next(plp_iter, &tid, &p_pos, &n_plp)) != 0) {
                    if (p_pos == pos_begin) { // only pay attention to pileups in out target position
                        for (int i = 0; i < n_plp; i++)  // take each read in turn
                        {
                            char c;
                            const bam_pileup1_t* p = plp + i;
                            // for the exact signification of "is_del" and "indel", see https://www.biostars.org/p/104301/#104311
                            const bool skip = (params.skip_indel && (p->indel != 0));
                            if ((! skip) and p->is_del != 1)
                                // based on samtools/bam.h
                                c = seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)];
                            else if (p->is_del == 1)
                                c = '-';
                            else
                                c = '*';
                                // NOTE historically, fil used to ignore SNVs nearby insertions/deletions, leading to disagreeing with B2W and causing bugs.
                                // Now there's a switch for that. Either both always keep SNVs (default behaviour and old b2w behaviour) or both ignore them when nearby insertions (old fil's behaviour)
                            if (bam_is_rev(p->b) == 0)  // forward read
                                st_freq_all[Strand::f]++;
                            else
                                st_freq_all[Strand::r]++;  // reverse read
                            if (c == var) {     // base is variant nucleotide ?
                                if (bam_is_rev(p->b) == 0) {
                                    st_freq_Mut[Strand::f]++;  // forward read
                                } else
                                    st_freq_Mut[Strand::r]++;  // reverse read
                            }
                        }
                        //break;
                    }
                }
            }

            /*
             * Phase 3:
             *
             * Do stats on the counts - is there a strand bias ?
             */
            
            // init sane defaults
            double pval = 1.;

            if (st_freq_all[Strand::f] == 0 and st_freq_all[Strand::r] == 0) {
                // historically this could happen due to B2W and FIL counting indels in different way.
                std::cerr << "Bug: No reads found at position " << pos <<" ?!?" << std::endl;
            } else if (st_freq_Mut[Strand::f] == 0 and st_freq_Mut[Strand::r] == 0) {
                std::cerr << "Critical: No mutations found at position " << pos <<" ?!?" << std::endl;
            } else {
                double mean = (double)st_freq_all[Strand::f] /
                            (st_freq_all[Strand::f] +
                            st_freq_all[Strand::r]);  // forward read ratio, all reads at this position
                double fr =
                    (double)st_freq_Mut[Strand::f] /
                    (st_freq_Mut[Strand::f] +
                    st_freq_Mut[Strand::r]);  // forward read ratio, only reads with variant at this position
                double alpha;
                double beta;
                double sigma = params.sig;
                double prMut;
                int m = st_freq_Mut[Strand::f] + st_freq_Mut[Strand::r];  // total reads with variant
                int k;
                double tail = 0.0;
                alpha = mean / sigma;  // alpha and beta for beta binomial
                beta = (1 - mean) / sigma;
                
                if (alpha < 1E-6 or beta < 1E-6) {
                    // ...at that point it's not even bias, it's a monopoly ! ;-)
                    std::cerr << "Warning: All reads at position " << pos <<" in the same " << ((alpha < 1E-6) ? "reverse" : "forward") << " orientation ?" << std::endl;
                } else {
                    if (fr < mean) {
                        for (k = 0; k <= st_freq_Mut[Strand::f]; k++) {
                            tail += gsl_sf_exp((
                                gsl_sf_lngamma((double)m + 1) - gsl_sf_lngamma((double)k + 1) -
                                gsl_sf_lngamma((double)m - (double)k + 1) + gsl_sf_lngamma(1 / sigma) +
                                gsl_sf_lngamma((double)k + (mean * (1 / sigma))) +
                                gsl_sf_lngamma((double)m + ((1 - mean) / sigma) - (double)k) -
                                gsl_sf_lngamma(mean * (1 / sigma)) - gsl_sf_lngamma((1 - mean) / sigma) -
                                gsl_sf_lngamma((double)m + (1 / sigma))));  // calculate cumulative distribution
                        }
                    } else {
                        for (k = st_freq_Mut[Strand::f]; k <= m; k++) {
                            tail += gsl_sf_exp(
                                (gsl_sf_lngamma((double)m + 1) - gsl_sf_lngamma((double)k + 1) -
                                gsl_sf_lngamma((double)m - (double)k + 1) + gsl_sf_lngamma(1 / sigma) +
                                gsl_sf_lngamma((double)k + (mean * (1 / sigma))) +
                                gsl_sf_lngamma((double)m + ((1 - mean) / sigma) - (double)k) -
                                gsl_sf_lngamma(mean * (1 / sigma)) - gsl_sf_lngamma((1 - mean) / sigma) -
                                gsl_sf_lngamma((double)m + (1 / sigma))));  // the other tail, if required
                        }
                    }
                    prMut = tail * 2;  // two sided test
                    if (prMut > 1) prMut = 1;
                    pval = prMut;  // p value
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

                    << st_freq_Mut[Strand::f] << '\t' << st_freq_Mut[Strand::r] << '\t'
                    << st_freq_all[Strand::f] << '\t' << st_freq_all[Strand::r] << '\t' << pval << '\n';
        }
    }
    snpsOut.close();
    fl.close();
    bam_plp_destroy(plp_iter);
    bam_hdr_destroy(header);
    hts_idx_destroy(idx);
    hts_close(inFile);
}
