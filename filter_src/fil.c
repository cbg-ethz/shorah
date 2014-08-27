/******************************
fil.c
Kerensa McElroy
ETH Zurich
adapted from samtools/calDep.c
 ******************************/
#include <stdio.h>
#include "sam.h"
#include "faidx.h"
#include <map>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <string>
#include <gsl/gsl_sf.h>
using namespace std;

typedef struct {
    int pos, pos1, SNP, forwM, revM, forwT, revT;
    double sig, pval;
    char nuc;
    samfile_t *in;
} tmpstruct_t;

//callback for bam_fetch()
static int fetch_func(const bam1_t *b, void *data)
{
    bam_plbuf_t *buf = (bam_plbuf_t*)data;
    bam_plbuf_push(b, buf); //push all reads covering variant to pileup buffer
    return 0;
}

//callback for bam_plbuf_init()
static int pileup_func(uint32_t tid, uint32_t pos, int n,
                       const bam_pileup1_t *pl, void *data)
{
    tmpstruct_t *tmp = (tmpstruct_t*)data;
    map<char, int> st_freq_all;
    map<char, int> st_freq_Mut;
    st_freq_all['f']=0; //store frequencies of forward and reverse reads
    st_freq_all['r']=0;
    st_freq_Mut['f']=0;
    st_freq_Mut['r']=0;
    if ((int)pos >= tmp->pos && (int)pos < tmp->pos+1)
    {
        for (int i=0; i<n; i++)           //take each read in turn
        {
            char c;
            const bam_pileup1_t *p = pl + i;
            if (p->indel==0 and p->is_del!=1)
                c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)]; //get base for this read
            else if (p->is_del==1)
                c = '-';
            else
                c = '*';
            if (bam1_strand(p->b)==0)    //forward read
                st_freq_all['f']++;
            else
                st_freq_all['r']++;      //reverse read
            if (c==tmp->nuc){            //base is variant
                if (bam1_strand(p->b)==0){
                    st_freq_Mut['f']++;  //forward read
                    }
                else
                    st_freq_Mut['r']++;  //reverse read
                }
         }
        double mean=(double)st_freq_all['f']/(st_freq_all['f']+st_freq_all['r']); //forward read ratio, all reads at this position
        double fr=(double)st_freq_Mut['f']/(st_freq_Mut['f']+st_freq_Mut['r']);   //forward read ratio, only reads with variant at this position
        double alpha;
        double beta;
        double sigma = tmp->sig;
        double prMut;
        int m = st_freq_Mut['f'] + st_freq_Mut['r']; //total reads with variant
        int k;
        double tail = 0.0;
        alpha = mean / sigma;   //alpha and beta for beta binomial
        beta = (1 - mean) / sigma;
		if (alpha < 1E-6 or beta < 1E-6){
	        fprintf(stderr, "All reads in the same orientation?\n");
	        return 11;
		}
        if (fr < mean) {
            for (k = 0; k <= st_freq_Mut['f']; k++){
                tail += gsl_sf_exp(
                        (gsl_sf_lngamma((double)m + 1) -
                         gsl_sf_lngamma((double)k + 1) -
                         gsl_sf_lngamma((double)m - (double)k + 1) +
                         gsl_sf_lngamma(1 / sigma) +
                         gsl_sf_lngamma((double)k + (mean * (1 / sigma))) +
                         gsl_sf_lngamma((double)m + ((1 - mean) / sigma) -
                                        (double)k) -
                         gsl_sf_lngamma(mean * (1 / sigma)) -
                         gsl_sf_lngamma((1 - mean) / sigma) -
                         gsl_sf_lngamma((double)m + (1 / sigma))
                         )
                                   );   //calculate cumulative distribution

                }
            }
        else {
            for (k = st_freq_Mut['f']; k <= m; k++) {
                tail += gsl_sf_exp(
                        (gsl_sf_lngamma((double)m + 1) -
                         gsl_sf_lngamma((double)k + 1) -
                         gsl_sf_lngamma((double)m - (double)k + 1) +
                         gsl_sf_lngamma(1 / sigma) +
                         gsl_sf_lngamma((double)k + (mean *(1 / sigma))) +
                         gsl_sf_lngamma((double)m + ((1 - mean) / sigma) -
                                        (double)k) -
                         gsl_sf_lngamma(mean * (1 / sigma)) -
                         gsl_sf_lngamma((1 - mean) / sigma) -
                         gsl_sf_lngamma((double)m + (1 / sigma))
                         )
                                   ); //the other tail, if required

                }
            }
        tmp->forwM=st_freq_Mut['f'];
        tmp->revM=st_freq_Mut['r'];
        tmp->forwT=st_freq_all['f'];
        tmp->revT=st_freq_all['r'];
        prMut=tail*2; //two sided test
        if (prMut>1)
            prMut=1;
        tmp->pval=prMut; //p value
      }
    return 0;
}

//main
int main(int argc, char *argv[])
{
    tmpstruct_t tmp;
    int c = 0;
    bam_index_t *idx;
	int amplicon = 0;
    while((c=getopt(argc, argv, "b:v:a"))!=EOF){
        switch(c){
			case 'b':
            	tmp.in = samopen(optarg, "rb", 0);
				idx = bam_index_load(optarg);
				break;
			case 'v':
            	tmp.sig = atof(optarg);
				break;
			case 'a':
				amplicon = 1;
        }
    }
    if (tmp.in==0){
        fprintf(stderr, "Failed to open BAM file %s\n", argv[1]);
        return 1;
        }
    if (idx==0){
        fprintf(stderr, "BAM indexing file is not available.\n");
        return 2;
    }
    int maxdepth;
    maxdepth = 100000;
    bam_plbuf_t *buf;
    buf = bam_plbuf_init(pileup_func, &tmp);
    bam_plp_set_maxcnt(buf->iter, maxdepth);
    FILE *fl = fopen("SNV.txt","rt");
    if (fl==NULL){
        fprintf(stderr, "Failed to open SNV file.\n");
        return 3;
        }
    char chr[100], ref, var, fr1[7], fr2[7], fr3[7], p1[7], p2[7], p3[7];
    int pos, name;
    char str_buf[200];
    char reg[200];
    char filename [50];
    ofstream snpsOut;
    sprintf(filename, "SNVs_%f.txt", tmp.sig);
    snpsOut.open(filename, ios::out);

	if (amplicon){
		// amplicon only has one window, parses a single frequency and posterior and uses fr1 and p1 only
		while (fgets(str_buf, 200, fl) !=NULL){
        	if (sscanf(str_buf, "%s %d %c %c %s %s", chr, &pos, &ref, &var, fr1, p1) == 6){
            	tmp.nuc = var;
				sprintf(reg, "%s:%d-%d", chr, pos, pos); 
				bam_parse_region(tmp.in->header, reg, &name, &tmp.pos, &tmp.pos1);
				bam_fetch(tmp.in->x.bam, idx, name, tmp.pos, tmp.pos1, buf, fetch_func);
				bam_plbuf_push(0, buf);
				bam_plbuf_reset(buf);
				snpsOut<<chr<<'\t'<<pos<<'\t'<<ref<<'\t'<<var<<'\t'<<fr1<<'\t'<<p1<<'\t'<<tmp.forwM<<'\t'<<tmp.revM<<'\t'<<tmp.forwT<<'\t'<<tmp.revT<<'\t'<<tmp.pval<<'\n';
			}
    	}
	}
	else{
		while (fgets(str_buf, 200, fl) !=NULL){
        	if (sscanf(str_buf, "%s %d %c %c %s %s %s %s %s %s", chr, &pos, &ref, &var, fr1, fr2, fr3, p1, p2, p3) == 10){
            	tmp.nuc = var;
				sprintf(reg, "%s:%d-%d", chr, pos, pos); 
				bam_parse_region(tmp.in->header, reg, &name, &tmp.pos, &tmp.pos1);
				bam_fetch(tmp.in->x.bam, idx, name, tmp.pos, tmp.pos1, buf, fetch_func);
				bam_plbuf_push(0, buf);
				bam_plbuf_reset(buf); 
				snpsOut<<chr<<'\t'<<pos<<'\t'<<ref<<'\t'<<var<<'\t'<<fr1<<'\t'<<fr2<<'\t'<<fr3<<'\t'<<p1<<'\t'<<p2<<'\t'<<p3<<'\t'<<tmp.forwM<<'\t'<<tmp.revM<<'\t'<<tmp.forwT<<'\t'<<tmp.revT<<'\t'<<tmp.pval<<'\n';
			}
    	}
	}
    snpsOut.close();
    fclose(fl);
    bam_index_destroy(idx);
    bam_plbuf_destroy(buf);
}
