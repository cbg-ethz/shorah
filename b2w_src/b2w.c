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
b2w.c
Kerensa McElroy
UNSW, Sydney Australia
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

using namespace std;

//data for fetch_func and pileup_func 
typedef struct {
	int b, e, win, inc, min_overlap, len, cov, max;
	samfile_t *in;
    uint32_t beg, end; //b, e store begin and end of entire region of interest
    bam_plbuf_t *buf; //beg, end store begin and end of current window
    faidx_t *fai;
    void *rd; //stores working read
    void *rLen; //stores coverage of read starts on reference
    ofstream outFile; //output file for current window
    ofstream reads;
} tmpstruct_t;

//callback one for bam_fetch() (for making windows)
static int fetch_func1(const bam1_t *b, void *data)
{
	tmpstruct_t *tmp = (tmpstruct_t*)data;
    int overlap = 0; //read overlap with window
    int wSize; //window size
    uint32_t Rstart, Rend; //read start and end

    Rstart = b->core.pos;
    Rend = bam_calend(&b->core, bam1_cigar(b));
    wSize = tmp->win;

    int *cv=(int*)tmp->rLen+Rstart; //tracks nr of reads starting at ref pos
    *cv += 1;

    if (*cv<tmp->max) { //limit coverage
        if (Rstart <= tmp->beg && Rend >= tmp->end) { //calculate overlap
        	overlap = tmp->end - tmp->beg + 1;
        }
        else if (Rstart >= tmp->beg && Rend <= tmp->end){
        	overlap = Rend - Rstart +1;
        }
        else if (Rstart <= tmp->beg && Rend <= tmp->end){
        	overlap = Rend - tmp->beg +1;
        }
        else if (Rstart >= tmp->beg && Rend >= tmp->end){
        	overlap = tmp->end - Rstart +1;
        }
        if (overlap >= tmp->min_overlap) {
            char rd[wSize+1]; //store blank read segment of 'N's
            for (int i=0; i<wSize; i++){
            	rd[i]='N';
            }
            rd[wSize]=0;
            tmp->rd=rd;
            bam_plbuf_push(b, tmp->buf); //push a read to pileup buffer     
            tmp->cov+=1; //increment window coverage
            bam_plbuf_push(0, tmp->buf); //but make sure only 1 read at a time
            bam_plbuf_reset(tmp->buf);
            tmp->outFile //write read header 
            << '>' 
            << bam1_qname(b) 
            << ' '
            << Rstart 
            << '\n';
            tmp->outFile //write read 
            << rd 
            << "\n";
        }
    }
    return 0;
}

//callback for bam_plbuf_init() (for making windows)
static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
{
	tmpstruct_t *tmp = (tmpstruct_t*)data;
    int i = pos - tmp->beg; //base position within window
    char *mm=(char*)tmp->rd + i; //location of base
    if (pos >= tmp->beg && pos <= tmp->end) { //make sure read position within window
        if (pl->is_del==1) { //remove deletions, replacing with reference
// char *rb;
// int td;
// td = pl->b->core.tid;
// rb = fai_fetch(tmp->fai, tmp->in->header->target_name[td], &tmp->len);
        	*mm='-';
        }
        else { //replace 'N' in tmp.read with current read base 
        	*mm=bam_nt16_rev_table[bam1_seqi(bam1_seq(pl->b),pl->qpos)];
    }
}
return 0;
}

//callback two (for printing reads)
static int fetch_func2(const bam1_t *b, void *data)
{
	tmpstruct_t *tmp = (tmpstruct_t*)data;
	uint32_t Rstart,Rend;
	Rstart=b->core.pos;
	Rend=bam_calend(&b->core, bam1_cigar(b));
	int *cv=(int*)tmp->rLen+Rstart;
	*cv+=1;
	int readLen=Rend-Rstart;
	if ((*cv<tmp->max) && (readLen>=tmp->min_overlap)) {
		char rd[readLen+1];
		for (int i=0;i<readLen;i++){
			rd[i]='N';
		}
		rd[readLen]=0;
		tmp->rd=rd;
		tmp->beg=Rstart;
		tmp->end=Rend;
		bam_plbuf_push(b, tmp->buf);
		bam_plbuf_push(0,tmp->buf);
		bam_plbuf_reset(tmp->buf);
		tmp->reads << bam1_qname(b) << "\t"
		<< tmp->b << "\t"
		<< tmp->e << "\t"
		<< Rstart+1 << "\t"
		<< Rend << "\t"
		<< rd << "\n";
	}
	return 0;
}

//main
int main(int argc, char *argv[])
{
    int c = 0; //for parsing command line arguments
    int ref = 0; //index for retrieving reference chromosome
    char filename [100]; //filename for window files

    ofstream covOut; //stores window coverages

    tmpstruct_t tmp; //data for callback functions

    char help_string[] = "\nUsage: b2w [options] <in.bam> <in.fasta> region\n\nOptions:\n\t-w: window length (INT)\n\t-i: increment (INT)\n\t-m: minimum overlap (INT)\n\t-x: max reads starting at a position (INT)\n\t-h: show this help\n\n";

    while((c = getopt(argc, argv, "w:i:m:x:h")) != EOF)
    {
    	switch(c)
    	{
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
    		fprintf(stderr, "%s", help_string);
    		return 1;
    	}
    }
    if (argc < 11 || argc > 12) {
    	fprintf(stderr, "%s", help_string);
    	return 1;
    }

    tmp.in = samopen(argv[optind], "rb", 0); //open bam file
    if (tmp.in == 0) {
    	fprintf(stderr, "Failed to open BAM file %s\n", argv[optind]);
    	return 2;
    }

    fai_build(argv[optind + 1]); //generate reference index
    tmp.fai = fai_load(argv[optind + 1]);

    int iBuild = bam_index_build(argv[optind]); //generate bam index
    bam_index_t *idx;
    idx = bam_index_load(argv[optind]); //load bam index
    if (idx == 0) {
    	fprintf(stderr, "BAM indexing file is not available.\n");
    	return 3;
    }

    int32_t n;
    uint32_t *ln;
    n = tmp.in->header->n_targets; //number of chromosomes
    ln = tmp.in->header->target_len; //chromosome lengths
    tmp.buf = bam_plbuf_init(pileup_func, &tmp); //initiate pileup buffer
    covOut.open("coverage.txt", ios::out); //open file to store window coverage
    tmp.reads.open("reads.fas", ios::out);
    if (argc == optind + 2) { //region not specified
        for (int i=0; i < n; i++) { //take each chromosome in turn
            int lnth = (int)ln[i]; //chromosome length
            tmp.b=0;
            tmp.len=lnth;
            tmp.e=lnth;
            int rLen[lnth]; //list for read start coverage
            for (int j=0;j<lnth;j++){
            	rLen[j]=0;
            }
            tmp.rLen=rLen;

            bam_fetch(tmp.in->x.bam, idx, i, tmp.b, tmp.e, &tmp, fetch_func2); //fetch and write reads
            tmp.beg = 0;
            tmp.end = tmp.win - 1;
            while (tmp.end <= ln[i] - 1) { //make windows over length of chromosome
            	tmp.cov = 0;
                for (int j=0; j < lnth; j++){ //restart read start coverage count for each window
                	rLen[j]=0;
                }
                tmp.rLen = rLen;
                sprintf(filename, "w-%s-%d-%d.reads.fas", //read window filename
                	tmp.in->header->target_name[i],
                	tmp.beg+1, tmp.end+1);
                tmp.outFile.open(filename, ios::out);    
                bam_fetch(tmp.in->x.bam, idx, i, tmp.beg, tmp.end, &tmp, fetch_func1); //fetch and write reads
                if (tmp.cov != 0){ //ignore empty windows
                    covOut << filename << "\t" //write out coverage information
                    << tmp.in->header->target_name[i] << "\t"
                    << tmp.beg+1 << "\t" 
                    << tmp.end+1 << "\t" 
                    << tmp.cov << "\n";
                }
                tmp.outFile.close();
                tmp.beg+=tmp.inc; //increment windows
                tmp.end+=tmp.inc;
                if (tmp.inc == 0){break;} // OZ: in amplicon mode, one window only
            }
        }
    }
    else { //parse specific region
    	bam_parse_region(tmp.in->header, argv[optind + 2], &ref, &tmp.b, &tmp.e);
    	if (ref < 0) {
    		fprintf(stderr, "Invalid region %s\n", argv[optind +2]);
    		return 4;
    	}
        tmp.b -= 3 * tmp.inc; //make sure start and end of region are covered by 3 windows
        tmp.e += 3 * tmp.inc;
        if (tmp.b < 0){
        	tmp.b = 0;
        }
        if ((uint32_t)tmp.e >= ln[ref]){
        	tmp.e = ln[ref] - 1;
        }
        int rLen[ln[ref]];
        tmp.len = ln[ref];

        for (uint32_t j=0; j<ln[ref]; j++){
        	rLen[j] = 0;
        } 
        tmp.rLen=rLen;
        bam_fetch(tmp.in->x.bam, idx, ref, tmp.b, tmp.e, &tmp, fetch_func2); //fetch and write reads
        tmp.beg = tmp.b;
        tmp.end = tmp.beg + tmp.win - 1;

        while (tmp.end < (uint32_t)tmp.e) { //make windows over region
        	tmp.cov = 0;
           for (uint32_t j=0; j<ln[ref]; j++){ //restart read start coverage count for each window
           	rLen[j] = 0;
           }
           tmp.rLen = rLen;
           sprintf(filename, "w-%s-%d-%d.reads.fas", //read window filename
           	tmp.in->header->target_name[ref],
           	tmp.beg + 1,
           	tmp.end + 1); 
           tmp.outFile.open(filename, ios::out);
           bam_fetch(tmp.in->x.bam, idx, ref, tmp.beg, tmp.end, &tmp, fetch_func1);  //fetch and write reads
           if (tmp.cov != 0) { //ignore empty windows
           	covOut << filename << "\t" 
           	<< tmp.in->header->target_name[ref] << "\t" 
           	<< tmp.beg+1<< "\t" 
           	<< tmp.end+1 << "\t"  
           	<< tmp.cov << "\n";
           }
           tmp.outFile.close();
           tmp.beg += tmp.inc; //increment windows
           tmp.end += tmp.inc;
		   if (tmp.inc == 0){break;} // OZ: in amplicon mode, one window only
		}
	}
	bam_plbuf_destroy(tmp.buf);
	bam_index_destroy(idx);
	samclose(tmp.in);
	covOut.close();
	tmp.reads.close();
	return 0;
}
