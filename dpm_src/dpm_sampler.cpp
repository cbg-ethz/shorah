/*
# Copyright 2007, 2008, 2009, 2010
# Niko Beerenwinkel,
# Arnab Bhattacharya,
# Nicholas Eriksson,
# Moritz Gerstung,
# Lukas Geyrhofer,
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <time.h>
#include <map>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <new>
#include <memory>
#include <sstream>
#include <algorithm>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
using namespace std;
#include "data_structures.h"
#include "dpm_sampler.h"

#define PROPHISTSIZE 100
int main(int argc, char** argv){
  
  unsigned int i, j, k, ll, K1=0, iter2, mesh, tot_untouch, new_proposed=0;
  int dk1, hapbases;
  int* p;
  cnode* tn;
  //rnode* rn;
  rnode* tr;
  ssret* samp_stat;
  double dt;
  time_t t1, t2;
  FILE* alphafile;
  FILE* iterfile;
  float alpha2;
  
  double_threshold_min = gsl_sf_log(DBL_MIN);
  double_threshold_max = gsl_sf_log(DBL_MAX);
  parsecommandline(argc, argv);
  
  string instr = filein;

  /// rename sampling file and debug file  
  int insize = instr.find('.');
  string statstr = instr.substr(0,insize).append(".dbg");
  string outstr = instr.substr(0,insize).append(".smp");
  string alphastr = instr.substr(0,insize).append(".alpha");
  string iterstr = instr.substr(0,insize).append(".iter");
  
  ofstream out_file(outstr.c_str());
  ofstream stat_file(statstr.c_str());
  
  char* astrchar = new char[strlen(alphastr.c_str())];
  strcpy(astrchar, alphastr.c_str());

  char* iterstrchar = new char[strlen(iterstr.c_str())];
  strcpy(iterstrchar, iterstr.c_str());
  
  stat_file << "# dna_code:\t";
  stat_file << i2dna_code;
  
  //  randseed = 1257501510;
  stat_file << "# randseed = " << randseed << "\n";
  
  // random number generator via gsl
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rg = gsl_rng_alloc(T);
  gsl_rng_set(rg, randseed);
  
  dist = (int*)calloc(2, sizeof(int));
  if (dist == NULL) exit(EXIT_FAILURE);
  res_dist = (int*)calloc(2, sizeof(int));
  if (res_dist == NULL) exit(EXIT_FAILURE);
  res = (ssret*)malloc(sizeof(ssret));
  if (res == NULL) exit(EXIT_FAILURE);
  
  cbase = (int*) malloc(B * sizeof(int)); // count base
  if (cbase == NULL) exit(EXIT_FAILURE);
  pbase = (double*) malloc(B * sizeof(double));
  if (pbase == NULL) exit(EXIT_FAILURE);
  log_pbase = (double*) malloc(B * sizeof(double)); 
  if (log_pbase == NULL) exit(EXIT_FAILURE);
  
  mesh = 1;//iter/100;
  
  read_data(filein, stat_file);
  ftable = (int**) calloc(J, sizeof(int*));
  ftable_sum = new int[J];

  //int tot_pos;

  for(j=0; j<J; j++){
    ftable[j] = new int[B];
    for(i=0; i<B; i++)
      ftable[j][i] = 0;
  }
 
  for (j=0; j<J; j++){
    ftable_sum[j] = 0;
    //ftable[j][i] = 0.0;
    for (k =0; k<n; k++){
      if(r[k][j] < B){
	ftable[j][r[k][j]] ++;
	ftable_sum[j]++;
      }
    }
    //for (i=0;i<B; i++){
    //  ftable[j][i] = ftable[j][i]/tot_pos;
    //}  
  }
  stat_file << "# Number of reads, n = "<<n<<endl<<"# Read length, J = "<<J<<endl<<"# J/10 + 1 = "<<J/10 + 1<<endl<<endl;

   // creads = (int**)calloc(n, sizeof(int*));

   readtable = (crnode**)calloc(n, sizeof(crnode*));

   for(i=0; i<n; i++){
     readtable[i] = (crnode*)calloc(1, sizeof(crnode));
     readtable[i]->creads = new int[J/10 + 1];
     readtable[i]->mindices = new int[LIMIT + 1];
   }

   for(i=0; i<n; i++){
     read_conversion(readtable[i], r[i], J);
     readtable[i]->mindex = 0;
     readtable[i]->weight = 1; // mapped = no
   }

   
   //multimap<string,int>::iterator it;
   //string stemp = ("");
   //int mapcount, mflag;
   /*
   for(j=0; j<J/10 + 1; j++)
     stemp += static_cast<ostringstream*>( &(ostringstream() << readtable[0]->creads[j]) )->str();

   readmap.insert ( pair<string,int>(stemp, 1) );
   readtable[0]->mindex = 0;

   
   for(i=1; i<n; i++){
     stemp = ("");
     for(j=0; j<J/10 + 1; j++)
       stemp += static_cast<ostringstream*>( &(ostringstream() << readtable[i]->creads[j]) )->str();
     mapcount = 0;
     mflag = 0;
     for(it = readmap.begin(); it != readmap.end(); it++, mapcount++){
       if(it->first.compare(stemp) == 0 && it->second < LIMIT){
	 it->second ++;
	 readtable[i]->mindex = mapcount;
	 mflag = 1;
	 break;
       }
     }
     if(mflag == 0){
       readmap.insert ( pair<string,int>(stemp, 1) );
       readtable[i]->mindex = mapcount;
       }
   }
   */

   int* hd;

   for(i=0; i<n; i++){
     if(readtable[i]->weight > 0){         // only if not mapped
       readtable[i]->mindex = i;          // then it is unique
       for(j=i+1; j<n; j++){
	 if(readtable[i]->weight < LIMIT){
	   hd = seq_distance_rr(readtable[i]->creads, readtable[j], J);
	   if(hd[0] == 0){
	     readtable[j]->mindex = i;
	     readtable[j]->weight = 0;     // not unique, mapped = yes
	     readtable[i]->weight++;
	   } 
	 }
       }
     }
   }
       
   /*  

   for(i=0; i<n; i++){
     hd = seq_distance_rr(readtable[i]->creads, readtable[readtable[i]->mindex], J); 
     cout<<" Read "<<i<<" is mapped to "<<readtable[i]->mindex<<" with distance = "<<hd[0]<<" and has weight = "<<readtable[i]->weight<<endl;
   }
   */
   //q = readmap.size();
   q = 0;

   for(i=0; i<n; i++){
     if(readtable[i]->weight > 0)
       q++;
   }

   stat_file << "# q = "<<q<<endl;

   readtable2 = (crnode**)calloc(q, sizeof(crnode*));

   for(i=0; i<q; i++){
     readtable2[i] = (crnode*)malloc(sizeof(crnode));
     //  readtable2[i]->creads = new int[J/10 + 1];
     readtable2[i]->weight = 0;
     readtable2[i]->mindex = 0;
     readtable2[i]->mindices = new int[LIMIT + 1];
   }
   
   int wcount;

   for(i=0; i<q; i++){
     for(wcount=0; wcount<LIMIT + 1; wcount++)
       readtable2[i]->mindices[wcount] = 0;
   }

   j = 0;
   unsigned int hope;
   int* itrack = new int[q];

   for(i=0; i<q; i++)
     itrack[i] = 0;

   for(i=0; i<n; i++){
     if(readtable[i]->weight > 0){
       readtable2[j]->creads = readtable[i]->creads;
       readtable2[j]->missing = readtable[i]->missing;
       readtable2[j]->mindex = i;
       readtable2[j]->weight = readtable[i]->weight;
       readtable2[j]->mindices[0] = i;
       itrack[j]++;
       j++;
     }
     else{
       for(hope=0; hope<j; hope++){
	 if(readtable2[hope]->mindex == readtable[i]->mindex){
	   readtable2[hope]->mindices[itrack[hope]] = i;
	   itrack[hope]++;
	 }
       }
     }
   }


  delete[] itrack;
  build_assignment(stat_file);

  stat_file << "# Assignment built"<<endl;

  stat_file << "# +++++++++++++++++++ BEFORE THE SAMPLING +++++++++++++++++++\n";
  //  print_stats(out_file, mxt, J);
  out_file << setw(6) << "#iter";
  out_file << setw(8) << "class";
  out_file << setw(8) << "untouch";
  out_file << setw(9) << "theta";
  out_file << setw(9) << "gamma\n";
  
  
  /*********************
    sampling procedure
  *********************/
  (void) time(&t1);
  
  int* temp;
  temp = new int[J/10 + 1];
  
  for(k=0; k<=iter; k++){
#ifdef DEBUG
    printf("-----------------------------------------------\n");
    printf("-----------> sampling the %ith time <-----------\n", k);
    printf("-----------------------------------------------\n");
#endif
    
    /// Modify alpha reading from an external file
    alphafile = fopen(astrchar, "r");
    if (alphafile != NULL){
      fscanf(alphafile, "%f", &alpha2);
      alpha = alpha2;
      fclose(alphafile);
    }

    /// Modify alpha reading from an external file
    iterfile = fopen(iterstrchar, "r");
    if (iterfile != NULL){
      fscanf(iterfile, "%ui", &iter2);
      if (iter != iter2 and iter2 > k+100){
	if (iter2 < iter){
	  HISTORY = iter2 - k - 1;
	  cerr << "HISTORY just changed to " << HISTORY << endl;
	}
	iter = iter2;
      }
      fclose(iterfile);
    }

    tot_untouch = 0;
    // sample classes for all reads
    for(i=0; i<q; i++){
      samp_stat = sample_class(i, k);
      tot_untouch  += samp_stat->untouched;
      new_proposed += samp_stat->proposed;
    }
    
    // sample haplotypes from reads
    tn = mxt;
    K1 = 0;
    dt = 0.0;
    dk1 = 0;
    hapbases = 0;
    while(tn != NULL){
      sample_hap(tn);
      conversion(temp, tn->h, J);
      tr = tn->rlist;
      while (tr != NULL){
	p = seq_distance_new(temp, readtable2[tr->ri], J);
	//p = seq_distance(tn->h, r[tr->ri], J);
	dt += p[1] * readtable2[tr->ri]->weight;
	//dt += p[1];
	//dt += p[0];
	tr = tr->next;
      }
      
      // compute distance of all reads to their haplotype
      for (ll=0; ll<q; ll++){
	p = seq_distance_new(temp, readtable2[ll], J);
	//p = seq_distance(tn->h, r[ll], J);
	tn->rd0[ll] = p[0];
	tn->rd1[ll] = p[1];
      }
      
      // compute distance of haplotypes to reference
      //p = seq_distance_new(conversion(tn->h, J), conversion(h, J), J);
      p = seq_distance(tn->h, h, J);
      dk1 += p[1];
      hapbases += p[1] + p[0];
      K1++;
      
      if(k == iter - HISTORY + 1){// starts recording
	if(record == 0){
	  //  create_history(k);
	  record = 1;
	}
      }
      
      if(record)
	record_conf(tn, HISTORY + k - iter - 1);
      

      tn = tn->next;
    }
    
    double b_alpha = 0.0, b_beta = 0.0;
    
    b_alpha =  dt + (eps1 * eps2 * totbases);
    b_beta = (totbases - dt) + eps2 * totbases * (1 - eps1);

    theta = gsl_ran_beta(rg, b_alpha, b_beta); 
    //theta = dt/(q * J) + gsl_ran_gaussian(rg, g_noise);
    //theta = dt/totbases + gsl_ran_gaussian(rg, g_noise);
    //theta = (totbases - dt)/totbases + gsl_ran_gaussian(rg, g_noise);
    gam = (double)dk1/hapbases;

    // HACK!!! theta=1 gives undesired behaviour; not elegant, but effective
    if(theta >= 1.0)
      theta = 0.9999;
    if(theta <= 0.0)
      theta = 0.0001;
    // sample_ref();    // sampling the reference every step gives strange behaviour...
    
#ifdef DEBUG
    cerr << "dt=" << dt <<  "\ttheta=" << theta<< "\tgamma=" << gamma << endl;
#endif
    //    out_file("iteration\t%i\t%i\t%i\t%f\t%f\n", k+1, count_classes(mxt), tot_untouch, theta, gam); 
    out_file << setw(6) << k+1;// << "\t";
    out_file << setw(8) << count_classes(mxt);// << "\t";
    out_file << setw(8) << tot_untouch;// << "\t";
    out_file << setw(10) << setprecision(6) << theta  << setw(10) << gam << endl;
  }
  
  if( remove( iterstrchar ) != 0 )
    stat_file << "# iter file was not created" << endl;
  else
    stat_file << "# iter file successfully deleted" << endl;
  
  if( remove( astrchar ) != 0 )
    stat_file << "# alpha file was not created" << endl;
  else
    stat_file << "# alpha file successfully deleted" << endl;
  
  delete []astrchar;
  delete []iterstrchar;
  
  (void) time(&t2);
  
  /*********************
      sampling ends
  *********************/
  
  stat_file << "# Survived sampling"<<endl;
  stat_file << "# sampling took " << difftime(t2, t1) << " seconds\n";
  
  //  write_assignment(k, new_proposed, mxt);
  
  stat_file << "\n\n# +++++++++++++++++++ AFTER THE SAMPLING ++++++++++++++++++++\n";
  print_stats(stat_file, mxt, J);
  stat_file << "\n# reference genome is:\n# ";
  for(j=0; j<J; j++)
    stat_file << i2dna_code[h[j]];
  stat_file << "\n\n\n";
  stat_file << "\n#gamma = " << gam << "\n";
  stat_file << "\n#theta = " << theta << "\n";
  stat_file << "\n#final number of components = " << K1 << "\n";
  stat_file << "\n#made "<< new_proposed << " new clusters\n";
  stat_file << "\n#number of haplotypes in history = " << haplotypecount << "\n";
  
  /******************
    write corrected
  ******************/
  
  /*
  fprintf(stderr, "# writing corrected\n");
  FILE* corr;
  if( (corr = fopen("corrected.tmp", "w")) == NULL ){
    printf("Impossible to write file corrected.tmp\n");
    exit(EXIT_FAILURE);
  }
  
  for(i=0;i<n;i++) {
    quality = setfinalhaplotype(i);
    j=0;
    while(id[i][j] != '\n') {
      fputc(id[i][j], corr);
      j++;
    }
    fprintf(corr, "|%f\n", quality);
    for(j=0; j<J; j++)fprintf(corr,"%c",i2dna_code[r[i][j]]);
    fprintf(corr,"\n");
  }
  
  if(ho_flag) {
    write_haplotype_frequencies(haplotype_output, ho_count);
  }
  */  
  
  write_posterior_files(instr);
  
  cleanup();
  (void) time(&t1);
  stat_file << "# after sampling took " << difftime(t1, t2) << " seconds";
  /*
  tn = mxt;
  cnode* tn2;
  int difference = 0, similarity = 0, alt_diff = 0, alt_sim = 0, it;
  
  while(tn != NULL){
    tr = tn->rlist;
    cout<<tn->ci<<endl;
    while(tr != NULL){
      for(it=0; it < readtable2[tr->ri]->weight; it++)
	cout<<readtable2[tr->ri]->mindices[it]<<"   ";
      cout<<"   "<<tr->ri<<"  weight = "<<readtable2[tr->ri]->weight<<"   "<<tn->rd0[tr->ri]<<"  ";
      // if(tn->rd0[tr->ri] >= 0){
      tn2 = mxt;
      while(tn2 != NULL){
	cout<<"distance to "<<tn2->ci<<" = "<<tn2->rd0[tr->ri]<<"   ";
	tn2 = tn2->next;
      }
	// }
      cout<<endl;
      difference += tn->rd0[tr->ri];
      similarity += tn->rd1[tr->ri];
      alt_diff += tn->rd0[tr->ri] * readtable2[tr->ri]->weight;
      alt_sim += tn->rd1[tr->ri] * readtable2[tr->ri]->weight;
      tr = tr->next;
    }
    cout<<" total weight = "<<weight(tn)<<endl;
    tn = tn->next;
  }

  cout<<"difference = "<<difference<<"\tsimilarity = "<<similarity<<endl;
  cout<<"alt_diff = "<<alt_diff<<"\talt_sim = "<<alt_sim<<endl;
  cout<<"\ttotbases = "<<totbases<<"\tq * J = "<<q * J<<endl;
  */
  delete[] temp;
  return 0;
}


void read_data(char* filename, std::ofstream& out_file)
{
  FILE* data;
  char c;
  unsigned int seq=0, i, j=0, k=0;
  out_file << "# reading data\n";
  if( (data = fopen(filename, "r")) == NULL){
    out_file << "Impossible to open file " << filename << "\n";
    exit(EXIT_FAILURE);
  }

  for (;;) // first sweep to count
    {
      c = fgetc(data);
      
      if (c == '>'){
	n++;
	seq = 0;
      }
      if (c == '\n')
	seq = 1;

      if (isdna(c) && seq) {
	totsites++;
	if (d2i(c) != B)
	  totbases++;
      }
      
      if (c == EOF)
	break;
    }
  
  if (n == 0) {
    cout << "Error: No data found!\n";
    exit(EXIT_FAILURE);
  }

  J = totsites/n; // assuming fixed length reads
  
  fclose(data);

  if(n == 1) {
    out_file << "Nothing to do, have only one read... (n == 1)\n";
    exit(0);
  }
  
  out_file << "# totbases = "<< totbases << "\ttotsites = " << totsites << "\n";
  
  ///allocate memory for r[i][j], id[i], read2hap[i]
  r = (short unsigned int**)calloc(n, sizeof(unsigned short int*));
  if (r == NULL) exit(EXIT_FAILURE);
  id = (char**)calloc(n, sizeof(char*));
  if (id == NULL) exit(EXIT_FAILURE);
  
  /// This defines an array of maps
  /// to store the haplotypes every read has been assigned to
  pair <sup_map*, ptrdiff_t> read2hap_pair = get_temporary_buffer<sup_map>(n);
  new (read2hap_pair.first) sup_map[n];
  read2hap = read2hap_pair.first;
  
  for (i = 0; i < n; i++){
    r[i]  = (unsigned short int*)calloc(J, sizeof(unsigned short int));
    id[i] = (char*)calloc(100, sizeof(char));
  }
  
  if( (data = fopen(filename, "r")) == NULL){
    printf("Impossible to REopen file %s\n", filename);
    exit(EXIT_FAILURE);
  }
  
  i = -1;
  for (;;)  {// second sweep to read
      c = fgetc(data);
      
      if (c == '>'){
	i++;
	seq = 0;
	k = 0;
      }
      if (seq == 0 and c != '\n'){
	id[i][k] = c;
	k++;
      }
      if (seq == 0 && c == '\n'){
	seq = 1;
	j=0;
      }
      
      if (isdna(c) && seq){
	r[i][j] = d2i(c);
	j++;
      }
      if (c == EOF)
	break;
  }
  
  fclose(data);

}

void read_conversion(crnode* b, unsigned short int* a, int seq_length){
  int i, j, temp;
//   crnode* b = (crnode*)malloc(sizeof(crnode));

//   b->creads = new int[seq_length/10 + 1];
//   b->missing = 0;
//   b->mindices = new int[LIMIT + 1];
 
  for(i=0; i<seq_length/10 + 1; i++)
    b->creads[i] = 0;

  for(i=0; i<seq_length/10; i++){
    for(j=10; j>0; j--){
      temp = (int) pow(8.0, (double)j-1);
      b->creads[i] += temp * a[10 * i + 10 - j];
      b->missing += (a[10 * i + 10 - j] == B);
    }
  }

  for(j=10; j>(10 - (seq_length%10)); j--){
    temp = (int) pow(8.0, (double)j-1);
    b->creads[i] += temp * a[10 * i + 10 -j];
    b->missing += (a[10 * i + 10 - j] == B);
  }
  
  return;
}


void conversion(int* b, unsigned short int* a, int seq_length){
  int i, j, temp;
  //b = new int[seq_length/10 + 1];

  for(i=0; i<seq_length/10 + 1; i++) b[i] = 0;

  for(i=0; i<seq_length/10; i++){
    for(j=10; j>0; j--){
      temp = (int) pow(8.0, (double)j-1);
      b[i] += temp * a[10 * i + 10 - j];
    }
  }

  for(j=10; j>(10 - (seq_length%10)); j--){
    temp = (int) pow(8.0, (double)j-1);
    b[i] += temp * a[10 * i + 10 -j];
  }

  return;// b;
}

int weight(cnode* wn){
  int w = 0;
  rnode* rn = NULL;
  
  rn = wn->rlist;
      
  while(rn != NULL){
    w += readtable2[rn->ri]->weight;
    rn = rn->next;
  }

  return w;
}

int weight_shift(cnode* wn, unsigned int i){
  int w = 0;
  rnode* rn = NULL;
  
  rn = wn->rlist;
      
  while(rn != NULL){
    if(rn->ri != i)
      w += readtable2[rn->ri]->weight;
    rn = rn->next;
  }

  return w;
}

int weight_final(const cnode* wn){
  int w = 0;
  rnode* rn = NULL;
  
  rn = wn->rlist;
      
  while(rn != NULL){
    w += readtable2[rn->ri]->weight;
    rn = rn->next;
  }

  return w;
}

void build_assignment(std::ofstream& out_file){
  
  unsigned int i, m, ll;
  unsigned int ci;
  cnode* cn;
  rnode* rn;
  int* p;
  double* p_k;
  //double* p_q;
  gsl_ran_discrete_t* g;
  
  out_file << "# building assignment\n";
  h = (short unsigned int*)calloc(J, sizeof(unsigned int));
  c_ptr = (cnode**)calloc(n, sizeof(cnode*));

  MAX_K = n + 1;
  P = (double*) calloc((MAX_K), sizeof(double));
  log_P = (double*) calloc((MAX_K), sizeof(double));
  cl_ptr = (cnode**) calloc((MAX_K), sizeof(cnode*));
  
  // if avgNK is set (>0), calculate K
  if(avgNK > 0.0) {
    K = (int)(1.*n/avgNK);
  }
  // otherwise K is already set as commandline option
  // K == 0 can happen with very few reads...

  // no empty clusters at the beginning 
  if ((n < K) || (K == 0))
    K = n;
  
  lowest_free = K;
  p_k = (double*)malloc(K*sizeof(double));
  //p_q = (double*)malloc(K*sizeof(double));  

  // make the class list of K (or q?) initial components 

  out_file << "# Initial number of clusters, K = "<<K<<endl;

  mxt = NULL;
  for(i=0; i<K; i++) {
    add_comp(&mxt, &hst, i, 0, NULL, NULL, NULL, NULL, J, 0, record);
    p_k[i] = 1.;
    //p_q[i] = 1;
  }
  cn = mxt;
  g = gsl_ran_discrete_preproc(K, p_k);
  
  for(m=0; m<q; m++){
    ci = gsl_ran_discrete(rg, g);
    cn = mxt;
    while(cn != NULL) {
      if(ci == cn->ci){
	add_read(&cn->rlist, m);
	cn->size++;
	break;
      }
      cn=cn->next;
    }
  }
  gsl_ran_discrete_free(g);
  free(p_k);
  //free(p_q);

  // there is at least the first cluster...
  cn = mxt;
  
#ifdef DEBUG
  printf("cluster %d has size %d\n",cn->ci,cn->size);
#endif
  
  cn->h = (unsigned short int *) malloc(J * sizeof(unsigned short int));
  cn->rd0 = (unsigned short int*) malloc(q * sizeof(unsigned short int));
  cn->rd1 = (unsigned short int*) malloc(q * sizeof(unsigned short int));
  
  while(cn->size == 0) {
    // or maybe not ...
    remove_comp(&mxt);
    cn = mxt;
    cn->h = (unsigned short int *) malloc(J * sizeof(unsigned short int));
    cn->rd0 = (unsigned short int*) malloc(q * sizeof(unsigned short int));
    cn->rd1 = (unsigned short int*) malloc(q * sizeof(unsigned short int));
  }
  sample_hap(cn);
  int* temp;
  temp = new int[J/10+1];
  conversion(temp, cn->h, J);
  for (ll=0; ll < q; ll++){
    p = seq_distance_new(temp, readtable2[ll], J);
    cn->rd0[ll] = p[0];
    cn->rd1[ll] = p[1];
  }
  
  // have to go through the list with ->next, because dont want to loose connection between the elements of the list
  // needed for remove_comp, which returns the ->next element of the deleted cluster at the position of the given argument
  while(cn->next != NULL){
    // have to allocate mem anyway, otherwise would get error, when freeing memory, which happens when removing cluster
    cn->next->h = (unsigned short int *) malloc(J * sizeof(unsigned short int));
    cn->next->rd0 = (unsigned short int*) malloc(q * sizeof(unsigned short int));
    cn->next->rd1 = (unsigned short int*) malloc(q * sizeof(unsigned short int));
    
#ifdef DEBUG
    printf("cluster %d has size %d\n", cn->next->ci, cn->next->size);
#endif

    if(cn->next->size == 0) {
#ifdef DEBUG
      printf("removing some node... p=%p\n",cn->next);
#endif
      remove_comp(&cn->next);
    }else{
      sample_hap(cn->next);
      conversion(temp, cn->next->h, J);
      for (ll=0; ll < q; ll++){
	p = seq_distance_new(temp, readtable2[ll], J);
	cn->next->rd0[ll] = p[0];
	cn->next->rd1[ll] = p[1];
      }
      cn = cn->next;
    }
  }
  delete[] temp;
  
  // define the predecessor of each read
  cn = mxt;
  
  rn = cn->rlist;
  while (rn != NULL){
    i = rn->ri;
    c_ptr[i] = NULL;
    rn = rn->next;
  }
  
  while (cn->next != NULL){
    
    rn = cn->next->rlist;
    while (rn != NULL){
      i = rn->ri;
      c_ptr[i] = cn;
      rn = rn->next;
    }
    
    cn = cn->next;
  }
  
  /*  
  for(i=0; i<n; i++)
    printf("predecessor of %i is %p\n", i, c_ptr[i]);
  */
  unsigned short int rr;
  
  int bmax = 0;
  unsigned short int bi;

  /* sample reference genome */
  for(i=0; i<J; i++){
    for(rr=0;rr<B;rr++)
      cbase[rr] = 0;
    
    for (rr=0; rr<q; rr++) {
      if(r[readtable2[rr]->mindex][i] < B) {
     	//cbase[r[rr][i]]++;
	cbase[r[readtable2[rr]->mindex][i]] += readtable2[rr]->weight;
	//bmax += 1;
	bmax += readtable2[rr]->weight;
      }
    }
    if(bmax > 0) { 
      bmax = cbase[0];
      bi = 0;
      for(rr=1; rr<B; rr++) {
	if(cbase[rr] > bmax) {
	  bmax = cbase[rr];
	  bi = rr;
	}
      }
      h[i] = bi;
    }
    else{
      h[i] = B;
    }
  }
  
  return;
}


double sample_ref() {
  cnode* cn;
  unsigned int i,j;
  double b1,b2;
  unsigned int K1 = 0;
  gsl_ran_discrete_t* g;
  ofstream err_file("error_ref.log");
  double max_log_pbase;
  int max_cbase;
  unsigned int base_id, dk1=0;
  unsigned int countbases = 0;

  for(j=0;j<J;j++) {

    K1 = 0;
    for(i=0;i<B;i++) {
      cbase[i] = 0;
      pbase[i] = 0.0;
    }

    cn = mxt;
    
    while(cn != NULL) {
      if(cn->size > 0) {
	if(cn->h[j] < B) {
	  cbase[cn->h[j]]++;
	  K1++;
	}
      }
      cn=cn->next;
    }

    countbases += K1;

    if(gam < .2)gam = 0.2;

    b1 = gam;
    b2 = (1.0-gam)/((double)B - 1.0);

    if(b1 == 0.0) {
      err_file << "# There is something wrong with GAMMA ( == 0.0), sample ref\n";
      print_stats(err_file, mxt, J);
      err_file << "# reference:\n# ";
      for(i=0;i<J;i++) {
	err_file << i2dna_code[h[i]];
      }
      err_file << "\n";
      exit(EXIT_FAILURE);
    }

    if(K1 > 0) {
      if(b1 != 1.0) {
	
	for(i=0; i<B; i++){
	  log_pbase[i] = cbase[i] * gsl_sf_log(b1) + (K1 - cbase[i]) * gsl_sf_log(b2);
	}
	max_log_pbase = log_pbase[0];
	base_id = 0;
	
	for(i=1; i<B; i++) {
	  if(log_pbase[i] > max_log_pbase) {
	    max_log_pbase = log_pbase[i];
	    base_id = i;
	  }
	}
	
	for(i=0; i<B; i++) {
	  log_pbase[i] -= max_log_pbase;
	  if(i == base_id) {
	    pbase[i] = 1.0;
	  }
	  else{
	    if(log_pbase[i] < double_threshold_min) {
	      pbase[i] = 0.0;
	    }
	    else{
	      pbase[i] = gsl_sf_exp(log_pbase[i]);
	    }
	  }
	}
	
	g = gsl_ran_discrete_preproc(B, pbase);
	h[j] = gsl_ran_discrete(rg, g);
	gsl_ran_discrete_free(g);
      }
      else{ // gamma == 1.0
	max_cbase = cbase[0];
	for(i=1; i<B; i++) {
	  if(cbase[i] >= max_cbase) {
	    max_cbase = cbase[i];
	  }
	}
	for(i=0; i< B; i++) {
	  if(cbase[i] == max_cbase) {
	    pbase[i] = 1.0;
	  }
	  else {
	    pbase[i] = 0.0;
	  }
	}
	g = gsl_ran_discrete_preproc(B, pbase);
	h[j] = gsl_ran_discrete(rg, g);
	gsl_ran_discrete_free(g);
      }
    }
    else{ // K1 == 0, that is all N's
      h[j] = B;
    }
    dk1 += cbase[h[j]];
  }
  return (double)dk1/(double)countbases;
}

void sample_hap(cnode* cn){
  
  rnode* tr;
  unsigned int i, j, tot_reads;
  gsl_ran_discrete_t* g;
  double b1, b2;

  double max_log_pbase;
  int max_cbase;
  unsigned int base_id;

#ifdef DEBUG
  printf("Haplotype %i is\n", cn->ci);
#endif

  for (j=0; j<J; j++){
    
    tot_reads = 0;
    for(i=0; i<B; i++) {
      cbase[i] = 0;
      pbase[i] = 0.0;
    }

    tr = cn->rlist;
    while (tr != NULL){ // correct for missing data
      if (r[readtable2[tr->ri]->mindex][j] < B) {
	cbase[r[readtable2[tr->ri]->mindex][j]] += readtable2[tr->ri]->weight;
	//cbase[r[readtable2[tr->ri]->mindex][j]] ++;
	//cbase[r[tr->ri][j]]++;
	tot_reads += readtable2[tr->ri]->weight;
	//tot_reads++;
      }
      tr = tr->next;
    }

    b1 = theta;
    b2 = (1. - theta)/((double)B-1);
    
    if(b1 == 0.0) {
      printf("# There is something wrong with THETA ( == 0.0)\n");
      exit(EXIT_FAILURE);
    }

    if(b1 != 1.0) {	
      for(i=0; i<B; i++){
	if(tot_reads > 0) { //base is not N: sample from reads in the cluster
	  log_pbase[i] = cbase[i] * gsl_sf_log(b1) + (tot_reads - cbase[i]) * gsl_sf_log(b2);
	}
	else //base is N: sample from all reads
	  log_pbase[i] = ftable[j][i] * gsl_sf_log(b1) + (ftable_sum[j] - ftable[j][i]) * gsl_sf_log(b2);	      
      }
      max_log_pbase = log_pbase[0];
      base_id = 0;
      for(i=1; i<B; i++) {
	if(log_pbase[i] > max_log_pbase) {
	  max_log_pbase = log_pbase[i];
	  base_id = i;
	}
      }
      for(i=0; i<B; i++) {
	if(i == base_id) {
	  pbase[i] = 1.0;
	}else{
	  log_pbase[i] -= max_log_pbase;
	  if(log_pbase[i] < double_threshold_min) {
	    pbase[i] = 0.0;
	  }else{
	    //cout<<"log_pbase["<<i<<"] = "<<log_pbase[i]<<endl;
	    pbase[i] = gsl_sf_exp(log_pbase[i]);
	  }
	}
      }
	
      g = gsl_ran_discrete_preproc(B, pbase);
      cn->h[j] = gsl_ran_discrete(rg, g);
      gsl_ran_discrete_free(g);
    }
    else{ // theta == 1.0
      if(tot_reads>0){//base not N: sample from reads in cluster.
	max_cbase = cbase[0];
	base_id = 0;
	for(i=1; i<B; i++) {
	  if(cbase[i] >= max_cbase) {
	    max_cbase = cbase[i];
	    base_id = i;
	  }
	}
	cn->h[j] = base_id;
      }
      else{//missing data: sample from all reads.
	max_cbase = ftable[j][0];
	base_id = 0;
	for(i=1; i<B; i++) {
	  if(ftable[j][i] >= max_cbase) {
	    max_cbase = ftable[j][i];
	    base_id = i;
	  }
	}
	cn->h[j] = base_id;
      }
    }
    
    /*
  else{ // tot_reads == 0
      //cn->h[j] = h[j]; // equals to the reference
      //cn->h[j] = B; // equals to N
      for(i=0; i<B; i++)
	pbase[i] = ftable[j][i];
      g = gsl_ran_discrete_preproc(B, pbase);
      cn->h[j] = gsl_ran_discrete(rg, g);
      gsl_ran_discrete_free(g);


    }
    */
#ifdef DEBUG
    printf("%i ", cn->h[j]);
#endif
  }
  
#ifdef DEBUG
  printf("\n");
#endif

  return;  
}


void check_size(const cnode* cst, unsigned int n){
  unsigned int this_n = 0;
  ofstream err_file("error_size.log");
  while(cst != NULL){
    this_n += cst->size;
    cst = cst->next;
  }
  
  if(this_n != n){
    print_stats(err_file, mxt, J);
    printf("!!!! STOP !!!!\n");
    printf("size is now %i\n", this_n);
    exit(EXIT_FAILURE);
  }
  
  return;
}

int count_classes(const cnode* cst){
  unsigned int ck = 0;
  
  while(cst != NULL){
    ck++;
    cst = cst->next;
  }

  return ck;
}

int isdna (char c){
  int ans = (c == 'A' || c == 'C' || c == 'G' || c == 'T' ||
	     c == 'a' || c == 'c' || c == 'g' || c == 't' ||
	     c == 'N' || c == 'n' || c == '-');
    return ans;
}

unsigned int d2i (char c){

  if(c == 'A' || c == 'a')
    return 0;
  if(c == 'C' || c == 'c')
    return 1;
  if(c == 'G' || c == 'g')
    return 2;
  if(c == 'T' || c == 't')
    return 3;
  if(c == '-')
    return 4;
  if(c == 'N' || c == 'n') // this stands for missing data
    return B;
  printf("%c not a DNA base", c);
  exit(EXIT_FAILURE);
  return 0;
}

int* seq_distance_new(int* A, crnode* B, int seq_length){
  int *X, i;
  X = new int[seq_length/10 + 1];

  dist[0] = 0;
  dist[1] = 0;

  for(i=0; i<seq_length/10 + 1; i++, A++){
    X[i] = *A ^ B->creads[i];
    X[i] = (X[i] & one_int) | ((X[i] & two_int) >> 1) | ((X[i] & four_int) >> 2);

    //Count the number of set bits (Knuth's algorithm)
    while(X[i]){
      dist[0] ++;
      X[i] &= X[i] - 1;
    }
  }
  dist[1] = seq_length - dist[0];
  dist[0] = dist[0] - B->missing;
  delete[] X;
  return dist;
}

int* seq_distance_rr(int* A, crnode* B, int seq_length){
  int *X, i;
  X = new int[seq_length/10 + 1];

  dist[0] = 0;
  dist[1] = 0;

  for(i=0; i<seq_length/10 + 1; i++, A++){
    X[i] = *A ^ B->creads[i];
    X[i] = (X[i] & one_int) | ((X[i] & two_int) >> 1) | ((X[i] & four_int) >> 2);

    //Count the number of set bits (Knuth's algorithm)
    while(X[i]){
      dist[0] ++;
      X[i] &= X[i] - 1;
    }
  }
  //ns = B->missing;
  //dist[1] = seq_length - dist[0] - N->missing;
  //dist[0] = dist[0] - ns;
  delete[] X;
  return dist;
}

int* seq_distance(unsigned short int* a, unsigned short int* b, int seq_length){
  int i, ns=0;
  //int d=sizeof(unsigned short int);
  dist[0] = 0;
  dist[1] = 0;
  for (i=0; i<seq_length; i++, a++, b++){
    if (*a != B && *b != B){
      dist[0] += (*a != *b);
    }else{
      ns++;
    }
  }
  dist[1] = seq_length - dist[0] - ns;
  return dist;
}

ssret* sample_class(unsigned int i, unsigned int step){
  /*****************************************
   the core of the Dirichlet process sampling 
  ******************************************/
  unsigned int dist, nodist, removed = 0, sz, ll;
  unsigned int tw;
  int* p;
  //  int local_ci;
#ifdef DEBUG
  unsigned int j;
#endif
  double b1, b2; //, pow1, pow2;
  cnode* to_class;
  cnode* from_class;
  size_t st = 0, this_class;
  gsl_ran_discrete_t* g;
  cnode* cn;
  rnode* rn = NULL;
  double min_log_P, max_log_P, delta_log;
  //  unsigned int class_id;
  int read_came_from_current;
  int* temp;
  temp = new int[J/10+1];
#ifdef DEBUG
  for(removed=0; removed < q; removed++)
    printf("c_ptr[%d]=%p\n", removed, c_ptr[removed]);
  removed = 0;
  printf("---------------------------------------------------\n");
  printf("-----------> sampling class for %ith read <----------\n", i);
  printf("---------------------------------------------------\n");
  //  printf("------------ c_ptr = %p with size = %d\n", c_ptr[i], c_ptr[i]->size);
  cout << "----------- mxt here is "<< mxt << "---------------------" << endl;
  //  printf("------------- NOW STATS----------\n");
  //  print_stats(out_file, mxt, J);
#endif

  // if the class is populated only by i, remove it

  // if it's in the head
  if(c_ptr[i] == NULL && mxt->size == 1){
    
    rn = mxt->next->rlist;
    
    while(rn != NULL){
      c_ptr[rn->ri] = NULL;
      rn = rn->next;
    }
    
    remove_comp(&mxt);
    removed = 1;
#ifdef DEBUG
    printf("----------- REMOVED SIZE 1 NODE ----------\n");
#endif
  }
  
  // if it's not in the head
  if(c_ptr[i] != NULL && c_ptr[i]->next->size == 1){
    
    // if i not in the last node update the following predecessor
    if(c_ptr[i]->next->next != NULL){
      
      rn = c_ptr[i]->next->next->rlist;
      while(rn != NULL){
	c_ptr[rn->ri] = c_ptr[i];
	rn = rn->next;
      }
      
    }
    
    remove_comp(&(c_ptr[i]->next));
    removed = 1;
#ifdef DEBUG
    printf("----------- REMOVED SIZE 1 NODE ----------\n");
#endif
  }
  
  /// run through the classes to count them
  cn = mxt;
  st = 0;
  while (cn->next != NULL){
    st++;
    cn = cn->next;
  }
  st++;
  //  log_P = (double*) realloc(log_P, st*sizeof(double));
  //  P = (double*) realloc(P, st*sizeof(double));
  
  /**********************************************************
   run through the populated classes to assign a probability
  ***********************************************************/
  
  cn = mxt;
  
  if( cn == NULL){
    printf("sampling classes, no classes found\n");
    exit(EXIT_FAILURE);
  }
  
  st = 0;
  
  b1 = theta;
  b2 = (1. - theta)/((double)B - 1.);

  if(theta == 0.0) {
    cerr << "# There is something wrong wih THETA! ( == 0.0 )\n";
    exit(EXIT_FAILURE);
  }
  if(gam == 0.0) {
    cerr << "# There is something wrong with GAMMA! ( == 0.0 )\n" ;
    exit(EXIT_FAILURE);
  }

  while (cn != NULL){
    if(cn->size > 0) {
      sz = cn->size;
      read_came_from_current = 0;
      if(removed == 0){
	if(c_ptr[i] != NULL && cn == c_ptr[i]->next) {
	  sz--;
	  read_came_from_current = 1;
	}
	if(c_ptr[i] == NULL && cn == mxt) {
	  sz--;
	  read_came_from_current = 1;
	}
      }
      
      if(sz == 0){
	cerr << "THIS IS ZERO!!!\n";
	exit(EXIT_FAILURE);
      }
      cl_ptr[st] = cn;

      tw = 0;

      tw = weight_shift(cn, i); //total weight of class - weight of read i if present in the class
      if (tw == 0){
	cout << cn << "cannot be zero "<< i<< endl;
	exit(EXIT_FAILURE);
      }
      
      dist = cn->rd0[i]; //seq_distance(cn->h, r[i], J);
      nodist = cn->rd1[i];
      
      if(b1 != 1.0) { // theta < 1.0
	log_P[st] = gsl_sf_log((double)tw);
	log_P[st] += nodist * gsl_sf_log(b1);
	log_P[st] += dist * gsl_sf_log(b2);
	P[st] = 1.0; // all probabilities, which should change afterwards, set to 1
      }
      
      else{ // theta == 1.0, not needed
	if (dist != 0){
	  log_P[st] = double_threshold_min - 1.0;
	  P[st] = 0.0;
	}
	else{
	  log_P[st] = gsl_sf_log((double)tw);
	  P[st] = 1.0; // same as above, P != 0 later
	}	
      }
    }
    else{
      cerr << "------********* CN->SIZE = 0 **********----------\n";
      P[st] = 0.0; // all prob, which shouldn't change, set to 0
    }
    
    st++;
    cn = cn->next;
  }
  // plus the newly instantiated one
  conversion(temp, h, J);  
  p = seq_distance_new(temp, readtable2[i], J);
  delete[] temp;
  dist = p[0];
  nodist = p[1];
  
  b1 = (theta * gam) + (1. - gam)*(1. - theta)/((double)B - 1.);
  b2 = (theta + gam + B*(1. - gam * theta) - 2.)/(gsl_sf_pow_int((double)B - 1., 2));
  
  if((theta == 1.0) && (gam == 1.0)) {
    if (dist == 0){
      log_P[st] = gsl_sf_log(alpha);
      P[st] = 1.0;
    }
    else{
      log_P[st] = double_threshold_min - 1.0;
      P[st] = 0.0;
    }
  }
  else {
    log_P[st] = gsl_sf_log(alpha) + gsl_sf_log((double)(readtable2[i]->weight)) +  nodist * gsl_sf_log(b1) + dist * gsl_sf_log(b2);
    P[st] = 1.0;
  }
  cl_ptr[st] = NULL;
  
  max_log_P = *max_element(log_P, log_P+st); //! renormalization
  min_log_P = *min_element(log_P, log_P+st); //! renormalization

  if (max_log_P >= 0)
    delta_log = -max_log_P;//0.5 * (min_log_P + max_log_P);
  else
    delta_log = max_log_P;//-0.5 * (min_log_P + max_log_P);

  
//   class_id = st;  
//   for(ll=0; ll<st; ll++) {
//     if((max_log_P == log_P[ll]) && (P[ll] > 0.0))
//       class_id = ll;
//   }

  for(ll=0; ll<=st; ll++) {
    if(P[ll] > 0.0) {      
      log_P[ll] += delta_log;
      if(log_P[ll] < double_threshold_min)
	P[ll] = DBL_MIN;
      else if (log_P[ll] > double_threshold_max)
	P[ll] = DBL_MAX;
      else {
        P[ll] = gsl_sf_exp(log_P[ll]);
      }
    }// else P[i] = 0, from above
  }
  
#ifdef DEBUG
  for (j=0; j<=st; j++)
    printf("with P[%i] = %e to class %p\n", j, P[j], cl_ptr[j]);
#endif
  
  g = gsl_ran_discrete_preproc(st+1, P);
  this_class = gsl_ran_discrete(rg, g);
  gsl_ran_discrete_free(g);
  
#ifdef DEBUG
  printf ("extracted class is = %lu\n", this_class);
#endif
  
  from_class = (c_ptr[i] != NULL) ? c_ptr[i]->next : mxt;
  to_class = cl_ptr[this_class];


  /***************************************
    move the read to the extracted class
  ****************************************/

  if( removed == 0 ){
    
    if ( from_class == to_class ){
#ifdef DEBUG
      printf ("from %p to itself\n", from_class);
#endif
      
      //      free(P);
      //      free(cl_ptr);
      
      
      res->untouched = 1;
      res->proposed = 0;
      res->to_class = to_class;
      return res;
    }
    
    else if ( to_class != NULL){
#ifdef DEBUG
      printf("moving the read from %p to %p\n", from_class, to_class);
#endif
      
      c_ptr[i] = c_ptr[ to_class->rlist->ri ]; // predecessor is taken from the already present read
      move_read(i, &from_class, &to_class);
      (from_class->size)--;
      (to_class->size)++;
      
      //      free(P);
      //      free(cl_ptr);
      
      
      res->untouched = 0;
      res->proposed = 0;
      res->to_class = to_class;
      return res;
    }
    
    else if (to_class == NULL){
#ifdef DEBUG    
      printf("moving %i to a new class from %p\n", index, from_class);
#endif
      
      remove_read( search_read(&from_class->rlist, i) );
      (from_class->size)--;

      add_comp(&mxt, &hst, lowest_free, 1, NULL, NULL, NULL, NULL, J, step, record);
      
      lowest_free++;
      
      add_read(&mxt->rlist, i);
      
      mxt->h = (unsigned short int*) malloc(J * sizeof(unsigned short int));
      sample_hap(mxt);
      temp = new int[J/10 + 1];
      conversion(temp, mxt->h, J);
   
      mxt->rd0 = (unsigned short int*) malloc(q * sizeof(unsigned short int));
      mxt->rd1 = (unsigned short int*) malloc(q * sizeof(unsigned short int));
      for (ll=0; ll < q; ll++){
	p = seq_distance_new(temp, readtable2[ll], J);
	//p = seq_distance(mxt->h, r[ll], J);
	mxt->rd0[ll] = p[0];
	mxt->rd1[ll] = p[1];
      }
      delete[] temp;
      c_ptr[i] = NULL;
      
      rn = mxt->next->rlist;
      
      if(rn == NULL){
	printf("STOP!!! There should be something\n");
	exit(21);
      }
      
      while (rn != NULL){
	c_ptr[rn->ri] = mxt;
	rn = rn->next;
      }
      
      //      free(P);
      //      free(cl_ptr);
      
      
      res->untouched = 0;
      res->proposed = 1;
      res->to_class = mxt;
      return res;
    }
    
  }
  
  else if (removed == 1){
#ifdef DEBUG
    printf("moving having removed\n");
#endif
    if (to_class != NULL){
      add_read(&to_class->rlist, i);

      (to_class->size)++;
      c_ptr[i] = c_ptr[ to_class->rlist->next->ri ];
      
      
      res->untouched = 0;
      res->proposed = 0;
      res->to_class = to_class;
    }
    else if (to_class == NULL){
      add_comp(&mxt, &hst, lowest_free, 1, NULL, NULL, NULL, NULL, J, step, record);
      lowest_free++;
      res->to_class = mxt;
      c_ptr[i] = NULL;
    
      cn = mxt->next;
      rn = cn->rlist;
      while (rn != NULL){
	c_ptr[rn->ri] = mxt;
	rn = rn->next;
      }
      
      add_read(&mxt->rlist, i);
    
      mxt->h = (unsigned short int*) malloc(J * sizeof(unsigned short int));
      sample_hap(mxt);
      temp = new int[J/10 + 1];
      conversion(temp, mxt->h, J);
      
      mxt->rd0 = (unsigned short int*) malloc(q * sizeof(unsigned short int));
      mxt->rd1 = (unsigned short int*) malloc(q * sizeof(unsigned short int));
      for (ll=0; ll < q; ll++){
	p = seq_distance_new(temp, readtable2[ll], J);
	//p = seq_distance(mxt->h, r[ll], J);
	mxt->rd0[ll] = p[0];
	mxt->rd1[ll] = p[1];
      }
      delete[] temp;
      
      res->untouched = 0;
      res->proposed = 1;
      
    }
    
    return res;
  }
  
  fprintf(stderr, "WARNING!!! I DIDN'T SAMPLE READ %d\n", i);
  // delete[] temp;

  return res; 
}


void write_assignment (unsigned int it, unsigned int new_proposed, const cnode* tn) {
  
  rnode* rn;
  unsigned int j = 0, k;
  
  fprintf(assign, "# %d iterations\n", it-1);
  fprintf(assign, "\n#gamma = %f\n", gam);
  fprintf(assign, "\n#made %i new clusters\n\n", new_proposed);
  
  while(tn != NULL){
    
    rn = tn->rlist;
    while(rn != NULL){
      k = 0;
      while(id[rn->ri][k] != '\n'){
	fputc(id[rn->ri][k], assign);
	k++;
      }
      fprintf(assign, " -> %i\n", j);
      rn = rn->next;
    }
    j++;
    tn = tn->next;
  }
  
  //  fprintf(assign, "\n");
  fprintf(assign, "\n#final number of components = %i\n\n", j);
}

void create_history(unsigned int step){
  // copy haplotypes from mxt to create the history
  
  cnode* tn;

  tn = mxt;
  while(tn != NULL){
    add_hap(&hst, NULL, tn->ci, J, step);
    tn->hh = hst;
    memcpy(hst->h, tn->h, J*sizeof(unsigned short int));
    tn = tn->next;
  }
  
  return;
}

void i2dna_string(char* ch, short unsigned int* h, int J){
  // char* ca;
  //ca = (char*) calloc(J, sizeof(char*));
  for (int i=0; i<J; i++){
    ch[i] = i2dna_code[h[i]];
  }
  return;// (string) ca;
}

void record_conf(cnode* tn, unsigned int step){
  //! record configuration of a single node at a single step
  
  rnode* tr;
  char* ca;
  ca = (char*) calloc(J, sizeof(char*));
  i2dna_string(ca, tn->h, J);
  string h = (string) ca;
  free(ca);
  freq_map::iterator freq_iter;
  ++support[h]; //! support is updated

  int tw = weight(tn);
  
  freq_iter = freq.find(h);
  if(freq_iter != freq.end()){ //! freq is updated if haplotype is present
    //freq_iter->second[step] = tn->size;
    freq_iter->second[step] = tw;
  }
  else{
    freq[h] = (int*)calloc(HISTORY, sizeof(unsigned int));
    freq[h][step] = tw; //! freq is updated if haplotype is NOT present
  }

  tr = tn->rlist;
  while (tr != NULL){
    read2hap[tr->ri][h]++;
    tr = tr->next;
  }
  h.clear();
  return;
}

void old_record_conf(cnode* tn, unsigned int step){
  /**
     record configuration of a single node at a single step
   */
  int* p;
  rnode* rn;
  
  if(tn->step < step) { // if the node has been created in a previous step
    //p = seq_distance_new(conversion(tn->h, J), conversion(tn->hh->h, J), J);
    p = seq_distance(tn->h, tn->hh->h, J);
    if(p[0] != 0) { // if the haplotype has changed
      add_hap(&hst, tn->hh, tn->ci, J, step);
      memcpy(hst->h, tn->h, J*sizeof(unsigned short int)); // copy into history 
      tn->hh = hst; // and then update
    }
  }else{ // update the haplotype
    memcpy(tn->hh->h, tn->h, J*sizeof(unsigned short int));
  }

  rn = tn->rlist;
  while(rn != NULL){
    //    ass_hist[rn->ri][step - iter + HISTORY - 1] = tn->hh;  
    rn = rn->next;
  }
  
  return;
}



void cleanup() {
  unsigned int i;
  freq_map::iterator fit;
  sup_map::iterator sit;

  free(P);
  free(log_P);
  free(cl_ptr);
  free(c_ptr);
  free(pbase);
  free(cbase);
  free(log_pbase);
  free(h);
  
  for(i=0; i<n; i++){
    delete[] readtable[i]->creads;
    delete[] readtable[i]->mindices;

    free(readtable[i]);
    free(r[i]);
    free(id[i]);
  }
  free(readtable);
  free(r);
  free(id);

  for (i=0; i<J; i++)
    delete[] ftable[i];
  free(ftable);
  
   for(i=0; i<q; i++){
     delete[] readtable2[i]->mindices;
//     readtable2[i]->creads = NULL;
//     delete[] readtable2[i]->creads;
     free(readtable2[i]);
   }
   free(readtable2);
  

  support.clear();
  
  for ( fit=freq.begin() ; fit != freq.end(); fit++ )
    free((*fit).second);
  freq.clear();
  
  for(i=0; i<n; i++)
    read2hap[i].clear();
  delete read2hap;

  delete[] ftable_sum;

  // should maybe cleanup all nodes ...

}


int compare(const void* a, const void* b) {
  unsigned int i=0;
  unsigned short int *h1, *h2;
  
  h1 = *(unsigned short int**)a;
  h2 = *(unsigned short int**)b;

  for(i=0; i<J; i++) {
    if (*h1 != *h2) {
      return *h1 - *h2;
    }
    h1++;
    h2++;
  }
  return 0;
}

int compare_hnss_seq(const void* a, const void* b) {
  unsigned int i=0;
  
  hnode_single h1 = **(hnode_single**)a;
  hnode_single h2 = **(hnode_single**)b;

  for(i=0; i<J; i++) {
    if (h1.h[i] != h2.h[i]) {
      return h1.h[i] - h2.h[i];
    }
  }
  return 0;
}

int compare_hnss_count(const void* a, const void* b) {
  hnode_single h1 = **(hnode_single**)a;
  hnode_single h2 = **(hnode_single**)b;
  
  return (h2.count - h1.count);
}


struct invcomp{
  //! used in multimap to store keys in descending order
  bool operator() (int i1, int i2) const
  {return i1>i2;}
};

void write_posterior_files(string instr){
  /** All the posterior information is written here 
   */
  //  cout << "HISTORY = " << HISTORY << endl;
  int insize = instr.rfind('.');
  string corstr = instr.substr(0,insize).append("-cor.fas");
  string supstr = instr.substr(0,insize).append("-support.fas");
  string freqstr = instr.substr(0,insize).append("-freq.csv");
  string str2;
  int i=0, rcount=0, wcount=0;
  float mean_freq, supp_fract;
  //  float supp_thresh = 0.5;
  
  ofstream cor_file(corstr.c_str());
  ofstream supp_file(supstr.c_str());
  ofstream freq_file(freqstr.c_str());
  
  sup_map::iterator si;
  multimap<int, string, invcomp> rev_sup; //! use a multimap to sort by value the support map
  multimap<int, string, invcomp>::iterator ri;
  
  for (si = support.begin(); si != support.end(); ++si)
    rev_sup.insert(pair<int, string>(si->second, si->first));
  
  // header of freq_file
  freq_file << setfill('0');
  freq_file << "#haplotype\tsupport";
  for (unsigned int k=0; k<HISTORY; k++)
    freq_file << "\treads";
  freq_file << "\n";
  
  for (ri = rev_sup.begin(); ri != rev_sup.end(); ++ri){
    supp_fract = (float)(ri->first) / (float)HISTORY;
    
    if (supp_fract >= 0.01){  
      mean_freq = 0;
      for (unsigned int k=0; k<HISTORY; k++)
	mean_freq += freq[ri->second][k];
      mean_freq /= HISTORY;    
      
      supp_file << ">hap_" << i << "|" << "posterior=" << supp_fract << " ave_reads=" << mean_freq << "\n";
      supp_file <<  ri->second << "\n";

      freq_file << ">hap_" << setw(5) << i << "\t" << ri->first;
      for (unsigned int k=0; k<HISTORY; k++)
	freq_file << "\t" << freq[ri->second][k];
      freq_file << "\n";
      ++i;
    }
    
  }
  
  for(unsigned int readi = 0; readi < q; readi++){
    rev_sup.clear(); // reuse the multimap
    for (si = read2hap[readi].begin(); si != read2hap[readi].end(); ++si)
      rev_sup.insert(pair<int, string>(si->second, si->first));
    
    ri = rev_sup.begin();
    for(wcount = 0; wcount < readtable2[readi]->weight; wcount++){
      cor_file << id[readtable2[readi]->mindices[wcount]] << "|posterior=" << float(ri->first)/HISTORY << "\n";
      cor_file << ri->second << "\n";
    }
    rcount++;
  }
  
  return;
}

double setfinalhaplotype(unsigned int i)  {
  unsigned int j,k;
  int hap; //last haplotype
  int *p;
  double quality;
  
  // sort lexicographically the haplotypes
  if(i==0) {
    haplotypes = (unsigned short int**)calloc(HISTORY, sizeof(unsigned short int*));
    for(k=0; k<HISTORY; k++) {
      haplotypes[k] = (unsigned short int*)calloc(J, sizeof(unsigned short int));
    }
    ch = (int*)calloc(HISTORY,sizeof(int));
  }
  
  //  for(k=0; k<HISTORY; k++)
  //    memcpy(haplotypes[k], ass_hist[i][k]->h, J*sizeof(unsigned short int));
  
  
  qsort(haplotypes, HISTORY, sizeof(unsigned short int*), compare);
  // qsort returns haplotype (pointers to the ascending sorted haplotypes)

  
  hap = 0;
  ch[0] = 1;
  
  for(k = 1; k < HISTORY; k++) {
    //p = seq_distance_new(conversion(haplotypes[k-1], J), conversion(haplotypes[k], J), J);
    p = seq_distance(haplotypes[k-1], haplotypes[k], J);
    
    if (p[0] == 0 ) {
      ch[hap] += 1;
      
    }else{
      hap = k;
      ch[hap] = 1;
    }
    
  }

  int max = 0;
  for(k=0; k<HISTORY; k++) {
    if (ch[k] >= max){
      max = ch[k];
      hap = k;
    }
  }
  
  quality = ch[hap]/(double)HISTORY;

  for(j=0;j<J;j++)
    if(r[i][j] < B)
      r[i][j] = haplotypes[hap][j];
  
  if(i == n-1) {
    for(k=0;k<HISTORY;k++) {
      free(haplotypes[k]);
    }
    free(haplotypes);
    free(ch);
  }
  
  return quality;

}


void write_haplotype_frequencies(char* filename, unsigned int hcount) {
  FILE* fp;
  unsigned int i, j, hap, running_sum=0;
  hnode_single** all_haplo;
  int* p;
  
  /*
  int k;
  for(i=0; i<n; i++) {
    printf("read_%d\n", i);
    for(j=0; j<HISTORY; j++) {
      printf("step_%d\n", j);
      for(k=0; k<J; k++) {
	printf("%c",i2dna_code[ass_hist[i][j]->h[k]]);
      }
      printf("\n");
    }
      printf("\n\n");
  }
  */
  
  all_haplo = (hnode_single**)malloc(n*HISTORY*sizeof(hnode_single*));
  for(i=0; i<n; i++) {
    for(j=0; j<HISTORY; j++) {
      all_haplo[i*HISTORY+j] = (hnode_single*)malloc(sizeof(hnode_single));
      all_haplo[i*HISTORY+j]->h = (unsigned short int*)malloc(J*sizeof(unsigned short int));
      //      all_haplo[i*HISTORY+j]->hi = ass_hist[i][j]->hi;
      //      all_haplo[i*HISTORY+j]->step = ass_hist[i][j]->step;
      all_haplo[i*HISTORY+j]->count = 1;
      //      memcpy(all_haplo[i*HISTORY+j]->h, ass_hist[i][j]->h, J*sizeof(unsigned short int));
    }
  }
  printf("# n=%d, J=%d, HISTORY=%d\n", n, J, HISTORY);
  printf("# n * J * HISTORY=%d\n", n*J*HISTORY);
  /*
  printf("# The number of bytes in a short int is %ld.\n", sizeof(short int));
  printf("# The number of bytes in a unsigned short int is %ld.\n", sizeof(unsigned short int));
  printf("# The number of bytes in ass_hist is %ld.\n", sizeof(ass_hist[1][1]));
  printf("# The number of bytes in all_haplo is %ld.\n", sizeof(hnode_single));
  */
  qsort(all_haplo, n*HISTORY, sizeof(hnode_single*), compare_hnss_seq);
  
  hap = 0;
  for(i=1; i<n*HISTORY; i++) {
    //p = seq_distance_new(conversion(all_haplo[i]->h, J), conversion(all_haplo[i-1]->h, J), J);
    p = seq_distance(all_haplo[i]->h, all_haplo[i-1]->h, J);
    if(p[0] == 0) {
      all_haplo[hap]->count++;
      all_haplo[i]->count = 0;
    }else{
      hap = i;
    }
  }
  
  qsort(all_haplo, n*HISTORY, sizeof(hnode_single**), compare_hnss_count);
  if(hcount <= 0) {
    hcount = n*HISTORY;
  }
  
  hap=0;
  fp = fopen(filename,"w");
  for(i=0; i<n*HISTORY; i++) {
    if(all_haplo[i]->count > 0) {
      running_sum += all_haplo[i]->count;
      fprintf(fp, ">haplotype_%04d|%.9f\n", hap, (double)(1.0*all_haplo[i]->count)/(1.0 * n * HISTORY));
      for(j=0; j<J; j++) {
	fprintf(fp,"%c", i2dna_code[all_haplo[i]->h[j]]);
      }
      fprintf(fp,"\n");
      hap++;
      if(hap >= hcount)
	break;
    }
  }
  fclose(fp);
  if (1.0*fabs(running_sum - (n*HISTORY))/(n*HISTORY) > 0.001){
    printf("WATCH OUT, running_sum=%d and not %d\n", running_sum, n*HISTORY);
    }
  for(i=0; i<n*HISTORY; i++) {
    free(all_haplo[i]->h);
  }
  free(all_haplo);

}

void print_stats(std::ofstream& out_file, const cnode* cn, unsigned int J){
  
  unsigned int i, j;
  rnode* p;
  out_file << "# -------- PRINTING STATS ---------\n";
  if( cn == NULL)
    out_file << "# no components in this list\n";
  
  while (cn != NULL){
    
    i=0;
    p = cn->rlist;
    
    if(p == NULL)
      out_file << "# no reads in this component\n";
    
#ifdef DEBUG
    if(p != NULL)
      out_file << "# reads in the list:\t";
#endif

    while(p != NULL){
      i++;
      p = p->next;
    }
    
    out_file << "\n# component " << cn->ci << " at " << (void*)cn << " has " << i << " reads size= " << weight_final(cn) <<"\n";
    out_file << "# haplotype is:\n# >h" << cn->ci <<"\n# ";
    for(j=0; j<J; j++)
      out_file << i2dna_code[cn->h[j]];
    out_file << "\n\n";

    cn = cn->next;
    
  }
  out_file << "# -------- STATS PRINTED ----------\n"; 
}

int parsecommandline(int argc, char** argv){
// get command line parameters
  
  char c;
  while((c = getopt(argc, argv,"i:j:a:t:o:m:K:k:R:c:Hh")) != -1){
    
    switch (c){
    case 'i':
      filein = optarg;
      break;
    case 'j':
      iter = atoi(optarg);
      break;
    case 'a':
      alpha = atof(optarg);
      break;
    case 't':
      HISTORY = atoi(optarg);
      break;
    case 'K':
      if(avgNK == 0.0) {
	K = atoi(optarg);
      }else{
	printf("can't use -k and -K at same time.\n");
	exit(1);
      }
      break;
    case 'k':
      if(K == 0) {
	avgNK = atof(optarg);
      }else{
	printf("can't use -k and -K at same time.\n");
	exit(1);
      }
      break;
    case 'R':
      randseed = atoi(optarg);
      break;
    default :
      fprintf(stdout, "%s [options]\n"
	      "\n"
	      "  files\n"
	      "\t-i <input data file>\n"
	      "  parameters\n"
	      "\t-j <sampling iterations>\n"
	      "\t-a <alpha>\n"
	      "\t-K <startvalue for number of clusters> not compat. with -k\n"
	      "\t-k <avg. number of reads in each startcluster> not compat. with -K\n"
	      "\t-t <history time>\n"
	      "\t-R <randomseed>\n"
	      "-----------------------------------------------------\n"
	      "\t-h\t\t this help!\n", argv[0]);
      exit(-1);
    }
    
  }
  if(HISTORY > iter)HISTORY=iter;
  if(randseed == 0) randseed = time(NULL);
  if( (K==0) && (avgNK <= 0.0))avgNK = default_avgNK;
  return 0;
} // ends parsecommandline
