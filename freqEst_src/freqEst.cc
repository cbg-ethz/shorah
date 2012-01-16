/*  freqEst.cc - runs the EM algorithm to estimate haplotype frequencies
 *  given reads and reconstructed haplotypes
 
 # Copyright 2007, 2008, 2009
 # Niko Beerenwinkel,
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

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>

using namespace std;

struct Read {
	unsigned int pos;
	unsigned int len;
	string seq;
};

// needed to sort, by Osvaldo
typedef struct comp_class {
  string* hap;
  double freq;
} com_pair;

bool comp_func (com_pair* pi, com_pair* pj) { return (pi->freq > pj->freq); }

//#define DEBUG
//#define SEEPROB

int CORRECTPARAM = 1;  // which definition of pr(r | h) do we use
int SEED=1; // random seed
int ITERMAX = 5000; // maximum iterations in EM
double ACCURACY = .000001; // convergence for EM
double KILLP = 0; // kill all probabilities below this
int RUNS= 5;
int N, M;
int unmatchingReads = 0;

inline int readMatchesSomeHaplotype (Read & r, vector<string> & haplotypes) {
	for (vector<string>::iterator  g = haplotypes.begin(); g != haplotypes.end(); ++g)
		if (r.seq == (*g).substr(r.pos,r.len)) 
			return 1;
	return 0;
}


void getReads (istream & infile, vector<Read>& ans, vector<string> & haplotypes) {
	/*
	 * get the reads from a file/stdin
	 * FIXME: this should be much more robust
	 */

	Read r;
	string t;
	string tmp;
	string throwaway;

	// line of file is "pos seq"
	while(!getline(infile,tmp,' ').eof()) {
		if (tmp[0] == '#') {
			getline(infile,throwaway);
			cerr << "Comment : " << tmp << " " << throwaway << endl;
			continue;       
		}
		r.pos = atoi(tmp.c_str());
		if ((r.pos < 0)) { 
			cerr << "Error: " << r.pos << " is an invalid start position\n";
			exit(1);
		}

		getline(infile,r.seq);
		r.len = r.seq.length();
		//cout << r.pos << " " << r.seq << endl;

		if (readMatchesSomeHaplotype(r, haplotypes)) {
			ans.push_back(r);
		} else {
			unmatchingReads++;
		}
	}
}


void getHaplotypes (istream & in, vector<string> & genotypesFinal) {
	/* read in the list of candidate haplotypes */
	string tmp;
	genotypesFinal.clear();
	vector<string> genotypes;
	//int i = 0;

	// line of file is "pos seq"
	while(!getline(in,tmp).eof()) {
		if (tmp[0] == '#') {
			cerr << "Comment : " << tmp << endl;
			continue;       
		}
		genotypes.push_back(tmp);
		//cout << i++ << " =" << tmp << ".\n";
	}
	if (genotypes.size() == 0) {
		cout << "Error, no haplotypes\n";
		exit(1);
	}
	// now genotypes looks like
	// PQITSLX
	// .QITSSX
	// PQITSS.
	// So we build a consensus haplotype and fill in the .'s with the consensus
	int seqLen = genotypes[0].size();
	char consChar = ' ';
	string cons;
	for (int pos = 0; pos < seqLen; ++pos) {
		map<char, int> col;
		for(vector<string>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
			col[(*g)[pos]]++;
		}
		int max = -1;
		// don't want to fill in with a '.'
		col['.'] = -100;
		for(map<char, int>::iterator base = col.begin(); base != col.end(); ++base) {
			if (base->second > max) {
				consChar = base->first;
				max = base->second;
			}
		}
		cons = cons + consChar;
	}
	cout << "consensus haplotype\n" << cons << endl;

	vector<int> correctionCounts (genotypes.size(), 0);
	int tmpCount;
	for(vector<string>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
		tmpCount = 0;
		for (int pos = 0; pos < seqLen; ++pos) {
			if ((*g)[pos] == '.') {
				tmpCount++;
				(*g)[pos] = cons[pos];
#ifdef DEBUG
				cout << "correcting . to " << cons[pos] << " in sequence " << *g << endl;
#endif
			}
		}
		correctionCounts.push_back(tmpCount);
	}
	
	//for(vector<int>::iterator a = correctionCounts.begin(); a != correctionCounts.end(); ++a) {
	//	cout << *a << " ";
	//}
	//cout << endl;

	// but the genotypes might not be unique right now
	map<string, int> uniq;
	for(vector<string>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
		uniq[*g]++;
	}
	for(map<string, int>::iterator x = uniq.begin(); x != uniq.end(); ++x) {
#ifdef DEBUG
		cout << "Haplotype was repeated " <<  x->second << " times\n";
#endif
		genotypesFinal.push_back(x->first);
	}
	
}

inline int match (string geno, Read r) {
	return (r.seq == geno.substr(r.pos,r.len));
}

inline int readsEqual(Read & r, Read & s) {
	if (r.pos == s.pos) {
		if (r.seq == s.seq) {
			return 1;
		}
	}
	return 0;
}

vector<int> countReads (vector<Read> & origReads, vector<Read> & reads) {
	// take orig reads and do
	// %hash = (rd -> number appearing in origReads
	// u = values(hash)
	// reads = keys(hash)
	//
	int tmp;
	for (vector<Read>::iterator r = origReads.begin(); r != origReads.end(); ++r) {
		tmp = 1;
		for (vector<Read>::iterator s = origReads.begin(); s != r; ++s) {
			if (readsEqual(*r, *s)) {
				tmp = 0;
			}
		}
		if (tmp) {
			reads.push_back(*r);
		}
	}

	vector<int> u(reads.size(), 0);

	//for(vector<Read>::iterator r = reads.begin(); r != reads.end(); ++r) {
	for (unsigned int i = 0 ; i < reads.size(); ++i) {
		for (vector<Read>::iterator s = origReads.begin(); s != origReads.end(); ++s) {
			if (readsEqual(reads.at(i), *s)) {
				u[i]++;
			}
		}
	}
#ifdef DEBUG
	cout << "Started with " << origReads.size() << " reads.\n";
	cout << reads.size() << " are unique\n";
	cout << "Counts:\n";
	for (unsigned int i = 0; i < u.size(); ++i) {
		cout << u[i] << " ";
	} 
	cout << endl;
#endif
	return u;
}

void Estep (vector<double> const & p, vector<vector<double> >& U, vector<vector<double> > const & Z, vector<int> const & u) {
	// given p, fills U with expected frequencies
	int i,j;
	double ProbY;

#ifdef DEBUG
	cout << "Estep input: " << endl;
	for (i = 0; i < N; ++i) {
		cout << p[i] << " ";
	}
	cout << endl;
#endif

	for (i = 0; i < M; ++i) {
		ProbY = 0;
		for (j = 0; j < N; ++j) {
			ProbY += Z[i][j] * p[j];
	//		cout << "ProbY = " << ProbY << endl;
		}
		for (j = 0; j < N; ++j) {
	//		cout << "U[i][j] = " << u[i] << " " << Z[i][j] << " " << p[j] << " " << ProbY << endl;
			U[i][j] = u[i]* ((Z[i][j] * p[j]) / ProbY);
		}
	}
//	cout << "Estep output: " << endl;
//	for (i = 0; i < M; ++i) {
//		for (j = 0; j < N; ++j) {
//			cout << U[i][j] << " ";
//		} 
//		cout << endl;
//	}

}

void Mstep (vector<double> & p, vector<vector<double> > const & U) {
	vector<double> v(N,0);
	double m = 0;
	int i,j;

	for (j = 0; j < N; ++j) {
		//cout << "." <<  v[j] << ".\n";
		for (i = 0; i < M; ++i) {
		//	cout << U[i][j] << " \n";
			v[j] += U[i][j];
		}
		m += v[j];
	}

	for (j = 0; j < N; ++j) {
		p[j] = v[j] / m;
	}

#ifdef DEBUG
	for (j = 0; j < N; ++j) {
		cout << p[j] << " ";
	}
	cout << endl;
#endif

}

double logLike (vector<double> & p,vector<vector<double> > const & Z, vector<int> const & u) {
	int i,j;

	double ell = 0;
	double Prob_Y;
	for (i= 0; i < M; i++) {
		Prob_Y = 0;
		for (j= 0; j < N; j++) {
			Prob_Y += Z[i][j] * p[j];
		}
		if (Prob_Y > 0) {
			ell += (u[i] * log(Prob_Y));
		}
	}
	return ell;
}

void round(vector<double> & p) {
	for (vector<double>::iterator i = p.begin(); i != p.end(); ++i) {
		if ((*i) < KILLP) 
			*i = 0;
	}
}

double EM (vector<double> & newP, vector<vector<double> > const & Z, vector<int> const & u) {
	double sum = 0;
	double newEll = 0;
	vector<double> p(N,0);
	vector<vector<double> > U(M,vector<double>(N,0));
	double ell = 0; 
	int iter = 0;
	int j;

	for (j = 0; j < N; ++j) {
		p[j] = rand();
		sum += p[j];
	}
	for (j = 0; j < N; ++j) {
		p[j] = p[j] / sum;
	}

#ifdef DEBUG
	for (j = 0; j < N; ++j) {
		cout << p[j] << " ";
	}
	cout << endl;
#endif

	while (((iter <= 2) || (abs(ell - newEll) > ACCURACY)) && (iter < ITERMAX)) {
		if (iter > 0) {
			round(newP);
			p = newP;
			ell = newEll;
		}

		Estep(p, U, Z, u); //  fills U
		Mstep(newP,U); // fills p
		newEll = logLike(newP,Z,u);


		// Print out some stats
		if (iter%100 == 50 || iter == 1) {

#ifdef SEEPROB
			printf("%4d\t%f\t", iter,newEll);

			for (j = 0; j < N; ++j) {
				if (p[j] > 0.05) {
					cout << j << ":" << p[j] << " ";
				}
			}
			cout << endl;
#else
			printf("%4d\t%f\n", iter,newEll);
#endif


		}
		//printf("%.3f %.3f %.3f ", newP[0], newP[1], newP[2]);
		//printf("%.3f %.3f %.3f ", newP[3], newP[4], newP[5]);
		//printf("%.3f %.3f %.3f\n", newP[6], newP[7], newP[8]);
		iter++;
	}
	return newEll;
}
void help (void) {
	cout << "Usage: freqEst -f basename [-p precision -i maxiter -r runs -h -k kill -?]\n";
	cout << "Expects basename.read and basename.geno\n";
	cout << "Outputs to basename.popl\n";
	exit(1);
}

int main (int argc, char **argv) {
	int c;
	vector<string> haplotypes;
	vector<Read> origReads, reads;
	string basename;
	string filename;
	string genotypeFilename, outfilename;
	int i,j;

	/* parse options
	 * -f file (expect file.read file.geno)
	 * -h                           (help)
	 */
	while (1) {
		c = getopt (argc, argv, "f:h?p:i:r:k:s:");
		if (c == -1)
			break;
		switch (c) {
			case 'k':
				KILLP = atof(optarg);
				break;
			case 's':
				SEED = atoi(optarg);
				break;
			case 'f':
				basename = optarg;
				break;
			case 'p':
				ACCURACY = atof(optarg);
				break;
			case 'i':
				ITERMAX = atoi(optarg);
				break;
			case 'r':
				RUNS = atoi(optarg);
				break;
			case '?':
				help();
				break;
			case 'h':
				help();
				break;
			default:
				help();
		}
	}
	if (basename == "") {
		help();
	}
	filename = basename + ".read";
	genotypeFilename = basename + ".geno";
	outfilename = basename + ".popl";

	srand(SEED);

	ifstream infile (filename.c_str());
	ifstream ginfile(genotypeFilename.c_str());
	if (ginfile.fail() || infile.fail()) {
		cout << "ERROR -- can't read " << basename << ".geno or .read\n";
		exit(1);
	}
	ofstream outfile(outfilename.c_str());
	//outfile.open (outfilename.c_str());

	getHaplotypes(ginfile,haplotypes);

	getReads(infile,origReads, haplotypes);


	// uniq reads
	vector<int> u;
	u = countReads(origReads, reads);

	M = reads.size();
	N = haplotypes.size();

	cout << N << " haplotypes\n" << M << " reads\n";
	cout << "Threw out " << unmatchingReads << " unexplained reads\n";

	// count how many reads a haplotype is compat with
	vector<int> K(N,0);

	// the M x N matrix of 0/1 incidences
	vector<vector<int> > A(M,vector<int>(N,0));

	for (i = 0; i < M; ++i) {
		for (j = 0; j < N; ++j) {
			if (match(haplotypes[j],reads[i])) {
				A[i][j] = 1;
				K[j]++;
			}
		}
	}

	//my @Z;  # Z[i][j] == Prob(Y = y_i | X = x_j)
	vector<vector<double> > Z (M,vector<double>(N,0));
	for (i = 0; i < M; ++i) {
		for (j = 0; j < N; ++j) {
			if (A[i][j] == 1) {
				if (CORRECTPARAM) {
					Z[i][j] = 1.0 / 1000;
				} else {
					Z[i][j] = 1.0 / K[j];
				}
				//				cout << "Z[i][j] = "<< Z[i][j] << endl;
			}
		}
	}

	vector<double> prob(N,0);
	vector<double> bestProb(N,0);
	double logL;
	double bestL = -1e100;

	for (int run = 0; run < RUNS; ++run) {
		cout << "run " << run << endl;
		logL = EM(prob,Z,u);
		if (logL > bestL) {
			bestProb = prob;
			bestL = logL;
		}
	}
	
	// sort the haplotypes by frequency before printing
	// added by Osvaldo
	string s6 (10, 'x');
	com_pair** freq_hap;
	freq_hap = (com_pair**) calloc(N, sizeof(com_pair*));
	for (int i = 0; i < N; ++i) {
	  freq_hap[i] = (com_pair*)calloc(1, sizeof(com_pair));
	  freq_hap[i]->hap = &(haplotypes[i]);
	  freq_hap[i]->freq = bestProb[i];
	}

	sort(freq_hap, freq_hap+N, comp_func);
	
	for (int i = 0; i < N; ++i) {
	  if (freq_hap[i]->freq > 0)  {
	    //outfile << ">HAP" << i << "_" << bestProb[i] <<  "\n" << haplotypes[i] << endl;
	    outfile << ">HAP" << i << "_" <<  freq_hap[i]->freq <<  "\n" << *freq_hap[i]->hap << endl;
	    cout << i << "\t" << freq_hap[i]->freq <<  "\t" << *freq_hap[i]->hap << endl;
	    //cout << i << " " << freq_hap[i]->freq << endl;
	  }
	}
	// The output is now sorted
	return 0;
}
