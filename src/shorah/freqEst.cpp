/*  freqEst.cpp - runs the EM algorithm to estimate haplotype frequencies
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

#include <unistd.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

struct Read
{
    unsigned int pos;
    unsigned int len;
    std::string seq;
};

// needed to sort, by Osvaldo
typedef struct comp_class
{
    std::string* hap;
    double freq;
} com_pair;

bool comp_func(com_pair* pi, com_pair* pj) { return (pi->freq > pj->freq); }

//#define SEEPROB

int CORRECTPARAM = 1;       // which definition of pr(r | h) do we use
int SEED = 1;               // random seed
int ITERMAX = 5000;         // maximum iterations in EM
double ACCURACY = .000001;  // convergence for EM
double KILLP = 0;           // kill all probabilities below this
int RUNS = 5;
int N, M;
int unmatchingReads = 0;

static inline int readMatchesSomeHaplotype(Read& r, std::vector<std::string>& haplotypes)
{
    for (std::vector<std::string>::iterator g = haplotypes.begin(); g != haplotypes.end(); ++g)
        if (r.seq == (*g).substr(r.pos, r.len)) return 1;
    return 0;
}

void getReads(std::istream& infile, std::vector<Read>& ans, std::vector<std::string>& haplotypes)
{
    /*
     * get the reads from a file/stdin
     * FIXME: this should be much more robust
     */

    Read r;
    std::string t;
    std::string tmp;
    std::string throwaway;

    // line of file is "pos seq"
    while (!std::getline(infile, tmp, ' ').eof()) {
        if (tmp[0] == '#') {
            std::getline(infile, throwaway);
            std::cerr << "Comment : " << tmp << " " << throwaway << std::endl;
            continue;
        }
        r.pos = atoi(tmp.c_str());
        if ((r.pos < 0)) {
            std::cerr << "Error: " << r.pos << " is an invalid start position\n";
            exit(1);
        }

        std::getline(infile, r.seq);
        r.len = r.seq.length();
        // std::cout << r.pos << " " << r.seq << std::endl;

        if (readMatchesSomeHaplotype(r, haplotypes)) {
            ans.push_back(r);
        } else {
            unmatchingReads++;
        }
    }
}

void getHaplotypes(std::istream& in, std::vector<std::string>& genotypesFinal)
{
    /* read in the list of candidate haplotypes */
    std::string tmp;
    genotypesFinal.clear();
    std::vector<std::string> genotypes;
    // int i = 0;

    // line of file is "pos seq"
    while (!std::getline(in, tmp).eof()) {
        if (tmp[0] == '#') {
            std::cerr << "Comment : " << tmp << std::endl;
            continue;
        }
        genotypes.push_back(tmp);
        // std::cout << i++ << " =" << tmp << ".\n";
    }
    if (genotypes.size() == 0) {
        std::cout << "Error, no haplotypes\n";
        exit(1);
    }
    // now genotypes looks like
    // PQITSLX
    // .QITSSX
    // PQITSS.
    // So we build a consensus haplotype and fill in the .'s with the consensus
    int seqLen = genotypes[0].size();
    char consChar = ' ';
    std::string cons;
    for (int pos = 0; pos < seqLen; ++pos) {
        std::map<char, int> col;
        for (std::vector<std::string>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
            col[(*g)[pos]]++;
        }
        int max = -1;
        // don't want to fill in with a '.'
        col['.'] = -100;
        for (std::map<char, int>::iterator base = col.begin(); base != col.end(); ++base) {
            if (base->second > max) {
                consChar = base->first;
                max = base->second;
            }
        }
        cons = cons + consChar;
    }
    std::cout << "consensus haplotype\n" << cons << std::endl;

    std::vector<int> correctionCounts(genotypes.size(), 0);
    int tmpCount;
    for (std::vector<std::string>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
        tmpCount = 0;
        for (int pos = 0; pos < seqLen; ++pos) {
            if ((*g)[pos] == '.') {
                tmpCount++;
                (*g)[pos] = cons[pos];
#ifndef NDEBUG
                std::cout << "correcting . to " << cons[pos] << " in sequence " << *g << std::endl;
#endif
            }
        }
        correctionCounts.push_back(tmpCount);
    }

    // for(std::vector<int>::iterator a = correctionCounts.begin(); a != correctionCounts.end();
    // ++a) {
    //	std::cout << *a << " ";
    //}
    // std::cout << std::endl;

    // but the genotypes might not be unique right now
    std::map<std::string, int> uniq;
    for (std::vector<std::string>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
        uniq[*g]++;
    }
    for (std::map<std::string, int>::iterator x = uniq.begin(); x != uniq.end(); ++x) {
#ifndef NDEBUG
        std::cout << "Haplotype was repeated " << x->second << " times\n";
#endif
        genotypesFinal.push_back(x->first);
    }
}

inline int match(std::string geno, Read r) { return (r.seq == geno.substr(r.pos, r.len)); }

inline int readsEqual(Read& r, Read& s)
{
    if (r.pos == s.pos) {
        if (r.seq == s.seq) {
            return 1;
        }
    }
    return 0;
}

std::vector<int> countReads(std::vector<Read>& origReads, std::vector<Read>& reads)
{
    // take orig reads and do
    // %hash = (rd -> number appearing in origReads
    // u = values(hash)
    // reads = keys(hash)
    //
    int tmp;
    for (std::vector<Read>::iterator r = origReads.begin(); r != origReads.end(); ++r) {
        tmp = 1;
        for (std::vector<Read>::iterator s = origReads.begin(); s != r; ++s) {
            if (readsEqual(*r, *s)) {
                tmp = 0;
            }
        }
        if (tmp) {
            reads.push_back(*r);
        }
    }

    std::vector<int> u(reads.size(), 0);

    // for(std::vector<Read>::iterator r = reads.begin(); r != reads.end(); ++r) {
    for (unsigned int i = 0; i < reads.size(); ++i) {
        for (std::vector<Read>::iterator s = origReads.begin(); s != origReads.end(); ++s) {
            if (readsEqual(reads.at(i), *s)) {
                u[i]++;
            }
        }
    }
#ifndef NDEBUG
    std::cout << "Started with " << origReads.size() << " reads.\n";
    std::cout << reads.size() << " are unique\n";
    std::cout << "Counts:\n";
    for (unsigned int i = 0; i < u.size(); ++i) {
        std::cout << u[i] << " ";
    }
    std::cout << std::endl;
#endif
    return u;
}

void Estep(std::vector<double> const& p, std::vector<std::vector<double> >& U,
           std::vector<std::vector<double> > const& Z, std::vector<int> const& u)
{
    // given p, fills U with expected frequencies
    int i, j;
    double ProbY;

#ifndef NDEBUG
    std::cout << "Estep input: " << std::endl;
    for (i = 0; i < N; ++i) {
        std::cout << p[i] << " ";
    }
    std::cout << std::endl;
#endif

    for (i = 0; i < M; ++i) {
        ProbY = 0;
        for (j = 0; j < N; ++j) {
            ProbY += Z[i][j] * p[j];
            //		std::cout << "ProbY = " << ProbY << std::endl;
        }
        for (j = 0; j < N; ++j) {
            //		std::cout << "U[i][j] = " << u[i] << " " << Z[i][j] << " " << p[j] << " " <<
            // ProbY
            //<<
            // std::endl;
            U[i][j] = u[i] * ((Z[i][j] * p[j]) / ProbY);
        }
    }
    //	std::cout << "Estep output: " << std::endl;
    //	for (i = 0; i < M; ++i) {
    //		for (j = 0; j < N; ++j) {
    //			std::cout << U[i][j] << " ";
    //		}
    //		std::cout << std::endl;
    //	}
}

void Mstep(std::vector<double>& p, std::vector<std::vector<double> > const& U)
{
    std::vector<double> v(N, 0);
    double m = 0;
    int i, j;

    for (j = 0; j < N; ++j) {
        // std::cout << "." <<  v[j] << ".\n";
        for (i = 0; i < M; ++i) {
            //	std::cout << U[i][j] << " \n";
            v[j] += U[i][j];
        }
        m += v[j];
    }

    for (j = 0; j < N; ++j) {
        p[j] = v[j] / m;
    }

#ifndef NDEBUG
    for (j = 0; j < N; ++j) {
        std::cout << p[j] << " ";
    }
    std::cout << std::endl;
#endif
}

double logLike(std::vector<double>& p, std::vector<std::vector<double> > const& Z,
               std::vector<int> const& u)
{
    int i, j;

    double ell = 0;
    double Prob_Y;
    for (i = 0; i < M; i++) {
        Prob_Y = 0;
        for (j = 0; j < N; j++) {
            Prob_Y += Z[i][j] * p[j];
        }
        if (Prob_Y > 0) {
            ell += (u[i] * log(Prob_Y));
        }
    }
    return ell;
}

void round(std::vector<double>& p)
{
    for (std::vector<double>::iterator i = p.begin(); i != p.end(); ++i) {
        if ((*i) < KILLP) *i = 0;
    }
}

double EM(std::vector<double>& newP, std::vector<std::vector<double> > const& Z,
          std::vector<int> const& u)
{
    double sum = 0;
    double newEll = 0;
    std::vector<double> p(N, 0);
    std::vector<std::vector<double> > U(M, std::vector<double>(N, 0));
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

#ifndef NDEBUG
    for (j = 0; j < N; ++j) {
        std::cout << p[j] << " ";
    }
    std::cout << std::endl;
#endif

    while (((iter <= 2) || (std::abs(ell - newEll) > ACCURACY)) && (iter < ITERMAX)) {
        if (iter > 0) {
            round(newP);
            p = newP;
            ell = newEll;
        }

        Estep(p, U, Z, u);  //  fills U
        Mstep(newP, U);     // fills p
        newEll = logLike(newP, Z, u);

        // Print out some stats
        if (iter % 100 == 50 || iter == 1) {

#ifdef SEEPROB
            printf("%4d\t%f\t", iter, newEll);

            for (j = 0; j < N; ++j) {
                if (p[j] > 0.05) {
                    std::cout << j << ":" << p[j] << " ";
                }
            }
            std::cout << std::endl;
#else
            printf("%4d\t%f\n", iter, newEll);
#endif
        }
        // printf("%.3f %.3f %.3f ", newP[0], newP[1], newP[2]);
        // printf("%.3f %.3f %.3f ", newP[3], newP[4], newP[5]);
        // printf("%.3f %.3f %.3f\n", newP[6], newP[7], newP[8]);
        iter++;
    }
    return newEll;
}
void help(void)
{
    std::cout << "Usage: freqEst -f basename [-p precision -i maxiter -r runs -h -k kill -?]\n";
    std::cout << "Expects basename.read and basename.geno\n";
    std::cout << "Outputs to basename.popl\n";
    exit(1);
}

int main(int argc, char** argv)
{
    int c;
    std::vector<std::string> haplotypes;
    std::vector<Read> origReads, reads;
    std::string basename;
    std::string filename;
    std::string genotypeFilename, outfilename;
    int i, j;

    /* parse options
     * -f file (expect file.read file.geno)
     * -h                           (help)
     */
    while (1) {
        c = getopt(argc, argv, "f:h?p:i:r:k:s:");
        if (c == -1) break;
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

    std::ifstream infile(filename.c_str());
    std::ifstream ginfile(genotypeFilename.c_str());
    if (ginfile.fail() || infile.fail()) {
        std::cout << "ERROR -- can't read " << basename << ".geno or .read\n";
        exit(1);
    }
    std::ofstream outfile(outfilename.c_str());
    // outfile.open (outfilename.c_str());

    getHaplotypes(ginfile, haplotypes);

    getReads(infile, origReads, haplotypes);

    // uniq reads
    std::vector<int> u;
    u = countReads(origReads, reads);

    M = reads.size();
    N = haplotypes.size();

    std::cout << N << " haplotypes\n" << M << " reads\n";
    std::cout << "Threw out " << unmatchingReads << " unexplained reads\n";

    // count how many reads a haplotype is compat with
    std::vector<int> K(N, 0);

    // the M x N matrix of 0/1 incidences
    std::vector<std::vector<int> > A(M, std::vector<int>(N, 0));

    for (i = 0; i < M; ++i) {
        for (j = 0; j < N; ++j) {
            if (match(haplotypes[j], reads[i])) {
                A[i][j] = 1;
                K[j]++;
            }
        }
    }

    // my @Z;  # Z[i][j] == Prob(Y = y_i | X = x_j)
    std::vector<std::vector<double> > Z(M, std::vector<double>(N, 0));
    for (i = 0; i < M; ++i) {
        for (j = 0; j < N; ++j) {
            if (A[i][j] == 1) {
                if (CORRECTPARAM) {
                    Z[i][j] = 1.0 / 1000;
                } else {
                    Z[i][j] = 1.0 / K[j];
                }
                //				std::cout << "Z[i][j] = "<< Z[i][j] << std::endl;
            }
        }
    }

    std::vector<double> prob(N, 0);
    std::vector<double> bestProb(N, 0);
    double logL;
    double bestL = -1e100;

    for (int run = 0; run < RUNS; ++run) {
        std::cout << "run " << run << std::endl;
        logL = EM(prob, Z, u);
        if (logL > bestL) {
            bestProb = prob;
            bestL = logL;
        }
    }

    // sort the haplotypes by frequency before printing
    // added by Osvaldo
    std::string s6(10, 'x');
    com_pair** freq_hap;
    freq_hap = (com_pair**)calloc(N, sizeof(com_pair*));
    for (int i = 0; i < N; ++i) {
        freq_hap[i] = (com_pair*)calloc(1, sizeof(com_pair));
        freq_hap[i]->hap = &(haplotypes[i]);
        freq_hap[i]->freq = bestProb[i];
    }

    std::sort(freq_hap, freq_hap + N, comp_func);

    for (int i = 0; i < N; ++i) {
        if (freq_hap[i]->freq > 0) {
            // outfile << ">HAP" << i << "_" << bestProb[i] <<  "\n" << haplotypes[i] << std::endl;
            outfile << ">HAP" << i << "_" << freq_hap[i]->freq << "\n"
                    << *freq_hap[i]->hap << std::endl;
            std::cout << i << "\t" << freq_hap[i]->freq << "\t" << *freq_hap[i]->hap << std::endl;
            // std::cout << i << " " << freq_hap[i]->freq << std::endl;
        }
    }
    // The output is now sorted
    return 0;
}
