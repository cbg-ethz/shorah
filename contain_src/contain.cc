/* 
 * contain.cc -- eliminates redundant reads
 *
 * usage: contain -f basename
 * input: basename.read
 * output: basename.rest
 * both in the format "start_pos read_string"
 *
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

using namespace std;

struct Read {
	unsigned int pos;
	unsigned int len;
	unsigned int end;
	string seq;
};

void help (void) {
	cout << "Usage: contain -f basename\n";
	cout << "Expects basename.read, outputs basename.rest\n";
	exit(1);
}

void getReads (istream & infile, vector<Read>& ans) {
	/*
	 * get the reads from a file/stdin
	 * FIXME: this should be much more robust
	   FIXME: die if the file isn't readable!!!
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
		r.end = r.pos + r.len - 1;

		if (r.len > 0) 
			ans.push_back(r);
	}
}

int contain (vector<Read> & reads, int i, int j) {
	// return i if read[i] contained in read[j]
	//        j if read[j] ............ read[i]
	//        -1 else
	//
	int tmp, overlap;

	if ((reads[i].pos > reads[j].end) || (reads[j].pos > reads[i].end)) {
		//cout << "no overlap\n";
		return -1;
	}

	// switch so that i starts first
	if ((reads[i].pos > reads[j].pos) || 
			((reads[i].pos == reads[j].pos) && (reads[i].len < reads[j].len))) {
		tmp = j;
		j = i;
		i = tmp;
	}
	overlap = reads[i].pos + reads[i].len - reads[j].pos;
	//cout << "read " << i << " at pos " << reads[i].pos << " length " << reads[i].len << "\n";
	//cout << "read " << j << " at pos " << reads[j].pos << " length " << reads[j].len << "\n";
	//cout << "overlap == " << overlap << "\n";
	if (overlap < reads[j].len) {
	//	cout << "not enough overlap\n";
		return -1;
	}
	if (overlap > reads[j].len) {
		overlap = reads[j].len;
	}
	if (reads[i].seq.substr(reads[j].pos - reads[i].pos, overlap) == 
			reads[j].seq.substr(0, overlap)) {
	//	cout << "Agree\n";
		return j;
	} else {
	//	cout << "Disagree\n";
		return -1;
	}
}

int main (int argc, char **argv) {
	int c;
	vector<Read> reads;
	string basename;
	string filename;
	string outfilename;
	int i,j;

	/* parse options
	 * -f file (expect file.read file.geno)
	 * -h                           (help)
	 */
	while (1) {
		c = getopt (argc, argv, "f:h?");
		if (c == -1)
			break;
		switch (c) {
			case 'f':
				basename = optarg;
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
	outfilename = basename + ".rest";

	ifstream infile (filename.c_str());
	ofstream outfile(outfilename.c_str());
	map<int,int> delReads;

	// pull all reads into file
	getReads(infile, reads);
	int idx;
	int finalReads = 0;
	int numDeleted = 0;
	cout << "readfile read" << endl;
	// do double loop over all reads.
	// contain returns the index to delete: either i, j, or -1
	// should be able to speed this up by skipping i or j if
	// delReads[i] or j is 1
	// but have to be careful
	//
	for (i = 0; i < reads.size(); ++i) {
		// if reads[i] already deleted, then we will have already
		// deleted everything reads[i] makes redundant
		if (delReads[i]) {
			continue;
		}
		if (i % 10 == 0) 
			cout << "finished " << i << " reads, " << numDeleted << " deleted\n";

		for (j = i + 1; j < reads.size(); ++j) {
			if (delReads[j])
				continue;
			idx = contain(reads,i,j);

			if (idx == i) {
				delReads[idx] = 1;
				// skip to next iteration of i loop
				j = reads.size();
				numDeleted++;
			} else if (idx == j) {
				delReads[idx] = 1;
				numDeleted++;
			}
			/*
			if (idx != -1) {
				cout << "del " << idx << " (from " << i << " " << j  << ")\n";
			}
			*/
		}
	}
	for (i = 0; i < reads.size(); ++i) {
		if (delReads[i] != 1) {
			finalReads++;
			outfile << reads[i].pos << " " << reads[i].seq << "\n";
		}
	}
	cout << "Started with " << reads.size() << " reads.\n";
	cout << "Ended with " << finalReads << " reads.\n";
}

