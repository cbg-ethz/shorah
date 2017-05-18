/*
 * contain.cpp -- eliminates redundant reads
 *
 * usage: contain -f basename
 * input: basename.read
 * output: basename.rest
 * both in the format "start_pos read_string"
 *
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

#include <unistd.h>
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
    unsigned int end;
    std::string seq;
};

void help(void)
{
    std::cout << "Usage: contain -f basename\n";
    std::cout << "Expects basename.read, outputs basename.rest\n";
    exit(1);
}

void getReads(std::istream& infile, std::vector<Read>& ans)
{
    /*
     * get the reads from a file/stdin
     * FIXME: this should be much more robust
       FIXME: die if the file isn't readable!!!
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
        r.end = r.pos + r.len - 1;

        if (r.len > 0) ans.push_back(r);
    }
}

int contain(std::vector<Read>& reads, int i, int j)
{
    // return i if read[i] contained in read[j]
    //        j if read[j] ............ read[i]
    //        -1 else
    //
    int tmp, overlap;

    if ((reads[i].pos > reads[j].end) || (reads[j].pos > reads[i].end)) {
        // std::cout << "no overlap\n";
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
    // std::cout << "read " << i << " at pos " << reads[i].pos << " length " << reads[i].len <<
    // "\n";
    // std::cout << "read " << j << " at pos " << reads[j].pos << " length " << reads[j].len <<
    // "\n";
    // std::cout << "overlap == " << overlap << "\n";
    if (overlap < reads[j].len) {
        //	std::cout << "not enough overlap\n";
        return -1;
    }
    if (overlap > reads[j].len) {
        overlap = reads[j].len;
    }
    if (reads[i].seq.substr(reads[j].pos - reads[i].pos, overlap) ==
        reads[j].seq.substr(0, overlap)) {
        //	std::cout << "Agree\n";
        return j;
    } else {
        //	std::cout << "Disagree\n";
        return -1;
    }
}

int main(int argc, char** argv)
{
    int c;
    std::vector<Read> reads;
    std::string basename;
    std::string filename;
    std::string outfilename;
    int i, j;

    /* parse options
     * -f file (expect file.read file.geno)
     * -h                           (help)
     */
    while (1) {
        c = getopt(argc, argv, "f:h?");
        if (c == -1) break;
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

    std::ifstream infile(filename.c_str());
    std::ofstream outfile(outfilename.c_str());
    std::map<int, int> delReads;

    // pull all reads into file
    getReads(infile, reads);
    int idx;
    int finalReads = 0;
    int numDeleted = 0;
    std::cout << "readfile read" << std::endl;
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
        if (i % 10 == 0) std::cout << "finished " << i << " reads, " << numDeleted << " deleted\n";

        for (j = i + 1; j < reads.size(); ++j) {
            if (delReads[j]) continue;
            idx = contain(reads, i, j);

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
                std::cout << "del " << idx << " (from " << i << " " << j  << ")\n";
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
    std::cout << "Started with " << reads.size() << " reads.\n";
    std::cout << "Ended with " << finalReads << " reads.\n";
}
