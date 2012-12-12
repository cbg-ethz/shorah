#!/usr/bin/env python

# Copyright 2007-2012
# Niko Beerenwinkel,
# Nicholas Eriksson,
# Moritz Gerstung,
# Lukas Geyrhofer,
# Kerensa McElroy,
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
''' amplian.py is the program that performs the analysis in amplicon mode.
    It creates a MSA of the reads and performs error correction with a single
    run of diri_sampler. Then, it performs SNV discovery by calling the
    program snv.py.
    '''
import os
import os.path
import sys

import logging
import logging.handlers

# Make a global logging object.
amplog = logging.getLogger("AmplianLog")
# set logging level
amplog.setLevel(logging.DEBUG)
# This handler writes everything to a file.
LOG_FILENAME = './amplian.log'
h = logging.handlers.RotatingFileHandler(LOG_FILENAME, 'w',
                                         maxBytes=100000, backupCount=5)
f = logging.Formatter("%(levelname)s %(asctime)s %(funcName)s\
                      %(lineno)d %(message)s")
h.setFormatter(f)
amplog.addHandler(h)

win_min_ext = 0.95

dn = os.path.dirname(__file__)


def run_child(exe_name, arg_string):
    '''use subrocess to run an external program with arguments'''
    import subprocess

    if not arg_string.startswith(' '):
        arg_string = ' ' + arg_string

    amplog.debug(exe_name + arg_string)

    try:
        retcode = subprocess.call(exe_name + arg_string, shell=True)
        if retcode > 0:
            amplog.error(exe_name + arg_string)
            amplog.error("Child %s terminated by signal" % exe_name, retcode)
        else:
            amplog.debug("Child %s returned %i" % (exe_name, retcode))
    except OSError, ee:
        amplog.error("Execution of %s failed:" % exe_name, ee)

    return retcode


def run_diagnostics(window_file, reads):
    '''Performs some basic diagnostics on the quality of the MC sampling
    '''
    import warnings
    import csv

    smp_file = window_file.split('.')[0] + '.smp'
    dbg_file = window_file.split('.')[0] + '.dbg'

    with open(dbg_file) as l:
        lines = l.readlines()
    for line in lines:
        if line.startswith('# q ='):
            q = int(line.split('=')[1])
        if line.startswith('#made'):
            new_clusters = int(line.split()[1])

    if new_clusters < reads:
        warnings.warn('Few clusters created, maybe alpha is too low')

    clusters = []
    untouched = []
    theta = []
    gamma = []
    with open(smp_file, 'rb') as smp_reader:
        lines = smp_reader.readlines()
        fields = lines[0].split()
        for line in lines[1:]:
            it, cl, unt = map(int, line.split()[:3])
            the, gam = map(float, line.split()[3:])
            clusters.append(cl)
            untouched.append(unt)
            theta.append(the)
            gamma.append(gam)

    untouched_hst = untouched[-2000:]
    unt_mean = float(sum(untouched_hst)) / len(untouched_hst)
    #unt_msg = '%f %% of untouched objects <should be around 90-95%>' %
        
    print >> sys.stderr, unt_mean


def main(in_bam='', in_fasta='', min_overlap=0.95, max_coverage=50000,
         alpha=0.1):
    '''
    Performs the amplicon analysis, running diri_sampler
    and analyzing the result
    '''
    from Bio import SeqIO

    ref_seq = list(SeqIO.parse(in_fasta, 'fasta'))[0]
    ref_name = ref_seq.id
    ref_length = len(ref_seq)

    # output the reads, aligned to the amplicon
    b2w_exe = os.path.join(dn, 'b2w')
    b2w_args = ' -i 0 -w %d -m %d -x %d %s %s' % \
        (ref_length, int(min_overlap * ref_length),
         max_coverage, in_bam, in_fasta)
    ret_b2w = run_child(b2w_exe, b2w_args)

    # run diri_sampler on the aligned reads
    win_file = 'w-%s-1-%d.reads.fas' % (ref_name, ref_length)
    h = list(open('coverage.txt'))[0]
    n_reads = int(h.split()[-1])
    assert os.path.exists(win_file), 'window file not found'
    diri_exe = os.path.join(dn, 'diri_sampler')
    iterations = min(20000, n_reads * 10)
    diri_args = '-i %s -j %d -a 0.1 -t 2000' % (win_file, iterations)
    ret_diri = run_child(diri_exe, diri_args)

    # diagnostics on the convergence
    run_diagnostics(win_file, n_reads)


if __name__ == "__main__":

    import optparse
    # parse command line
    optparser = optparse.OptionParser()
    opts = main.func_defaults  # set the defaults (see http://bit.ly/2hCTQl)

    optparser.add_option("-b", "--bam",
                         help="file with aligned reads in .bam format",
                         default=opts[0], type="string", dest="in_bam")

    optparser.add_option("-f", "--fasta",
                         help="reference genome in fasta format",
                         default=opts[1], type="string", dest="in_fasta")

    optparser.add_option("-m", "--min_overlap",
                         help="fraction of read overlap to be included",
                         default=opts[2], type="float", dest="min_overlap")

    optparser.add_option("-x", "--maxcov",
                         help="approximate max coverage allowed <%default>",
                         default=opts[3], type="int", dest="max_coverage")

    optparser.add_option("-a", "--alpha",
                         help="alpha in dpm sampling <%default>",
                         default=opts[4], type="float", dest="alpha")

    (options, args) = optparser.parse_args()

    supported_formats = {
        'bam': 'aligned reads',
        'fasta': 'reference genome'
    }
    amplog.info(' '.join(sys.argv))
    # check the input file is in supported format
    try:
        tmp_filename = os.path.split(options.in_bam)[1]
        [in_stem, in_format] = [tmp_filename.split('.')[0],
                                tmp_filename.split('.')[-1]]
        t = supported_formats[in_format]
    except IndexError:
        amplog.error('The input file must be filestem.format')
        print 'The input file must be filestem.format'
        print 'Supported formats are'
        for sf in supported_formats.iteritems():
            print sf[0], ':', sf[1]
        sys.exit()
    except KeyError:
        amplog.error('format unknown')
        print in_format, 'format unknown'
        print 'Supported formats are'
        for sf in supported_formats.iteritems():
            print sf[0], ':', sf[1]
        sys.exit()
    main(*args, **vars(options))
