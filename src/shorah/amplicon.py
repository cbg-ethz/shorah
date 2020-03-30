#!/usr/bin/env python3

# Copyright 2007-2018
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
"""amplian.py is the program that performs the analysis in amplicon mode.
   It creates a MSA of the reads and performs error correction with a
   single run of diri_sampler. Then, it performs SNV discovery by
   calling the program snv.py.
   """
from __future__ import division
from __future__ import print_function
import os
import os.path
import sys
import shlex
import re
import shutil

import logging
import logging.handlers
from pkg_resources import resource_filename

from Bio import SeqIO

dn_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if __name__ == '__main__':
    if __package__ is None:
        os.sys.path.insert(1, dn_dir)
        mod = __import__('shorah')
        sys.modules["shorah"] = mod
        import shorah_snv
else:
    from . import shorah_snv

diri_exe = resource_filename(__name__, 'bin/diri_sampler')
b2w_exe = resource_filename(__name__, 'bin/b2w')

if not (os.path.exists(diri_exe) and os.path.exists(b2w_exe)):
    diri_exe = shutil.which('diri_sampler')
    b2w_exe = shutil.which('b2w')
    if not (diri_exe and b2w_exe):
        logging.error(
            'Executables b2w and diri_sampler not found, compile first.')
        sys.exit('Executables b2w and diri_sampler not found, compile first.')

win_min_ext = 0.95
cutoff_depth = 10000


def run_child(exe_name, arg_string):
    """use subrocess to run an external program with arguments"""
    import subprocess

    if not arg_string.startswith(' '):
        arg_string = ' ' + arg_string

    logging.debug(exe_name + arg_string)

    try:
        retcode = subprocess.call(exe_name + arg_string, shell=True)
        if retcode > 0:
            logging.error(exe_name + arg_string)
            logging.error("Child %s terminated by signal %s",
                          exe_name, retcode)
        else:
            logging.debug("Child %s returned %i", exe_name, retcode)
    except OSError as ee:
        logging.error("Execution of %s failed: %s", exe_name, ee)

    return retcode


def run_diagnostics(window_file, reads):
    """Performs some basic diagnostics on the quality of the MC sampling
    """
    import warnings

    # greedy re match to handle situation where '.reads' appears in the ID
    stem = re.match(r'^(?P<stem>.*).reads', window_file).group('stem')
    smp_file = stem + '.smp'
    dbg_file = stem + '.dbg'
    del stem

    with open(dbg_file) as l:
        lines = l.readlines()
    for line in lines:
        if line.startswith('# q ='):
            q = int(line.strip().split('=')[1])
        if line.startswith('#made'):
            new_clusters = int(line.split()[1])

    if new_clusters < q:
        warnings.warn('clusters created: %d, q: %d maybe alpha is too low'
                      % (new_clusters, q))

    clusters = []
    untouched = []
    theta = []
    gamma = []
    with open(smp_file, 'rb') as smp_reader:
        lines = smp_reader.readlines()
        # fields = lines[0].split()
        for line in lines[1:]:
            it, cl, unt = list(map(int, line.split()[:3]))
            the, gam = list(map(float, line.split()[3:]))
            clusters.append(cl)
            untouched.append(unt)
            theta.append(the)
            gamma.append(gam)
            del it

    logging.info('sample has %d reads', reads)
    untouched_hst = untouched[-2000:]
    unt_mean = sum(untouched_hst) / len(untouched_hst)
    unt_ratio = 100 * unt_mean / q
    unt_msg = '%3.1f %% of untouched objects <should be around 90-95%%>' % \
        unt_ratio
    logging.info(unt_msg)
    if unt_ratio < 90.0:
        warnings.warn(unt_msg)


def matchremove(matchobj):
    """Callback function used in mpileup manipulation"""
    match = matchobj.group(0)
    if match.startswith('+') or match.startswith('-'):
        c = int(match[1:])
    else:
        c = 0
    del c
    return ''


def shannon_entropy(bases):
    """Shannon entropy of a mpileup column: Pseudocount = 1 and
    returns 0 if length = 1"""
    import math
    if len(bases) == 1:
        return 0.0
    letters = ['A', 'C', 'G', 'T']
    counts = [float(bases.count(l)) if bases.count(l) else 1.0
              for l in letters]
    sc = sum(counts)
    return -sum([(c / sc) * math.log(c / sc) / math.log(2.0) for c in counts])


def plot_entropy(pos_ent, pos_coords, ave_ent, win_coords):
    """Plot entropies and window to a pdf file with matplotlib"""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logging.error('could not import matplotlib, no pdf produced')
        return
    high_start, high_stop = win_coords

    ent_start, ent_stop = pos_coords
    X = list(range(ent_start, ent_stop))

    fig, ax1 = plt.subplots()
    # ax1.plot(X, pos_ent, '-', color='#E69F00', alpha=0.8)
    ax1.vlines(X, 0, pos_ent[ent_start:ent_stop], color='#E69F00', alpha=0.8)
    ax1.set_xlabel('position on reference')
    # Make the y-axis label and tick labels match the line color.
    ax1.set_ylabel('entropy per position', color='#E69F00')
    for tl in ax1.get_yticklabels():
        tl.set_color('#E69F00')

    ax2 = ax1.twinx()
    ax2.plot(X, ave_ent[ent_start:ent_stop], ':', lw=2., color='#56B4E9',
             alpha=1.0)
    ax2.set_ylabel('window average', color='#56B4E9')
    for tl in ax2.get_yticklabels():
        tl.set_color('#56B4E9')

    ax2.axvspan(high_start, high_stop, color='#CC79A7', alpha=0.4,
                label='highest entropy window')

    plt.title('chosen entropy window is %d-%d' % (high_start, high_stop))
    plt.savefig('entropy.pdf')
    del fig


def highest_entropy(bam_file, fasta_file, ent_sel='relative'):
    """Parse reads to have their length distribution and compute the
    trimmed mean read length"""

    import warnings

    read_len = []
    run_child('samtools', ' view %s | cut -f 10 > rl.txt' % bam_file)
    for l in open('rl.txt'):
        read_len.append(len(l.strip()))
    os.remove('rl.txt')
    read_len = sorted(read_len)
    n_reads = len(read_len)
    logging.info('n_reads: %d', n_reads)
    # max_len = max(read_len)
    trimmed_mean = sum([read_len[i] for i in range(int(0.1 * n_reads),
                                                   int(0.9 * n_reads))])
    trimmed_mean /= (0.8 * n_reads)
    trimmed_mean = int(round(trimmed_mean, 0))
    logging.info('trimmed_mean: %d', trimmed_mean)
    # Build the mpileup and compute the entropy per position
    ref_seq = list(SeqIO.parse(fasta_file, 'fasta'))[0]
    entropy = [None] * (len(ref_seq) + 1)
    run_child('samtools', 'mpileup  -f %s -d %d %s > sample.mpu' %
              (fasta_file, cutoff_depth, bam_file))
    for l in open('sample.mpu'):
        pos, refbase, depth = int(l.split()[1]), l.split()[2], \
            int(l.split()[3])
        if depth == 0:
            continue
        readbase = l.split()[4]
        # remove read start, insertions, deletions
        column = re.sub(r'\^.', matchremove, readbase)
        column = re.sub(r'-[0-9]+', matchremove, column)
        column = re.sub(r'\+[0-9]+', matchremove, column)
        column = column.replace('$', '')
        column = column.upper()
        # still not perfect control over mpileup format
        if abs(depth - len(column)) > 10:
            warnings.warn('mpileup column not fully parsed')
            # print(refbase, depth, len(column))
            # print(readbase)
            # print(column.replace(',', '').replace('.', ''))
            # sys.exit()
        column = column.replace(',', refbase).replace('.', refbase)
        entropy[pos] = shannon_entropy(column)

    # identifies the start and stop of the high entropy region
    start, stop = None, None
    for i, e in enumerate(entropy):
        if start is None and e is not None and e > 0.0:
            start = i
        if start and e is None:
            stop = i
            break
        stop = i
    logging.info('start: %d, stop: %d', start, stop)

    # mean entropy
    ent_mean = [None] * (len(ref_seq) + 1)
    delta = trimmed_mean // 2  # used to center the moving window
    for i in range(start, stop - trimmed_mean):
        ent_mean[i + delta] = sum(entropy[i:i + trimmed_mean]) / trimmed_mean

    # max entropy per position, excluding the first and last 10 positions
    max_ent_per_pos = -1.0
    for i in range(start + 10, stop - 10):
        if entropy[i] > max_ent_per_pos:
            max_ent_per_pos = entropy[i]
            highest_ent_pos = i
    logging.info('highest entropy found at position %d', highest_ent_pos)

    # the window is chosen as the absolute max mean_entropy or as the
    # max mean entropy covering the position with max entropy
    max_ent = -1.0
    high_ent_start = -1

    if ent_sel == 'absolute':
        rsta = start
        rsto = stop - trimmed_mean
    elif ent_sel == 'relative':
        rsta = max(start, highest_ent_pos - trimmed_mean + 1)
        rsto = min(stop - trimmed_mean, highest_ent_pos + 1)

    for i in range(rsta, rsto):
        ent_mean[i + delta] = sum(entropy[i:i + trimmed_mean]) / trimmed_mean
        if ent_mean[i + delta] >= max_ent:
            max_ent = ent_mean[i + delta]
            high_ent_start = i
    high_ent_stop = high_ent_start + trimmed_mean

    # print entropy file
    eh = open('entropy.csv', 'w')
    eh.write('pos,entropy,mean_entropy,high_entropy_window\n')
    for i in range(start, stop):
        se = str(round(entropy[i], 4))
        try:
            sm = str(round(ent_mean[i], 4))
        except IndexError:
            sm = 'NA'
        except TypeError:
            sm = 'NA'
        if i >= high_ent_start and i <= high_ent_stop:
            sh = '1'
        else:
            sh = '0'
        ltw = '%d,%s,%s,%s\n' % (i, se, sm, sh)
        eh.write(ltw)
    eh.close()

    # plot entropy; requires matplotlib
    plot_entropy(entropy, (start, stop), ent_mean,
                 (high_ent_start, high_ent_stop))

    return high_ent_start, high_ent_stop


# def main(in_bam, in_fasta, min_overlap=0.95, max_coverage=50000,
#         alpha=0.5, s=0.01, region='', diversity=False):
def main(args):
    """
    Performs the amplicon analysis, running diri_sampler
    and analyzing the result
    """
    in_bam = args.b
    in_fasta = args.f
    region = args.r
    max_coverage = args.max_coverage
    alpha = args.a
    seed = args.seed
    sigma = args.sigma
    cov_thrd = args.cov_thrd
    diversity = args.diversity
    min_overlap = args.min_overlap
    ignore_indels = args.ignore_indels

    logging.info(' '.join(sys.argv))
    # info on reference and region if given, or discover high entropy one
    ref_seq = list(SeqIO.parse(in_fasta, 'fasta'))[0]
    ref_name = ref_seq.id
    if region:
        # handles situation where ':' or '-' appears in the ID
        reg_bound = re.search(r':(?P<start>\d+)-(?P<stop>\d+)$', region)
        reg_start, reg_stop = int(reg_bound.group(
            'start')), int(reg_bound.group('stop'))
        ref_length = reg_stop - reg_start + 1
        del reg_bound
    elif region == '' and diversity:
        reg_start, reg_stop = highest_entropy(in_bam, in_fasta)
        ref_length = reg_stop - reg_start + 1
        region = '%s:%d-%d' % (ref_seq.id, reg_start, reg_stop)
    elif region == '' and not diversity:
        reg_start = 1
        ref_length = len(ref_seq)
        reg_stop = ref_length

    logging.info('analysing region from %d to %d', reg_start, reg_stop)

    # output the reads, aligned to the amplicon
    d = ' -d' if ignore_indels else ''

    b2w_args = ' -i 0 -w %d -m %d -x %d -c %i%s %s %s %s' % (ref_length, int(
        min_overlap * ref_length), max_coverage, cov_thrd, d, in_bam, in_fasta, region)
    ret_b2w = run_child(shlex.quote(b2w_exe), b2w_args)
    logging.debug('b2w returned %d', ret_b2w)

    # run diri_sampler on the aligned reads
    win_file = 'w-%s-%u-%u.reads.fas' % (ref_name, reg_start, reg_stop)
    # TODO clean ref_name of special caracters - that would be an alternative to processing everything with regex down the line
    # BUG the solution currently used by ShoRAH can still fail when path '/'
    # (or on windows '\\' and ':') characters are present in the ref_seq.id
    h = list(open('coverage.txt'))[0]
    n_reads = int(h.split()[-1])
    assert os.path.exists(win_file), 'window file %s not found' % win_file

    iterations = min(30000, n_reads * 20)
    diri_args = '-i %s -j %d -a %f -t 2000' % (win_file, iterations, alpha)
    ret_diri = run_child(shlex.quote(diri_exe), diri_args)
    logging.debug('diri_sampler returned %d', ret_diri)

    # diagnostics on the convergence
    run_diagnostics(win_file, n_reads)

    # run snv.py to parse single nucleotide variants
    args.increment = 1
    shorah_snv.main(args)
