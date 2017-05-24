# Copyright 2007-2012
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

'''dec.py draws windows on the bam file, cuts them into separate fasta files,
    calls diri_sampler to correct them, and merges the correction into a
    file of corrected reads
'''

from __future__ import division
from __future__ import print_function
import numpy as np
import os
import pipes
import sys

import logging
import logging.handlers

# Make a global logging object.
declog = logging.getLogger(__name__)

#################################################
# a common user should not edit above this line #
#################################################
# parameters not controlled by command line options
fasta_length = 80   # controls line length in fasta files
win_min_ext = 0.85  # if read covers at least win_min_ext fraction of
# the window, fill it with Ns
hist_fraction = 0.20  # fraction that goes into the history
min_quality = 0.9  # quality under which discard the correction
min_x_thresh = 10  # threshold of X in the correction
init_K = 20  # initial number of clusters in diri_sampler

#################################################
# a common user should not edit below this line #
#################################################

# Dictionary storing by read ID (key) all subsequences which are part of
# different windows
correction = {}
# Dictionary storing by read ID (key) the posterior for each of window
quality = {}

count = {}
count['A'] = 0
count['C'] = 0
count['G'] = 0
count['T'] = 0
count['X'] = 0
count['-'] = 0
clusters = [[]]
untouched = [[]]


def gzip_file(f_name):
    '''Gzip a file and return the name of the gzipped, removing the original
    '''
    import gzip
    f_in = open(f_name, 'rb')
    f_out = gzip.open(f_name + '.gz', 'wb')
    f_out.writelines(f_in)
    f_out.close()
    f_in.close()
    os.remove(f_in.name)

    return f_out.name


def parse_aligned_reads(reads_file):
    """Parse reads from a file with aligned reads
    """
    out_reads = {}

    if not os.path.isfile(reads_file):
        declog.error('There should be a file here: ' + reads_file)
        sys.exit('There should be a file here: ' + reads_file)
    else:
        declog.info('Using file ' + reads_file + ' of aligned reads')

    handle = open(reads_file)
    declog.debug('Parsing aligned reads')

    for h in handle:
        name, start, stop, mstart, mstop, this_m = h.rstrip().split('\t')
        if this_m == '':
            declog.warning('parsing empty read: %s' % h)
        out_reads[name] = [None, None, None, None, []]
        # start of the region of interest (0-based indexing)
        out_reads[name][0] = int(start)
        # end of the region of interest (0-based inde
        out_reads[name][1] = int(stop)
        # start index of the aligned read w.r.t reference (1-based indexing)
        out_reads[name][2] = int(mstart)
        # end index of the aligned read w.r.t reference (1-based indexing)
        out_reads[name][3] = int(mstop)
        out_reads[name][4] = this_m

    return out_reads


def windows(run_settings):
    """run b2w to make windows from bam
    """
    import subprocess
    bam, fasta, w, i, m, x, reg = run_settings
    dn = sys.path[0]
    my_prog = os.path.join(dn, 'b2w')
    my_arg = ' -w %i -i %i -m %i -x %i %s %s %s' % \
        (w, i, m, x, bam, fasta, reg)

    try:
        retcode = subprocess.call(my_prog + my_arg, shell=True)
        if retcode > 0:
            declog.error('%s %s' % (my_prog, my_arg))
            declog.error('b2w returned %i' % retcode)
        else:
            declog.debug('Finished making windows')
            declog.debug('b2w returned %i' % retcode)
    except OSError as ee:
        declog.error('Execution of b2w failed: %s' % ee)
    return retcode


def run_dpm(run_setting):
    """run the dirichlet process clustering
    """

    import shutil
    import subprocess

    filein, j, a, seed = run_setting

    # if cor.fas.gz exists, skip
    stem = filein.split('.reads')[0]
    corgz = 'corrected/%s.reads-cor.fas.gz' % stem
    if os.path.exists(corgz):
        declog.debug('file %s already analysed, skipping' % filein)
        return

    # if already run before, extract the read file
    fstgz = 'raw_reads/%s.reads.fas.gz' % stem
    if os.path.exists(filein):
        pass
    elif os.path.exists(fstgz):
        shutil.move(fstgz, './')
        subprocess.check_call(["gunzip", "%s-reads.gz" % stem])

    dn = sys.path[0]
    my_prog = os.path.join(dn, 'diri_sampler')
    my_arg = ' -i %s -j %i -t %i -a %f -K %d -R %d' % \
        (pipes.quote(filein), j, int(j * hist_fraction), a, init_K, seed)

    try:
        os.remove('./corrected.tmp')
        # os.remove('./assignment.tmp')
    except OSError:
        pass
    declog.debug(my_prog + my_arg)
    # runs the gibbs sampler for the dirichlet process mixture
    try:
        retcode = subprocess.call(my_prog + my_arg, shell=True)
        if retcode < 0:
            declog.error('%s %s' % (my_prog, my_arg))
            declog.error('Child %s terminated by SIG %d' % (my_prog, -retcode))
        else:
            declog.debug('run %s finished' % my_arg)
            declog.debug('Child %s returned %i' % (my_prog, retcode))
    except OSError as ee:
        declog.error('Execution of %s failed: %s' % (my_prog, ee))

    return


def correct_reads(chr_c, wstart, wend):
    ''' Parses corrected reads (in fasta format) and correct the reads
    '''
    # out_reads[match_rec.id][0] = qstart
    # out_reads[match_rec.id][1] = qstop
    # out_reads[match_rec.id][2] = mstart
    # out_reads[match_rec.id][3] = mstop
    # out_reads[match_rec.id][4] = Sequence...

    from Bio import SeqIO
    import gzip
    try:
        if os.path.exists('corrected/w-%s-%s-%s.reads-cor.fas.gz' %
                          (chr_c, wstart, wend)):
            cor_file = 'corrected/w-%s-%s-%s.reads-cor.fas.gz' % \
                (chr_c, wstart, wend)
            handle = gzip.open(cor_file)
        else:
            cor_file = 'w-%s-%s-%s.reads-cor.fas' % (chr_c, wstart, wend)
            handle = open(cor_file, 'r')

        for seq_record in SeqIO.parse(handle, 'fasta'):
            assert '\0' not in str(seq_record.seq), 'binary file!!!'
            read_id = seq_record.id
            try:
                correction[read_id][wstart] = list(str(seq_record.seq))
                quality[read_id][wstart] = \
                    float(seq_record.description.split('|')[1].split('=')[1])
                assert quality[read_id][wstart] <= 1.0, \
                    'try: quality must be < 1, %s' % cor_file
            except KeyError:
                correction[read_id] = {}
                quality[read_id] = {}
                correction[read_id][wstart] = list(str(seq_record.seq))
                quality[read_id][wstart] = \
                    float(seq_record.description.split('|')[1].split('=')[1])
                assert quality[read_id][wstart] <= 1.0, \
                    'except: quality must be < 1, %s' % cor_file
        handle.close()
        return
    except IOError:
        declog.warning('No reads in window %s?' % wstart)
        return


def get_prop(filename):
    """fetch the number of proposed clusters from .dbg file
    """
    import gzip

    if os.path.exists(filename):
        h = open(filename)
    elif os.path.exists(filename + '.gz'):
        h = gzip.open(filename + '.gz')
    elif os.path.exists('debug/' + filename):
        h = open('debug/' + filename)
    elif os.path.exists('debug/' + filename + '.gz'):
        h = gzip.open('debug/' + filename + '.gz')
    else:
        return 'not found'

    for l in h:
        if l.startswith('#made'):
            prop = int(l.split()[1])
            break
    h.close()
    return prop


def base_break(baselist):
    """break the tie if different corrections are found
    """
    import random

    for c1 in count:
        count[c1] = 0
    for c in baselist:
        if c.upper() != 'N':
            count[c.upper()] += 1

    maxm = 0
    out = []
    for b in count:
        if count[b] >= maxm:
            maxm = count[b]
    for b in count:
        if count[b] == maxm:
            out.append(b)

    rc = random.choice(out)

    return rc


def win_to_run(alpha_w, seed):
    '''returns windows to run on diri_sampler
    '''

    rn_list = []
    try:
        file1 = open('coverage.txt')
    except IOError:
        sys.exit('Coverage file generated by b2w not found.')

    for f1 in file1:
        winFile, chr1, beg, end, cov = f1.rstrip().split('\t')
        j = min(300000, int(cov) * 15)
        rn_list.append((winFile, j, alpha_w, seed))

    del(end)
    del(beg, chr1)
    return rn_list


def merge_corrected_reads(aligned_read):

    ID = aligned_read[0]
    seq = aligned_read[1][4]
    corrected_read = correction.get(ID)
    merged_corrected_read = []

    if corrected_read is not None:
        # rlen: length of the orginal read
        rlen = len(seq)
        # rstart: start index of the original read with respect to reference
        rstart = aligned_read[1][2]
        # posterior: fraction of times a given read was assigned to current
        #            cluster among those iterations that were recorded
        posterior = quality.get(ID)

        # kcr: extract start index of the aligned reads in all windows
        # vcr: extract sequence of the aligned reads in all windows
        kcr = np.array(list(corrected_read.keys()), dtype=int)
        vcr = np.array(list(corrected_read.values()))
        vcr_len = [v.size for v in vcr]

        for rpos in range(rlen):
            tp = rstart + rpos - kcr
            mask = np.logical_and(
                np.logical_and(np.greater_equal(tp, 0), np.less(tp, vcr_len)),
                np.greater(list(posterior.values()), min_quality)
            )
            if np.sum(mask > 0):
                # data structure needs this
                idx = np.argwhere(mask)
                this = [vcr[k][tp[k]]
                        for k in idx[0]]  # this is unlikely to be optimal
                if len(this) > 1:
                    corrected_base = base_break(this)
                else:
                    corrected_base = this[0]
            else:
                corrected_base = 'X'
            merged_corrected_read.append(corrected_base)

    return(ID, merged_corrected_read)


def main(in_bam, in_fasta, win_length=201, win_shifts=3, region='',
         max_coverage=10000, alpha=0.1, keep_files=True, seed=None):
    '''
    Performs the error correction analysis, running diri_sampler
    and analyzing the result
    '''
    from multiprocessing import Pool, cpu_count
    import glob
    import shutil
    import time

    import shorah_snv

    # set logging level
    declog.setLevel(logging.DEBUG)
    # This handler writes everything to a file.
    LOG_FILENAME = './dec.log'
    hl = logging.handlers.RotatingFileHandler(LOG_FILENAME, 'w',
                                              maxBytes=100000, backupCount=5)
    f = logging.Formatter("%(levelname)s %(asctime)s\
                          %(funcName)s %(lineno)d %(message)s")
    hl.setFormatter(f)
    declog.addHandler(hl)
    declog.info(' '.join(sys.argv))

    # check options
    if win_length % win_shifts != 0:
        sys.exit('Window size must be divisible by win_shifts')
    if win_min_ext < 1 / win_shifts:
        declog.warning('Some bases might not be covered by any window')
    if max_coverage / win_length < 1:
        sys.exit('Please increase max_coverage')
    if not os.path.isfile(in_bam):
        sys.exit("File '%s' not found" % in_bam)
    if not os.path.isfile(in_fasta):
        sys.exit("File '%s' not found" % in_fasta)
    if seed is None:
        seed = np.random.randint(100, size=1)

    incr = win_length // win_shifts
    max_c = max_coverage // win_length
    keep_all_files = keep_files

    # run b2w
    retcode = windows((in_bam, in_fasta, win_length, incr,
                       win_min_ext * win_length, max_c, region))
    if retcode is not 0:
        sys.exit('b2w run not successful')

    aligned_reads = parse_aligned_reads('reads.fas')
    r = list(aligned_reads.keys())[0]
    gen_length = aligned_reads[r][1] - aligned_reads[r][0]

    if win_length > gen_length:
        sys.exit('The window size must be smaller than the genome region')

    declog.info('%s reads are being considered' % len(aligned_reads))

    ############################################
    # Now the windows and the error correction #
    ############################################

    runlist = win_to_run(alpha, seed)
    declog.info('will run on %d windows' % len(runlist))
    # run diri_sampler on all available processors but one
    max_proc = max(cpu_count() - 1, 1)
    pool = Pool(processes=max_proc)
    pool.map(run_dpm, runlist)
    pool.close()
    pool.join()

    # prepare directories
    if keep_all_files:
        for sd_name in ['debug', 'sampling', 'freq', 'support',
                        'corrected', 'raw_reads']:
            try:
                os.mkdir(sd_name)
            except OSError:
                pass

    # parse corrected reads
    proposed = {}
    for i in runlist:
        winFile, j, a, s = i
        del(a)  # in future alpha might be different on each window
        parts = winFile.split('.')[0].split('-')
        chrom = '-'.join(parts[1:-2])
        beg = int(parts[-2])
        end = int(parts[-1])
        declog.info('reading windows for start position %s' % beg)
        # correct reads populates correction and quality, globally defined
        correct_reads(chrom, beg, end)
        stem = 'w-%s-%s-%s' % (chrom, beg, end)
        declog.info('this is window %s' % stem)
        dbg_file = stem + '.dbg'
        # if os.path.exists(dbg_file):
        proposed[beg] = (get_prop(dbg_file), j)
        declog.info('there were %s proposed' % str(proposed[beg][0]))

    # (re)move intermediate files
    if not keep_all_files:
        declog.info('removing intermediate files')
        tr_files = glob.glob('./w*reads.fas')
        tr_files.extend(glob.glob('./*.smp'))
        tr_files.extend(glob.glob('./w*.dbg'))
        for trf in tr_files:
            os.remove(trf)

        tr_files = glob.glob('./w*reads-cor.fas')
        tr_files.extend(glob.glob('./w*reads-freq.csv'))
        tr_files.extend(glob.glob('./w*reads-support.fas'))
        for trf in tr_files:
            if os.stat(trf).st_size == 0:
                os.remove(trf)
    else:

        for dbg_file in glob.glob('./w*dbg'):
            if os.stat(dbg_file).st_size > 0:
                gzf = gzip_file(dbg_file)
                try:
                    os.remove('debug/%s' % gzf)
                except OSError:
                    pass
                shutil.move(gzf, 'debug/')
            else:
                os.remove(dbg_file)

        for smp_file in glob.glob('./w*smp'):
            if os.stat(smp_file).st_size > 0:
                gzf = gzip_file(smp_file)
                try:
                    os.remove('sampling/%s' % gzf)
                except OSError:
                    pass
                shutil.move(gzf, 'sampling/')
            else:
                os.remove(smp_file)

        for cor_file in glob.glob('./w*reads-cor.fas'):
            if os.stat(cor_file).st_size > 0:
                gzf = gzip_file(cor_file)
                try:
                    os.remove('corrected/%s' % gzf)
                except OSError:
                    pass
                shutil.move(gzf, 'corrected/')
            else:
                os.remove(cor_file)

        for sup_file in glob.glob('./w*reads-support.fas'):
            if os.stat(sup_file).st_size > 0:
                gzf = gzip_file(sup_file)
                try:
                    os.remove('support/%s' % gzf)
                except OSError:
                    pass
                shutil.move(gzf, 'support/')
            else:
                os.remove(sup_file)

        for freq_file in glob.glob('./w*reads-freq.csv'):
            if os.stat(freq_file).st_size > 0:
                gzf = gzip_file(freq_file)
                try:
                    os.remove('freq/%s' % gzf)
                except OSError:
                    pass
                shutil.move(gzf, 'freq/')
            else:
                os.remove(freq_file)

        for raw_file in glob.glob('./w*reads.fas'):
            if os.stat(raw_file).st_size > 0:
                gzf = gzip_file(raw_file)
                try:
                    os.remove('raw_reads/%s' % gzf)
                except OSError:
                    pass
                shutil.move(gzf, 'raw_reads/')
            else:
                os.remove(raw_file)

    ############################################
    ##      Print the corrected reads         ##
    ##
    ## correction[read_id][wstart] = sequence ##
    ## quality[read_id][wstart] = posterior   ##
    # ##########################################
    declog.info('Merging windows of corrected reads')
    # Multi-threaded version
    params = list(aligned_reads.items())
    pool = Pool(processes=max_proc)
    to_correct = pool.map(merge_corrected_reads, params)
    pool.close()
    pool.join()

    declog.info('All corrected reads have been merged')

    ccx = {}
    cin_stem = '.'.join(os.path.split(in_bam)[1].split('.')[:-1])
    fch = open('%s.cor.fas' % cin_stem, 'w')
    declog.debug('writing to file %s.cor.fas' % cin_stem)
    for ID, seq_list in to_correct:
        cor_read = ''.join(seq_list)
        init_x = len(cor_read.lstrip('-')) - len(cor_read.lstrip('-X'))
        fin_x = len(cor_read.rstrip('-')) - len(cor_read.rstrip('-X'))
        cx = seq_list.count('X') - init_x - fin_x
        ccx[cx] = ccx.get(cx, 0) + 1
        if cx <= min_x_thresh and cor_read.lstrip('-X') != '':
            fch.write('>%s %d\n' % (ID, aligned_reads[ID][2] + init_x
                                    - aligned_reads[ID][0]))
            cc = 0
            for c in cor_read.lstrip('-X'):
                if c != 'X':
                    fch.write(str(c))
                    fch.flush()
                    cc = cc + 1
                    if cc % fasta_length == 0:
                        fch.write('\n')

            if cc % fasta_length != 0:
                fch.write('\n')
    print(ccx)
    fch.close()

    # write proposed_per_step to file
    ph = open('proposed.dat', 'w')
    ph.write('#base\tproposed_per_step\n')
    for kp in sorted(proposed):
        if proposed[kp] != 'not found':
            ph.write('%s\t%f\n' %
                     (kp, proposed[kp][0] / proposed[kp][1]))
    ph.close()

    declog.info('running snv.py')
    shorah_snv.main(reference=in_fasta, bam_file=in_bam,
                    increment=win_length // win_shifts, max_coverage=max_coverage)

    # tidy snvs
    try:
        os.mkdir('snv')
    except OSError:
        os.rename('snv', 'snv_before_%d' % int(time.time()))
        os.mkdir('snv')
    for snv_file in glob.glob('./SNV*'):
        shutil.move(snv_file, 'snv/')

    declog.info('dec.py ends')
