#!/usr/bin/env python3

# Copyright 2007-2018
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

"""dec.py draws windows on the bam file, cuts them into separate fasta files,
    calls diri_sampler to correct them, and merges the correction into a
    file of corrected reads
"""
import os
import pipes
import sys
import logging
import re
import shutil
import numpy as np

import libshorah

dn_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if __name__ == '__main__':
    if __package__ is None:
        os.sys.path.insert(1, dn_dir)
        mod = __import__('shorah')
        sys.modules["shorah"] = mod
        import shorah_snv
        import b2w
        import tiling
else:
    from . import shorah_snv
    from . import b2w
    from . import tiling

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

count = {
    'A': 0,
    'C': 0,
    'G': 0,
    'T': 0,
    'X': 0,
    '-': 0
}
clusters = [[]]
untouched = [[]]


def gzip_file(f_name):
    """Gzip a file and return the name of the gzipped, removing the original
    """
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
        logging.error('There should be a file here: %s', reads_file)
        sys.exit('There should be a file here: ' + reads_file)
    else:
        logging.info('Using file %s of aligned reads', reads_file)

    handle = open(reads_file)
    logging.debug('Parsing aligned reads')

    for h in handle:
        name, start, stop, mstart, mstop, this_m = h.rstrip().split('\t')
        if this_m == '':
            logging.warning('parsing empty read: %s', h)
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


def b2w_logging(run_settings):
    """run b2w to make windows from bam
    """
    bam, fasta, w, i, m, x, c, reg, ignore_indels = run_settings
    d = ' -d' if ignore_indels else ''
    my_arg = '-w %i -i %i -m %i -x %i -c %i%s %s %s %s' % \
        (w, i, m, x, c, d, bam, fasta, reg)
    logging.debug(f'To run standalone: python3 b2w.py {my_arg}')


def run_dpm(run_setting):
    """run the dirichlet process clustering
    """
    import subprocess

    filein, j, a, seed = run_setting

    # if cor.fas.gz exists, skip
    # greedy re match to handle situation where '.reads' appears in the ID
    stem = re.match(r'^(?P<stem>.*).reads', filein).group('stem')
    corgz = 'corrected/%s.reads-cor.fas.gz' % stem
    if os.path.exists(corgz):
        logging.debug('file %s already analysed, skipping', filein)
        return

    # if already run before, extract the read file
    fstgz = 'raw_reads/%s.reads.fas.gz' % stem
    if os.path.exists(filein):
        pass
    elif os.path.exists(fstgz):
        shutil.move(fstgz, './')
        subprocess.check_call(["gunzip", "%s-reads.gz" % stem])

    # dn = sys.path[0]
    #my_prog = shlex.quote(diri_exe)  # os.path.join(dn, 'diri_sampler')
    #my_arg = ' -i %s -j %i -t %i -a %f -K %d -R %d' % \
    #    (pipes.quote(filein), j, int(j * hist_fraction), a, init_K, seed)

    # TODO integration
    logging.debug('Exec dpm_sampler')
    try:
        os.remove('./corrected.tmp')
        # os.remove('./assignment.tmp')
    except OSError:
        pass


    # runs the gibbs sampler for the dirichlet process mixture
    try:
        logging.debug(f"{filein} {j} {a} {int(j * hist_fraction)} {init_K} {seed}")
        retcode = libshorah.exec_dpm_sampler(
            pipes.quote(filein),
            j,
            a,
            int(j * hist_fraction),
            K_cluster_start=init_K,
            R_seed=seed
        )
        if retcode == 0:
            logging.debug(f'{filein} - Run finished successfully.')
        else:
            logging.error(f'{filein} - Run failed with return code %i.', retcode)
    except Exception as e:
        logging.error(f'{filein} - Run failed: {e}')


    return


def correct_reads(chr_c, wstart, wend):
    """ Parses corrected reads (in fasta format) and correct the reads
    """
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
            handle = gzip.open(
                cor_file, 'rb' if sys.version_info < (3, 0) else 'rt')
        else:
            cor_file = 'w-%s-%s-%s.reads-cor.fas' % (chr_c, wstart, wend)
            handle = open(cor_file, 'r')

        for seq_record in SeqIO.parse(handle, 'fasta'):
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
        logging.warning('No reads in window %s?', wstart)
        return


def get_prop(filename):
    """fetch the number of proposed clusters from .dbg file
    """
    import gzip

    if os.path.exists(filename):
        h = open(filename)
    elif os.path.exists(filename + '.gz'):
        h = gzip.open(filename + '.gz',
                      'rb' if sys.version_info < (3, 0) else 'rt')
    elif os.path.exists('debug/' + filename):
        h = open('debug/' + filename)
    elif os.path.exists('debug/' + filename + '.gz'):
        h = gzip.open('debug/' + filename + '.gz',
                      'rb' if sys.version_info < (3, 0) else 'rt')
    else:
        return 'not found'

    prop = 'not found'
    for l in h:
        if l.startswith('#made'):
            prop = int(l.split()[1])
            break
    h.close()
    return prop


def base_break(baselist):
    """Break the tie if different corrections are found."""
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
    """Return windows to run on diri_sampler."""

    rn_list = []
    try:
        file1 = open('coverage.txt')
    except IOError:
        sys.exit('Coverage file generated by b2w not found.')

    for f1 in file1:
        winFile, chr1, beg, end, cov = f1.rstrip().split('\t')
        j = min(300000, int(cov) * 20)
        rn_list.append((winFile, j, alpha_w, seed))

    del end
    del(beg, chr1)
    return rn_list


def merge_corrected_reads(aligned_read):
    if aligned_read is None:
        print("empty window found", file=sys.stderr)
        return (None, [])

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
        vcr = np.array([np.array(v) for v in corrected_read.values()], dtype=object) # FIXED dtype
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


# def main(in_bam, in_fasta, win_length=201, win_shifts=3, region='',
#         max_coverage=10000, alpha=0.1, keep_files=True, seed=None):
def main(args):
    """
    Performs the error correction analysis, running diri_sampler
    and analyzing the result
    """
    from multiprocessing import Pool, cpu_count
    import glob
    import math
    import time
    import pysam

    in_bam = args.b
    in_fasta = args.f
    win_length = args.w # TODO remove this var
    win_shifts = args.win_shifts # TODO remove this var
    region = args.r
    max_coverage = args.max_coverage
    alpha = args.a
    cov_thrd = args.cov_thrd
    keep_files = args.keep_files
    seed = args.seed
    ignore_indels = args.ignore_indels
    maxthreads = args.maxthreads
    path_insert_file = args.path_insert_file

    logging.info(' '.join(sys.argv))

    # check options
    if win_length % win_shifts != 0:
        sys.exit('Window size must be divisible by win_shifts')
    if win_min_ext < 1 / win_shifts:
        logging.warning('Some bases might not be covered by any window')
    if max_coverage / win_length < 1:
        sys.exit('Please increase max_coverage')
    if not os.path.isfile(in_bam):
        sys.exit("File '%s' not found" % in_bam)
    if not os.path.isfile(in_fasta):
        sys.exit("File '%s' not found" % in_fasta)
    if seed is None:
        seed = np.random.randint(100, size=1)

    incr = win_length // win_shifts
    keep_all_files = keep_files

    # run b2w

    logging.info('starting b2w')
    try:
        if ignore_indels == True:
            raise NotImplementedError('This argument was deprecated.')
        b2w_logging((in_bam, in_fasta, win_length, incr, win_min_ext *
                       win_length, max_coverage, cov_thrd, region, ignore_indels))

        if path_insert_file == None and region == "": # special case if no region defined
            samfile = pysam.AlignmentFile(
                in_bam,
                "r", # auto-detect bam/cram (rc)
                reference_filename=in_fasta,
                threads=1
            )
            if samfile.nreferences != 1:
                raise NotImplementedError("There are multiple references in this alignment file.")
            strategy = tiling.EquispacedTilingStrategy(
                f"{samfile.references[0]}:1-{samfile.lengths[0]}",
                win_length,
                incr,
                False,
                True
            )
        elif path_insert_file == None:
            strategy = tiling.EquispacedTilingStrategy(region, win_length, incr, True)
        else:
            strategy = tiling.PrimerTilingStrategy(path_insert_file)
            if region != "":
                logging.warn(f"region is set to {region} but is not used with this tiling strategy")

        logging.info(f"Using tiling strategy: {type(strategy).__name__}")

        b2w.build_windows(
            in_bam,
            strategy,
            math.floor(win_min_ext * win_length),
            max_coverage,
            cov_thrd,
            in_fasta
        )
        logging.info('finished b2w')

    except Exception as e:
        logging.debug(e)
        sys.exit('b2w run not successful')

    aligned_reads = parse_aligned_reads('reads.fas')
    if len(aligned_reads) == 0:
        msg = 'No reads found in the requested region %s' % region
        logging.debug(msg)
        print(msg, file=sys.stderr)
        logging.info('shotgun run ends with no processing')
        # technically it is a success: we did produce what we were asked for.
        # it just happens that we were asked to produce nothing and thus got nothing to output
        sys.exit(0)

    r = list(aligned_reads.keys())[0]
    gen_length = aligned_reads[r][1] - aligned_reads[r][0]

    if win_length > gen_length:
        sys.exit('The window size must be smaller than the genome region')

    logging.info('%s reads are being considered', len(aligned_reads))

    ############################################
    # Now the windows and the error correction #
    ############################################

    runlist = win_to_run(alpha, seed)
    logging.info('will run on %d windows', len(runlist))
    # run diri_sampler on all available processors but one
    max_proc = max(cpu_count() - 1, 1)
    if maxthreads:
        max_proc = min(max_proc, maxthreads)
    logging.info('CPU(s) count %u, max thread limit %u, will run %u parallel dpm_sampler', cpu_count(), maxthreads, max_proc)
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
        del a  # in future alpha might be different on each window
        del s
        # greedy re match to handle situation where '.' or '-' appears in the
        # ID
        parts = re.match(
            r'^w-(?P<chrom>.*)-(?P<beg>\d+)-(?P<end>\d+).reads', winFile)
        chrom = parts.group('chrom')
        beg = int(parts.group('beg'))
        end = int(parts.group('end'))
        del parts
        logging.info('reading windows for start position %s', beg)
        # correct reads populates correction and quality, globally defined
        correct_reads(chrom, beg, end)
        stem = 'w-%s-%u-%u' % (chrom, beg, end)
        logging.info('this is window %s', stem)
        dbg_file = stem + '.dbg'
        # if os.path.exists(dbg_file):
        proposed[beg] = (get_prop(dbg_file), j)
        logging.info('there were %s proposed', str(proposed[beg][0]))

    # (re)move intermediate files
    if not keep_all_files:
        logging.info('removing intermediate files')
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

        for raw_file in glob.glob('./w*reads.fas') + glob.glob('./w*ref.fas'):
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
    logging.info('Merging windows of corrected reads')
    # Multi-threaded version
    params = list(aligned_reads.items())
    pool = Pool(processes=max_proc)
    to_correct = pool.map(merge_corrected_reads, params)
    pool.close()
    pool.join()

    logging.info('All corrected reads have been merged')

    ccx = {}
    # handle case where bamfile has no dots in name
    cin_stem = re.sub(r'\.[^.]+$', r'', os.path.split(in_bam)[1])
    logging.debug('writing to file %s.cor.fas', cin_stem)
    with open('%s.cor.fas' % cin_stem, 'w') as fch:
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

    # write proposed_per_step to file
    ph = open('proposed.dat', 'w')
    ph.write('#base\tproposed_per_step\n')
    for kp in sorted(proposed):
        if proposed[kp] != 'not found' and 'not found' not in proposed[kp]:
            ph.write('%s\t%f\n' %
                     (kp, proposed[kp][0] / proposed[kp][1]))
    ph.close()

    logging.info('running snv.py')
    args.increment = win_length // win_shifts # TODO remove dependency on these vars
    shorah_snv.main(args)

    # tidy snvs
    try:
        os.mkdir('snv')
    except OSError:
        os.rename('snv', 'snv_before_%d' % int(time.time()))
        os.mkdir('snv')

    for snv_file in glob.glob('./raw_snv*') + glob.glob('./SNV*'):
        shutil.move(snv_file, 'snv/')

    logging.info('shotgun run ends')
