#!/usr/bin/env python3

# Copyright 2007-2018
# Niko Beerenwinkel,
# Nicholas Eriksson,
# Moritz Gerstung,
# Lukas Geyrhofer,
# Osvaldo Zagordi,
# Kerensa McElroy,
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


'''
    ------------
    Output:
    a file of raw snvs, parsed from the directory support,
    and a directory containing snvs resulting from strand
    bias tests with different sigma values

    Output Files:
    - raw_snv_collected.tsv: contains all SNVs from all windows. Each in a seperate row.
    - raw_snv.tsv: contains all SNVs
    - SNV.tsv: contains only SNVs that are covered by at least {min_windows_coverage}
               windows
    - SNVs_*.tsv, strand bias filter applied
    - raw_snvs_*.tsv: strand bias filter applied
    - SNVs_*_final.vcf: final file
    - SNVs_*_final.csv: file file
    ------------
'''

import glob
import gzip
import os
import sys
import warnings
from collections import namedtuple
from dataclasses import dataclass

import logging

import libshorah

SNP_id = namedtuple('SNP_id', ['pos', 'var'])
@dataclass
class SNV:
    chrom: str
    pos: int
    ref: str
    var: str
    freq: float = 0.0
    support: float = 0.0

standard_header_row = ['Chromosome', 'Pos', 'Ref', 'Var', 'Frq', 'Pst']

def deletion_length(seq):
    """Determines the length of the deletion. Note that a sequence migth have
       more than one deletion
       seq: substring of the reconstructed haplotype
    """
    count = 0
    for c in seq:
        if c == '-':
            count += 1
        else:
            break
    return count


def parseWindow(line, ref1, threshold=0.9):
    """SNVs from individual support files, getSNV will build
        the consensus SNVs
        It returns a dictionary called snp with the following structure
        key:   pos.allele (position on the reference file and mutated base)
        value: reference name, position, reference_base, mutated base,
               average number of reads, posterior times average n of reads
    """
    from Bio import SeqIO
    from re import search

    snp = {}
    reads = 0.0
    winFile, chrom, beg, end, cov = line.rstrip().split('\t')
    del([winFile, cov])
    filename = 'w-%s-%s-%s.reads-support.fas' % (chrom, beg, end)

    # take cares of locations/format of support file
    if os.path.exists(filename):
        pass
    elif os.path.exists('support/' + filename):
        filename = 'support/' + filename
    elif os.path.exists('support/' + filename + '.gz'):
        filename = 'support/' + filename + '.gz'
    elif os.path.exists(filename + '.gz'):
        filename = filename + '.gz'

    try:
        if filename.endswith('.gz'):
            window = gzip.open(
                filename, 'rb' if sys.version_info < (3, 0) else 'rt')
        else:
            window = open(filename, 'r')
    except IOError:
        logging.error('File not found')
        return snp

    beg = int(beg)
    end = int(end)
    refSlice = ref1[chrom][beg - 1:end]
    max_snv = -1
    # sequences in support file exceeding the posterior threshold
    for s in SeqIO.parse(window, 'fasta'):
        seq = str(s.seq).upper()
        match_obj = search('posterior=(.*)\s*ave_reads=(.*)', s.description)
        post, av = float(match_obj.group(1)), float(match_obj.group(2))
        if post > 1.0:
            warnings.warn('posterior = %4.3f > 1' % post)
            logging.warning('posterior = %4.3f > 1' % post)
        if post >= threshold:
            reads += av
            pos = beg
            tot_snv = 0
            aux_del = -1
            for idx, v in enumerate(refSlice):  # iterate on the reference
                if v != seq[idx]:  # SNV detected, save it
                    if seq[idx] == '-':
                        # Avoid counting multiple times a long deletion in the
                        # same haplotype
                        if idx > aux_del:
                            tot_snv += 1
                            # Check for gap characters and get the deletion
                            # length
                            del_len = deletion_length(seq[idx:])
                            aux_del = idx + del_len
                            snp_id = SNP_id(pos=pos, var=seq[idx:aux_del])

                            if snp_id in snp:
                                # Aggregate counts for long deletions which
                                # are observed in multiple haplotypes
                                snp[snp_id].freq += av
                                snp[snp_id].support += post * av
                            else:
                                # Comply with the convention to report deletion
                                # in VCF format. Position correspond to the
                                # preceding position w.r.t. the reference
                                # without a deletion
                                pos_prev = pos - 1
                                reference_seq = ref1[chrom][
                                    (pos_prev - 1):(pos_prev + del_len)]
                                snp[snp_id] = SNV(
                                    chrom, pos_prev, reference_seq,
                                    reference_seq[0], av, post * av)
                    else:
                        tot_snv += 1
                        snp_id = SNP_id(pos=pos, var=seq[idx])
                        if snp_id in snp:
                            snp[snp_id].freq += av
                            snp[snp_id].support += post * av
                        else:
                            snp[snp_id] = SNV(
                                chrom, pos, v, seq[idx], av, post * av)
                pos += 1
            if tot_snv > max_snv:
                max_snv = tot_snv

    logging.info('max number of snvs per sequence found: %d', max_snv)
    # normalize
    for k, v in snp.items():
        v.support /= v.freq
        v.freq /= reads

    window.close() # TODO

    return snp


def add_SNV_to_dict(all_dict, add_key, add_val):
    if add_key in all_dict.keys():
        all_dict[add_key].append(add_val)
    else:
        all_dict.update({add_key: [add_val]})

    return all_dict


def getSNV(ref, window_thresh=0.9):
    """Parses SNV from all windows and output the dictionary with all the
    information.

    Writes file: collected_snv.tsv

    Note: SNVs that were discovered by several windows will occurr in several
    rows of the file (one row = one window).

    Returns a dict containing all SNVs, already grouped together per SNV_id.
    """

    all_snp = {}
    # cycle over all windows reported in coverage.txt
    with open('coverage.txt') as cov_file, open('raw_snv_collected.tsv', 'w') as f_collect:
        f_collect.write('\t'.join(standard_header_row) + '\n')
        for i in cov_file:
            snp = parseWindow(i, ref, window_thresh)
            #all_snp = join_snp_dict(all_snp, snp)
            for SNV_id, val in sorted(snp.items()):
                all_snp = add_SNV_to_dict(all_snp, SNV_id, val)
                f_collect.write('\t'.join(map(str, [val.chrom, val.pos, val.ref, val.var, val.freq, val.support])) + "\n")

    return all_snp

def writeRaw(all_snp, min_windows_coverage):
    """Write the SNVs that were collected for each window into
        - raw_snv.tsv: contains all SNVs
        - SNV.tsv: contains only SNVs that are covered by at least {min_windows_coverage}
                   windows
    """
    header_row =  ['Chromosome', 'Pos', 'Ref', 'Var']

    import numpy as np
    max_number_window_covering_SNV = np.max([len(val) for _, val in sorted(all_snp.items())])

    header_row = header_row + ['Frq'+str(k+1) for k in range(max_number_window_covering_SNV)]
    header_row = header_row + ['Pst'+str(k+1) for k in range(max_number_window_covering_SNV)]

    with open('raw_snv.tsv', 'w') as f_raw_snv, open('SNV.tsv', 'w') as f_SNV:
        f_raw_snv.write('\t'.join(header_row) + '\n')
        f_SNV.write('\t'.join(header_row) + '\n')

        for SNV_id, val in sorted(all_snp.items()):

            write_line = [val[0].chrom, val[0].pos, val[0].ref, val[0].var]
            freq_list = [single_val.freq for single_val in val]
            support_list = [single_val.support for single_val in val]

            number_window_covering_SNV = len(freq_list)

            if number_window_covering_SNV < max_number_window_covering_SNV:
                filler = (max_number_window_covering_SNV - len(freq_list)) * ['*']
                freq_list = freq_list + filler
                support_list = support_list + filler

            write_line = write_line + freq_list + support_list

            f_raw_snv.write('\t'.join(map(str, write_line)) + '\n')

            if number_window_covering_SNV >= min_windows_coverage:
                f_SNV.write('\t'.join(map(str, write_line)) + '\n')


def sb_filter(in_bam, file_to_append, out_file_prefix, sigma, amplimode="",
              drop_indels="", max_coverage=100000): # TODO max_coverage is 10 times higher than in Cpp
    """run strand bias filter calling the external program 'fil'
    """

    logging.debug('Running fil')
    logging.debug(f"{in_bam} {file_to_append} {out_file_prefix} {sigma} {max_coverage}")
    retcode = libshorah.fil(
        in_bam,
        file_to_append,
        out_file_prefix,
        sigma,
        max_coverage,
        False if amplimode == "" else True,
        False if drop_indels == "" else True
    )
    return retcode


# is last column of final output file
def BH(p_vals, n):
    """performs Benjamini Hochberg procedure, returning q-vals'
       you can also see http://bit.ly/QkTflz
    """
    # p_vals contains the p-value and the index where it has been
    # found, necessary to assign the correct q-value
    q_vals_l = []
    prev_bh = 0
    for i, p in enumerate(p_vals):
        # Sometimes this correction can give values greater than 1,
        # so we set those values at 1
        bh = p[0] * n / (i + 1)
        bh = min(bh, 1)
        # To preserve monotonicity in the values, we take the
        # maximum of the previous value or this one, so that we
        # don't yield a value less than the previous.
        bh = max(bh, prev_bh)
        prev_bh = bh
        q_vals_l.append((bh, p[1]))
    return q_vals_l


def main(args):
    '''main code
    '''
    from Bio import SeqIO
    from math import log10
    import csv
    import inspect
    from datetime import date

    reference = args.f
    bam_file = args.b
    sigma = args.sigma
    increment = args.increment
    max_coverage = args.max_coverage
    ignore_indels = args.ignore_indels
    posterior_thresh = args.posterior_thresh

    logging.info(str(inspect.getfullargspec(main)))
    ref_m = dict([[s.id, str(s.seq).upper()]
                  for s in SeqIO.parse(reference, 'fasta')])

    # snpD_m is the file with the 'consensus' SNVs (from different windows)
    logging.debug('now parsing SNVs')
    all_SNVs = getSNV(ref_m, posterior_thresh)
    writeRaw(all_SNVs, min_windows_coverage=1) # min_windows_coverage = 2 usually 

    with open('raw_snv.tsv') as f_raw_snv:
        windows_header_row = f_raw_snv.readline().split('\t')
        windows_header_row[-1] = windows_header_row[-1].split('\n')[0]

    d = ' -d' if ignore_indels else ''

    a = ' -a' if increment == 1 else '' # TODO when is increment == 1 (amplimode)

    # run strand bias filter
    #retcode_m = sb_filter(bam_file, "raw_snv.tsv", "raw_snvs_", sigma, amplimode=a, drop_indels=d,
    #                      max_coverage=max_coverage)
    retcode_n = sb_filter(bam_file, "SNV.tsv", "SNVs_", sigma, amplimode=a, drop_indels=d,
                          max_coverage=max_coverage)

    if retcode_n != 0 :
        logging.error('sb_filter exited with error %d', retcode_m)
        sys.exit()

    # parse the p values from raw_snv_* file
    snpFile = glob.glob('SNVs_*.tsv')[0]  # takes the first file only
    logging.debug(f"For BH - selected file: {snpFile}")

    write_list = []
    d = {}

    with open(snpFile) as f:
        for line_no, line in enumerate(f):
            parts = line.rstrip().split('\t')
            write_list.append(parts)
            idx = parts[0] + parts[1] + parts[2] + parts[3]
            if idx in d:
                d[idx][1].append(line_no)
            else:
                d[idx] = (float(parts[-1]), [line_no])

    p_vals = list(d.values())

    # sort p values, correct with Benjamini Hochberg and append to output
    p_vals.sort()
    q_vals = BH(p_vals, len(p_vals))

    for q, indices in q_vals:
        for i in indices:
            write_list[i].append(q)

    # Write ShoRAH csv output file
    if 'csv' in args.format:
        csv_file = '{}_final.csv'.format(os.path.splitext(snpFile)[0])
        header_row = windows_header_row + ['Fvar', 'Rvar', 'Ftot', 'Rtot', 'Pval', 'Qval']
        with open(csv_file, 'w') as cf:
            writer = csv.writer(cf)
            writer.writerow(header_row)
            # only print when q >= 5%
            for wl in write_list:
                if wl[-1] >= 0.05:
                    writer.writerow(wl)

    if 'vcf' in args.format:
        VCF_file = f'{os.path.splitext(snpFile)[0]}_final.vcf'
        VCF_meta = [
            '##fileformat=VCFv4.2',
            f'##fileDate={date.today():%Y%m%d}',
            f'##source=ShoRAH_{args.version}',
            f'##reference={args.f}'
        ]
        for ref_name, ref_seq in ref_m.items():
            VCF_meta.append(f'##contig=<ID={ref_name},length={len(ref_seq)}>',)
        VCF_meta.extend([
            '##INFO=<ID=Fvar,Number=1,Type=Integer,Description="Number of forward reads with variant">',
            '##INFO=<ID=Rvar,Number=1,Type=Integer,Description="Number of reverse reads with variant">',
            '##INFO=<ID=Ftot,Number=1,Type=Integer,Description="Total number of forward reads">',
            '##INFO=<ID=Rtot,Number=1,Type=Integer,Description="Total number of reverse reads">',
            '##INFO=<ID=Pval,Number=1,Type=Float,Description="P-value for strand bias">',
            '##INFO=<ID=Qval,Number=1,Type=Float,Description="Q-value for strand bias">'
        ])
        VCF_meta.extend([
                '##INFO=<ID=FreqX,Number=1,Type=Float,Description="Frequency of the variant in window X">',
                '##INFO=<ID=PostX,Number=1,Type=Float,Description="Posterior probability of the variant in window x">',
            ])

        with open(VCF_file, 'w') as vcf:
            vcf.write('\n'.join(VCF_meta))
            # VCFv4.2 HEADER line
            vcf.write('\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')
            # Iterate over single SNV lines and write them to output file
            for wl in write_list:
                # only print when q >= 5%
                if wl[-1] >= 0.05:
                    info = f'Fvar={wl[-6]};Rvar={wl[-5]};Ftot={wl[-4]};' \
                        f'Rtot={wl[-3]};Pval={wl[-2]};Qval={wl[-1]}'
                    if increment == 1:
                        post_avg = min([1, float(wl[5])])
                        info = f'Freq={wl[4]};Post={wl[5]};' + info
                    else:
                        freq_str = ';'.join([f'Freq{i+1}={j}'
                            for i, j in enumerate(wl[4:7]) if j != '*'])
                        post_str = ';'.join([f'Post{i+1}={j}'
                            for i, j in enumerate(wl[7:10]) if j != '*'])
                        info = f'{freq_str};{post_str};{info}'.replace('-', '0')
                        post_all = []
                        for freq, post in zip(wl[4:7], wl[7:10]):
                            if freq == '*':
                                pass
                            elif freq == '-':
                                post_all.append(0)
                            else:
                                post_all.append(min([1, float(post)]))
                        # Calculate posterior average
                        post_avg = sum(post_all) / len(post_all)
                    # Calculate a Phred quality score where the base calling
                    # error probabilities is set to (1 - posterior avg).
                    # Maximum is set to 100.
                    try:
                        qual_norm = -10 * log10(1 - post_avg)
                    except ValueError:
                        qual_norm = 100

                    vcf.write(f'\n{wl[0]}\t{wl[1]}\t.\t{wl[2]}\t{wl[3]}'
                              f'\t{qual_norm}\tPASS\t{info}')
