#!/usr/bin/env python

import pysam
from dataclasses import dataclass
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal


@dataclass
class StrandCounts:
    forward: int = 0
    reverse: int = 0


class AlignedRead():

    def __init__(self, read):
        self.read = read

    def is_contained(self, start, end):

        cond = (self.read.reference_start <= start and
                self.read.reference_end >= end)
        return True if cond else False


def read_counts(chrm, start, end, bamfile):
    read_count = StrandCounts()

    # load alignment file
    with pysam.AlignmentFile(bamfile, 'rb') as alnfile:
        for read in alnfile.fetch(reference=chrm, start=start, end=end):
            aligned_read = AlignedRead(read)
            if aligned_read.is_contained(start, end):
                if read.is_reverse:
                    read_count.reverse += 1
                else:
                    read_count.forward += 1

    return read_count.forward, read_count.reverse


def count_matching_variants(chrm, start, length, variant, is_del, bamfile):
    variant_count = StrandCounts()

    with pysam.AlignmentFile(bamfile, 'rb') as alnfile:
        # The 'region' parameter expects positions using 1-based indexing
        region = f'{chrm}:{start + 1}-{start + 1}'
        for column in alnfile.pileup(region=region,
                                     truncate=True, max_depth=100000,
                                     ignore_overlaps=False):
            # Pileup returns a matrix containing all reads which cover a
            # specified region. However, it returns all positions covered by at
            # least one read - including positions outside the region of
            # interest, unless 'tuncate=True'
            for read in column.pileups:
                if is_del:
                    # Handle deletions
                    if -read.indel == length:
                        if read.alignment.is_reverse:
                            variant_count.reverse += 1
                        else:
                            variant_count.forward += 1
                else:
                    # Handle subsitutions
                    if not read.is_del and not read.is_refskip:
                        base = read.alignment.query_sequence[
                            read.query_position]
                        if base == variant:
                            if read.alignment.is_reverse:
                                variant_count.reverse += 1
                            else:
                                variant_count.forward += 1

    return variant_count.forward, variant_count.reverse


def main(bamfile, snvsfile, outfile):
    # Parse SNV and for each SNV count number of read covering the target
    # region and the number of reads supporting the SNV
    df_snvs = pd.read_csv(snvsfile, sep="\t", header=0, compression=None)
    # convert to 0-based
    df_snvs["Pos"] -= 1
    # Deletion length
    aux_len = df_snvs["Ref"].str.len() - 1
    del_mask = aux_len > 0
    df_snvs["Del_len"] = np.nan
    df_snvs.loc[del_mask, "Del_len"] = aux_len[del_mask]
    df_snvs["Is_del"] = del_mask

    df_snvs[["Variant_forward", "Variant_reverse"]] = df_snvs.apply(
        lambda x: count_matching_variants(
            x["Chromosome"], x["Pos"], x["Del_len"], x["Var"], x["Is_del"],
            bamfile), axis=1, result_type="expand")

    # Temporarily change to start position for the deletions to count the
    # number of reads that cover the long deletion in full length
    df_snvs.loc[del_mask, "Pos"] += 1
    df_snvs.loc[del_mask, "Pos_end"] = df_snvs["Pos"] + aux_len
    df_snvs.loc[~del_mask, "Pos_end"] = df_snvs["Pos"] + 1
    df_snvs[["Depth_forward", "Depth_reverse"]] = df_snvs.apply(
        lambda x: read_counts(
            x["Chromosome"], x["Pos"], x["Pos_end"], bamfile), axis=1,
        result_type="expand")

    # Revert operations for the reporting position
    df_snvs.loc[~del_mask, "Pos"] += 1

    # Clean-up for comparison
    df_snvs = df_snvs.drop(columns=["Del_len", "Is_del", "Pos_end"])

    # Load output produced by ShoRAH for comparison
    df_out = pd.read_csv(outfile, sep="\t", header=None, compression=None)
    df_out = df_out.iloc[:, :-1]

    df_out.columns = df_snvs.columns
    assert_frame_equal(df_snvs, df_out)


if __name__ == '__main__':
    # Input data
    bamfile = "test_aln.cram"
    snvsfile = "SNV.txt"
    outfile = "SNVs_0.010000.txt"

    main(bamfile=bamfile, snvsfile=snvsfile, outfile=outfile)
