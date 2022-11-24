import pysam
from typing import Optional
from shorah.tiling import TilingStrategy, EquispacedTilingStrategy
import numpy as np
import math

def _write_to_file(lines, file_name):
    with open(file_name, "w") as f:
        f.writelines("%s\n" % l for l in lines)

def _calc_location_maximum_reads(samfile, reference_name, maximum_reads):
    budget = dict()
    indel_map = set() # TODO quick fix because pileup created duplicates; why?

    for pileupcolumn in samfile.pileup(
        reference_name,
        max_depth=1_000_000, # TODO big enough?
        stepper="nofilter",
        multiple_iterators=False,
        ignore_overlaps=False,
        ignore_orphans=False,
        min_base_quality=0,):
        budget[pileupcolumn.reference_pos] = min(
            pileupcolumn.nsegments,
            maximum_reads-1 # minus 1 because because maximum_reads is exclusive
        )

        for pileupread in pileupcolumn.pileups: # TODO generalize
            if pileupread.indel > 0 or pileupread.is_del:
                indel_map.add((
                    pileupread.alignment.query_name,
                    pileupread.alignment.reference_start,
                    pileupcolumn.reference_pos,
                    pileupread.indel,
                    pileupread.is_del
                ))

    # ascending reference_pos are necessary for later steps
    indel_map = sorted(indel_map, key=lambda tup: tup[2])

    return budget, indel_map

def _run_one_window(samfile, window_start, reference_name, window_length,
        minimum_overlap, permitted_reads_per_location, counter,
        exact_conformance_fix_0_1_basing_in_reads, indel_map):

    arr = []
    arr_read_summary = []
    arr_read_qualities_summary = []

    iter = samfile.fetch(
        reference_name,
        window_start, # 0 based
        window_start + window_length # arg exclusive as per pysam convention
    )

    for read in iter:

        if (read.reference_start is None) or (read.reference_end is None):
            continue
        first_aligned_pos = read.reference_start # this is 0-based
        last_aligned_pos = read.reference_end - 1 #reference_end is exclusive


        if permitted_reads_per_location[first_aligned_pos] == 0:
            continue
        else:
            permitted_reads_per_location[first_aligned_pos] -= 1

        # start_cut_out is 0-based while window_start is 1-based
        start_cut_out = window_start - first_aligned_pos
        end_cut_out = start_cut_out + window_length

        s = slice(max(0, start_cut_out), end_cut_out)
        full_read = list(read.query_sequence)
        full_qualities = list(read.query_qualities)

        for ct_idx, ct in enumerate(read.cigartuples):
            if ct[0] in [0,1,2,7,8]: # 0 = BAM_, 1 = BAM_CINS, 2 = BAM_CDEL, 7 = BAM_CEQUAL, 8 = BAM_CDIFF
                pass
            elif ct[0] == 4: # 4 = BAM_CSOFT_CLIP
                for _ in range(ct[1]):
                    k = 0 if ct_idx == 0 else len(full_read)-1
                    full_read.pop(k)
                    full_qualities.pop(k)
                if ct_idx != 0 and ct_idx != len(read.cigartuples)-1:
                    raise ValueError("Soft clipping only possible on the edges of a read.")
            else:
                raise NotImplementedError("CIGAR op code found that is not implemented:", ct[0])

        extended_window_mode = False # TODO

        for i in indel_map:
            if i[0] == read.query_name and i[1] == first_aligned_pos:
                if i[4] == 1: # if del
                    full_read.insert(i[2] - first_aligned_pos, "-")
                    full_qualities.insert(i[2] - first_aligned_pos, "2")
                if i[4] == 0 and not extended_window_mode:
                    assert i[3] > 0
                    for _ in range(i[3]):
                        full_read.pop(i[2] + 1 - first_aligned_pos)
                        full_qualities.pop(i[2] + 1 - first_aligned_pos)
                if i[3] != 0 and i[4] == 1:
                    # TODO implement
                    raise NotImplementedError("Deletions larger than 1 not implemented.")

            # TODO new
            # if extended_window_mode:
            #     del_exists_already = False
            #     if i[0] != read.query_name and first_aligned_pos <= i[2] <= last_aligned_pos and i[4] == 1:
            #         for j in indel_map:
            #             if j[0] == read.query_name and j[1] == first_aligned_pos and j[4] == 1 and i[2] == j[2]:
            #                 del_exists_already = True
            #                 break
            #         if not del_exists_already:
            #             full_read.insert(i[2] - first_aligned_pos, "-")
            #             change_in_reference_space += 1



        full_read = ("".join(full_read))

        if (first_aligned_pos < window_start + 1 + window_length - minimum_overlap
                and last_aligned_pos >= window_start + minimum_overlap - 2 # TODO justify 2
                and len(full_read) >= minimum_overlap):

            cut_out_read = full_read[s]
            cut_out_qualities = full_qualities[s]

            k = window_start + window_length - 1 - last_aligned_pos
            if k > 0:
                cut_out_read = cut_out_read + k * "N"
                cut_out_qualities = cut_out_qualities + k * [2]
                # Phred scores have a minimal value of 2, where an “N” is inserted
                # https://www.zymoresearch.com/blogs/blog/what-are-phred-scores
            if start_cut_out < 0:
                cut_out_read = -start_cut_out * "N" + cut_out_read
                cut_out_qualities = -start_cut_out * [2] + cut_out_qualities

            assert len(cut_out_read) == window_length, (
                "read unequal window size", read.query_name, first_aligned_pos, cut_out_read, len(cut_out_read)
            )
            assert len(cut_out_qualities) == window_length, (
                "quality unequal window size"
            )

            if exact_conformance_fix_0_1_basing_in_reads == False:
                # first_aligned_pos is 0-based
                arr_line = f'>{read.query_name} {first_aligned_pos}\n{cut_out_read}'
            else:
                arr_line = f'>{read.query_name} {first_aligned_pos+1}\n{cut_out_read}'

            arr.append(arr_line)
            arr_read_qualities_summary.append(np.asarray(cut_out_qualities))

        if read.reference_start >= counter and len(full_read) >= minimum_overlap:
            arr_read_summary.append(
                (read.query_name, read.reference_start + 1, read.reference_end, full_read)
            )

    counter = window_start + window_length

    return arr, arr_read_qualities_summary, arr_read_summary, counter


def build_windows(alignment_file: str, tiling_strategy: TilingStrategy,
    win_min_ext: float, maximum_reads: int, minimum_reads: int,
    reference_filename: str,
    exact_conformance_fix_0_1_basing_in_reads: Optional[bool] = False) -> None:
    """Summarizes reads aligned to reference into windows.
    Three products are created:
    #. Multiple FASTA files (one for each window position)
    #. A coverage file that lists all files in (1)
    #. A FASTA file that lists all reads used in (1)
        .. caution::
            ``reads.fas`` does not comply with the FASTA format.
    Args:
        alignment_file: Path to the alignment file in CRAM format.
        tiling_strategy: A strategy on how the genome is partitioned.
        win_min_ext: Minimum percentage of bases to overlap between reference
            and read to be considered in a window. The rest (i.e.
            non-overlapping part) will be filled with Ns.
        maximum_reads: Upper (exclusive) limit of reads allowed to start at the
            same position in the reference genome. Serves to reduce
            computational load.
        minimum_reads: Lower (exclusive) limit of reads allowed in a window.
            Serves to omit windows with low coverage.
        reference_filename: Path to a FASTA file of the reference sequence.
        exact_conformance_fix_0_1_basing_in_reads: Fixes an incorrect 0-basing
            of reads in the window file in the old C++ version. 1-basing is
            applied everywhere now. Set this flag to `False` only for exact
            conformance with the old version (in tests).
    """
    assert 0 <= win_min_ext <= 1

    pysam.index(alignment_file)
    samfile = pysam.AlignmentFile(
        alignment_file,
        "r", # auto-detect bam/cram (rc)
        reference_filename=reference_filename,
        threads=1
    )
    reffile = pysam.FastaFile(reference_filename)

    cov_arr = []
    reads = open("reads.fas", "w")
    counter = 0
    reference_name = tiling_strategy.get_reference_name()
    tiling = tiling_strategy.get_window_tilings()
    region_end = tiling_strategy.get_region_end()

    permitted_reads_per_location, indel_map = _calc_location_maximum_reads(
        samfile,
        reference_name,
        maximum_reads
    )

    for idx, (window_start, window_length) in enumerate(tiling):
        arr, arr_read_qualities_summary, arr_read_summary, counter = _run_one_window(
            samfile,
            window_start - 1, # make 0 based
            reference_name,
            window_length,
            math.floor(win_min_ext * window_length),
            dict(permitted_reads_per_location), # copys dict ("pass by value")
            counter,
            exact_conformance_fix_0_1_basing_in_reads,
            indel_map
        )

        window_end = window_start + window_length - 1
        file_name = f'w-{reference_name}-{window_start}-{window_end}'

        # TODO solution for backward conformance
        if len(tiling) > 1:
            end_extended_by_a_window = region_end + (tiling[1][0]-tiling[0][0])*3
        else:
            end_extended_by_a_window = region_end + window_length*3
        for read in arr_read_summary:
            if idx == len(tiling) - 1 and read[1] > end_extended_by_a_window:
                continue
            # TODO reads.fas not FASTA conform, +-0/1 mixed
            # TODO global end does not really make sense, only for conformance
            # read name, global start, global end, read start, read end, read
            reads.write(
                f'{read[0]}\t{tiling[0][0]-1}\t{end_extended_by_a_window}\t{read[1]}\t{read[2]}\t{read[3]}\n'
            )

        if (idx != len(tiling) - 1 # except last
            and len(arr) > 0) or len(tiling)==1: # suppress output if window empty

            _write_to_file(arr, file_name + '.reads.fas')
            with open(file_name + '.qualities.npy', 'wb') as f:
                np.save(f, np.asarray(arr_read_qualities_summary, dtype=np.int64), allow_pickle=True)

            _write_to_file([
                f'>{reference_name} {window_start}\n' +
                reffile.fetch(reference=reference_name, start=window_start-1, end=window_end)
            ], file_name + '.ref.fas')

            if len(arr) > minimum_reads:
                line = (
                    f'{file_name}.reads.fas\t{reference_name}\t{window_start}\t'
                    f'{window_end}\t{len(arr)}'
                )
                cov_arr.append(line)

    samfile.close()
    reads.close()

    _write_to_file(cov_arr, "coverage.txt")


if __name__ == "__main__":
    import argparse

    # Naming as in original C++ version
    parser = argparse.ArgumentParser(description='b2w')
    parser.add_argument('-w', '--window_length', nargs=1, type=int,
        help='window length', required=True)
    parser.add_argument('-i', '--incr', nargs=1, type=int, help='increment',
        required=True)
    parser.add_argument('-m', nargs=1, type=float, help='minimum overlap in percent',
        required=True)
    parser.add_argument('-x', nargs=1, type=int,
        help='max reads starting at a position', required=True)
    parser.add_argument('-c', nargs=1, type=int,
        help='coverage threshold. Omit windows with low coverage.',
        required=True)

    parser.add_argument('-d', nargs='?',
        help='drop SNVs that are adjacent to insertions/deletions (alternate behaviour).',
        const=True)

    parser.add_argument('alignment_file', metavar='ALG', type=str)
    parser.add_argument('reference_filename', metavar='REF', nargs='?', default=None, type=str)
    parser.add_argument('region', metavar='REG', type=str)


    args = parser.parse_args()

    if args.d != None:
        raise NotImplementedError('This argument was deprecated.')

    eqsts = EquispacedTilingStrategy(args.region, args.window_length[0], args.incr[0])

    build_windows(
        alignment_file = args.alignment_file,
        tiling_strategy = eqsts,
        minimum_overlap = args.m[0],
        maximum_reads = args.x[0], # 1e4 / window_length, TODO why divide?
        minimum_reads = args.c[0],
        reference_filename = args.reference_filename
    )
