import pysam
from typing import Optional
from shorah.tiling import TilingStrategy, EquispacedTilingStrategy
import numpy as np
import math

def _write_to_file(lines, file_name):
    with open(file_name, "w") as f:
        f.writelines("%s\n" % l for l in lines)

def _calc_via_pileup(samfile, reference_name, maximum_reads):
    budget = dict()
    max_ins_at_pos = dict()
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

        max_at_this_pos = 0
        for pileupread in pileupcolumn.pileups:
            if pileupread.indel > 0 or pileupread.is_del:
                indel_map.add((
                    pileupread.alignment.query_name, # TODO is unique?
                    pileupread.alignment.reference_start, # TODO is unique?
                    pileupcolumn.reference_pos,
                    pileupread.indel,
                    pileupread.is_del
                ))
            if pileupread.indel > max_at_this_pos:
                max_at_this_pos = pileupread.indel

        if max_at_this_pos > 0:
            max_ins_at_pos[pileupcolumn.reference_pos] = max_at_this_pos

    # ascending reference_pos are necessary for later steps
    indel_map = sorted(indel_map, key=lambda tup: tup[2])

    print(max_ins_at_pos)

    return budget, indel_map, max_ins_at_pos

def _build_one_full_read(full_read: list[str], full_qualities: list[int]|list[str],
    read_query_name: str|None, first_aligned_pos, last_aligned_pos,
    window_start, indel_map, max_ins_at_pos,
    extended_window_mode) -> tuple[str, list[int]]:

    all_inserts = dict()
    own_inserts = set()

    change_in_reference_space_ins = 0

    for name, start, ref_pos, indel_len, is_del in indel_map:
        if name == read_query_name and start == first_aligned_pos:
            if is_del == 1: # if del
                if indel_len != 0:
                    raise NotImplementedError("Deletions larger than 1 not expected.")
                full_read.insert(ref_pos - first_aligned_pos + change_in_reference_space_ins, "-")
                full_qualities.insert(ref_pos - first_aligned_pos + change_in_reference_space_ins, "2")
                continue

            elif is_del == 0 and not extended_window_mode:
                assert indel_len > 0
                for _ in range(indel_len):
                    full_read.pop(ref_pos + 1 - first_aligned_pos)
                    full_qualities.pop(start + 1 - first_aligned_pos)
                continue

            elif is_del == 0 and extended_window_mode and ref_pos >= window_start:
                own_inserts.add((ref_pos, indel_len))
                change_in_reference_space_ins += indel_len
                all_inserts[ref_pos] = max_ins_at_pos[ref_pos]

            else:
                pass


        if (extended_window_mode and
            (name != read_query_name or start != first_aligned_pos) and
            window_start <= ref_pos < last_aligned_pos and is_del == 0 and
            first_aligned_pos <= ref_pos): # TODO edge values left
            all_inserts[ref_pos] = max_ins_at_pos[ref_pos]

    if extended_window_mode:
        change_in_reference_space = 0
        own_inserts_pos = []
        own_inserts_len = []
        if len(own_inserts) != 0:
            [own_inserts_pos, own_inserts_len] = [list(t) for t in zip(*own_inserts)]

        for pos in sorted(all_inserts):
            n = all_inserts[pos]
            if (pos, n) in own_inserts:
                change_in_reference_space += n
                continue

            L = max_ins_at_pos[pos]
            in_idx = pos + 1 - first_aligned_pos + change_in_reference_space
            if pos in own_inserts_pos:
                k = own_inserts_len[own_inserts_pos.index(pos)]
                L -= k
                in_idx += k
            for _ in range(L):
                full_read.insert(in_idx, "-")
                full_qualities.insert(in_idx, "2")

            change_in_reference_space += max_ins_at_pos[pos]

    full_read = ("".join(full_read))

    return full_read, full_qualities # TODO return same data type twice


def _run_one_window(samfile, window_start, reference_name, window_length,
        minimum_overlap, permitted_reads_per_location, counter,
        exact_conformance_fix_0_1_basing_in_reads, indel_map, max_ins_at_pos,
        extended_window_mode):

    arr = []
    arr_read_summary = []
    arr_read_qualities_summary = []

    iter = samfile.fetch(
        reference_name,
        window_start, # 0 based
        window_start + window_length # arg exclusive as per pysam convention
    )

    original_window_length = window_length
    if extended_window_mode:
        for pos, val in max_ins_at_pos.items():
            if window_start <= pos < window_start + original_window_length:
                window_length += val

    for read in iter:

        if (read.reference_start is None) or (read.reference_end is None):
            continue
        first_aligned_pos = read.reference_start # this is 0-based
        last_aligned_pos = read.reference_end - 1 #reference_end is exclusive


        if permitted_reads_per_location[first_aligned_pos] == 0:
            continue
        else:
            permitted_reads_per_location[first_aligned_pos] -= 1

        full_read = list(read.query_sequence)
        full_qualities = list(read.query_qualities)

        for ct_idx, ct in enumerate(read.cigartuples):
            if ct[0] in [0,1,2,7,8]: # 0 = BAM_CMATCH, 1 = BAM_CINS, 2 = BAM_CDEL, 7 = BAM_CEQUAL, 8 = BAM_CDIFF
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

        full_read, full_qualities = _build_one_full_read(full_read, full_qualities,
            read.query_name, first_aligned_pos, last_aligned_pos, window_start,
            indel_map, max_ins_at_pos, extended_window_mode)

        if (first_aligned_pos < window_start + 1 + window_length - minimum_overlap
                and last_aligned_pos >= window_start + minimum_overlap - 2 # TODO justify 2
                and len(full_read) >= minimum_overlap):

            num_inserts_right_of_read = 0
            num_inserts_left_of_read = 0
            if extended_window_mode:
                for pos, val in max_ins_at_pos.items():
                    if last_aligned_pos <= pos < window_start + original_window_length:
                        num_inserts_right_of_read += val
                    if window_start <= pos < first_aligned_pos:
                        num_inserts_left_of_read += val # TODO no tests

            start_cut_out = window_start - first_aligned_pos
            end_cut_out = start_cut_out + window_length - num_inserts_left_of_read
            s = slice(max(0, start_cut_out), end_cut_out)

            cut_out_read = full_read[s]
            cut_out_qualities = full_qualities[s]

            k = (window_start + original_window_length - 1 - last_aligned_pos
                + num_inserts_right_of_read)

            if k > 0:
                cut_out_read = cut_out_read + k * "N"
                cut_out_qualities = cut_out_qualities + k * [2]
                # Phred scores have a minimal value of 2, where an “N” is inserted
                # https://www.zymoresearch.com/blogs/blog/what-are-phred-scores
            if start_cut_out < 0:
                cut_out_read = (-start_cut_out + num_inserts_left_of_read) * "N" + cut_out_read
                cut_out_qualities = (-start_cut_out + num_inserts_left_of_read) * [2] + cut_out_qualities

            assert len(cut_out_read) == window_length, (
                "read unequal window size",
                read.query_name, first_aligned_pos, cut_out_read, window_start, window_length
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
    exact_conformance_fix_0_1_basing_in_reads: Optional[bool] = False,
    extended_window_mode: Optional[bool] = False) -> None:
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
        extended_window_mode: Mode where inserts are not deleted but kept. The
            windows are instead extended.
    """
    assert 0 <= win_min_ext <= 1
    extended_window_mode = True # TODO

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

    permitted_reads_per_location, indel_map, max_ins_at_pos = _calc_via_pileup(
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
            indel_map,
            max_ins_at_pos,
            extended_window_mode
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
            and len(arr) > 0) or len(tiling) == 1: # suppress output if window empty

            _write_to_file(arr, file_name + '.reads.fas')
            with open(file_name + '.qualities.npy', 'wb') as f:
                np.save(f, np.asarray(arr_read_qualities_summary, dtype=np.int64), allow_pickle=True)

            ref = reffile.fetch(reference=reference_name, start=window_start-1, end=window_end)
            _write_to_file([
                f'>{reference_name} {window_start}\n' + ref
            ], file_name + '.ref.fas')

            if extended_window_mode:
                print("HERE")
                # _write_to_file([
                #     f'>{reference_name} {window_start}\n' +
                #     _build_one_full_read(
                #         list(ref), list(ref), None,
                #         window_start, window_start + window_length - 1, window_start,
                #         indel_map, max_ins_at_pos, extended_window_mode)[0]
                # ], file_name + '.extended-ref.fas')

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
