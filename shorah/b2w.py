import pysam
from typing import Optional
from shorah.tiling import TilingStrategy, EquispacedTilingStrategy

def _write_to_file(lines, file_name):
    with open(file_name, "w") as f:
        f.writelines("%s\n" % l for l in lines)

def _calc_location_maximum_reads(samfile, reference_name, maximum_reads):
    budget = dict()
    for pileupcolumn in samfile.pileup(reference_name, multiple_iterators=False):
        budget[pileupcolumn.reference_pos] = min(
            pileupcolumn.nsegments, 
            maximum_reads-1 # minus 1 because because maximum_reads is exclusive
        )

    return budget

def _run_one_window(samfile, window_start, reference_name, window_length, 
        minimum_overlap, permitted_reads_per_location, counter):

    arr = []
    arr_read_summary = []

    iter = samfile.fetch(
        reference_name, 
        window_start, 
        window_start + window_length # arg exclusive as per pysam convention
    ) 

    for idx, read in enumerate(iter):

        first_aligned_pos = read.reference_start
        last_aligned_post = read.reference_end - 1 #reference_end is exclusive


        if permitted_reads_per_location[first_aligned_pos] == 0:
            continue
        else:
            permitted_reads_per_location[first_aligned_pos] = (
                permitted_reads_per_location[first_aligned_pos] - 1
            )

        # 0- vs 1-based correction
        start_cut_out = window_start - first_aligned_pos - 1

        end_cut_out = start_cut_out + window_length 

        s = slice(max(0, start_cut_out), end_cut_out)
        full_read = list(read.query_sequence)
        
        diff_counter = 0
        for idx, pair in enumerate(read.get_aligned_pairs()):
            if pair[0] == None:
                full_read.insert(idx - diff_counter, "-")
            if pair[1] == None:
                full_read.pop(idx - diff_counter)
                diff_counter = diff_counter + 1
        
        full_read = ("".join(full_read))
        
        if (first_aligned_pos < window_start + window_length - minimum_overlap
                and last_aligned_post >= window_start + minimum_overlap - 3 # TODO justify 3
                and len(full_read) >= minimum_overlap): 

            cut_out_read = full_read[s]

            # TODO justify 2
            k = (window_start + window_length) - last_aligned_post - 2 
            if k > 0:
                cut_out_read = cut_out_read + k * "N" 
            if start_cut_out < 0:
                cut_out_read = -start_cut_out * "N" + cut_out_read

            assert len(cut_out_read) == window_length, (
                "read unequal window size"
            )

            arr.append(
                f'>{read.query_name} {first_aligned_pos}\n{cut_out_read}'
            )

        if read.reference_start >= counter and len(full_read) >= minimum_overlap:
            arr_read_summary.append(
                (read.query_name, read.reference_start + 1, read.reference_end, full_read)
            )

    counter = window_start + window_length

    return arr, arr_read_summary, counter


def build_windows(alignment_file: str, tiling_strategy: TilingStrategy, 
    minimum_overlap: int, maximum_reads: int, minimum_reads: int, 
    reference_filename: Optional[str] = None) -> None:
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
        minimum_overlap: Minimum number of bases to overlap between reference
            and read to be considered in a window. The rest (i.e. 
            non-overlapping part) will be filled with Ns.
        maximum_reads: Upper (exclusive) limit of reads allowed to start at the
            same position in the reference genome. Serves to reduce 
            computational load.
        minimum_reads: Lower (exclusive) limit of reads allowed in a window.
            Serves to omit windows with low coverage.
        reference_filename: Path to a FASTA file of the reference sequence.
            Only necessary if this information is not included in the CRAM file.
    """

    pysam.index(alignment_file)
    samfile = pysam.AlignmentFile(
        alignment_file, 
        "r", # auto-detect bam/cram (rc)
        reference_filename=reference_filename,
        threads=1
    )

    cov_arr = []
    reads = open("reads.fas", "w")
    counter = 0
    reference_name = tiling_strategy.get_reference_name()
    tiling = tiling_strategy.get_window_tilings()
    region_end = tiling_strategy.get_region_end()

    permitted_reads_per_location = _calc_location_maximum_reads(
        samfile, 
        reference_name, 
        maximum_reads
    )

    for idx, (window_start, window_length) in enumerate(tiling):
        arr, arr_read_summary, counter = _run_one_window(
            samfile, 
            window_start, 
            reference_name, 
            window_length, 
            minimum_overlap,
            dict(permitted_reads_per_location), # copys dict ("pass by value")
            counter
        )
        window_end = window_start + window_length - 1
        file_name = f'w-{reference_name}-{window_start}-{window_end}.reads.fas'

        # TODO solution for backward conformance
        end_extended_by_a_window = region_end + (tiling[1][0]-tiling[0][0])*3
        for read in arr_read_summary:
            if idx == len(tiling) - 1 and read[1] > end_extended_by_a_window:
                continue
            # TODO reads.fas not FASTA conform, +-0/1 mixed
            # TODO global end does not really make sense, only for conformance
            # read name, global start, global end, read start, read end, read
            reads.write(
                f'{read[0]}\t{tiling[0][0]-1}\t{end_extended_by_a_window}\t{read[1]}\t{read[2]}\t{read[3]}\n'
            )
        
        if idx != len(tiling) - 1: # except last

            _write_to_file(arr, file_name) 

            if len(arr) > minimum_reads: 
                line = (
                    f'{file_name}\t{reference_name}\t{window_start}\t'
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
    parser.add_argument('-m', nargs=1, type=int, help='minimum overlap', 
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