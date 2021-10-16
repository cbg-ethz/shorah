import pysam

def _parse_region(region):
    tmp = region.split(":")
    reference_name = tmp[0]
    tmp = tmp[1].split("-")
    start = int(tmp[0]) 
    end = int(tmp[1]) 
    return reference_name, start, end # indexed 1 like samtools

def _write_to_file(lines, file_name):
    with open(file_name, "w") as f:
        f.writelines("%s\n" % l for l in lines)

def _run_one_window(samfile, region_start, reference_name, window_length, 
        minimum_overlap, maximum_reads, counter):

    arr = []
    arr_read_summary = []

    iter = samfile.fetch(
        reference_name, 
        region_start, 
        region_start + window_length # arg exclusive as per pysam convention
    ) 

    for idx, read in enumerate(iter):
        # For loop limited by maximum_reads 
        # TODO might not be random 
        #if read_idx > maximum_reads: 
        #    break

        first_aligned_pos = read.reference_start
        last_aligned_post = read.reference_end - 1 #reference_end is exclusive

        # 0- vs 1-based correction
        start_cut_out = region_start - first_aligned_pos - 1

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
        
        if (first_aligned_pos < region_start + window_length - minimum_overlap  
                and last_aligned_post >= region_start + minimum_overlap - 4): 
                # TODO justify 4

            cut_out_read = full_read[s]

            # TODO justify 2
            k = (region_start + window_length) - last_aligned_post - 2 
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
            arr_read_summary.append( # TODO reads.fas not FASTA conform, +-0/1
                f'{read.query_name}\t2267\t3914\t{read.reference_start + 1}\t{read.reference_end}\t{full_read}'
            )
            counter = read.reference_start

    counter = counter + 1
    print(f'GLOBAL: {counter}')

    return arr, arr_read_summary, counter


def b2w(window_length: int, incr: int, minimum_overlap: int, maximum_reads: int, 
        minimum_reads: int) -> None:
    """Summarizes reads aligned to reference into windows. 

    Three products are created:

    #. Multiple FASTA files (one for each window position)
    #. A coverage file that lists all files in (1)
    #. A FASTA file that lists all reads used in (1) #TODO not really FASTA

    Args:
        window_length: Number of bases considered at once per loop.
        incr: Increment between each window.
        minimum_overlap: Minimum number of bases to overlap between reference
            and read to be considered in a window.
        maximum_reads: Upper (inclusive) limit of reads allowed in a window.
            Serves to reduce computational load.
        minimum_reads: Lower (inclusive) limit of reads allowed in a window.
            Serves to omit windows with low coverage.
        
    Returns:
        None.
    """
    alignment_file = "data/test_aln.cram" # TODO
    region = "HXB2:2469-3713"
    reference_name, start, end = _parse_region(region)

    pysam.index(alignment_file)
    samfile = pysam.AlignmentFile(alignment_file, "rc")

    window_positions = range(
        start - window_length, # TODO corrected start
        end + window_length, 
        incr 
    )

    print(window_positions)

    cov_arr = []
    arr_read_summary_all = []
    counter = 0
    for region_start in window_positions:
        arr, arr_read_summary, counter = _run_one_window(
            samfile, 
            region_start, 
            reference_name, 
            window_length, 
            minimum_overlap,
            maximum_reads,
            counter
        )
        region_end = region_start + window_length - 1
        file_name = f'w-{reference_name}-{region_start}-{region_end}.reads.fas'
        if len(arr) >= max(minimum_reads, 1):

            # TODO write to file earlier to free up memory
            _write_to_file(arr, file_name) 

            line = (
                f'{file_name}\t{reference_name}\t{region_start}\t'
                f'{region_end}\t{len(arr)}'
            )
            cov_arr.append(line)

            arr_read_summary_all.extend(arr_read_summary)
        
    samfile.close()

    _write_to_file(cov_arr, "coverage.txt")
    _write_to_file(arr_read_summary_all, "reads.fas")


def main():
    """
        -w: window length (INT)
        -i: increment (INT)
        -m: minimum overlap (INT)
        -x: max reads starting at a position (INT)
        -c: coverage threshold. Omit windows with low coverage (INT)

        -d: drop SNVs that are adjacent to insertions/deletions 
            (alternate behaviour)
        -h: show this help
    """

    #if args.d == True:
    #    raise NotImplementedError # TODO

    window_length = 201
    incr = window_length//3

    # if read covers at least win_min_ext fraction of the window, 
    # fill it with Ns
    win_min_ext = 0.85 

    b2w(
        window_length, 
        incr, 
        window_length * win_min_ext, 
        1e4 / window_length, # TODO why divide?
        0
    )

    print("main")