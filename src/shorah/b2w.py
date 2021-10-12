import pysam

def _write_to_file(lines, file_name):
    with open(file_name, "w") as f:
        f.writelines("%s\n" % l for l in lines)

def _run_one_window(samfile, region_start, reference_name, window_length, 
        minimum_overlap, maximum_reads):
    arr = []
    arr_read_summary = []

    iter = samfile.fetch(
        reference_name, 
        region_start, 
        region_start + window_length # arg exclusive as per pysam convention
    ) 

    for read in iter:
        # For loop limited by maximum_reads 
        # TODO might not be random 
        #if read_idx > maximum_reads: 
        #    break

        first_aligned_pos = read.reference_start
        last_aligned_post = read.reference_end - 1 #reference_end is exclusive
        if (first_aligned_pos < region_start + window_length - minimum_overlap  
                and last_aligned_post >= region_start + minimum_overlap - 4): 
                # TODO justify 4
            
            # 0- vs 1-based correction
            start_cut_out = region_start - first_aligned_pos - 1

            end_cut_out = start_cut_out + window_length 
            print(f"{start_cut_out} {end_cut_out}")
            s = slice(max(0, start_cut_out), end_cut_out)
            cut_out_read = list(read.query_sequence)
            
            diff_counter = 0
            for idx, pair in enumerate(read.get_aligned_pairs()):
                if pair[0] == None:
                    cut_out_read.insert(idx - diff_counter, "-")
                if pair[1] == None:
                    cut_out_read.pop(idx - diff_counter)
                    diff_counter = diff_counter + 1
            
            cut_out_read = ("".join(cut_out_read))[s]

            # TODO justify 2
            k = (region_start + window_length) - last_aligned_post - 2 
            if k > 0:
                print(k)
                cut_out_read = cut_out_read + k * "N" 
            if start_cut_out < 0:
                cut_out_read = -start_cut_out * "N" + cut_out_read

            assert len(cut_out_read) == window_length, (
                "read unequal window size"
            )

            arr.append(
                f'>{read.query_name} {first_aligned_pos}\n{cut_out_read}'
            )
            arr_read_summary.append((
                read.query_name, 
                read.reference_start + 1, # conversion to 1-based
                read.reference_end, # already 1-based
                read.query_sequence
            ))

    return arr, arr_read_summary


def b2w(window_length: int, incr: int, minimum_overlap: int, maximum_reads: int, 
        minimum_reads: int) -> None:
    """
    Creates three products:
    (1) Multiple FASTA files (one for each window position)
    (2) A coverage file that lists all files in (1)
    (3) A FASTA file that lists all reads used in (1) #TODO not really FASTA

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
    reference_name = "HXB2" # TODO

    pysam.index(alignment_file)
    samfile = pysam.AlignmentFile(alignment_file, "rc")

    window_positions = range(
        incr-10, # TODO corrected start
        samfile.get_reference_length(reference_name), 
        incr
    )

    cov_arr = []
    for region_start in window_positions:
        arr = _run_one_window(
            samfile, 
            region_start, 
            reference_name, 
            window_length, 
            minimum_overlap,
            maximum_reads
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
        
    samfile.close()

    _write_to_file(cov_arr, "coverage.txt")


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