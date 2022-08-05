#!/usr/bin/env python3

import sys
import os
import pickle
import logging

# my python-scripts
from . import preparation
from . import analyze_results
from . import cavi

logging.basicConfig(
    filename="shorah_inference.log", encoding="utf-8", level=logging.INFO
)


def gzip_file(f_name):
    """Gzip a file and return the name of the gzipped, removing the original"""
    import gzip

    f_in = open(f_name, "rb")
    f_out = gzip.open(f_name + ".gz", "wb")
    f_out.writelines(f_in)
    f_out.close()
    f_in.close()
    os.remove(f_in.name)

    return f_out.name


def main(
    freads_in,
    fref_in,
    fname_qualities,
    output_dir,
    n_starts,
    K,
    alpha0,
    alphabet="ACGT-",
):

    window_id = freads_in.split("/")[-1][:-4]  # freads_in is absolute path
    output_name = output_dir + window_id + "-"

    os.makedirs(output_dir, exist_ok=True)

    # Read in reads
    reference_binary, ref_id = preparation.load_reference_seq(fref_in, alphabet)

    reads_list, qualities = preparation.load_fasta_and_qualities(
        freads_in, fname_qualities, alphabet
    )
    reads_seq_binary, reads_weights = preparation.reads_list_to_array(reads_list)
    reads_log_error_proba = preparation.compute_reads_log_error_proba(
        qualities, reads_seq_binary, len(alphabet)
    )

    result_list = [
        cavi.run_cavi(
            K,
            alpha0,
            alphabet,
            reference_binary,
            reads_list,
            reads_seq_binary,
            reads_weights,
            reads_log_error_proba,
            0,
            output_name,
        )
    ]

    logging.info("reference " + fref_in)
    logging.info("reads " + freads_in)
    logging.info("lenght of sequences " + str(reads_list[0].seq_binary.shape[0]))
    logging.info("number of reads " + str(len(reads_list)))

    # Find best run
    sort_elbo = [
        (idx, state_run[1]["elbo"]) for idx, state_run in enumerate(result_list)
    ]
    # sort list of tuple by ELBO value
    sort_elbo.sort(key=lambda x: x[1], reverse=True)

    best_run_idx = sort_elbo[0][0]
    best_run_elbo = sort_elbo[0][1]
    logging.info("Maximal ELBO " + str(best_run_elbo) + "in run " + str(best_run_idx))

    sorted_results = [result_list[tuple_idx_elbo[0]] for tuple_idx_elbo in sort_elbo]
    exit_meassage = sorted_results[0][1]["exit_message"]
    logging.info("CAVI termination " + str(exit_meassage))

    with open(output_name + "all_results.pkl", "wb") as f:
        pickle.dump(sorted_results,f)

    # TODO: Would be nicer to use json dump.
    #import json
    #with open(output_name + "all_results.pkl", "wb") as fp:
    #    json.dump(sorted_results, fp)

    logging.info(
        "Results dicts of all runs written to " + output_name + "all_results.pkl"
    )

    state_curr_dict = result_list[best_run_idx][0]
    summary = analyze_results.summarize_results(
        state_curr_dict,
        alphabet,
        reads_seq_binary,
        reads_weights,
        reads_list,
        reads_log_error_proba,
        reference_binary,
    )
    state_curr_dict.update(summary)

    # write output like in original shorah
    analyze_results.haplotypes_to_fasta(state_curr_dict, output_name + "support.fas")
    analyze_results.correct_reads(state_curr_dict, output_name + "cor.fas")

    # f_best_run = open(output_name+'best_run.txt','w')
    # f_best_run.write(str(max_idx))
    # f_best_run.close()

    # clean up Files
    os.makedirs(output_dir + "inference/", exist_ok=True)

    import glob
    import shutil

    inference_files = glob.glob("./w*results*.pkl")

    for inf_file in inference_files:
        if os.stat(inf_file).st_size > 0:
            gzf = gzip_file(inf_file)
            try:
                os.remove("inference/%s" % gzf)
            except OSError:
                pass
            shutil.move(gzf, "inference/")
        else:
            try:
                os.remove(inf_file)
            except OSError:
                pass

    logging.info("Files cleaned up.")


if __name__ == "__main__":
    main(
        sys.argv[1],
        sys.argv[2],
        sys.argv[3],
        int(sys.argv[4]),
        int(sys.argv[5]),
        float(sys.argv[6]),
        sys.argv[7],
    )
# freads_in, fref_in, output_dir, n_starts, K, alpha0, alphabet = 'ACGT-'
