#!/usr/bin/env python3

import pandas as pd
import skbio
import numpy as np
import sys
import os
import pickle
from timeit import default_timer as timer
import logging
logging.basicConfig(filename='shorah_inference.log', encoding='utf-8', level=logging.INFO)

# my python-scripts
from . import preparation
from . import update_eqs as update_eqs
from . import initialization
from . import analyze_results
from . import elbo_eqs
from . import cavi

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

def main(freads_in, fref_in, fname_qualities, output_dir, n_starts, K, alpha0, alphabet = 'ACGT-'):
    start_time = timer()

    window_id = freads_in.split('/')[-1][:-4] # freads_in is absolute path
    window=[int(window_id.split('-')[2])-1,int(window_id.split('-')[3].split('.')[0])]
    dict_runtime={'window_id': window_id}

    output_name = output_dir+window_id+'-'

    if os.path.exists(output_dir)==False: # Check whether the specified path exists or not
        os.makedirs(output_dir) # Create a new directory because it does not exist

    # Read in reads
    reference_seq, ref_id = preparation.load_reference_seq(fref_in)
    reference_binary = preparation.reference2binary(reference_seq, alphabet)
    if freads_in.endswith('fasta') or freads_in.endswith('fas'):
        reads_list= preparation.load_fasta2reads_list(freads_in, alphabet)
    elif freads_in.endswith('bam'):
        reads_list= preparation.load_bam2reads_list(freads_in, alphabet)

    reads_seq_binary, reads_weights = preparation.reads_list_to_array(reads_list)
    qualities = preparation.get_average_qualities(fname_qualities, reads_list)
    reads_log_error_proba = preparation.get_reads_log_error_proba(qualities, reads_seq_binary, len(alphabet))

    end_time_init = timer()
    dict_runtime.update({'time_preparation': end_time_init-start_time})
    dict_runtime.update({'n_starts': n_starts})

    result_list = cavi.multistart_cavi(K, alpha0, alphabet, reference_binary, reference_seq, reads_list, reads_seq_binary, reads_weights, reads_log_error_proba, n_starts, output_name)

    end_time_optimization = timer()
    dict_runtime.update({'time_optimization': end_time_optimization-end_time_init})

    logging.info('reference '+ fref_in)
    logging.info('reads '+ freads_in)
    logging.info('lenght of sequences '+ str(reads_list[0].seq_binary.shape[0]))
    logging.info('number of reads '+ str(len(reads_list)))

    # Find best run
    sort_elbo = [(idx,state_run[1]['elbo']) for idx, state_run in enumerate(result_list)]
    sort_elbo.sort(key=lambda x:x[1], reverse=True) # sort list of tuple by ELBO value

    max_idx = sort_elbo[0][0]
    max_elbo = sort_elbo[0][1]
    sort_results = [result_list[tuple_idx_elbo[0]] for tuple_idx_elbo in sort_elbo]


    f2 = open(output_name+'all_results.pkl',"wb")
    pickle.dump(sort_results,f2)
    f2.close()
    logging.info("Results dicts of all runs written to " + output_name+'all_results.pkl')

    state_curr_dict = result_list[max_idx][1]
    logging.info('Maximal ELBO '+str(max_elbo) + 'in run '+ str(max_idx))
    # TODO maybe write this to logging
    # outfile = analyze_results.write_info2file(state_curr_dict, outfile, alphabet, reads_list, reads_seq_binary, reads_weights,reference_binary, reference_seq)

    # write output like in original shorah
    analyze_results.haplotypes_to_fasta(state_curr_dict, output_name+'support.fas')
    analyze_results.correct_reads(state_curr_dict, output_name+'cor.fas')

    #f_best_run = open(output_dir+'debug/best_run.txt','w')
    f_best_run = open(output_name+'best_run.txt','w')
    f_best_run.write(str(max_idx))
    f_best_run.close()


    end_time_total = timer()
    dict_runtime.update({'time_total': end_time_total-start_time})

    f = open(output_name+'runtime.pkl',"wb")
    pickle.dump(dict_runtime,f)
    f.close()

    # clean up Files
    if os.path.exists(output_dir+'inference/')==False:
        os.makedirs(output_dir+'inference/')

    import glob
    import shutil

    inference_files = glob.glob('./w*best_run.txt') + glob.glob('./w*history_run*.csv') + glob.glob('./w*runtime.pkl') + glob.glob('./w*results*.pkl')

    for inf_file in inference_files:
        if os.stat(inf_file).st_size > 0:
            gzf = gzip_file(inf_file)
            try:
                os.remove('inference/%s' % gzf)
            except OSError:
                pass
            shutil.move(gzf, 'inference/')
        else:
            os.remove(inf_file)

    logging.info('Files cleaned up.')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6]), sys.argv[7])
#freads_in, fref_in, output_dir, n_starts, K, alpha0, alphabet = 'ACGT-'
