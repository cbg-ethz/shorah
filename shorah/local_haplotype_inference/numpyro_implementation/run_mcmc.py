#!/usr/bin/env python3
import numpyro
import numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS, DiscreteHMCGibbs, Predictive
from jax import random
import jax
import os
import sys
import jax.numpy as jnp
import logging
logging.basicConfig(filename='numpyro_run_mcmc.log', encoding='utf-8', level=logging.INFO)

numpyro.set_platform('cpu')
#print('jax-version ',jax.__version__) #0.2.3
#print('numpyro-version ', numpyro.__version__) #0.4.1
#print('jax.config.FLAGS.jax_backend_target ', jax.config.FLAGS.jax_backend_target) #local
#print('jax.lib.xla_bridge.get_backend().platform ',jax.lib.xla_bridge.get_backend().platform) #gpu

from . import preprocessing
from . import models_NumPyro

# post processing
import scipy
import numpy as np

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def convert_h_to_seq(h, alphabet):
    '''
    Convert numeric representation to DNA-seq representation.
    '''
    seq = [alphabet[int(h[k])] for k in range(h.shape[0])]
    return ''.join(seq)

def corrected_reads_to_fasta(posterior_samples, model, input_data, rng_key, list_read_ids, fname_output_corr, alphabet,last_x_samples=100):

    # TODO: Do we want to cut posterior samples?
    reads = input_data[1]

    # only consider the last {last_x_samples} samples to summarize results
    cluster_assignments = posterior_samples['cluster_assignments'][-last_x_samples:,:,:]
    haplotypes = posterior_samples['haplotypes'][-last_x_samples:,:,:]

    average_cluster_assignment = scipy.stats.mode(cluster_assignments, axis=0)[0][0]
    average_haplotypes = scipy.stats.mode(haplotypes, axis=0)[0][0]

    # computes empirical posterior from posterior_samples p(theta | data)
    posterior_predictive = Predictive(model, posterior_samples)
    posterior_predictions = posterior_predictive(rng_key, input_data=input_data)['obs'][-last_x_samples:,:,:]
    posterior = [(posterior_predictions[:,n,:]==reads[n][:]).all(-1).sum() / last_x_samples for n in range(len(list_read_ids))]

    records = []
    for n ,read_id in enumerate(list_read_ids):
        hap_seq = convert_h_to_seq(average_haplotypes[average_cluster_assignment[n]][0], alphabet)
        header = '|posterior=' + str(posterior[n])
        records.append(SeqRecord(Seq(hap_seq), id=read_id, description=header))

    SeqIO.write(records, fname_output_corr, "fasta")

def haplotypes_to_fasta(posterior_samples, model, input_data, rng_key, fname_output_corr, alphabet, last_x_samples=100):
    genome_length = posterior_samples['haplotypes'].shape[2]
    alphabet_length = len(alphabet)

    # only consider the last {last_x_samples} samples to summarize results
    cluster_assignments = posterior_samples['cluster_assignments'][-last_x_samples:,:,:]
    haplotypes = posterior_samples['haplotypes'][-last_x_samples:,:,:]

    # number of reads assigned to each haplo
    hap_ids, ave_reads = np.unique(cluster_assignments, return_counts=True)

    # collapse haplotypes into unique set
    average_haplotypes = scipy.stats.mode(haplotypes, axis=0)[0][0]
    unique_haplotypes = np.unique(average_haplotypes, axis=0)
    idx_collapsed_haplotypes = [np.where(np.all(average_haplotypes==unique_hap,axis=1)) for unique_hap in unique_haplotypes]
    idx_collapsed_haplotypes = [np.intersect1d(idx_col,hap_ids) for idx_col in idx_collapsed_haplotypes]
    # map idx_collapsed_haplotypes to hap_ids order
    outer = []
    for idx_collapsed in idx_collapsed_haplotypes:
        outer.append([np.argwhere(hap_ids==i)[0][0] for i in idx_collapsed])
    mapped_idx_collapsed_haplotypes = outer

    # collapse number of average reads to unique haplotypes
    ave_reads = np.array([np.sum(ave_reads[idx_hap]) for idx_hap in mapped_idx_collapsed_haplotypes])
    ave_reads = ave_reads / last_x_samples

    # empirical posteriors of haplotypes (collapsed like in ave_reads computation)
    posterior_predictive = Predictive(model,
                                      posterior_samples,
                                      return_sites=['haplotypes']
                                     )
    posterior_predictions = posterior_predictive(rng_key,
                                                 input_data=input_data
                                                )['haplotypes'][-last_x_samples:,:,:]
    posterior = [(posterior_predictions[:,hap_ids[idx_hap],:]==average_haplotypes[hap_ids[idx_hap]][:]).all(-1).sum() / last_x_samples for idx_hap in mapped_idx_collapsed_haplotypes]

    # write to fasta
    records = []
    for k in range(unique_haplotypes.shape[0]):
        head = ' | posterior='+str(posterior[k])+' ave_reads='+str(ave_reads[k])
        seq = convert_h_to_seq(unique_haplotypes[k], alphabet)
        records.append(SeqRecord(Seq(seq), id='haplotype'+str(k), description=head))
    SeqIO.write(records, fname_output_corr, "fasta")


def main(freads_in, fref_in, output_dir, alpha0, cluster_num, max_num_samples, str_model, alphabet, thres_ess_max, thres_rhat):
    # alphabet = 'ACGT-'

    if os.path.exists(output_dir)==False: # Check whether the specified path exists or not
        os.makedirs(output_dir)

    window_id = freads_in.split('/')[-1][:-4] # freads_in is absolute path
    window=[int(window_id.split('-')[2])-1,int(window_id.split('-')[3].split('.')[0])]
    output_name = output_dir+window_id+'-'

    logging.info('Input reads: ' + str(freads_in))
    logging.info('Inference model: '+ str(str_model))
    logging.info('Start main ')

    if str_model == "infiniteSBP":
        model = models_NumPyro.model_infiniteSBP
    elif str_model == "finiteDPM":
        model = models_NumPyro.model_finiteDPM
    elif str_model == "finiteDPM_extended":
        model = models_NumPyro.model_finiteDPM_extended
    elif str_model == "infiniteSBP_extended":
        model = models_NumPyro.model_infiniteSBP_extended
    else:
        raise ValueError("Input model is not defined.")

    alphabet_length = len(alphabet) # size alphabet

    reference, ref_id = preprocessing.fasta2ref(fref_in, alphabet)
    reads, list_read_ids = preprocessing.fasta2reads(freads_in, alphabet)

    # reference shorter than the windows
    if len(reference) < window[1]-window[0]:
        missing_pos = window[1]-window[0] - len(reference)
        # shorten reads to length of reference --> we don't want nan in ref
        reads = reads[:, :-missing_pos]

    # all positions in reads that are != N
    is_observed = jnp.invert(jnp.isnan(reads))

    if str_model in ["infiniteSBP"]:
        input_data = reference, reads, alphabet_length, cluster_num, is_observed
    elif str_model in ["finiteDPM"]:
        input_data = reference, reads, alphabet_length, cluster_num, is_observed, alpha0
    elif str_model in ["finiteDPM_extended", "infiniteSBP_extended"]:
        pos_del_ref =  jnp.argwhere(reference==4)
        if str_model in ["finiteDPM_extended"]:
            input_data = reference, reads, alphabet_length, cluster_num, is_observed, pos_del_ref, alpha0
        elif str_model in ["infiniteSBP_extended"]:
            input_data = reference, reads, alphabet_length, cluster_num, is_observed, pos_del_ref

    #numpyro.set_host_device_count(n_cores)
    rng_key = jax.random.PRNGKey(0)
    num_samples = max_num_samples
    num_warmup = 10 #int(num_samples / 2)
    logging.info('Start sampling ')
    # Run NUTS
    kernel = NUTS(model)
    mcmc = MCMC(
        DiscreteHMCGibbs(kernel),
        num_warmup=num_warmup,
        num_samples=num_samples,
        num_chains=1,
        progress_bar=True,
        chain_method="parallel"
    )
    mcmc.run(rng_key, input_data)
    logging.info('Finished sampling.')
    #logging.info(mcmc.print_summary())

    posterior_samples = mcmc.get_samples()

    haplotypes_to_fasta(posterior_samples, model, input_data, rng_key, output_name+'support.fas' , alphabet, last_x_samples=100)
    logging.info('haplo to fasta done')
    corrected_reads_to_fasta(posterior_samples, model, input_data, rng_key, list_read_ids, output_name+'cor.fas', alphabet,last_x_samples=100)
    logging.info('reads to fasta')

    """
    Check convergence of the chain every x samples.
    1. Check convergence by looking at Gelman Rubin and ESS
    Gelman Rubin: https://arxiv.org/pdf/1812.09384.pdf
    \hat{R} <1.1 or 1.01
    2. Continue running the chain using  post_warmup_state
    Example:
    mcmc = MCMC(NUTS(model), num_warmup=100, num_samples=100)
    mcmc.run(random.PRNGKey(0))
    first_100_samples = mcmc.get_samples()
    mcmc.post_warmup_state = mcmc.last_state
    mcmc.run(mcmc.post_warmup_state.rng_key)  # or mcmc.run(random.PRNGKey(1))
    second_100_samples = mcmc.get_samples()
    """


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8] )
