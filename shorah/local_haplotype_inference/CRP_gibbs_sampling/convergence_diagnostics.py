import pandas as pd
import numpy as np
import pypesto
from pypesto.sample.auto_correlation import autocorrelation_sokal


def compute_effective_sample_size(chain):
    """
    Calculate the effective sample size of the MCMC chains.
    ess: Estimate of the effective sample size of the MCMC chains.

    Adapted from: pyPESTO/pypesto/sample/diagnostics.py
    """
    # compute autocorrelation (An estimate of the integrated autocorrelation time of the MCMC chain.)
    # Estimate the integrated autocorrelation time of a MCMC chain using Sokal's adaptive truncated
    # periodogram estimator.
    autocorrelation = autocorrelation_sokal(chain)

    # Get length of the converged chain
    N = chain.shape[0]

    # Calculate effective sample size
    ess = N / (1 + autocorrelation)

    return ess

def compute_ess_max(chain):

    N=chain.shape[0] # number of samples
    burn_in_freq = [0.1, 0.2, 0.3, 0.4]
    ess_max=0
    ess_max_chain=1
    for freq in burn_in_freq:
        burn_in_idx = int(freq*N)
        reduced_chain = chain[burn_in_idx:]
        chain_ess = compute_effective_sample_size(reduced_chain)
        if np.min(chain_ess) > ess_max:
            ess_max=np.min(chain_ess)
            ess_max_chain= chain_ess
    return ess_max_chain

def check_convergence_ess_max(df_history,proportion_threshold=0.5):
    """
    Chain passed ess_max convergence test if ess_max >= N*proportion_threshold with N being the number of samples in the chain.
    """
    relevant_cols = ['theta','gamma'] + [col for col in df_history.columns if col.startswith('c_')]
    chain = df_history[relevant_cols].to_numpy()
    N=chain.shape[0]
    ess_max_chain=compute_ess_max(chain)

    if np.min(ess_max_chain)>N*proportion_threshold:
        return True, np.min(ess_max_chain)
    else:
        return False, np.min(ess_max_chain)


from pypesto.sample.geweke_test import calculate_zscore, burn_in_by_sequential_geweke

def check_convergence_geweke(df_history, current_iter):
    """
    Based on the burn-in index computed using sequential geweke convergence is decided
    (converged iff burn-in index < current_iter).
    """
    relevant_cols = ['theta','gamma'] + [col for col in df_history.columns if col.startswith('c_')]
    chain = df_history[relevant_cols].to_numpy()
    try:
        burn_in_idx=burn_in_by_sequential_geweke(chain)
    except Exception:
        return "geweke failed to compute"

    if burn_in_idx<current_iter:
        return True
    else:
        return False

import tensorflow_probability as tfp
tfd = tfp.distributions
import tensorflow as tf
tf.config.experimental.enable_tensor_float_32_execution(False)

def compute_gelman_rubin(list_chains):
    n_samples = np.min([chain.shape[0] for chain in list_chains])
    corr_list_chains = [chain[:n_samples, :]for chain in list_chains]
    n_parameters = list_chains[0].shape[1]
    n_chains = len(list_chains)

    chain_merge = np.asarray(corr_list_chains)
    chain_merge = chain_merge.reshape((n_samples, n_chains, n_parameters))
    chain_merge = np.asarray(chain_merge).astype('float32')

    rhat = tfp.mcmc.diagnostic.potential_scale_reduction(chain_merge, independent_chain_ndims=n_chains)

    return rhat

def check_convergence_gelman_rubin(list_df_history, threshold_rhat):
    relevant_cols = ['theta','gamma'] + [col for col in list_df_history[0].columns if col.startswith('c_')]

    list_chains =[]
    for df_temp in list_df_history:
        list_chains.append(df_temp[relevant_cols].to_numpy())

    rhat = compute_gelman_rubin(list_chains)

    if rhat < threshold_rhat:
        return True, rhat
    else:
        return False, rhat
