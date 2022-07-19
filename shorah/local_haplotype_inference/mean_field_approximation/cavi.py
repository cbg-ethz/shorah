import numpy as np
import pickle
import pandas as pd
import multiprocessing as mp
from scipy.special import digamma
from scipy.stats._multivariate import _lnB as lnB
from scipy.special import betaln

# my python scripts
from . import initialization
from . import update_eqs
from . import elbo_eqs
from . import analyze_results

"""
Parallelizing with Pool following:
https://www.machinelearningplus.com/python/parallel-processing-python/
"""

results = []


def collect_result(result):
    global results
    results.append(result)


def multistart_cavi(
    n_cluster,
    alpha0,
    alphabet,
    reference_binary,
    reads_list,
    reads_seq_binary,
    reads_weights,
    reads_log_error_proba,
    n_starts,
    output_dir,
):

    pool = mp.Pool(mp.cpu_count())
    for start in range(n_starts):
        pool.apply_async(
            run_cavi,
            args=(
                n_cluster,
                alpha0,
                alphabet,
                reference_binary,
                reads_list,
                reads_seq_binary,
                reads_weights,
                reads_log_error_proba,
                start,
                output_dir,
            ),
            callback=collect_result,
        )

    pool.close()
    pool.join()

    return results


def run_cavi(
    n_cluster,
    alpha0,
    alphabet,
    reference_binary,
    reads_list,
    reads_seq_binary,
    reads_weights,
    reads_log_error_proba,
    start_id,
    output_dir,
):
    """
    Runs cavi (coordinate ascent variational inference).
    """
    dict_result = {
        "run_id": start_id,
        "n_reads": len(reads_list),
        "n_cluster": n_cluster,
        "alpha0": alpha0,
        "alphabet": alphabet,
    }

    history_alpha = []
    history_mean_log_pi = []
    history_mean_log_gamma = []
    history_mean_haplo = []
    history_mean_cluster = []
    history_elbo = []

    state_init_dict = initialization.draw_init_state(
        n_cluster, alpha0, alphabet, reads_list, reference_binary
    )
    state_init_dict.update(
        {
            "lnB_alpha0": lnB(state_init_dict["alpha"]),
            "betaln_a0_b0": betaln(
                state_init_dict["gamma_a"], state_init_dict["gamma_b"]
            ),
        }
    )

    history_alpha = [state_init_dict['alpha']]
    history_mean_log_pi = [state_init_dict['mean_log_pi']]
    history_mean_log_gamma = [state_init_dict['mean_log_gamma']]
    history_mean_cluster = [state_init_dict['mean_cluster']]
    history_elbo = []

    # Iteratively update mean values
    iter = 0
    message = ""
    converged = False
    elbo = 0
    state_curr_dict = state_init_dict
    while converged == False:

        if iter <= 1:
            digamma_alpha_sum = digamma(state_curr_dict["alpha"].sum(axis=0))
            digamma_a_b_sum = digamma(
                state_curr_dict["gamma_a"] + state_curr_dict["gamma_b"]
            )
            state_curr_dict.update({"digamma_alpha_sum": digamma_alpha_sum})
            state_curr_dict.update({"digamma_a_b_sum": digamma_a_b_sum})

        state_curr_dict = update_eqs.update(
            reads_seq_binary,
            reads_weights,
            reference_binary,
            reads_log_error_proba,
            state_init_dict,
            state_curr_dict,
        )
        elbo = elbo_eqs.compute_elbo(
            reads_weights,
            reference_binary,
            reads_log_error_proba,
            state_init_dict,
            state_curr_dict,
        )

        history_elbo.append(elbo)
        history_mean_log_pi.append(state_curr_dict["mean_log_pi"])
        history_mean_log_gamma.append(state_curr_dict["mean_log_gamma"])
        history_mean_cluster.append(state_curr_dict["mean_cluster"])

        if iter > 1:
            if (history_elbo[-2] > elbo) and np.abs(elbo - history_elbo[-2]) > 1e-08:
                exit_message = "Error: ELBO is decreasing."
                break
            elif np.abs(elbo - history_elbo[-2]) < 1e-03:
                converged = True
                exit_message = "ELBO converged."

        state_curr_dict.update({"elbo": elbo})

        iter += 1
    # End: While-loop

    state_curr_dict.update({"elbo": elbo})

    dict_result.update(
        {
            "exit_message": exit_message,
            "n_iterations": iter,
            "converged": converged,
            "elbo": elbo,
            "history_elbo": history_elbo,
            "history_alpha": history_alpha,
            "history_mean_log_pi": history_mean_log_pi,
            "history_mean_log_gamma": history_mean_log_gamma,
            "history_mean_cluster": history_mean_cluster,
        }
    )

    result = (state_curr_dict, dict_result)

    return result
