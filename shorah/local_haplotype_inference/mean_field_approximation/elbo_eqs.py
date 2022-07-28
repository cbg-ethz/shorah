from scipy.special import betaln
import numpy as np
from scipy.stats._multivariate import _lnB as lnB


def compute_elbo(
    reads_weights, reference_binary, reads_log_error_proba, state_init, state_curr
):

    alpha0 = state_init["alpha"]
    a = state_init["gamma_a"]
    b = state_init["gamma_b"]

    lnB_alpha0 = state_init["lnB_alpha0"]
    betaln_a0_b0 = state_init["betaln_a0_b0"]

    mean_z = state_curr["mean_cluster"]
    mean_h = state_curr["mean_haplo"]

    alpha_updated = state_curr["alpha"]
    mean_log_pi = state_curr["mean_log_pi"]

    a_updated = state_curr["gamma_a"]
    b_updated = state_curr["gamma_b"]
    mean_log_gamma = state_curr["mean_log_gamma"]

    elbo = elbo_data(reads_weights, mean_z, mean_h, reads_log_error_proba)
    elbo += elbo_pi(alpha0, lnB_alpha0, alpha_updated, mean_log_pi)
    elbo += elbo_cluster(mean_z, mean_log_pi, reads_weights)
    elbo += elbo_haplo(reference_binary, mean_h, mean_log_gamma)
    elbo += elbo_gamma(a, b, betaln_a0_b0, mean_log_gamma, a_updated, b_updated)

    return elbo


def elbo_data(reads_weights, mean_cluster, mean_haplo, reads_log_error_proba):

    haplo_error_rate_part = np.einsum("NLB,KLB->NK", reads_log_error_proba, mean_haplo)
    mean_cluster_weight = np.einsum("N,NK->NK", reads_weights, mean_cluster)
    final = np.einsum("NK,NK->", mean_cluster_weight, haplo_error_rate_part)

    return final


def elbo_cluster(mean_cluster, mean_log_pi, reads_weights):

    mean_cluster_weight = np.einsum("N,NK->NK", reads_weights, mean_cluster)
    p_part = np.matmul(mean_cluster_weight, mean_log_pi).sum(axis=0)

    null_pos = mean_cluster > 0
    q_part = np.einsum("NK,NK->", np.log(mean_cluster, where=null_pos), mean_cluster)

    return p_part - q_part


def elbo_pi(alpha0, lnB_alpha0, updated_alpha, mean_log_pi):
    # p_part = (-1)*lnB(alpha0) + np.multiply(alpha0-1,mean_log_pi).sum(axis=0)
    p_part = (-1) * lnB_alpha0 + np.multiply(alpha0 - 1, mean_log_pi).sum(axis=0)
    q_part = (-1) * lnB(updated_alpha) + np.multiply(
        updated_alpha - 1, mean_log_pi
    ).sum(axis=0)
    return p_part - q_part


def elbo_haplo(reference_table, mean_haplo, mean_log_gamma):

    B = mean_haplo.shape[2]

    b1 = mean_log_gamma[0]
    b2 = mean_log_gamma[1] - np.log(B - 1)

    haplo_p = b1 * np.einsum("KLB,LB->", mean_haplo, reference_table)
    haplo_p += b2 * np.einsum("KLB,LB->", mean_haplo, (1 - reference_table))

    null_pos = mean_haplo > 0
    haplo_q = np.einsum("KLB,KLB->", np.log(mean_haplo, where=null_pos), mean_haplo)

    return haplo_p - haplo_q


def elbo_gamma(a0, b0, betaln_a0_b0, mean_log_gamma, a_up, b_up):
    p_part = (a0 - 1) * mean_log_gamma[0] + (b0 - 1) * mean_log_gamma[1] - betaln_a0_b0
    q_part = (
        (a_up - 1) * mean_log_gamma[0]
        + (b_up - 1) * mean_log_gamma[1]
        - betaln(a_up, b_up)
    )
    return p_part - q_part
