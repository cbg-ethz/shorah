import numpy as np
from scipy.special import digamma

from . import update_eqs as update_eqs


def draw_init_state(K, alpha0, alphabet, reads_list, reference_binary):

    genome_length = reads_list[0].seq_binary.shape[0]  # length of seq
    n_reads = len(reads_list)  # number of reads
    # fixed parameters
    alpha_temp = alpha0 * np.ones(
        K
    )  # concentration parameter for Dirichlet prior of components

    # initialization of mean values
    digamma_alpha_sum = digamma(alpha_temp.sum(axis=0))
    mean_log_pi = update_eqs.get_mean_log_pi(alpha_temp, digamma_alpha_sum)

    matches = 10
    mismatch = 2

    k = np.random.uniform(low=0.5, high=1.0, size=4)
    a, b = matches * k[0], mismatch * k[1]
    c, d = matches * k[2], mismatch * k[3]

    theta0 = np.random.beta(c, d)
    gamma0 = np.random.beta(a, b)

    mean_log_gamma = np.log(gamma0), np.log(1 - gamma0)
    mean_log_theta = np.log(theta0), np.log(1 - theta0)

    mean_h = init_mean_haplo(K, genome_length, alphabet, mean_log_gamma, reference_binary)

    mean_z = np.ones((n_reads, K)) / K
    for n in range(n_reads):
        mean_z[n] = np.random.dirichlet((alpha_temp) * 100)

    state_init_dict = dict(
        {
            "alpha": alpha_temp,
            "mean_log_pi": mean_log_pi,
            "theta_c": c,
            "theta_d": d,
            "mean_log_theta": mean_log_theta,
            "gamma_a": a,
            "gamma_b": b,
            "mean_log_gamma": mean_log_gamma,
            "mean_haplo": mean_h,
            "mean_cluster": mean_z,
        }
    )

    return state_init_dict


def count_mis_and_matches_wrt_ref(reads_list, reference_table):
    matches = 0
    mismatch = 0
    totbase = 0
    for n in range(len(reads_list)):  # iterate over reads
        matches += reads_list[n].weight * (
            np.multiply(reference_table, reads_list[n].seq_binary)
            .sum(axis=0)
            .sum(axis=0)
        )
        totbase += reads_list[n].weight * reads_list[n].n_non_N

    mismatch = totbase - matches

    return matches, mismatch


def init_mean_haplo(K, L, alphabet, mean_log_gamma, reference_table):
    base_true = np.exp(mean_log_gamma[0])
    base_false = (1.0 - np.exp(mean_log_gamma[1])) / (len(alphabet) - 1)
    mean_haplo = np.ones((K, L, len(alphabet)))
    for k in range(K):
        for position in range(L):
            for base in range(len(alphabet)):
                mean_haplo[k][position][base] = (
                    base_true ** reference_table[position][base]
                ) * (base_false ** (1 - reference_table[position][base]))
    return mean_haplo
