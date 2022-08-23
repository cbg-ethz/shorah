import numpy as np
from scipy.special import digamma
from . import quality_scores_update_eqs as update_eqs


def draw_init_state(n_clusters, alpha0, alphabet, reads_list, reference_binary, qualities):

    genome_length = reads_list[0].seq_binary.shape[0]
    n_reads = len(reads_list)
    size_alphabet = len(alphabet)

    # concentration parameter for Dirichlet prior of components
    alpha_temp = alpha0 * np.ones(n_clusters)

    # initialization of mean values
    digamma_alpha_sum = digamma(alpha_temp.sum(axis=0))
    mean_log_pi = update_eqs.get_mean_log_pi(alpha_temp, digamma_alpha_sum)

    matches = 10
    mismatch = 2

    k = np.random.uniform(low=0.5, high=1.0, size=2)
    a, b = matches * k[0], mismatch * k[1]

    gamma0 = np.random.beta(a, b)
    mean_log_gamma = np.log(gamma0), np.log(1 - gamma0)
    mean_h = init_mean_haplo(
        n_clusters, genome_length, size_alphabet, mean_log_gamma, reference_binary
    )

    mean_z = init_mean_cluster(n_clusters, n_reads, alpha0)

    state_init_dict = dict(
        {
            "alpha": alpha_temp,
            "mean_log_pi": mean_log_pi,
            "gamma_a": a,
            "gamma_b": b,
            "mean_log_gamma": mean_log_gamma,
            "mean_haplo": mean_h,
            "mean_cluster": mean_z,
        }
    )

    if not qualities:
        # draw the error params
        k = np.random.uniform(low=0.5, high=1.0, size=2)
        a, b = matches * k[0], mismatch * k[1]

        theta0 = np.random.beta(a, b)
        mean_log_theta = np.log(theta0), np.log(1-theta0)
        state_init_dict['mean_log_theta'] = mean_log_theta
        state_init_dict['theta_c'] = a
        state_init_dict['theta_d'] = b


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


def init_mean_cluster(n_clusters, n_reads, alpha0):

    mean_z = np.random.dirichlet(np.ones(n_clusters) * alpha0, size=n_reads)

    if np.any(np.isnan(mean_z)):
        alpha_new = alpha0 * 10
        mean_z = init_mean_cluster(n_clusters, n_reads, alpha_new)

    return mean_z


def init_mean_haplo(
    n_clusters, genome_length, size_alphabet, mean_log_gamma, reference_table
):
    base_true = np.exp(mean_log_gamma[0]) * np.ones(
        (n_clusters, genome_length, size_alphabet)
    )

    base_false = (
        (1.0 - np.exp(mean_log_gamma[1]))
        / (size_alphabet - 1)
        * np.ones((n_clusters, genome_length, size_alphabet))
    )

    return np.multiply(
        np.power(base_true, reference_table), np.power(base_false, 1 - reference_table)
    )
