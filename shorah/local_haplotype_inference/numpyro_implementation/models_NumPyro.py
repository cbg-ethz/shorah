import numpyro
import numpyro.distributions as dist

import jax
import jax.numpy as jnp

from . import helper

custom_put_along_axis = helper.custom_put_along_axis

def model_finiteDPM_extended(input_data):
    '''
    This model is only for alphabet={A,C,G,T,-} as deletion mutations and errors are treated seperatly
    (wrt to subsitutions).
    '''
    reference, read_data, alphabet_length, max_cluster_num_guess, is_observed, pos_del_ref, alpha0 = input_data

    # parameters
    read_count = read_data.shape[0]
    genome_length = read_data.shape[1]
    alphabet_length = alphabet_length
    cluster_num = max_cluster_num_guess

    # define rates
    mutation_rate = numpyro.sample('mutation_rate', dist.Dirichlet(jnp.ones(4)/4))
    error_rate = numpyro.sample('error_rate', dist.Dirichlet(jnp.ones(4)/4))

    # create matrix of rates
    mutation_rate_matrix = jnp.full(
        (genome_length, alphabet_length), mutation_rate[1] / (alphabet_length - 2)
    ) # subsitutions
    mutation_rate_matrix = mutation_rate_matrix.at[:,4].set(mutation_rate[2]) # deletion in haplotype

    mutation_rate_matrix = mutation_rate_matrix.at[pos_del_ref,:].set(mutation_rate[3] / (alphabet_length - 1)) # deletion in reference

    mutation_rate_matrix = custom_put_along_axis(
        mutation_rate_matrix, reference.reshape(genome_length, 1), mutation_rate[0], axis=1
    ) # matches
    mutation_rate_matrix = mutation_rate_matrix / mutation_rate_matrix.sum(axis=1)[:,jnp.newaxis] #normalize


    alpha0 = alpha0 # alpha must be bigger than zero
    alpha = alpha0 * jnp.ones(cluster_num) / cluster_num
    cluster_weights = numpyro.sample("cluster_weights", dist.Dirichlet(jnp.ones(cluster_num)))

    genome_axis = numpyro.plate("genome_axis", genome_length, dim=-1)
    with numpyro.plate("haplotype_axis", cluster_num, dim=-2):

        with genome_axis:
            haplotypes = numpyro.sample("haplotypes", dist.Categorical(mutation_rate_matrix))  # cluster centers

    with numpyro.plate("read_axis", read_count, dim=-2):
        cluster_assignments = numpyro.sample("cluster_assignments", dist.Categorical(cluster_weights))
        with genome_axis:
            error_rate_matrix = jnp.full(
                (read_count, genome_length, alphabet_length),
                error_rate[1] / (alphabet_length - 2),
            ) # substitutions
            error_rate_matrix = error_rate_matrix.at[:,:,4].set(error_rate[2]) # deletion in read

            error_rate_del_in_haplo = jnp.full(
                (read_count, genome_length, alphabet_length), error_rate[3] / (alphabet_length - 1))

            error_rate_matrix = jnp.where(
                haplotypes[cluster_assignments].reshape(read_count, genome_length, 1) == 4,
                error_rate_del_in_haplo,
                error_rate_matrix
            ) # deletion in haplotype

            error_rate_matrix = custom_put_along_axis(
                error_rate_matrix,
                haplotypes[cluster_assignments].reshape(read_count, genome_length, 1),
                error_rate[0],
                axis=2,
            ) # matches

            obs=numpyro.sample("obs", dist.Categorical(error_rate_matrix).mask(is_observed), obs=read_data)

def model_finiteDPM(input_data):
    reference, read_data, alphabet_length, max_cluster_num_guess, is_observed, alpha0 = input_data

    # parameters
    read_count = read_data.shape[0]
    genome_length = read_data.shape[1]
    alphabet_length = alphabet_length
    cluster_num = max_cluster_num_guess

    # define rates
    mutation_rate = numpyro.sample('mutation_rate', dist.Beta(1, 1))
    error_rate = numpyro.sample('error_rate', dist.Beta(1, 1))

    # create matrix of rates
    mutation_rate_matrix = jnp.full(
        (genome_length, alphabet_length), (1 - mutation_rate) / (alphabet_length - 1)
    )
    mutation_rate_matrix = custom_put_along_axis(
        mutation_rate_matrix, reference.reshape(genome_length, 1), mutation_rate, axis=1
    )

    alpha = alpha0 * jnp.ones(cluster_num) / cluster_num
    cluster_weights = numpyro.sample("cluster_weights", dist.Dirichlet(jnp.ones(cluster_num)))

    genome_axis = numpyro.plate("genome_axis", genome_length, dim=-1)
    with numpyro.plate("haplotype_axis", cluster_num, dim=-2):

        with genome_axis:
            haplotypes = numpyro.sample("haplotypes", dist.Categorical(mutation_rate_matrix))  # cluster centers

    with numpyro.plate("read_axis", read_count, dim=-2):
        cluster_assignments = numpyro.sample("cluster_assignments", dist.Categorical(cluster_weights))

        with genome_axis:

            error_rate_matrix = jnp.full(
                (read_count, genome_length, alphabet_length), (1 - error_rate) / (alphabet_length - 1)
            )
            error_rate_matrix = custom_put_along_axis(
                error_rate_matrix, haplotypes[cluster_assignments].reshape(read_count, genome_length, 1), error_rate, axis=2
            )

            obs=numpyro.sample("obs", dist.Categorical(error_rate_matrix).mask(is_observed), obs=read_data)


def mix_weights(beta):
    beta1m_cumprod = (1 - beta).cumprod(-1)
    return jnp.pad(beta, (0, 1), constant_values=1) * jnp.pad(
        beta1m_cumprod, (1, 0), constant_values=1
    )

def model_infiniteSBP(input_data):
    reference, read_data, alphabet_length, max_cluster_num_guess, is_observed = input_data

    max_cluster_num_guess = max_cluster_num_guess  # guess for maximum number of clusters

    # parameters
    read_count = read_data.shape[0]
    genome_length = reference.shape[0]
    alphabet_length = alphabet_length #jnp.unique(reference).shape[0]

    # define rates
    mutation_rate = numpyro.sample("mutation_rate", dist.Beta(1, 1))
    error_rate = numpyro.sample("error_rate", dist.Beta(1, 1))

    # create matrix of rates
    mutation_rate_matrix = jnp.full(
        (genome_length, alphabet_length), (1 - mutation_rate) / (alphabet_length - 1)
    )
    mutation_rate_matrix = custom_put_along_axis(
        mutation_rate_matrix, reference.reshape(genome_length, 1), mutation_rate, axis=1
    )

    with numpyro.plate("beta_plate", max_cluster_num_guess - 1):
        beta = numpyro.sample("beta", dist.Beta(1, 1.1))

    cluster_weights = numpyro.sample(
        "cluster_weights", dist.Dirichlet(mix_weights(beta))
    )

    genome_axis = numpyro.plate("genome_axis", genome_length, dim=-1)
    with numpyro.plate("haplotype_axis", max_cluster_num_guess, dim=-2):
        with genome_axis:
            haplotypes = numpyro.sample(
                "haplotypes", dist.Categorical(mutation_rate_matrix)
            )  # cluster centers

    with numpyro.plate("read_axis", read_count, dim=-2):
        cluster_assignments = numpyro.sample(
            "cluster_assignments", dist.Categorical(cluster_weights)
        )

        with genome_axis:
            error_rate_matrix = jnp.full(
                (read_count, genome_length, alphabet_length),
                (1 - error_rate) / (alphabet_length - 1),
            )
            error_rate_matrix = custom_put_along_axis(
                error_rate_matrix,
                haplotypes[cluster_assignments].reshape(read_count, genome_length, 1),
                error_rate,
                axis=2,
            )
            obs=numpyro.sample("obs", dist.Categorical(error_rate_matrix).mask(is_observed), obs=read_data)


def model_infiniteSBP_extended(input_data):
    '''
    This model is only for alphabet={A,C,G,T,-} as deletion mutations and errors are treated seperatly
    (wrt to subsitutions).
    '''
    reference, read_data, alphabet_length, max_cluster_num_guess, is_observed, pos_del_ref = input_data

    # parameters
    read_count = read_data.shape[0]
    genome_length = read_data.shape[1]
    alphabet_length = alphabet_length
    cluster_num = max_cluster_num_guess

    # define rates
    mutation_rate = numpyro.sample('mutation_rate', dist.Dirichlet(jnp.ones(4)/4))
    error_rate = numpyro.sample('error_rate', dist.Dirichlet(jnp.ones(4)/4))

    # create matrix of rates
    mutation_rate_matrix = jnp.full(
        (genome_length, alphabet_length), mutation_rate[1] / (alphabet_length - 2)
    ) # subsitutions
    mutation_rate_matrix = mutation_rate_matrix.at[:,4].set(mutation_rate[2]) # deletion in haplotype

    mutation_rate_matrix = mutation_rate_matrix.at[pos_del_ref,:].set(mutation_rate[3] / (alphabet_length - 1)) # deletion in reference

    mutation_rate_matrix = custom_put_along_axis(
        mutation_rate_matrix, reference.reshape(genome_length, 1), mutation_rate[0], axis=1
    ) # matches
    mutation_rate_matrix = mutation_rate_matrix / mutation_rate_matrix.sum(axis=1)[:,jnp.newaxis] #normalize

    with numpyro.plate("beta_plate", max_cluster_num_guess - 1):
        beta = numpyro.sample("beta", dist.Beta(1, 1.1))

    cluster_weights = numpyro.sample(
        "cluster_weights", dist.Dirichlet(mix_weights(beta))
    )

    genome_axis = numpyro.plate("genome_axis", genome_length, dim=-1)
    with numpyro.plate("haplotype_axis", max_cluster_num_guess, dim=-2):
        with genome_axis:
            haplotypes = numpyro.sample(
                "haplotypes", dist.Categorical(mutation_rate_matrix)
            )  # cluster centers

    with numpyro.plate("read_axis", read_count, dim=-2):
        cluster_assignments = numpyro.sample(
            "cluster_assignments", dist.Categorical(cluster_weights)
        )

        with genome_axis:
            error_rate_matrix = jnp.full(
                (read_count, genome_length, alphabet_length),
                error_rate[1] / (alphabet_length - 2),
            ) # substitutions
            error_rate_matrix = error_rate_matrix.at[:,:,4].set(error_rate[2]) # deletion in read

            error_rate_del_in_haplo = jnp.full(
                (read_count, genome_length, alphabet_length), error_rate[3] / (alphabet_length - 1))

            error_rate_matrix = jnp.where(
                haplotypes[cluster_assignments].reshape(read_count, genome_length, 1) == 4,
                error_rate_del_in_haplo,
                error_rate_matrix
            ) # deletion in haplotype

            error_rate_matrix = custom_put_along_axis(
                error_rate_matrix,
                haplotypes[cluster_assignments].reshape(read_count, genome_length, 1),
                error_rate[0],
                axis=2,
            ) # matches

            obs=numpyro.sample("obs", dist.Categorical(error_rate_matrix).mask(is_observed), obs=read_data)
