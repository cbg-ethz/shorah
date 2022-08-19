from scipy.special import digamma
from scipy.special import betaln
import numpy as np

def update(reads_seq_binary, reads_weights, reference_binary, reads_log_error_proba, state_init, state_curr):

    #sssstart_time = timer()
    alpha0 = state_init['alpha']
    a = state_init['gamma_a']
    b = state_init['gamma_b']

    mean_z = state_curr['mean_cluster']
    mean_h = state_curr['mean_haplo']

    #alpha_updated = state_curr['alpha']
    mean_log_pi = state_curr['mean_log_pi']
    digamma_alpha_sum = state_curr['digamma_alpha_sum']
    mean_log_gamma = state_curr['mean_log_gamma']
    digamma_a_b_sum=state_curr['digamma_a_b_sum']
    mean_h = update_mean_haplo(reads_weights,reference_binary, reads_log_error_proba, mean_z,mean_log_gamma)
    mean_z = update_mean_cluster(mean_log_pi,mean_h,reads_log_error_proba)
    alpha_updated = update_alpha(alpha0, mean_z,reads_weights)
    mean_log_pi = get_mean_log_pi(alpha_updated, digamma_alpha_sum)
    a_updated,b_updated = update_a_and_b(reference_binary,mean_h,a,b)
    mean_log_gamma =get_mean_log_beta_dist(a_updated,b_updated,digamma_a_b_sum)

    state_curr_dict_new = dict({'alpha': alpha_updated,
                            'mean_log_pi': mean_log_pi,
                            'gamma_a': a_updated,
                            'gamma_b': b_updated,
                            'digamma_alpha_sum':digamma_alpha_sum,
                            'digamma_a_b_sum':digamma_a_b_sum,
                            'mean_log_gamma': mean_log_gamma,
                            'mean_haplo': mean_h,
                            'mean_cluster': mean_z
                            })

    return state_curr_dict_new

def get_mean_log_pi(alpha, digamma_alpha_sum):
    """
    Note that the digamma function can be inefficient.
    """
    #digamma_alpha_sum = digamma(alpha.sum(axis=0))
    mean_log_pi = digamma(alpha)-digamma_alpha_sum
    return mean_log_pi

def get_mean_log_beta_dist(a,b,digamma_sum):
    # I tested this one and it seemed not to make the difference
    #digamma_sum = digamma(a+b)
    mean_log_gamma = digamma(a)-digamma_sum
    mean_log_gamma_inv = digamma(b)-digamma_sum
    return mean_log_gamma, mean_log_gamma_inv

def update_mean_cluster(mean_log_pi,mean_haplo,reads_log_error_proba):

    haplo_error_rate_part = np.einsum('NLB,KLB->NK', reads_log_error_proba, mean_haplo)

    temp_c = haplo_error_rate_part
    temp_c[:] += mean_log_pi

    del haplo_error_rate_part

    max_z = np.max(temp_c, axis=1)
    max_z = max_z[:, np.newaxis]
    mean_z= np.exp(temp_c-max_z)
    c_normalize = mean_z.sum(axis=1)
    c_normalize = c_normalize[:, np.newaxis]
    mean_z = mean_z/c_normalize

    return mean_z

def update_mean_haplo(reads_weights,reference_table, reads_log_error_proba, mean_cluster,mean_log_gamma):

    B=reference_table.shape[1] # size of alphabet
    ref_part = reference_table*mean_log_gamma[0]+(1-reference_table)*(mean_log_gamma[1]-np.log(B-1))

    mean_cluster_weight = np.einsum('N,NK->NK',reads_weights,mean_cluster)
    cluster_assignment_part = np.einsum('NLB,NK->KLB', reads_log_error_proba, mean_cluster_weight)

    log_mean_haplo = cluster_assignment_part # shape: (K,L,B)
    log_mean_haplo[:] += ref_part

    del ref_part
    del cluster_assignment_part
    del mean_cluster_weight

    max_hap = np.max(log_mean_haplo, axis=2) # shape: (K,L)
    max_hap = max_hap[:, :, np.newaxis]
    mean_haplo = np.exp(log_mean_haplo-max_hap)
    c_normalize=mean_haplo.sum(axis=2)
    c_normalize = c_normalize[:,:, np.newaxis]
    mean_haplo= mean_haplo/c_normalize

    return mean_haplo

def update_a_and_b(reference_table,mean_haplo,a,b):
    # update a and b for mutation rate gamma

    up_a = a.copy()
    up_a += np.einsum('KLB,LB->',mean_haplo,reference_table)

    up_b = b.copy()
    up_b += np.einsum('KLB,LB->',mean_haplo,(1-reference_table))

    return up_a,up_b

def update_alpha(alpha, mean_cluster,reads_weights):

    temp_alpha = alpha.copy()
    temp_mean_cluster = mean_cluster.copy()
    mean_cluster_weight = np.einsum('N,NK->NK',reads_weights,temp_mean_cluster)
    temp_alpha+=mean_cluster_weight.sum(axis=0)

    return temp_alpha
