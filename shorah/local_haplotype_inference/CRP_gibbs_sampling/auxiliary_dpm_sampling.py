import pandas as pd
import skbio
from skbio import Sequence
from skbio.sequence.distance import hamming
import numpy as np
import random


class State:
    """
    The state of the sampled variables.
    """
    def __init__(self, cluster_list, theta, gamma, window_size, alphabet ,alpha):
        self.cluster_list = cluster_list
        self.theta = theta
        self.gamma = gamma
        self.B = len(alphabet)
        self.alpha = alpha
        self.alphabet = alphabet
        self.window_size = window_size
        self.cluster_assignments = None
        self.sampled_haplotypes = None

    def get_cluster_assignments(self, reads_list):
        return [temp_read.metadata['cluster'] for temp_read in reads_list]

    def get_sampled_haplotypes(self):
        return [temp_cluster.haplotype_seq for temp_cluster in self.cluster_list]

    def to_dict(self, reads_list):
        temp_dict = {
                    'theta': self.theta,
                    'gamma': self.gamma,
                    'alpha': self.alpha,
                    'alphabet': str(self.alphabet),
                    'n_cluster': len(self.cluster_list)#,
                    #'haplotypes': str(self.get_sampled_haplotypes())
                    }
        for idx_read , assigned_cluster in enumerate(self.get_cluster_assignments(reads_list)):
            temp_dict.update({'c_'+str(idx_read): assigned_cluster})
        return temp_dict


class DPM_cluster:
    """
    Class for clusters that are sampled, includes the corresponding sampled haplotype.
    """
    def __init__(self, reads_idx_list,reference_seq):
        self.reads_idx_list = reads_idx_list # list of indexes of reads in cluster
        self.haplotype_seq = None
        self.distance2reads = [distance(read, self.haplotype_seq) for read in reads_list] if self.haplotype_seq is not None else None #distance between haplo and all reads
        self.matches2reads = [matches(read, self.haplotype_seq) for idx_read in reads_list] if self.haplotype_seq is not None else None  # all reads
        self.matches2reference = matches(reference_seq, self.haplotype_seq) if self.haplotype_seq is not None else None
        self.distance2reference = distance(reference_seq, self.haplotype_seq) if self.haplotype_seq is not None else None

    def update_distances(self, reads_list,reference_seq):
        if self.haplotype_seq is not None:
            self.distance2reads = [distance(read, self.haplotype_seq) for read in reads_list]
            self.matches2reads = [matches(read, self.haplotype_seq) for read in reads_list]
            self.matches2reference = matches(reference_seq, self.haplotype_seq)
            self.distance2reference = distance(reference_seq, self.haplotype_seq)

def sample_haplotype(reads_idx_list, reads_list, state, total_bases_counts, totbases):
    #'----- Haplotype sampling ----'

    b1= np.log(state.theta) # theta = rate that base is drawn without error
    b2= np.log((1-state.theta)/(state.B-1))

    haplotype_sequence = []

    alphabet = state.alphabet
    theta = state.theta

    for j in range(state.window_size):

        # count how often a specific base occures in reads from cluster
        bases_counts = np.zeros(state.B) #[A,C,G,T,-], corresponds to cbase in ShoRAH
        log_pbase=np.zeros(state.B)
        exp_pbase=np.zeros(state.B)
        pbase =np.zeros(state.B)
        #alphabet=['A','C','G','T','-'] # we do not sample N in the haplotype

        totreads=0
        for reads_idx in reads_idx_list:
            if str(reads_list[reads_idx].values[j])[2] != 'N':
                totreads+=reads_list[reads_idx].metadata['weight']
                base_idx=alphabet.index(str(reads_list[reads_idx].values[j])[2])
                bases_counts[base_idx]+=reads_list[reads_idx].metadata['weight']

        if totreads==0: # there is only Ns in the cluster reads sample from all reads
            for i in range(len(bases_counts)):
                log_pbase[i]=total_bases_counts[j][i]*b1+(totbases - total_bases_counts[j][i])*b2 # eq9
        else:
            for i in range(len(bases_counts)):
                log_pbase[i]=bases_counts[i]*b1 +(totreads - bases_counts[i])*b2 # eq9

        max_log_pbase=np.max(log_pbase)
        base_id = np.argmax(log_pbase)

        for i in range(len(bases_counts)):
            if i == base_id:
                pbase[i]=1.0 # = np.exp(log_pbase[i]-max_log_pbase)
                log_pbase[i]-=max_log_pbase
            else:
                log_pbase[i]-=max_log_pbase
                pbase[i]=np.exp(log_pbase[i])

        if theta != 1.0:
            x = one_shot_discrete_B(pbase)[0]
            haplotype_sequence.append(x)
        else:
            haplotype_sequence.append(alphabet[np.argmax(bases_counts)]) # don't sample just take the majority one

    return Sequence(''.join(haplotype_sequence))


def sample_class(idx_read,reads_list,reference_seq,state,  total_bases_counts, totbases):
    list_of_clusters = state.cluster_list
    theta = state.theta
    gamma = state.gamma
    B = state.B
    alpha = state.alpha

    # run through the populated classes to assign a probability
    log_P_list =[]
    P_list = []
    for count, temp_cluster in enumerate(list_of_clusters):
        # exclude cluster where the read is assigned to if it only contains this read
        if idx_read in temp_cluster.reads_idx_list:
            old_cluster=count
        if not (temp_cluster.reads_idx_list == [idx_read]):
            x1, x2 = log_P(theta, B, idx_read, temp_cluster, reads_list)
            log_P_list.append(x1)
            P_list.append(x2)

    # probablity of creating a new cluster
    x1, x2 = log_P_reference(theta, gamma, alpha, B, idx_read, reads_list)
    log_P_list.append(x1)
    P_list.append(x2)

    #renormalization
    max_log_P = np.max(log_P_list)

    if max_log_P>=0:
        delta_log=-max_log_P
    else:
        delta_log=max_log_P

    for ll in range(len(log_P_list)):
        if P_list[ll]>0:
            log_P_list[ll]+=delta_log
            # in sample_haplotype we use the - and from the theory I would also tend to -
            # but the clustering for seq_verysimple.fasta works with both.
            # --> acutally it does not change what we do
            P_list[ll]=np.exp(log_P_list[ll])

    this_cluster = one_shot_discrete(len(log_P_list),P_list)[0]
    #print('extracted class is ' + str(this_cluster))

    #update the clusters
    if this_cluster == len(log_P_list)-1:
        #print('read ' ,idx_read, ' is assigned to new cluster (old_cluster=)', old_cluster)
        list_of_clusters.append(DPM_cluster([idx_read],reference_seq))
        list_of_clusters[-1].haplotype_seq = sample_haplotype([idx_read], reads_list, state, total_bases_counts, totbases)
        list_of_clusters[-1].update_distances(reads_list,reference_seq) # new to avoid this error: --> seemed to solve the problem
        """
         File "/Users/lfuhrmann/Documents/Projects/ShoRAH2.0/Code/ShoRAH2.0/auxiliary_dpm_sampling.py", line 132, in sample_class
    x1, x2 = log_P(theta, B, idx_read, temp_cluster, reads_list)
  File "/Users/lfuhrmann/Documents/Projects/ShoRAH2.0/Code/ShoRAH2.0/auxiliary_dpm_sampling.py", line 262, in log_P
    log_P= np.log(tw)+cluster.matches2reads[reads_idx]*b1 + cluster.distance2reads[reads_idx]*b2
TypeError: 'NoneType' object is not subscriptable
        """
        list_of_clusters[old_cluster].reads_idx_list.remove(idx_read) # remove read from old cluster

    else:
        #print('read ',idx_read, ' is assigned to cluster ', this_cluster, '(old_cluster=)', old_cluster)
        if idx_read not in list_of_clusters[this_cluster].reads_idx_list:
            list_of_clusters[this_cluster].reads_idx_list.append(idx_read)
            list_of_clusters[old_cluster].reads_idx_list.remove(idx_read)

    reads_list[idx_read].metadata.update({'cluster': this_cluster})

    #if old_cluster empty
    if list_of_clusters[old_cluster].reads_idx_list ==[]:
        list_of_clusters.remove(list_of_clusters[old_cluster])

    state.cluster_list = list_of_clusters
    return list_of_clusters


### ----------------------------------------------------
### ---- Functions used in each sampling iteration -----
### ----------------------------------------------------

def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]

def remove_at(list_positions, s):
    for i in list_positions:
        s = s[:i] + s[i+1:]
    return s

def distance(read_a, read_b):
    '''
    Positions where either read_a or read_b display an N are not accounted for.
    '''
    N_positions = find(str(read_a),'N')+find(str(read_b),'N')

    read_a = remove_at(N_positions, str(read_a))
    read_b = remove_at(N_positions, str(read_b))

    return hamming(Sequence(str(read_a)), Sequence(str(read_b)))*len(read_a)

def matches(read_a, read_b):
    '''
    Positions where either read_a or read_b display an N are not accounted for.
    '''
    N_positions = find(str(read_a),'N')+find(str(read_b),'N')

    read_a = remove_at(N_positions, str(read_a))
    read_b = remove_at(N_positions, str(read_b))

    return len(read_a) - hamming(Sequence(str(read_a)), Sequence(str(read_b)))*len(read_a)

def weight_shift(reads_idx, curr_cluster,reads_list):
    '''
    Compute weight of cluster (number of reads_list in cluster)
    excluding weight of read with index reads_list_idx.
    '''
    w = weight_of_cluster(curr_cluster, reads_list)
    if reads_idx in curr_cluster.reads_idx_list:
        w-=reads_list[reads_idx].metadata['weight']
    return w

def weight_of_cluster(curr_cluster, reads_list):
    w=0
    for idx in curr_cluster.reads_idx_list:
        w+=reads_list[idx].metadata['weight']
    return w


def one_shot_discrete_B(weights):
    return random.choices(population=['A','C','G','T','-'],weights=weights)

def one_shot_discrete(n, P):
    '''
    Simple naive implementation of a discrete distribution, using cumulated weights
    Boost implementation:
    https://www.boost.org/doc/libs/1_70_0/libs/compute/doc/html/boost/compute/discrete_distribution.html

    n: numbers to choose from
    P: weights
    '''
    return random.choices(population=range(n),weights=P)

def log_P(theta, B, reads_idx, cluster, reads_list):
    """
    reads_idx: index of read in reads_list (list)
    cluster: cluster for which the log(P) is computed
    """
    b1=np.log(theta)
    b2=np.log((1-theta)/(B-1))

    # weighted clusre size, excluding read i if it is present in the current cluster
    # tw = total weight of cluster - weigth of read i if present in the cluster
    tw = weight_shift(reads_idx,cluster, reads_list)
    if tw<0: print('ERROR with tw, value of tw=', tw)

    if theta!=1:
        log_P= np.log(tw)+cluster.matches2reads[reads_idx]*b1 + cluster.distance2reads[reads_idx]*b2
        P=1.0
    else:
        if distance2reads_list[reads_idx]!=0:
            log_P = double_threshold_min - 1 #???
            P=0.0
        else:
            #TODO log_P = ???? line 1479
            P=1.0

    return log_P, P

def log_P_reference(theta, gamma, alpha, B, read_idx, reads_list):
    '''
    computation of the probablity to assign read read_idx to a new cluster
    '''
    b1=np.log(theta*gamma+(1-gamma)*(1-theta)/(B-1))
    b2=np.log((theta+gamma+B*(1-gamma*theta)-2)/((B-1)**2))

    if theta!=1 or gamma!=1:
        log_P=np.log(alpha)+np.log(reads_list[read_idx].metadata['weight'])+reads_list[read_idx].metadata['matches2reference']*b1+reads_list[read_idx].metadata['distance2reference']*b2
        P=1.0
    if (theta==1 and gamma==1):
        if reads_list[read_idx].distance2reference ==0:
            lop_P = np.log(alpha)
            P=1.0
        else:
            log_P= double_threshold_min - 1 #???
            P=0.0

    return log_P, P

def get_idx_cluster_of_read(read_index):
    for idx_cluster, cluster in enumerate(list_of_clusters):
        if read_index in cluster.reads_idx_list:
            return idx_cluster

    print('read has not been assigned to a cluster')
