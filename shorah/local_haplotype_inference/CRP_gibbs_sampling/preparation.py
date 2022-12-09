import pandas as pd
import skbio
from skbio import Sequence
from skbio.sequence.distance import hamming
import numpy as np
import random

from . import auxiliary_dpm_sampling as aux_dpm


### ---------------------------------------
### ---------------------------------------
### ---------------------------------------

def history_dataframe(n_reads):
    import pandas as pd
    columns =  ['n_iter', 'alpha', 'alphabet', 'n_cluster', 'haplotypes']
    columns += ['c_'+str(n) for n in range(n_reads)]
    columns += ['gamma', 'theta']
    # convergence diagnostics
    columns += ['ess_c_'+str(n) for n in range(n_reads)]
    columns += ['ess_gamma', 'ess_theta']
    columns += ['geweke_c_'+str(n) for n in range(n_reads)]
    columns += ['geweke_gamma', 'geweke_theta']

    return pd.DataFrame(columns=columns)

### ---------------------------------------
### ----- Functions for initializing ----
### ---------------------------------------

def init_clusters(reads_list, reference_seq):
    N = len(reads_list)
    read_idx_list = list(range(N))
    list_of_clusters = []
    import random
    n_clusters = int((1+N*0.2)+(0.5*N-(1+N*0.2))*random.uniform(0, 1))
    random.shuffle(read_idx_list)
    for read_idxs in [read_idx_list[i::n_clusters] for i in range(n_clusters)]:
        list_of_clusters.append(aux_dpm.DPM_cluster(read_idxs,reference_seq))
    return list_of_clusters


### ---------------------------------------
### ----- Functions for preparation ----
### ---------------------------------------

def load_reference_seq(reference_file):
    for seq in skbio.io.read(reference_file, format='fasta'):
         return seq

# This function is used for data/simple/seq.fasta and data/simple/seq_verysimple.fasta
def load_fasta2reads_list(reads_file):
    '''
    reads_list fasta with sequences and writes unique reads_list in reads_list_idx_list,
    adds also meta information
    weight = the number of reads_list that are the same
    '''
    reads_list = [seq for seq in skbio.io.read(reads_file, format='fasta')]
    for temp_read in reads_list:
        temp_read.metadata.update({'weight': 1.0})
        temp_read.metadata.update({'identical_reads': [temp_read.metadata['id']]})

    # test for unique reads_list
    reads_list = unique_reads_list(reads_list)

    return reads_list

def unique_reads_list(reads_list):
    # test for unique reads_list
    for i, temp_read in enumerate(reads_list):
        if temp_read.metadata['weight']>0.0:
            for j in range(i+1,len(reads_list)):
                hd = hamming(temp_read, reads_list[j])
                if hd==0:
                    temp_read.metadata['weight']+=1
                    temp_read.metadata['identical_reads'].append(reads_list[j].metadata['id'])
                    reads_list[j].metadata['weight']-=1

    # keep only unique reads_list
    reads_list=[read for read in reads_list if read.metadata['weight']>0]
    return reads_list

def load_bam2reads_list(bam_file):
    import pysam
    samfile = pysam.AlignmentFile("data/amplicon_test/ampli_sorted.bam", "rb")
    reads_list = []
    for read in samfile:
        read_dict = read.to_dict()
        seq_obj = Sequence(read_dict['seq'])
        seq_obj.metadata['id']=read_dict['name']
        seq_obj.metadata.update({'weight': 1.0})
        seq_obj.metadata.update({'identical_reads': [read_dict['name']]})
        reads_list.append(seq_obj)
    return reads_list

def count_totbases(reads_list):
    # total number of bases over all reads_list, do not count 'N'
    totbases=0
    for read in reads_list:
        for base in iter(read):
            if base!='N':
                totbases+=read.metadata['weight']
    return totbases

def base_counts_per_position(reads_list, window_size, alphabet):
    '''
    output of this function is used in sample_haplotype()
    '''
    total_bases_counts = np.zeros(shape=(window_size,5)) #[A,C,G,T,-], corresponds to ftable in ShoRAH
    for j in range(window_size):
        for read in reads_list:
            if str(read.values[j])[2] != 'N':
                base_idx=alphabet.index(str(read.values[j])[2])
                total_bases_counts[base_idx]+=read.metadata['weight']
    return total_bases_counts
