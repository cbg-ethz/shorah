import skbio
import numpy as np
from skbio.sequence.distance import hamming
from skbio import Sequence

class Read:
    def __init__(self, seq_string, seq_id):
        self.seq_string=seq_string
        self.weight =1
        self.id = seq_id
        self.seq_binary=[]
        self.identical_reads = []
        self.n_non_N = len(seq_string) - seq_string.count('N')

    def seq2binary(self, alphabet):
        length_seq = len(self.seq_string)
        seq_table = np.zeros((length_seq, len(alphabet)))
        for base_position, base in enumerate(str(self.seq_string)):
            if alphabet.find(base)>=0:
                seq_table[base_position][alphabet.index(base)]=1

        self.seq_binary = seq_table


def get_average_qualities(fname_qualities, reads_list):

    with open(fname_qualities, 'rb') as f:
        qualities = np.load(f, allow_pickle = True)

    # get average of qualties scores for same reads
    unique_qualities = np.full((len(reads_list), int(qualities.shape[1])),1)

    # get average of qualties scores for same reads
    for i, temp_read in enumerate(reads_list):
        unique_qualities[i]+= qualities[i]
        if len(temp_read.identical_reads) > 0:
            for j in temp_read.idx_identical_reads:
                unique_qualities[i]+= qualities[j]
            unique_qualities[i] = unique_qualities[i]/temp_read.weight

    return unique_qualities

def get_reads_log_error_proba(qualities, reads_seq_binary, size_alphabet):
    """
    \log \theta_{n,l}^{r_{nl}^i}  \left( \frac{1-\theta_{n,l}}{B-1}\right)^{(1-r_{nl}^i)} \right)
    Do I need the weights? No this can be integrated later in update_eqs.py
    """
    reads_log_error_proba = 1 - 10**(-qualities/10)
    reads_log_error_proba = reads_log_error_proba[:,:,np.newaxis]
    reads_log_error_proba = np.tile(reads_log_error_proba, (1,1,size_alphabet))

    reads_log_inv_error_proba = (10**(-qualities/10)) / (size_alphabet-1)
    # this was written in the old version I dont know why there is the 1- () in the denominator
    # reads_log_inv_error_proba = (10**(-qualities/10)) / (1 - (size_alphabet-1))
    reads_log_inv_error_proba = reads_log_inv_error_proba[:,:,np.newaxis]
    reads_log_inv_error_proba = np.tile(reads_log_inv_error_proba, (1,1,size_alphabet))

    reads_log_error_proba = np.power(reads_log_error_proba, reads_seq_binary)
    reads_log_error_proba += np.power(reads_log_inv_error_proba, (1-reads_seq_binary))
    reads_log_error_proba = np.log(reads_log_error_proba) # dimension: NxLxB

    # if reads_list[n].seq_binary[l].sum(axis=0)=0 then "N" at position l then position l is ignored
    # there will be a zero in the row n,l
    # dimension: NL
    all_N_pos=reads_seq_binary.sum(axis=2)>0
    all_N_pos= all_N_pos[:,:,np.newaxis]
    all_N_pos = np.tile(all_N_pos, (1,1,size_alphabet))

    # write zero where there is an "N"  in the position
    reads_log_error_proba[~all_N_pos]=0

    return reads_log_error_proba # dimension: NxLxB


def reads_list_to_array(reads_list):

    reads_binary=[reads_list[n].seq_binary for n in range(len(reads_list))]
    reads_binary_array = np.asarray(reads_binary)

    reads_weights = [reads_list[n].weight for n in range(len(reads_list))]
    reads_weights_array = np.asarray(reads_weights)

    """
    print('number of reads ', len(reads_list))
    print('lenght of sequence ', reads_list[0].seq_binary.shape[0])
    print('number of bases ', reads_list[0].seq_binary.shape[1])
    print('shape of reads_seq_binary ', reads_binary_array.shape)
    print(' shape weights ', reads_weights_array.shape)

    Output of those prints:
    number of reads  158
    lenght of sequence  140
    number of bases  5
    shape of reads_seq_binary  (158, 140, 5)
    shape weights  (158,)
    """

    return reads_binary_array, reads_weights_array

def load_fasta2reads_list(reads_fasta_file, alphabet):
    # go through each sequence in fasta file
    reads_list =[]
    for seq in skbio.io.read(reads_fasta_file, format='fasta'):
        reads_list.append(Read(str(seq),seq.metadata['id']))
        reads_list[-1].seq2binary(alphabet)
    # unique reads_list
    #reads_list = unique_reads_list(reads_list)

    return reads_list

def load_bam2reads_list(reads_fasta_file, alphabet):
    # go through each sequence in fasta file
    reads_list =[]
    for seq in skbio.io.read(reads_fasta_file, format='fasta'):
        reads_list.append(Read(str(seq),seq.metadata['id']))
        reads_list[-1].seq2binary(alphabet)
    # unique reads_list
    reads_list = unique_reads_list(reads_list)

    return reads_list


def load_bam2reads_list(bam_file, alphabet):
    import pysam
    samfile = pysam.AlignmentFile(bam_file, "rb")
    reads_list = []
    for read in samfile:
        read_dict = read.to_dict()
        #print(read_dict['flag'])
        #break
        reads_list.append(Read(str(read_dict['seq']),read_dict['name']))
        reads_list[-1].seq2binary(alphabet)
    # unique reads_list
    reads_list = unique_reads_list(reads_list)

    return reads_list


def unique_reads_list(reads_list):
    # test for unique reads_list
    for i, temp_read in enumerate(reads_list):
        if temp_read.weight>0.0:
            for j in range(i+1,len(reads_list)):
                hd = hamming(Sequence(temp_read.seq_string), Sequence(reads_list[j].seq_string))
                if hd==0:
                    temp_read.weight+=1
                    temp_read.identical_reads.append(reads_list[j].id)
                    reads_list[j].weight-=1

    # keep only unique reads_list
    reads_list=[read for read in reads_list if read.weight>0]
    return reads_list

def load_reference_seq(reference_file):
    for seq in skbio.io.read(reference_file, format='fasta'):
         return seq, seq.metadata['id']

def reference2binary(reference_seq, alphabet):
    length_seq = len(reference_seq)
    reference_table = np.zeros((length_seq, len(alphabet)))
    for base_position, base in enumerate(str(reference_seq)):
        reference_table[base_position][alphabet.index(base)]=1
    return reference_table
