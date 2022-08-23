import numpy as np
from Bio import SeqIO
from scipy.spatial.distance import hamming


class Read:
    def __init__(self, seq_string, seq_id):
        self.seq_string = seq_string
        self.weight = 1
        self.phred_quality_score = None
        self.id = seq_id
        self.seq_binary = []
        self.identical_reads = []
        self.idx_identical_reads = []
        self.n_non_N = len(seq_string) - seq_string.count("N")

    def seq2binary(self, alphabet):
        length_seq = len(self.seq_string)
        seq_table = np.zeros((length_seq, len(alphabet)))
        for base_position, base in enumerate(str(self.seq_string)):
            if alphabet.find(base) >= 0:
                seq_table[base_position][alphabet.index(base)] = 1

        self.seq_binary = seq_table


def compute_reads_log_error_matrix(
    theta, one_minus_theta, reads_seq_binary, size_alphabet
):
    log_theta = np.log(theta)
    log_theta = log_theta[:, :, np.newaxis]
    log_theta = np.tile(log_theta, (1, 1, size_alphabet))

    log_one_minus_theta = np.log(one_minus_theta)
    log_one_minus_theta = log_one_minus_theta[:, :, np.newaxis]
    log_one_minus_theta = np.tile(log_one_minus_theta, (1, 1, size_alphabet))

    final = np.einsum("NLB,NLB->NLB", log_theta, reads_seq_binary)
    final += np.einsum("NLB,NLB->NLB", log_one_minus_theta, 1 - reads_seq_binary)

    # if reads_list[n].seq_binary[l].sum(axis=0)=0 then "N" at position l then position l is ignored
    # there will be a zero in the row n,l
    # dimension: NxL
    all_N_pos = reads_seq_binary.sum(axis=2) > 0
    all_N_pos = all_N_pos[:, :, np.newaxis]
    all_N_pos = np.tile(all_N_pos, (1, 1, size_alphabet))
    # write zero where there is an "N"  in the position
    final[~all_N_pos] = 0

    return final  # dimension: NxLxB


def compute_reads_log_error_proba(qualities, reads_seq_binary, size_alphabet):
    """
    \log \theta_{n,l}^{r_{nl}^i}  \left( \frac{1-\theta_{n,l}}{B-1}\right)^{(1-r_{nl}^i)} \right)
    Do I need the weights? No this can be integrated later in update_eqs.py

    theta_{n,l} = probabiltiy that at positiion l in read n the base was called without error
    1 - theta_{n,l} = error probablity, e.g. probablity that at position l in read n the base was called erroneous.

    input-qualities:
    Q_{n,l} = confidence of the sequencer that base at position l in read n was called correctly.
    """
    theta = 1 - 10 ** (-qualities / 10)  # dimension: N X L
    one_minus_theta = (1 - theta) / (size_alphabet - 1)

    log_theta = np.log(theta)
    log_theta = log_theta[:, :, np.newaxis]
    log_theta = np.tile(log_theta, (1, 1, size_alphabet))

    log_one_minus_theta = np.log(one_minus_theta)
    log_one_minus_theta = log_one_minus_theta[:, :, np.newaxis]
    log_one_minus_theta = np.tile(log_one_minus_theta, (1, 1, size_alphabet))

    final = np.einsum("NLB,NLB->NLB", log_theta, reads_seq_binary)
    final += np.einsum("NLB,NLB->NLB", log_one_minus_theta, 1 - reads_seq_binary)

    # if reads_list[n].seq_binary[l].sum(axis=0)=0 then "N" at position l then position l is ignored
    # there will be a zero in the row n,l
    # dimension: NxL
    all_N_pos = reads_seq_binary.sum(axis=2) > 0
    all_N_pos = all_N_pos[:, :, np.newaxis]
    all_N_pos = np.tile(all_N_pos, (1, 1, size_alphabet))
    # write zero where there is an "N"  in the position
    final[~all_N_pos] = 0

    return final  # dimension: NxLxB


def reads_list_to_array(reads_list):

    reads_binary = [reads_list[n].seq_binary for n in range(len(reads_list))]
    reads_binary_array = np.asarray(reads_binary)

    reads_weights = [reads_list[n].weight for n in range(len(reads_list))]
    reads_weights_array = np.asarray(reads_weights)

    return reads_binary_array, reads_weights_array


def load_fasta_and_qualities(fname_fasta, fname_qualities, alphabet,unique_modus):

    with open(fname_qualities, "rb") as f:
        qualities = np.load(f, allow_pickle=True)

    # go through each sequence in fasta file
    reads_list = []
    for idx, seq in enumerate(SeqIO.parse(fname_fasta, "fasta")):
        reads_list.append(Read(seq.seq, seq.id))
        reads_list[-1].seq2binary(alphabet)
        reads_list[-1].phred_quality_score = qualities[idx]

    if unique_modus:
        reads_list = unique_reads_list_qualities(reads_list)
    qualities = get_qualities(reads_list)

    return reads_list, qualities


def unique_reads_list_qualities(reads_list):
    # test for unique reads_list
    for i, temp_read in enumerate(reads_list):
        if temp_read.weight > 0.0:
            for j in range(i + 1, len(reads_list)):
                hd = hamming(
                    temp_read.seq_string, reads_list[j].seq_string
                )
                if hd == 0:
                    temp_read.weight += 1
                    temp_read.identical_reads.append(reads_list[j].id)
                    temp_read.phred_quality_score = (
                        temp_read.phred_quality_score
                        + reads_list[j].phred_quality_score
                    ) / 2
                    reads_list[j].weight -= 1

    # keep only unique reads_list
    reads_list = [read for read in reads_list if read.weight > 0]
    return reads_list

def load_reads_fasta(fname_fasta, alphabet,unique_modus):
    # go through each sequence in fasta file
    reads_list = []
    for idx, seq in enumerate(SeqIO.parse(fname_fasta, "fasta")):
        reads_list.append(Read(seq.seq, seq.id))
        reads_list[-1].seq2binary(alphabet)

    if unique_modus:
        reads_list = unique_reads_list(reads_list)

def unique_reads_list(reads_list):
    # test for unique reads_list
    for i, temp_read in enumerate(reads_list):
        if temp_read.weight > 0.0:
            for j in range(i + 1, len(reads_list)):
                hd = hamming(
                        temp_read.seq_string, reads_list[j].seq_string
                )
                if hd == 0:
                    temp_read.weight += 1
                    temp_read.identical_reads.append(reads_list[j].id)
                    reads_list[j].weight -= 1

    # keep only unique reads_list
    reads_list = [read for read in reads_list if read.weight > 0]
    return reads_list


def check_qualities(qualities):
    return np.where(qualities == 0, 2, qualities)


def get_qualities(reads_list):
    qualities = [temp_read.phred_quality_score for temp_read in reads_list]
    qualities = np.asarray(qualities)
    return check_qualities(qualities)


def load_reference_seq(reference_file, alphabet):
    for seq in SeqIO.parse(reference_file, "fasta"):
        return reference2binary(seq.seq, alphabet), seq.id


def reference2binary(reference_seq, alphabet):
    length_seq = len(reference_seq)
    reference_table = np.zeros((length_seq, len(alphabet)))
    for base_position, base in enumerate(str(reference_seq)):
        reference_table[base_position][alphabet.index(base)] = 1
    return reference_table
