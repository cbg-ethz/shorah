import skbio
import numpy as np
from skbio.sequence.distance import hamming
from skbio import Sequence


class Read:
    def __init__(self, seq_string, seq_id):
        self.seq_string = seq_string
        self.weight = 1
        self.id = seq_id
        self.seq_binary = []
        self.identical_reads = []
        self.n_non_N = len(seq_string) - seq_string.count("N")

    def seq2binary(self, alphabet):
        length_seq = len(self.seq_string)
        seq_table = np.zeros((length_seq, len(alphabet)))
        for base_position, base in enumerate(str(self.seq_string)):
            if alphabet.find(base) >= 0:
                seq_table[base_position][alphabet.index(base)] = 1

        self.seq_binary = seq_table


def reads_list_to_array(reads_list):

    reads_binary = [reads_list[n].seq_binary for n in range(len(reads_list))]
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
    reads_list = []
    for seq in skbio.io.read(reads_fasta_file, format="fasta"):
        reads_list.append(Read(str(seq), seq.metadata["id"]))
        reads_list[-1].seq2binary(alphabet)
    # unique reads_list
    reads_list = unique_reads_list(reads_list)

    return reads_list


def load_bam2reads_list(reads_fasta_file, alphabet):
    # go through each sequence in fasta file
    reads_list = []
    for seq in skbio.io.read(reads_fasta_file, format="fasta"):
        reads_list.append(Read(str(seq), seq.metadata["id"]))
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
        # print(read_dict['flag'])
        # break
        reads_list.append(Read(str(read_dict["seq"]), read_dict["name"]))
        reads_list[-1].seq2binary(alphabet)
    # unique reads_list
    reads_list = unique_reads_list(reads_list)

    return reads_list


def unique_reads_list(reads_list):
    # test for unique reads_list
    for i, temp_read in enumerate(reads_list):
        if temp_read.weight > 0.0:
            for j in range(i + 1, len(reads_list)):
                hd = hamming(
                    Sequence(temp_read.seq_string), Sequence(reads_list[j].seq_string)
                )
                if hd == 0:
                    temp_read.weight += 1
                    temp_read.identical_reads.append(reads_list[j].id)
                    reads_list[j].weight -= 1

    # keep only unique reads_list
    reads_list = [read for read in reads_list if read.weight > 0]
    return reads_list


def load_reference_seq(reference_file):
    for seq in skbio.io.read(reference_file, format="fasta"):
        return seq, seq.metadata["id"]


def reference2binary(reference_seq, alphabet):
    length_seq = len(reference_seq)
    reference_table = np.zeros((length_seq, len(alphabet)))
    for base_position, base in enumerate(str(reference_seq)):
        reference_table[base_position][alphabet.index(base)] = 1
    return reference_table
