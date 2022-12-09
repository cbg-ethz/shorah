import numpy as np
from skbio.sequence.distance import hamming
from skbio import Sequence
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def correct_reads(state_curr_dict, fname_output_corr):
    # TODO compute posterior for each read
    haplo_list = [key for key in state_curr_dict.keys() if key.startswith('haplotype')]
    records = []
    for k in range(len(haplo_list)):
        for read_id in state_curr_dict['assignedReads'+str(k)]:
            records.append(SeqRecord(Seq(str(state_curr_dict['haplotype'+str(k)])), id= read_id, description='|posterior=1'))

        SeqIO.write(records, fname_output_corr, "fasta")


def haplotypes_to_fasta(state_curr_dict, output_dir):
    haplo_list = [key for key in state_curr_dict.keys() if key.startswith('haplotype')]
    records = []
    for k in range(len(haplo_list)):
        # ave_reads = assignedReads
        #head = 'hap_'+str(k)+' | posterior='+str(state_curr_dict['HapPosterior'+str(k)])+' ave_reads='+str(state_curr_dict['weight'+str(k)])
        head = 'hap_'+str(k)+' | posterior='+str(1.0)+' ave_reads='+str(state_curr_dict['weight'+str(k)])
        records.append(SeqRecord(Seq(str(state_curr_dict['haplotype'+str(k)])), id= 'haplotype'+str(k), description=head))

        SeqIO.write(records, output_dir, "fasta")
