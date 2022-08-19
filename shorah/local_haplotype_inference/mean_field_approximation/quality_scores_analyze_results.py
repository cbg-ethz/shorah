import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import quality_scores_update_eqs as update_eqs


def correct_reads(state_curr_dict, fname_output_corr):
    # TODO compute posterior for each read
    haplo_list = [key for key in state_curr_dict.keys() if key.startswith("haplotype")]
    records = []
    for k in range(len(haplo_list)):
        for read_id in state_curr_dict["assignedReads" + str(k)]:
            records.append(
                SeqRecord(
                    Seq(state_curr_dict["haplotype" + str(k)]),
                    id=read_id,
                    description="|posterior=1",
                )
            )

        SeqIO.write(records, fname_output_corr, "fasta")


def haplotypes_to_fasta(state_curr_dict, output_dir):
    haplo_list = [key for key in state_curr_dict.keys() if key.startswith("haplotype")]
    records = []
    for k in range(len(haplo_list)):
        # ave_reads = assignedReads
        ave_reads = state_curr_dict["weight" + str(k)]
        if ave_reads==0:
            # this haplotype will not be reported as there are no reads
            # supporting it.
            continue

        head = (
            "hap_"
            + str(k)
            + " | posterior="
            + str(state_curr_dict["approximatePosterior" + str(k)])
            + " ave_reads="
            + str(ave_reads)
        )
        records.append(
            SeqRecord(
                Seq(state_curr_dict["haplotype" + str(k)]),
                id="haplotype" + str(k),
                description=head,
            )
        )

        SeqIO.write(records, output_dir, "fasta")


def summarize_results(
    state_curr_dict,
    alphabet,
    reads_seq_binary,
    reads_weights,
    reads_list,
    reads_log_error_proba,
    reference_binary,
):

    mean_z = state_curr_dict["mean_cluster"]
    mean_h = state_curr_dict["mean_haplo"]
    mean_log_gamma = state_curr_dict["mean_log_gamma"]

    dict_summary = {}

    unique_haplo = get_unique_haplotypes(mean_h, alphabet)
    unique_cluster = merge_cluster_assignments(mean_z, unique_haplo)
    unique_mean_h = update_eqs.update_mean_haplo(
        reads_weights,
        reference_binary,
        reads_log_error_proba,
        unique_cluster,
        mean_log_gamma,
    )
    unique_haplo_posterior = compute_unique_haplo_posterior(
        unique_mean_h, unique_haplo, alphabet
    )

    for k in range(len(unique_haplo)):
        # write info for haplotype k
        dict_summary.update(
            {
                "haplotype" + str(k): unique_haplo[k][0],
                "approximatePosterior" + str(k): unique_haplo_posterior[k],
                "assignedUniqueReads"
                + str(k): get_assigned_unique_reads(unique_cluster, k, reads_list),
                "assignedReads"
                + str(k): get_assigned_all_reads(unique_cluster, k, reads_list),
                "weight" + str(k): get_cluster_weight(unique_cluster, k, reads_list),
            }
        )

    return dict_summary


def get_haplotype(mean_h_k, alphabet):
    haplotype_sequence = []
    for position in range(mean_h_k.shape[0]):
        haplotype_sequence.append(alphabet[np.argmax(mean_h_k[position])])

    return "".join(haplotype_sequence)


def compute_unique_haplo_posterior(unique_mean_h, unique_haplo_var, alphabet):
    posterior = []
    for idx_haplo, haplo in enumerate(unique_haplo_var):
        proba = 1
        for position, base in enumerate(haplo[0]):
            idx_base = alphabet.index(base)
            proba = proba * unique_mean_h[idx_haplo][position][idx_base]
        posterior.append(proba)
    return posterior


def get_assigned_unique_reads(mean_z, k, reads_list):
    reads = []
    for n in range(len(reads_list)):
        max_val = np.max(mean_z[n])

        if k in set([i for i in range(len(mean_z[n])) if mean_z[n][i] >= max_val]):
            reads.append(reads_list[n].id)
    return reads


def get_assigned_all_reads(mean_z, k, reads_list):
    reads = []
    for n in range(len(reads_list)):
        max_val = np.max(mean_z[n])

        if k in set([i for i in range(len(mean_z[n])) if mean_z[n][i] >= max_val]):
            reads.append(reads_list[n].id)
            reads = reads + reads_list[n].identical_reads
    return reads


def list_duplicates(seq):
    from collections import defaultdict

    tally = defaultdict(list)
    for i, item in enumerate(seq):
        tally[item].append(i)
    return ((key, locs) for key, locs in tally.items() if len(locs) > 0)


def get_unique_haplotypes(mean_h, alphabet):
    haplo_seqs = []
    for k in range(mean_h.shape[0]):
        haplo_seqs.append(get_haplotype(mean_h[k], alphabet))
    return sorted(list_duplicates(haplo_seqs))


def merge_cluster_assignments(mean_z, unique_haplo):
    new_cluster = np.zeros((mean_z.shape[0], len(unique_haplo)))
    for n in range(mean_z.shape[0]):
        for k in range(len(unique_haplo)):
            for i in unique_haplo[k][1]:
                new_cluster[n][k] += mean_z[n][i]

    return new_cluster


def get_cluster_weight(unique_cluster, k, reads_list):
    weight = 0
    for n in range(len(reads_list)):
        max_val = np.max(unique_cluster[n])

        if k in set(
            [
                i
                for i in range(len(unique_cluster[n]))
                if unique_cluster[n][i] >= max_val
            ]
        ):
            weight += reads_list[n].weight

    return weight
