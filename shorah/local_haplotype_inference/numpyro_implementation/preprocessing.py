import jax.numpy as jnp
import skbio


def seq_mapping(seq, alphabet):
    # Coding scheme
    # 0:A, 1:C, 2:G, 3:T 4:- (NOT YET:, 5:N)
    mapped = []
    for base in seq:
        idx = alphabet.find(base)
        if idx != -1:
            mapped.append(idx)
        else:
            mapped.append(jnp.nan)
    return jnp.array(mapped)


def fasta2ref(fref_in, alphabet):
    # Coding scheme
    # 0:A, 1:C, 2:G, 3:T 4:- (NOT YET:, 5:N)
    for seq in skbio.io.read(fref_in, format="fasta"):
        ref = seq_mapping(str(seq), alphabet)
    # print('reference ', str(seq))
    return ref, seq.metadata['id']


def fasta2reads(freads_in, alphabet):
    # Coding scheme
    # 0:A, 1:C, 2:G, 3:T 4:- (NOT YET:, 5:N)
    reads_mapped = []
    read_id = []
    # print('reads')
    for seq in skbio.io.read(freads_in, format="fasta"):
        # print(str(seq))
        read_id.append(seq.metadata['id'])
        reads_mapped.append(seq_mapping(str(seq), alphabet))
    return jnp.array(reads_mapped), read_id
