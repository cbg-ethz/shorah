def sequence_compare(seq_a, seq_b, outfile=None):
    len1 = len(seq_a)
    len2 = len(seq_b)
    mismatches = []
    for pos in range (0, min(len1, len2)) :
        if str(seq_a[pos]) != str(seq_b[pos]):
            mismatches.append('|')
        else:
            mismatches.append('-')
    if outfile==None:
        print(seq_a)
        print("".join(mismatches))
        print(seq_b)
    else:
        outfile.write(str(seq_a)+'\n')
        outfile.write("".join(mismatches)+'\n')
        outfile.write(str(seq_b)+'\n')
    return outfile

'''
Example:

sequence_compare(list_of_clusters[0].haplotype_seq, reads[1])
distance(list_of_clusters[0].haplotype_seq, reads[1])
matches(list_of_clusters[0].haplotype_seq, reads[1])

Out:
GTTGAAAAGGTGTTGAGGGCGGAGAAATGCAAGTTATTACAAAAGTTAACGTAACAAAGAATCTGGTAGGGGTGAGTT
-----|||-------------------------------------------------------------||-------
GTTGAGGGGGTGTTGAGGGCGGAGAAATGCAAGTTATTACAAAAGTTAACGTAACAAAGAATCTGGTAGAAGTGAGTT

73.0

'''

def print_cluster_overview(cluster_list, reads_list, outfile=None):
    for i,cluster in enumerate(cluster_list):
        read_name= [reads_list[idx].metadata['id'] for idx in cluster.reads_idx_list]
        if outfile==None:
            print('In cluster ',i, 'are reads: ', read_name)
        else:
            outfile.write('In cluster '+str(i)+ ' are reads: '+str(read_name)+'\n')
    return outfile



'''

#go through unique reads
idx_read = 0
read = reads[0]
# sample cluster for read

# run through the populated classes to assign a probability
log_P_list =[]
P_list = []
for temp_cluster in list_of_clusters:
    if not (temp_cluster.size ==1 and temp_cluster.reads_list == [idx_read]): # exclude cluster where the read is assigned to if it only contains this read
        x1, x2 = log_P(theta, B, idx_read, temp_cluster)
        log_P_list.append(x1)
        P_list.append(x2)

# probablity of creating a new cluster
x1, x2 = log_P_reference(theta, gamma, alpha, B, idx_read)
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
        P_list[ll]=np.exp(log_P_list[ll])
        print('with weight ' +str(P_list[ll])+' to cluster '+ str(ll))

this_cluster = one_shot_discrete(len(log_P_list),P_list)[0]

print('extracted class is ' + str(this_cluster))

#update the clusters
if this_cluster == len(log_P_list)-1:
    print('read ' ,idx_read, ' is assigned to new cluster')
    list_of_clusters.append(cluster([idx_read]))
else:
    print('read ',idx_read, ' is assigned to cluster ', this_cluster)
    list_of_clusters[this_cluster].reads_list.append(reads_idx)

n_cluster=len(list_of_clusters)
dt=0.0
dk1=0
hapbases=0
# sample haplotypes
for cluster_temp in list_of_clusters:
    # sample haplotypes
    sample_hap(cluster_temp)
    # update distances of all reads to the new haplotype (=cluster center)

    cluster_temp.distance2reads = [distance(read, cluster_temp.haplotype_seq) for read in reads]
    cluster_temp.matches2reads = [matches(read, cluster_temp.haplotype_seq) for read in reads]
    cluster_temp.matches2reference = matches(reference_seq, cluster_temp.haplotype_seq)
    cluster_temp.distance2reference = distance(reference_seq, cluster_temp.haplotype_seq)

    for i in cluster_temp.reads_list: dt+=cluster_temp.matches2reads[i]*reads[i].metadata['weight']

    dk1+=cluster_temp.matches2reference
    hapbases+=cluster_temp.matches2reference+cluster_temp.distance2reference

b_alpha=dt+(eps1*eps2*totbases)
b_beta=(totbases - dt) + eps2 * totbases * (1 - eps1)
theta = np.random.beta(a=b_alpha, b= b_beta,size=1)
theta=np.max(np.min(theta, 0.9999),0.0001)

gamma= dk1/hapbases
'''
