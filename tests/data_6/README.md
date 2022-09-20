### Sample files to test `shorah shotgun`

Use files in this directory to test shorah in shotgun mode. The reads data have been generated with V-pipe's benchmarking framework (simulated with parameters (seq_tech~illumina__seq_mode~shotgun__seq_mode_param~nan__read_length~90__genome_size~90__coverage~100__haplos~5@5@10@5@10@geom@0.75)

The reads stam from one single amplicon of length 90, meaning the reference is of length 90 and each read is of length 90bps.

To run ShoRAH's original Gibbs sampler use the following command:
```
poetry run shorah shotgun -f reference.fasta -b reads.shotgun.bam -w 90 --shorah
```

If only one window should cover the amplicon set in addtion the parameter `-s 1`.

To use the new inference method using the sequencing quality scores use:
```
poetry run shorah shotgun -f reference.fasta -b reads.shotgun.bam -w 90 --use_quality_scores --alpha 0.0001 --n_max_haplotypes 100 --n_mfa_starts 1 --conv_thres 0.0001
```
To use the new inference method learning the sequencing error parameter:
```
poetry run shorah shotgun -f reference.fasta -b reads.shotgun.bam -w 90 --learn_error_params --alpha 0.0001 --n_max_haplotypes 100 --n_mfa_starts 1 --conv_thres 0.0001
```

In the new inference method reads are filtered (and weighted respectively) such that only a set of unique reads are processed. This mode can be switch off by setting the parameter `--non-unique_modus`,  e.g.:
```
poetry run shorah shotgun -f reference.fasta -b reads.shotgun.bam -w 90 --learn_error_params --alpha 0.0001 --n_max_haplotypes 100 --n_mfa_starts 1 --conv_thres 0.0001 --non-unique_modus
