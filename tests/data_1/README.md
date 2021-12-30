### Sample files to test `shorah shotgun`

Use files in this directory to test shorah in shotgun mode.
The reads data comes from the 
[test-data](https://github.com/cbg-ethz/V-pipe/tree/master/testdata/2VM-sim/20170904/raw_data)
of [V-pipe](https://cbg-ethz.github.io/V-pipe/)
and has been processed with the pipeline using the `bwa` 
[option](https://github.com/cbg-ethz/V-pipe/wiki/options#aligner):

```ini
[general]
aligner = bwa
```

The sorted bam file has been further compressed with samtools for space saving:

    [user@host shotgun_test]$ samtools view -T test_ref.fasta -C -O cram,embed_ref,use_bzip2,use_lzma,level=9,seqs_per_slice=1000000 -o test_aln.cram V-pipe/work/samples/2VM-sim/20170904/alignments/REF_aln.bam

You can then run `shorah shotgun` as follows

    [user@host shotgun_test]$ shorah shotgun -b test_aln.cram -f test_ref.fasta

The output files will be `snv/SNVs_0.010000_final.vcf` and `snv/SNVs_0.010000_final.csv`.
