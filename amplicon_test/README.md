### Sample files to test `amplian.py`

Use files in this directory to test shorah in amplicon mode.
The sorted bam file has been obtained from the reads with
[`smalt`](http://www.sanger.ac.uk/resources/software/smalt/)
to align reads to the reference and samtools to index and
sort the alignment. The commands given are

    [user@host amplicon_test]$ smalt index -k 7 -s 2 ref reference.fasta
    [user@host amplicon_test]$ smalt map -f sam -o ampli.sam ref amplicon_reads.fastq
    [user@host amplicon_test]$ samtools view -S -b ampli.sam | samtools sort - ampli_sorted

You can then run `amplian` as follows
    [user@host amplicon_test]$ amplian.py -b ampli_sorted.bam -f reference.fasta

Some warning will be raised. The output files will be `SNV.txt` and `SNVs_0.010000_final.csv`.
Only SNVs exceeding the threshold of 5% are returned (see [paper](http://www.biomedcentral.com/1471-2164/14/501)).

#### Tips
You can inspect the alignment with `samtools tview`:

    [user@host amplicon_test]$ samtools index ampli_sorted.bam
    [user@host amplicon_test]$ samtools tview ampli_sorted.bam reference.fasta

Within tview, press ? for help.

----

`amplian` is still under development and thus should be considered in
beta more than the rest of ShoRAH (say, alpha).
