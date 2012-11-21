---
layout: slate_posts
title: How to prepare the input
permalink: input.html 
---
### First of all, align the reads
Since the mass adoption of next-generation sequencing, the bioinformatics
community has produced an enormous number of read mappers (aligners). An
incomplete list is [here](http://lh3lh3.users.sourceforge.net/NGSalign.shtml).

Whatever is your chosen aligner, it will most likely have the option to
output a file in [SAM format][samtools]. From the alignment in SAM format, you
need to create a sorted bam alignment.

### Use `samtools` to convert and sort the alignment
The most common set of tools to manipulate SAM alignments is [samtools]. This
is shipped with ShoRAH, so you can navigate to the `samtools` directory and
install it from there, or check if there is a more recent release, download
and install it.

In the following, we assume that a sam alignment file `my_reads.sam` has been
created by aligning reads to a reference file `reference.fasta`. The following
converts to sorted bam

    [user@host]$ samtools view -b -T reference.fasta my_reads.sam | samtools sort - my_reads_sorted

The output is a sorted bam file called `my_reads_sorted.bam`.

The above code is a shortcut to run conversion and sorting. The two commands
would read

    [user@host]$ samtools view -b -T reference.fasta -o my_reads.bam my_reads.sam
    [user@host]$ samtools sort my_reads.bam my_reads_sorted

[samtools]: http://samtools.sourceforge.net/ "samtools"

---

Go back [home](index.html)
