#!/usr/bin/env bash
#samtools faidx HIV_amplicon.fasta
samtools view -S -b amplicon_reads.sam | samtools sort - amplicon_reads_sorted
#samtools index amplicon_reads_sorted.bam
