#!/bin/bash

for r in 'alt' 'ref'; do
  wgsim -N 10000 -1 250 -2 250 -e 0.01 -S 349 "${r}.fasta" "${r}1.fq" "${r}2.fq" > "wginfo_${r}.txt"
done

rm reads1.fastq reads2.fastq

for r in 'alt' 'ref'; do
  cat ${r}1.fq >> reads1.fastq
  cat ${r}2.fq >> reads2.fastq
  rm ${r}1.fq ${r}2.fq
done

bwa index ref.fasta
bwa mem -t 3 ref.fasta reads1.fastq reads2.fastq | samtools view -u | samtools sort -o aln_sorted.bam
rm ref.fasta.*
samtools index aln_sorted.bam
