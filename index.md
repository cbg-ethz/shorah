---
layout: slate_mod
title: Shorah - Short Reads Assembly into Haplotypes
permalink: index.html 
---
What is ShoRAH
======
ShoRAH is an open source project for the analysis of next generation sequencing
data. It is designed to analyse genetically heterogeneous samples. Its tools
are written in different programming languages and provide error correction,
haplotype reconstruction and estimation of the frequency of the different
genetic variants present in a mixed sample.

---

The software suite ShoRAH (Short Reads Assembly into Haplotypes) consists of
several programs, the most imporant of which are:

`amplian.py`   - amplicon based analysis

`dec.py`       - local error correction based on diri_sampler

`diri_sampler` - Gibbs sampling for error correction via Dirichlet
process mixture

`contain`      - removal of redundant reads

`mm.py`        - maximum matching haplotype construction

`freqEst`      - EM algorithm for haplotype frequency

`snv.py`       - detects single nucleotide variants, taking strand bias into
account

`shorah.py`    - wrapper for everything

---

### Citation
If you use shorah, please cite the application note paper
on [BMC Bioinformatics](http://www.biomedcentral.com/1471-2105/12/119).

---

### Usage

#### Dependencies and installation
Please download and install:

- [Biopython](http://biopython.org/wiki/Download), following the online
  instructions.
- [GNU scientific library GSL](http://www.gnu.org/software/gsl/),
  installation is described in the included README and INSTALL files.
- [ncurses](http://www.gnu.org/software/ncurses/ncurses.html) is
  required by samtools. It is usually already included in Linux/Mac OS X.

Please note that these dependencies can be satisfied also using the package
manager of many operating system. For example
[MacPorts](http://www.macports.org/) on Mac OS X,
[yum](http://yum.baseurl.org/) on several linux installations and so on.


Type 'make' to build the C++ programs. This should be enough in most cases.
If your GSL installation is not standard, you might need to edit the relevant
lines in the `Makefile` (location `/opt/local/` is already included).

#### Run

The input is a sorted bam file. Analysis can be performed in local or global
mode.

##### [Local analysis](local.html)

The local analysis alone can be run invoking `dec.py` or `amplian.py` (program
for the amplicon mode). They work by cutting windows from the multiple sequence
alignment, invoking `diri_sampler` on the windows and calling `snv.py` for the
SNV calling.

##### [Global analysis](global.html)

The whole global reconstruction consists of the following steps:

1. error correction (*i.e.* local haplotype reconstruction);
2. SNV calling;
3. removal of redundant reads;
4. global haplotype reconstruction;
5. frequency estimation.

These can be run one after the other, or one can invoke `shorah.py`, that runs
the whole process from bam file to frequency estimation and SNV calling.

---

### Getting help
See the dedicated [page](help.html).

---

### Authors and Contributors
Niko Beerenwinkel

Arnab Bhattacharya

Nicholas Eriksson

Moritz Gerstung

Lukas Geyrhofer

Fabio Luciani

Kerensa McElroy

Osvaldo Zagordi
