shorah
======

Repo for the software suite ShoRAH (Short Reads Assembly into Haplotypes)

Full online documentation available [here]
(https://wiki-bsse.ethz.ch/display/ShoRAH/Documentation)

ShoRAH consists of several programs, the most imporant of which are

>* `dec.py`       - local error correction based on diri_sampler
>* `diri_sampler` - Gibbs sampling for error correction via Dirichlet
>process mixture
>* `contain`      - removal of redundant reads
>* `mm.py`        - maximum matching haplotype construction
>* `freqEst`      - EM algorithm for haplotype frequency
>* `snv.py`       - detects single nucleotide variants, taking strand bias into
>account
>* `shorah.py`    - wrapper for everything

GENERAL USAGE:
======

### Install

type 'make' to build the C++ programs. This should be enough in most cases.

### Run

The whole process can be run one step after the other, or one can
invoke `shorah.py`, that runs the whole process from bam file to
frequency estimation and SNV calling. The local analysis alone
can be run invoking	`dec.py` or directly `diri_sampler`.
