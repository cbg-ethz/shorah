What is ShoRAH?
===============
[![Build Status](https://travis-ci.org/cbg-ethz/shorah.svg?branch=master)](https://travis-ci.org/cbg-ethz/shorah)
[![Bioconda package](https://img.shields.io/conda/dn/bioconda/shorah.svg?label=Bioconda)](https://bioconda.github.io/recipes/shorah/README.html)
[![Docker container](https://quay.io/repository/biocontainers/shorah/status)](https://quay.io/repository/biocontainers/shorah)


ShoRAH is an open source project for the analysis of next generation sequencing
data. It is designed to analyse genetically heterogeneous samples. Its tools
are written in different programming languages and provide error correction,
haplotype reconstruction and estimation of the frequency of the different
genetic variants present in a mixed sample.

---

The software suite ShoRAH (Short Reads Assembly into Haplotypes) consists of
several programs, the most imporant of which are:

| Tool           | What it does                                                        |
| -------------- | ------------------------------------------------------------------- |
| `shorah`       | wrapper for everything                                              |
|`shorah shotgun`| shotgun sequencing analysis                                         |
| `shorah snv`   | detects single nucleotide variants, taking strand bias into account |
| `b2w`          | splitting shotgun sequencing .BAM into multiple overlapping windows |
| `diri_sampler` | Gibbs sampling for error correction via Dirichlet process mixture   |
| `fil`          | strand bias test                                                    |

### Dependencies
shorah requires the following pieces of software:

1. **Python 3**

2. **HTSlib** which is used to access bam/cram/sam and fasta files.

3. **Boost C++ library**, for random number generation.

### Installation
For installation miniconda is recommended: https://docs.conda.io/en/latest/miniconda.html
We recommend to install ShoRAH in a clean conda enviroment:
`conda create --name env_shorah python=3.9`

`conda activate env_shorah`

Next install, **HTSlib** and **Boost C++ library** using conda:

`conda install -c bioconda htslib`

`conda install -c conda-forge boost`

Then install this git repository:

`pip install git+https://github.com/LaraFuhrmann/shorah@master `

You might need to downgrade your pip version:
`pip install pip==21.3.1`

### Example
To test your installation, we recommend running the program on `tests/data_1`.

If the sequencing amplicon strategy is known, we recommend using the amplicon-mode of the program, which takes as input the `<smth>.insert.bed` - file:
`shorah shotgun -b test_aln.cram -f test_ref.fasta -z scheme.insert.bed --sampler use_quality_scores`

If the sequencing quality scores are not trustable, the sequencing error parameters can also be learned:
`shorah shotgun -b test_aln.cram -f test_ref.fasta -z scheme.insert.bed --sampler learn_error_params`.

If there is no information on the sequencing amplicon strategy available, run:
`shorah shotgun -b test_aln.cram -f test_ref.fasta --sampler use_quality_scores`


## Development/CI with Docker
The following command in the root directory will let you interact with the project locally through Docker.
```bash
docker run --rm -w="/usr/app" -it $(docker build -q .) bash
```
Then run `poetry install --only-root` to install `shorah` inside Docker.

This is the same setup as used in the CI at [`.github/workflows/test.yaml`](.github/workflows/test.yaml).
