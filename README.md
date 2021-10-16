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

More information [here](http://cbg-ethz.github.io/shorah).

---

The software suite ShoRAH (Short Reads Assembly into Haplotypes) consists of
several programs, the most imporant of which are:

| Tool           | What it does                                                        |
| -------------- | ------------------------------------------------------------------- |
| `shorah`       | wrapper for everything                                              |
|`shorah amplicon`| amplicon based analysis                                            |
|`shorah shotgun`| shotgun sequencing analysis                                         |
| `shorah snv`   | detects single nucleotide variants, taking strand bias into account |
| `b2w`          | splitting shotgun sequencing .BAM into multiple overlapping windows |
| `diri_sampler` | Gibbs sampling for error correction via Dirichlet process mixture   |
| `fil`          | strand bias test                                                    |

## Citation
If you use shorah, please cite the application note paper _Zagordi et al._ on
[BMC Bioinformatics](http://www.biomedcentral.com/1471-2105/12/119).

## General usage
### Pre-built packages
ShoRAH and its dependencies are all
[available in bioconda](https://bioconda.github.io/recipes/shorah/README.html).
We strongly advise you to install this package for a hassle-free experience.

	conda install shorah

### Dependencies
shorah requires the following pieces of software:

1. **Python 3** The required dependencies are:

   a) **Biopython**, and 
   b) **NumPy**.
   These packages can be downloaded using pip or conda

2. **HTSlib** which is used to access bam/cram/sam and fasta files.

3. **zlib**, which is used by HTSlib for compressing bam files

3. **pkg-config**, for discovering dependencies, which most Unix-like systems include

4. **Boost C++ library**, for random number generation

In addition, if you want to bootstrap the git version of shorah instead of using the provided tarballs,
you will need the GNU Autotools:

1. **Autoconf** 2.69

2. **Automake** 1.15

3. **m4**, which most Unix-like system include

### Installation
We strongly recommend you use one of the versioned tarballs from the releases page. ShoRAH uses Autoconf
and Automake, and these tarballs include all necessary scripts and files required for installation, whereas
the git tree only contains the bare minimum of files required for bootstrapping.

Further, we strongly recommend you use a virtualenv for python installation that shares the same directory
root as where you'd like to install shorah to. Not using a virtualenv means that the python dependencies will
not be located in the installation root, which will likely require you to specify `PYTHONPATH`, making the
installation more brittle.

Say for instance, you would like to install shorah to `/usr/local/shorah`. The first step consists of installing
the required python dependencies. Create a virtualenv:

	/opt/local/bin/virtualenv-3.6 /usr/local/shorah

where `/opt/local/bin/virtualenv-3.6` is the virtualenv command for python 3.6 on MacPorts. Now install
the python dependencies:

	/usr/local/shorah/bin/pip install Biopython numpy

Now call the `configure` script from the shorah tarball, taking care to specify the **absolute** path of the
python interpreter (or the relative one if it is in your `PATH`), as this gets inserted into the shebang line of all python scripts:

	./configure --prefix=/usr/local/shorah PYTHON=/usr/local/shorah/bin/python3.6

The configure script finds the dependencies using pkg-config. Once it completes, run:

	make -j4

where `4` specifies the number of compilation threads to use. Finally, after compilation, install using:

	make install

All the programs should now be located in `/usr/local/shorah/bin`.

#### Boostrapping from git
If you opted to clone the git repository instead of downloading a prepared tarball, you will need to bootstrap
the configure script:

	autoreconf -vif -I m4

After this, you can run the `configure` script as described previously.

#### Meson

For the developers who prefer this, it is alternatively possible to compile the C++ components with Meson and then install ShoRAH with pip in development mode using :

```bash
mkdir -p build
cd build
meson ../
ninja
cd ..
pip3 install -e .
```

#### Windows users
Since Windows 10, Microsoft provides the 
[Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/about) 
which enables a linux environment under Windows.

Users of older versions of Windows can can install and run `shorah` with [Cygwin](http://www.cygwin.com).
Please see the relevant paragraph on the
[documentation page](http://cbg-ethz.github.io/shorah/).

### Run
The input is a sorted bam file. Analysis can be performed in local mode.

**Note:** Currently in ShoRAH2 global haplotype reconstruction is disabled.
The last version of ShoRAH1 with global mode is
[v1.1.3](https://github.com/cbg-ethz/shorah/tree/v1.1.3)
and can easily be obtained from
[bioconda](https://bioconda.github.io/recipes/shorah/README.html).
If you wish to perform global reconstruction, we suggest that you consider the
[global haplotype options](https://github.com/cbg-ethz/V-pipe/wiki/options#haplotype_reconstruction)
available in [V-pipe](https://cbg-ethz.github.io/V-pipe/)

#### Local analysis
The local analysis alone can be run invoking `shorah shotgun` or `shorah amplicon` (program
for the amplicon mode). They work by cutting window from the multiple sequence
alignment, invoking `diri_sampler` on the windows and calling `shorah snv` for the
SNV calling. See the
[`README`](https://github.com/cbg-ethz/shorah/blob/master/examples/shotgun_test/README.md)
file in directory
[`shotgun_test`](https://github.com/cbg-ethz/shorah/blob/master/examples/shotgun_test/)
and the
[`README`](https://github.com/cbg-ethz/shorah/blob/master/examples/amplicon_test/README.md)
file in directory
[`amplicon_test`](https://github.com/cbg-ethz/shorah/blob/master/examples/amplicon_test/).


## Coding style
All changes to the C++ code in `src/cpp` should always be formatted according to the included `.clang-format` style by doing

	clang-format -style=file -i src/cpp/*.[ch]pp

in the root of the repository.

All changes to the python code in `src/shorah` should always be formatted conforming to the [PEP 8](https://www.python.org/dev/peps/pep-0008/) style guide. To this end, we advise to use [autopep8](https://pypi.python.org/pypi/autopep8).

## Development/CI with Docker
The following command will run the CI locally within Docker. 
```bash
docker run -w="/usr/app" -it $(docker build -q .) bash
```

## Contact

ShoRAH is maintained as part of the [V-pipe virus NGS pipeline](https://cbg-ethz.github.io/V-pipe/)
and you can [easily reach out the developers on its website](https://cbg-ethz.github.io/V-pipe/contact).

You can also report your problems in the [issue tracker on GitHub](https://github.com/cbg-ethz/shorah/issues).
