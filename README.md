What is ShoRAH?
===============
[![Build Status](https://travis-ci.org/ozagordi/shorah.svg?branch=master)](https://travis-ci.org/ozagordi/shorah)

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
| `amplian.py`   | amplicon based analysis                                             |
| `dec.py`       | local error correction based on diri_sampler                        |
| `diri_sampler` | Gibbs sampling for error correction via Dirichlet process mixture   |
| `contain`      | removal of redundant reads                                          |
| `mm.py`        | maximum matching haplotype construction                             |
| `freqEst`      | EM algorithm for haplotype frequency                                |
| `snv.py`       | detects single nucleotide variants, taking strand bias into account |
| `shorah.py`    | wrapper for everything                                              |

## Citation
If you use shorah, please cite the application note paper _Zagordi et al._ on
[BMC Bioinformatics](http://www.biomedcentral.com/1471-2105/12/119).

## General usage
### Dependencies
shorah requires the following pieces of software:

1. **Python 2**, which is generally available on most Unix-like system. The required dependencies are:

   a) **Biopython**, which can be downloaded using pip or anaconda

2. **Perl**, for some scripts

3. **zlib**, which is used by the bundled samtools for compressing bam files

4. **pkg-config**, for discovering dependencies, which most Unix-like systems include

5. **GNU scientific library**, for random number generation

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

	/opt/local/bin/virtualenv-2.7 /usr/local/shorah

where `/opt/local/bin/virtualenv-2.7` is the virtualenv command for python 2.7 on MacPorts. Now install
the python dependencies:

	/usr/local/shorah/bin/pip install Biopython

Now call the `configure` script from the shorah tarball, taking care to specify the **absolute** path of the
python interpreter (or the relative one if it is in your `PATH`), as this gets inserted into the shebang line of all python scripts:

	./configure --prefix=/usr/local/shorah PYTHON=/usr/local/shorah/bin/python2.7

The configure script finds the dependencies using pkg-config. Once it completes, run:

	make -j4

where `4` specifies the number of compilation threads to use. Finally, after compilation, install using:

	make install

All the programs should now be located in `/usr/local/shorah/bin`.

### Boostrapping from git
If you opted to clone the git repository instead of downloading a prepared tarball, you will need to bootstrap
the configure script:

	autoreconf -vif

After this, you can run the `configure` script as described previously.


#### Windows users
You can install and run `shorah` with [Cygwin](http://www.cygwin.com).
Please see the relevant paragraph on the
[documentation page](http://cbg-ethz.github.io/shorah/).

### Run
The input is a sorted bam file. Analysis can be performed in local or global
mode.

#### Local analysis
The local analysis alone can be run invoking `dec.py` or `amplian.py` (program
for the amplicon mode). They work by cutting window from the multiple sequence
alignment, invoking `diri_sampler` on the windows and calling `snv.py` for the
SNV calling. See the
[`README`](https://github.com/cbg-ethz/shorah/blob/master/examples/amplicon_test/README.md)
file in directory
[`amplicon_test`](https://github.com/cbg-ethz/shorah/blob/master/examples/amplicon_test/).

#### Global analysis
The whole global reconstruction consists of the following steps:

1. error correction (*i.e.* local haplotype reconstruction);
2. SNV calling;
3. removal of redundant reads;
4. global haplotype reconstruction;
5. frequency estimation.

These can be run one after the other, or one can invoke `shorah.py`, that runs
the whole process from bam file to frequency estimation and SNV calling.
