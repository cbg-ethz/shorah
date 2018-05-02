What is ShoRAH?
===============
[![Build Status](https://travis-ci.org/cbg-ethz/shorah.svg?branch=master)](https://travis-ci.org/cbg-ethz/shorah) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/shorah/README.html)
ShoRAH is an open source project for the analysis of next generation sequencing
data. It is designed to analyse genetically heterogeneous samples. Its tools
are written in different programming languages and provide error correction,
haplotype reconstruction and estimation of the frequency of the different
genetic variants present in a mixed sample.

More information [here](http://cbg-ethz.github.io/shorah).

---

ShoRAH (Short Reads Assembly into Haplotypes) is an open source project for the analysis of sequencing data.
designed to analyse genetically heterogeneous samples. It is written in Python (with some C++
code for performance reasons) and provides error correction, haplotype reconstruction and
estimation of the frequency of the different genetic variants present in a mixed sample.

ShoRAH 2 is an API breaking update that removes the option to run a _global haplotype reconstruction_.
While this was motivated by the emergence of better tools for this difficult problem, a welcome
consequence was the simplification of the codebase (_e.g._ by removing some Perl legacy
code), resulting in a simpler interface and easier development.


ShoRAH2 can be run in two different modes:

|     mode     |                         Use                               |   old name   |
|:------------:|:---------------------------------------------------------:|:------------:|
| `amplicon`   | useful on small regions entirely covered by a single read | `amplian.py` |
| `shotgun`    | local error correction of reads covering a long region    |   `dec.py `  |

See `shorah -h` for details and the [documentation](http://cbg-ethz.github.io/shorah).

## Citation
If you use shorah, please cite the application note paper _Zagordi et al._ on
[BMC Bioinformatics](http://www.biomedcentral.com/1471-2105/12/119).

## General usage

The easiest way to install Shorah is via [bioconda](http://bioconda.github.io/recipes/shorah/README.html).

Users who want to manually build find the relevant instructions here below.

### Dependencies
shorah requires the following pieces of software:

1. **Python 2 or Python 3**, backward compatibility is provided as some current Linux distributions and OS X systems
are still using 2.x as default. The required dependencies are:

   a) **Biopython**, and
   b) **NumPy**.
   These packages can be downloaded using pip or anaconda

2. **zlib**, which is used by the bundled samtools for compressing bam files

3. **pkg-config**, for discovering dependencies, which most Unix-like systems include

4. **GNU scientific library**, for random number generation

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

### Boostrapping from git
If you opted to clone the git repository instead of downloading a prepared tarball, you will need to bootstrap
the configure script:

	autoreconf -vif

After this, you can run the `configure` script as described previously.


### Run
The input is a sorted bam file. Analysis can be performed in shotgun or amplicon mode.


## Coding style
All changes to the C++ code in `src/cpp` should always be formatted according to the included `.clang-format` style by doing

	clang-format -style=file -i src/cpp/*.[ch]pp

in the root of the repository.

All changes to the python code in `src/shorah` should always be formatted conforming to the [PEP 8](https://www.python.org/dev/peps/pep-0008/) style guide. To this end, we advise to use [autopep8](https://pypi.python.org/pypi/autopep8).

## Collaborators

- Niko Beerenwinkel
- Arnab Bhattacharya
- Nicholas Eriksson
- Moritz Gerstung
- Lukas Geyrhofer
- Susana Posada Céspedes
- David Seifert
- Ivan Topolsky
- Osvaldo Zagordi

We are grateful to Manuel Holtgrewe (@holtgrewe) and George Kettleborough (@georgek) for helpful commits.
