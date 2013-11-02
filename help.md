---
layout: slate_mod
title: Getting help 
permalink: help.html 
---
### Get help
Despite our efforts to make the installation procedure as easy as possible,
the variety of computers, operating systems, distributions and libraries make
the installation outcome very difficult to predict. Should you have any
problem, you might

- check issues reported on this page,
- ask your local administrator,
- open an [issue](https://github.com/ozagordi/shorah/issues?state=open),
- [contact us](mailto:shorah@bsse.ethz.ch).

---

### Problems occurring during compilation
####GSL
The Makefile has been written to match typical location for the installation of
GSL on Linux and Mac OS X system. In case of different installation of the GSL,
please edit CFLAGS and XLIBS in the Makefile. The [online documentation of GSL][1]
provides a [sample program][2]. Try to compile this; it will give information
on the correct setting of the flags.

---

### Problems occurring during execution
#### Memory footprint
We work to keep the memory requirements at bay. Nevertheless, for large samples
the required RAM might exceed 4 Gb. For example, a run of `diri_sampler` on a
single window of 65000 reads from 454 required 4.1 Gb of RAM. If `dec.py` is
used, this will run several jobs at the same time (parallelization), thus
increasing the memory needs.

#### BAM files give me problems
We have heard from some users errors that might be due to colons `:` in the
reference names.

_reported by KMD_

#### `diri_sampler` and memory leak
This seems to be solved now (thanks to
[Manuel Holtgrewe](https://github.com/ozagordi/shorah/pull/1)).

It happened that `diri_sampler` failed to run with a cryptic message that reads
like: `memory corruption: 0x0000000000. . .`
Usually, changing the name of the input file to something shorter was enough.
For obscure reasons, it seemed to fail when the file name is 18 letter plus
the suffix `.far`.

---

Go back [home](index.html)
[1]: http://www.gnu.org/software/gsl/manual/html_node/index.html
[2]: http://www.gnu.org/software/gsl/manual/html_node/An-Example-Program.html
