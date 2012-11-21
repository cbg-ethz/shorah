---
layout: slate_posts
title: Getting help 
permalink: help.html 
---
Despite our efforts to make the installation procedure as easy as possible,
the variety of computers, operating systems, distributions and libraries make
the installation outcome very difficult to predict. Should you have any
problem, you might

- check issues reported on this page,
- ask your local administrator,
- [contact us](mailto:shorah@bsse.ethz.ch).

## Problems occurring during compilation
###GSL
The Makefile has been written to match typical location for the installation of
GSL on Linux and Mac OS X system. In case of different installation of the GSL,
please edit CFLAGS and XLIBS in the Makefile.


## Problems occurring during execution
### Memory footprint
We work to keep the memory requirements at bay. Nevertheless, for large samples
the required RAM might exceed 4 Gb. For example, a run of `diri_sampler` on a
single window of 65000 reads from 454 required 4.1 Gb of RAM. If `dec.py` is
used, this will run several jobs at the same time (parallelization), thus
increasing the memory needs.

### BAM files give me problems
We have heard from some users errors that might be due to colons `:` in the
reference names.

_reported by KMD_

### `diri_sampler` and memory leak
It happens that `diri_sampler` fails to run with a cryptic message that reads
like: `memory corruption: 0x0000000000. . .`
Usually, changing the name of the input file to something shorter is enough.
For obscure reasons, it seems to fail when the file name is 18 letter plus
the suffix `.far`. If you have the faintest idea why this happens, we would
like to know.

---

Go back [home](index.html)
