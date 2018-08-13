# install scripts in ${bindir}
bin_SCRIPTS = \
	shorah
CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = \
	.version \
	meson.build \
	setup.py \
	src/meson.build \
	src/cpp/meson.build
BUILT_SOURCES = .version

shorah: $(top_srcdir)/src/shorah/__main__.py Makefile
	$(SED) 's,#!/usr/bin/env python3,#!/usr/bin/env $(PYTHON),g' $(top_srcdir)/src/shorah/__main__.py > $@

# actual C/C++ programs doing the real computation
bin_PROGRAMS = diri_sampler b2w fil

diri_sampler_SOURCES = \
	src/cpp/data_structures.hpp \
	src/cpp/dpm_sampler.hpp \
	src/cpp/dpm_sampler.cpp
diri_sampler_CPPFLAGS = $(GSL_CFLAGS)
diri_sampler_LDADD = $(GSL_LIBS)

b2w_SOURCES = \
	src/cpp/b2w.cpp
b2w_CXXFLAGS = $(PTHREAD_CFLAGS)
b2w_CPPFLAGS = $(HTSLIB_CFLAGS)
b2w_LDADD = $(HTSLIB_LIBS) $(ZLIB_LIBS) $(PTHREAD_LIBS)

fil_SOURCES = \
	src/cpp/fil.cpp
fil_CXXFLAGS = $(PTHREAD_CFLAGS)
fil_CPPFLAGS = $(GSL_CFLAGS) $(HTSLIB_CFLAGS)
fil_LDADD = $(HTSLIB_LIBS) $(ZLIB_LIBS) $(GSL_LIBS) $(PTHREAD_LIBS)


# python modules
dist_pkgpython_PYTHON = \
	src/shorah/__main__.py \
	src/shorah/amplicon.py \
	src/shorah/cli.py \
	src/shorah/shorah_snv.py \
	src/shorah/shotgun.py
pkgpython_PYTHON = .version

# example files
dist_pkgdata_DATA = \
	examples/ref_genome.fasta \
	examples/sample_454.fasta

ampliconpkgdatadir = $(pkgdatadir)/amplicon_test
dist_ampliconpkgdata_DATA = \
	examples/amplicon_test/README.md \
	examples/amplicon_test/ampli_sorted.bam \
	examples/amplicon_test/amplicon_reads.fastq \
	examples/amplicon_test/reference.fasta

# version handling
$(top_srcdir)/.version:
	echo $(VERSION) > $@-t && mv $@-t $@
dist-hook:
	echo $(VERSION) > $(distdir)/.tarball-version
