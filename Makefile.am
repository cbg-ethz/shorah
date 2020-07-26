# install scripts in ${bindir}
bin_SCRIPTS = \
	shorah
CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = \
	.version \
	LICENSE \
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
diri_sampler_CPPFLAGS = $(BOOST_CPPFLAGS)

b2w_SOURCES = \
	src/cpp/b2w.cpp
b2w_CXXFLAGS = $(PTHREAD_CFLAGS)
b2w_CPPFLAGS = $(HTSLIB_CFLAGS)
b2w_LDADD = $(HTSLIB_LIBS) $(ZLIB_LIBS) $(PTHREAD_LIBS)

fil_SOURCES = \
	src/cpp/fil.cpp
fil_CXXFLAGS = $(PTHREAD_CFLAGS)
fil_CPPFLAGS =$(BOOST_CPPFLAGS) $(HTSLIB_CFLAGS)
fil_LDADD = $(HTSLIB_LIBS) $(ZLIB_LIBS) $(PTHREAD_LIBS)


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
	examples/amplicon_test/reference.fasta \
	examples/amplicon_test/amplicon_test

shotgunpkgdatadir = $(pkgdatadir)/shotgun_test
dist_shotgunpkgdata_DATA = \
	examples/shotgun_test/README.md \
	examples/shotgun_test/test_aln.cram \
	examples/shotgun_test/test_ref.fasta \
	examples/shotgun_test/test.csv \
	examples/shotgun_test/shotgun_test

# tests
dist_check_DATA = \
	examples/run_end2end_test \
	$(dist_ampliconpkgdata_DATA) \
	$(dist_shotgunpkgdata_DATA)

AM_TESTS_ENVIRONMENT = \
	export PYTHONPATH='$(realpath $(top_srcdir)/src)'; \
	export PATH="$(realpath $(builddir)/):$(PATH)";

TESTS = examples/amplicon_test/amplicon_test \
	examples/shotgun_test/shotgun_test \
	examples/shotgun_test/test.csv

examples/shotgun_test/test.log: examples/shotgun_test/shotgun_test.log

TEST_EXTENSIONS = .csv
CSV_LOG_COMPILER = diff
AM_CSV_LOG_FLAGS = -u $(builddir)/shotgun_test/snv/SNVs_0.010000_final.csv
# NOTE cheating by considering the data files as test (normally it should be wraping scripts, see next)
LOG_COMPILER = $(top_srcdir)/examples/run_end2end_test
AM_LOG_FLAGS = $(PYTHON)

# add test files to cleanup
clean-local:
	-rm -rf amplicon_test/ shotgun_test/

# version handling
$(top_srcdir)/.version:
	echo $(VERSION) > $@-t && mv $@-t $@
dist-hook:
	echo $(VERSION) > $(distdir)/.tarball-version
