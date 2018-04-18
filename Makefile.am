# install scripts in ${bindir}
bin_SCRIPTS = \
	shorah
CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = \
	meson.build \
	setup.py \
	src/meson.build \
	src/cpp/meson.build \
	src/samtools/COPYING \
	src/samtools/meson.build

shorah: $(top_srcdir)/src/python/__main__.py Makefile
	$(SED) 's,#!/usr/bin/env python3,#!/usr/bin/env $(PYTHON),g' $(top_srcdir)/src/python/__main__.py > $@

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
b2w_CPPFLAGS = -I$(top_srcdir)/src/samtools
b2w_LDADD = libbam.a $(ZLIB_LIBS) $(PTHREAD_LIBS)

fil_SOURCES = \
	src/cpp/fil.cpp
fil_CXXFLAGS = $(PTHREAD_CFLAGS)
fil_CPPFLAGS = $(GSL_CFLAGS) -I$(top_srcdir)/src/samtools
fil_LDADD = libbam.a $(ZLIB_LIBS) $(GSL_LIBS) $(PTHREAD_LIBS)

noinst_LIBRARIES = libbam.a
libbam_a_SOURCES = \
	src/samtools/bam.c \
	src/samtools/bam.h \
	src/samtools/bam_aux.c \
	src/samtools/bam_cat.c \
	src/samtools/bam_endian.h \
	src/samtools/bam_import.c \
	src/samtools/bam_index.c \
	src/samtools/bam_lpileup.c \
	src/samtools/bam_md.c \
	src/samtools/bam_pileup.c \
	src/samtools/bam_reheader.c \
	src/samtools/bam_sort.c \
	src/samtools/bedidx.c \
	src/samtools/bgzf.c \
	src/samtools/bgzf.h \
	src/samtools/faidx.c \
	src/samtools/faidx.h \
	src/samtools/kaln.h \
	src/samtools/khash.h \
	src/samtools/knetfile.c \
	src/samtools/knetfile.h \
	src/samtools/kprobaln.c \
	src/samtools/kprobaln.h \
	src/samtools/kseq.h \
	src/samtools/ksort.h \
	src/samtools/kstring.c \
	src/samtools/kstring.h \
	src/samtools/razf.c \
	src/samtools/razf.h \
	src/samtools/sam.c \
	src/samtools/sam.h \
	src/samtools/sam_header.c \
	src/samtools/sam_header.h
libbam_a_CPPFLAGS = $(ZLIB_CFLAGS) -D_USE_KNETFILE
libbam_a_CFLAGS = $(PTHREAD_CFLAGS)

# python modules
dist_pkgpython_PYTHON = \
	src/python/__main__.py \
	src/python/amplicon.py \
	src/python/cli.py \
	src/python/shorah_snv.py \
	src/python/shotgun.py

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
