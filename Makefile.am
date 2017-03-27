# install scripts in ${bindir}
bin_SCRIPTS = \
	src/scripts/amplian.py \
	src/scripts/dec.py \
	src/scripts/mm.py \
	src/scripts/shorah.py \
	src/scripts/snv.py \
	fas2read.pl
CLEANFILES = fas2read.pl
EXTRA_DIST = \
	src/scripts/amplian.py.in \
	src/scripts/dec.py.in \
	src/scripts/mm.py.in \
	src/scripts/shorah.py.in \
	src/scripts/snv.py.in \
	src/scripts/fas2read.pl.in \
	src/samtools/COPYING

fas2read.pl: src/scripts/fas2read.pl.in Makefile
	$(SED) 's,[@]perldir[@],$(perldir),g' $< > $@

# actual C/C++ programs doing the real computation
bin_PROGRAMS = diri_sampler contain freqEst b2w fil

diri_sampler_SOURCES = \
	src/shorah/data_structures.hpp \
	src/shorah/dpm_sampler.hpp \
	src/shorah/dpm_sampler.cpp
diri_sampler_CPPFLAGS = $(GSL_CFLAGS)
diri_sampler_LDADD = $(GSL_LIBS)

contain_SOURCES = \
	src/shorah/contain.cpp

freqEst_SOURCES = \
	src/shorah/freqEst.cpp

b2w_SOURCES = \
	src/shorah/b2w.cpp
b2w_CXXFLAGS = $(PTHREAD_CFLAGS)
b2w_CPPFLAGS = -I$(top_srcdir)/src/samtools
b2w_LDADD = libbam.a $(ZLIB_LIBS) $(PTHREAD_LIBS)

fil_SOURCES = \
	src/shorah/fil.cpp
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

# perl modules
dist_perl_DATA = \
	src/perl/DNA.pm \
	src/perl/NikoUtil.pm \
	src/perl/ReadUtil.pm \
	src/perl/Util.pm

# python modules
dist_pkgpython_DATA = \
	src/python/__init__.py \
	src/python/matching.py
