bam_sources = files([
  'bam.c',
  'bam_aux.c',
  'bam_cat.c',
  'bam_import.c',
  'bam_index.c',
  'bam_lpileup.c',
  'bam_md.c',
  'bam_pileup.c',
  'bam_reheader.c',
  'bam_sort.c',
  'bedidx.c',
  'bgzf.c',
  'faidx.c',
  'knetfile.c',
  'kprobaln.c',
  'kstring.c',
  'razf.c',
  'sam.c',
  'sam_header.c',
  'bam.h',
  'bam_endian.h',
  'bgzf.h',
  'faidx.h',
  'kaln.h',
  'khash.h',
  'knetfile.h',
  'kprobaln.h',
  'kseq.h',
  'ksort.h',
  'kstring.h',
  'razf.h',
  'sam.h',
  'sam_header.h'])

zlib_dep = dependency('zlib')

bam_library = static_library(
  'libbam',
  bam_sources,
  dependencies : zlib_dep,
  c_args : ['-D_USE_KNETFILE'],
  install : false)

bam_inc = include_directories('.')
