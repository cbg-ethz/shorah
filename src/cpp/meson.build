message('---------- And now the CPP part ----------')
extra_cpp_args = ['-I../src/samtools', '-DNDEBUG']
message('Is popcnt available?')
code = '''#include <nmmintrin.h>
int main() { return _mm_popcnt_u64(0); }
'''
cpp = meson.get_compiler('cpp')
success = cpp.links(code, args : '-mpopcnt', name : 'pocnt')
if success
    message('popcnt intrinsics are available')
    extra_cpp_args += ['-mpopcnt']
else
    message('popcnt is broken')
endif

gsl_dep = dependency('gsl')
thread_dep = dependency('threads')
zlib_dep = dependency('zlib')

bam_sources = files([
'../samtools/bam.c',
'../samtools/bam_aux.c',
'../samtools/bam_cat.c',
'../samtools/bam_import.c',
'../samtools/bam_index.c',
'../samtools/bam_lpileup.c',
'../samtools/bam_md.c',
'../samtools/bam_pileup.c',
'../samtools/bam_reheader.c',
'../samtools/bam_sort.c',
'../samtools/bedidx.c',
'../samtools/bgzf.c',
'../samtools/faidx.c',
'../samtools/knetfile.c',
'../samtools/kprobaln.c',
'../samtools/kstring.c',
'../samtools/razf.c',
'../samtools/sam.c',
'../samtools/sam_header.c',
'../samtools/bam.h',
'../samtools/bam_endian.h',
'../samtools/bgzf.h',
'../samtools/faidx.h',
'../samtools/kaln.h',
'../samtools/khash.h',
'../samtools/knetfile.h',
'../samtools/kprobaln.h',
'../samtools/kseq.h',
'../samtools/ksort.h',
'../samtools/kstring.h',
'../samtools/razf.h',
'../samtools/sam.h',
'../samtools/sam_header.h'
])

bam_library = static_library('libbam', bam_sources, dependencies: zlib_dep, c_args: ['-D_USE_KNETFILE'])

message('---------- Compiling DPM sampler ----------')
dpm_sources = files(['data_structures.hpp', 'dpm_sampler.hpp', 'dpm_sampler.cpp'])
dpm_exe = executable('diri_sampler', dpm_sources, dependencies: gsl_dep, cpp_args: extra_cpp_args)

message('---------- Compiling b2w ----------')
b2w_exe = executable('b2w', 'b2w.cpp', cpp_args: extra_cpp_args, link_with: bam_library)

message('---------- Compiling fil ----------')
fil_exe = executable('fil', 'fil.cpp', cpp_args: extra_cpp_args, link_with: bam_library, dependencies: gsl_dep)
