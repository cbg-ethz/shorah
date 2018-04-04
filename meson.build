project('shorah',
    ['c', 'cpp'],
    version : '1.1.1',
    default_options : [
        'buildtype=release',
        'default_library=shared',
        'warning_level=3'],
    license : 'GPL')

# info
as_version = meson.project_version()
ver_arr = as_version.split('.')
as_major_version = ver_arr[0]
as_minor_version = ver_arr[1]
as_micro_version = ver_arr[2]

message('configure python scripts')

conf_data = configuration_data()
conf_data.set('PYTHON', 'python3')
py_files = [['src/scripts/amplian.py.in', 'amplian.py'],
            ['src/scripts/dec.py.in', 'dec.py'],
            ['src/scripts/shorah.py.in', 'shorah.py'],
            ['src/scripts/snv.py.in', 'snv.py']
            ]
foreach p : py_files
    configure_file(input : p[0],
                   output : p[1],
                   configuration : conf_data)
endforeach


message('---------- And now the CPP part ----------')
extra_cpp_args = ['-I../src/samtools']
message('Is popcnt available?')
code = '''#include <nmmintrin.h>
int main() { return _mm_popcnt_u64(0); }
'''
cpp = meson.get_compiler('cpp')
success = cpp.links(code, args : '-mpopcnt', name : 'link check')
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
'src/samtools/bam.c',
'src/samtools/bam_aux.c',
'src/samtools/bam_cat.c',
'src/samtools/bam_import.c',
'src/samtools/bam_index.c',
'src/samtools/bam_lpileup.c',
'src/samtools/bam_md.c',
'src/samtools/bam_pileup.c',
'src/samtools/bam_reheader.c',
'src/samtools/bam_sort.c',
'src/samtools/bedidx.c',
'src/samtools/bgzf.c',
'src/samtools/faidx.c',
'src/samtools/knetfile.c',
'src/samtools/kprobaln.c',
'src/samtools/kstring.c',
'src/samtools/razf.c',
'src/samtools/sam.c',
'src/samtools/sam_header.c',
'src/samtools/bam.h',
'src/samtools/bam_endian.h',
'src/samtools/bgzf.h',
'src/samtools/faidx.h',
'src/samtools/kaln.h',
'src/samtools/khash.h',
'src/samtools/knetfile.h',
'src/samtools/kprobaln.h',
'src/samtools/kseq.h',
'src/samtools/ksort.h',
'src/samtools/kstring.h',
'src/samtools/razf.h',
'src/samtools/sam.h',
'src/samtools/sam_header.h'
])

bam_library = static_library('libbam.a', bam_sources, dependencies: zlib_dep, c_args: ['-D_USE_KNETFILE'])

message('---------- Compiling DPM sampler ----------')
dpm_sources = files(['src/shorah/data_structures.hpp', 'src/shorah/dpm_sampler.hpp', 'src/shorah/dpm_sampler.cpp'])
dpm_exe = executable('diri_sampler', dpm_sources, dependencies: gsl_dep, cpp_args: extra_cpp_args)

message('---------- Compiling b2w ----------')
b2w_exe = executable('b2w', 'src/shorah/b2w.cpp', cpp_args: extra_cpp_args, link_with: bam_library)
