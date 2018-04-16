message('Is popcnt available?')
code = '''
  #include <nmmintrin.h>
  int main() { return _mm_popcnt_u64(0); }'''
cpp = meson.get_compiler('cpp')
success = cpp.links(code, args : '-mpopcnt', name : 'pocnt')
if success
    message('popcnt intrinsics are available')
    extra_cpp_args = ['-mpopcnt']
else
    message('popcnt is broken')
    extra_cpp_args = []
endif

gsl_dep = dependency('gsl')
thread_dep = dependency('threads')

dpm_exe = executable(
  'diri_sampler', files([
    'data_structures.hpp',
    'dpm_sampler.hpp',
    'dpm_sampler.cpp']),
  dependencies : [ gsl_dep, thread_dep ],
  include_directories : bam_inc,
  cpp_args : extra_cpp_args,
  install : true)

b2w_exe = executable(
  'b2w', files([
    'b2w.cpp']),
  dependencies : thread_dep,
  include_directories : bam_inc,
  cpp_args : extra_cpp_args,
  link_with : bam_library,
  install : true)

fil_exe = executable(
  'fil', files([
    'fil.cpp']),
  dependencies : [ gsl_dep, thread_dep ],
  include_directories : bam_inc,
  cpp_args : extra_cpp_args,
  link_with : bam_library,
  install : true)
