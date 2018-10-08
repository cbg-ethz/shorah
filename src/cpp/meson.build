message('Is popcnt available?')
code = '''
  #include <nmmintrin.h>
  int main() { return _mm_popcnt_u64(0); }'''
cpp = meson.get_compiler('cpp')
success = cpp.links(code, args : '-mpopcnt', name : 'popcnt')
if success
    message('popcnt intrinsics are available')
    extra_cpp_args = ['-mpopcnt']
else
    message('popcnt is broken')
    extra_cpp_args = []
endif

boost_dep = dependency('boost', version : '>=1.56' ) # we only use headers : boost::random is templates-only
gsl_dep = dependency('gsl')
thread_dep = dependency('threads')
htslib_dep = dependency('htslib', version : '>=1.7')

dpm_exe = executable(
  'diri_sampler', files([
    'data_structures.hpp',
    'dpm_sampler.hpp',
    'dpm_sampler.cpp']),
  dependencies : [
    thread_dep,
    boost_dep],
  cpp_args : extra_cpp_args,
  install : true)

b2w_exe = executable(
  'b2w', files([
    'b2w.cpp']),
  dependencies : [
    thread_dep,
    htslib_dep],
  cpp_args : extra_cpp_args,
  install : true)

fil_exe = executable(
  'fil', files([
    'fil.cpp']),
  dependencies : [
    thread_dep,
    gsl_dep,
    htslib_dep],
  cpp_args : extra_cpp_args,
  install : true)
