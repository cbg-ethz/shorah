project(
  'shorah',
  'cpp',
  version : run_command('build-aux/git-version-gen', '@0@/.tarball-version'.format(meson.source_root()), check : true).stdout().strip(),
  default_options : [
    'buildtype=release',
    'cpp_std=c++11',
    'warning_level=3',
    'b_ndebug=if-release'],
  license : 'GPL-2+',
  meson_version : '>= 0.49.0')

meson.add_dist_script('sh', '-c', 'echo @0@>"$MESON_DIST_ROOT/.tarball-version"'.format(meson.project_version()))

# info
as_version = meson.project_version()
ver_arr = as_version.split('.')
as_major_version = ver_arr[0]
as_minor_version = ver_arr[1]
as_micro_version = ver_arr[2]

subdir('src')
