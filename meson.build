project(
  'shorah',
  'cpp',
  version : '1.99.1',
  default_options : [
    'buildtype=release',
    'cpp_std=c++11',
    'warning_level=3',
    'b_ndebug=if-release'],
  license : 'GPL-2+',
  meson_version : '>= 0.46.0')

# info
as_version = meson.project_version()
ver_arr = as_version.split('.')
as_major_version = ver_arr[0]
as_minor_version = ver_arr[1]
as_micro_version = ver_arr[2]

subdir('src')
