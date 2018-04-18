project('shorah',
  ['c', 'cpp'],
  version : '1.9.99',
  default_options : [
    'buildtype=release',
    'default_library=shared',
    'warning_level=3',
    'b_ndebug=if-release'],
  license : 'GPL-2+')

# info
as_version = meson.project_version()
ver_arr = as_version.split('.')
as_major_version = ver_arr[0]
as_minor_version = ver_arr[1]
as_micro_version = ver_arr[2]

subdir('src')
