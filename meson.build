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

# python 3
py3_mod = import('python3')
py3 = py3_mod.find_python()
py3_conf_data = configuration_data()
py3_conf_data.set('PYTHON', py3.path())

py_files = ['cli.py', 'shotgun.py', '__main__.py', 'shorah_snv.py', 'amplicon.py']

foreach p : py_files
  py_script = configure_file(
    input : 'src/python/' + p + '.in',
    output : p,
    configuration : py3_conf_data)
endforeach

subdir('src')
