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

message('---------- Configure python scripts ----------')
root_dir = meson.source_root()
conf_data = configuration_data()
conf_data.set('PYTHON', 'python3')

py_files = [['src/python/cli.py.in', 'cli.py'],
            ['src/python/shotgun.py.in', 'shotgun.py'],
            ['src/python/__main__.py.in', '__main__.py'],
            ['src/python/shorah_snv.py.in', 'shorah_snv.py'],
            ['src/python/amplicon.py.in', 'amplicon.py']
            ]
#            ['src/scripts/amplian.py.in', 'amplian.py'],
#            ['src/scripts/dec.py.in', 'dec.py'],

foreach p : py_files
    configure_file(input: p[0],
                   output: p[1],
                   configuration: conf_data)
endforeach

meson.add_install_script('shorah-install.sh')
subdir('src/cpp')
