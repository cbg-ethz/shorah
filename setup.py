#!/usr/bin/env python
import sys
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand

class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)


setup(
    use_scm_version=True,
    setup_requires=['setuptools_scm', 'setuptools_scm_git_archive'],
    install_requires=['setuptools_scm'],
    tests_require=['pytest', 'flake8', 'pep257'],
    cmdclass={'test': PyTest},
    name='ShoRAH',
    description='SHOrt Reads Assembly into Haplotypes',
    url='http://github.com/cbg-ethz/shorah',
    # author='Osvaldo Zagordi',
    # author_email='firstname.lastname@gmail.com',
    packages=find_packages('src'),  # include all packages under src
    package_dir={'': 'src'},  # tell setuptools packages are under src
    # package_data={'virmet': ['db/*', 'db/template.tex']},
    # scripts=['bin/virmet'],
    # data_files=[('Rscripts', ['Rscripts/covplot.R'])],
    entry_points={
        'console_scripts': ['shorah = shorah.cli:main']
    },
    license='GPL 3.0',
    long_description='''
    ShoRAH is designed to analyse genetically heterogeneous samples. It provides error correction,
    local haplotype reconstruction, and estimation of the frequency of the different genetic variants
    present in a mixed sample.
    ''',
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 5 - Production/Stable',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',

        # Pick your license as you wish (should match "license" above)
        'License :: GPL 3.0',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        # 'Programming Language :: Python :: 2',
        # 'Programming Language :: Python :: 2.6',
        # 'Programming Language :: Python :: 2.7',
        # 'Programming Language :: Python :: 3',
        # 'Programming Language :: Python :: 3.2',
        # 'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6']
)
