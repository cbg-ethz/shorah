#!/usr/bin/env python3

# Copyright 2007-2018
# Niko Beerenwinkel,
# Nicholas Eriksson,
# Moritz Gerstung,
# Lukas Geyrhofer,
# Kerensa McElroy,
# Osvaldo Zagordi,
# ETH Zurich

# This file is part of ShoRAH.
# ShoRAH is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ShoRAH is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with ShoRAH.  If not, see <http://www.gnu.org/licenses/>.
"""
Module that contains the command line app.

Why does this file exist, and why not put this in __main__?

  You might be tempted to import things from __main__ later, but that will cause
  problems: the code will get executed twice:

  - When you run `python -mminvar` python will execute
    ``__main__.py`` as a script. That means there won't be any
    ``minvar.__main__`` in ``sys.modules``.
  - When you import __main__ it will get executed again (as a module) because
    there's no ``minvar.__main__`` in ``sys.modules``.

  Also see (1) from http://click.pocoo.org/5/setuptools/#setuptools-integration
"""
import os
import sys

from pkg_resources import (get_distribution, DistributionNotFound)

try:
    __version__ = get_distribution('shorah').version
except DistributionNotFound:
    # probably installed using Autotools
    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), '.version'), 'r') as version_file:
        __version__ = version_file.read()

# manipulate path to import functions
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if __name__ == '__main__':
    if __package__ is None:
        os.sys.path.insert(1, parent_dir)
        mod = __import__('shorah')
        sys.modules["shorah"] = mod
        #from common import hcv_map, hiv_map, org_dict, wobbles
        #from stats import (genome_coverage, start_stop_coverage)
# else:
    #from .common import hcv_map, hiv_map, org_dict, wobbles
    #from .stats import (genome_coverage, start_stop_coverage)

# each subcommand takes one of these functions as default


def shotgun_run(args):
    """Default function for command line parser."""
    from shorah import shotgun
    shotgun.main(args)


def amplicon_run(args):
    """Default function for command line parser."""
    from shorah import amplicon
    amplicon.main(args)


def snv_run(args):
    from shorah import shorah_snv
    shorah_snv.main(args)


def main():
    """Parse command line, run default functions."""
    import argparse
    # parse command line
    # create the top-level parser
    version_parser = argparse.ArgumentParser(add_help=False)

    version_parser.add_argument(
        '-v', '--version', action='version', version=__version__)

    parent_parser = argparse.ArgumentParser(add_help=False)

    required = parent_parser.add_argument_group('required arguments')

    required.add_argument("-b", "--bam", metavar='BAM', required=True,
                          type=str, dest="b", help="sorted bam format alignment file")

    required.add_argument("-f", "--fasta", metavar='REF', required=True,
                          type=str, dest="f", help="reference genome in fasta format")

    parent_parser.add_argument("-a", "--alpha", metavar='FLOAT', required=False,
                               type=float, dest="a", default=0.1, help="alpha in dpm sampling")

    parent_parser.add_argument("-r", "--region", metavar='chrm:start-stop', required=False, type=str,
                               dest="r", default='', help="region in format 'chr:start-stop', e.g. 'chrm:1000-3000'")

    parent_parser.add_argument("-R", "--seed", metavar='INT', required=False,
                               type=int, dest="seed", default=None, help="set seed for reproducible results")

    parent_parser.add_argument("-x", "--maxcov", metavar='INT', required=False, type=int,
                               default=10000, dest="max_coverage", help="approximate max coverage allowed")

    parent_parser.add_argument("-S", "--sigma", metavar='FLOAT', default=0.01,
                               type=float, dest="sigma", help="sigma value to use when calling SNVs")

    parent_parser.add_argument("-I", "--ignore_indels", action="store_true", default=False, dest="ignore_indels",
                               help="ignore SNVs adjacent to insertions/deletions\n(legacy behaviour of 'fil', ignore this option if you don't understand)")

    parser = argparse.ArgumentParser(
        usage='%(prog)s <subcommand> [options]', epilog="Run `shorah subcommand -h` for more help", parents=[version_parser])

    subparsers = parser.add_subparsers(
        title='sub-commands', help='available sub-commands')

    # create the parser for command "shotgun"
    parser_shotgun = subparsers.add_parser(
        'shotgun', help='run local analysis in shotgun mode', parents=[parent_parser, version_parser])

    parser_shotgun.add_argument("-w", "--windowsize", metavar='INT',
                                required=False, type=int, dest="w", default=201, help="window size")

    parser_shotgun.add_argument("-s", "--winshifts", metavar='INT', required=False,
                                type=int, default=3, dest="win_shifts", help="number of window shifts")

    parser_shotgun.add_argument("-k", "--keep_files", required=False, action='store_true',
                                default=True, dest="keep_files", help="keep all intermediate files")

    parser_shotgun.set_defaults(func=shotgun_run)

    parser_amplicon = subparsers.add_parser(
        'amplicon', help='run local analysis in amplicon mode', parents=[parent_parser, version_parser])

    parser_amplicon.add_argument("-d", "--diversity", action="store_true", default=False,
                                 dest="diversity", help="detect the highest entropy region and run there")

    parser_amplicon.add_argument("-m", "--min_overlap", metavar='FLOAT', default=0.95,
                                 type=float, dest="min_overlap", help="fraction of read overlap to be included")

    parser_amplicon.set_defaults(func=amplicon_run)

    # create the parser for command "snv"
    parser_snv = subparsers.add_parser(
        'snv', help='run single-nucleotide-variant calling', parents=[parent_parser, version_parser])

    parser_snv.add_argument("-i", "--increment", metavar='INT', default=1, type=int, required=False,
                            dest="increment", help="value of increment to use when calling\nSNVs (1 used by amplian.py)")

    parser_snv.set_defaults(func=snv_run)

    # exit so that log file is not written
    if len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
        parser.print_help()
        sys.exit()

    # logging configuration
    import logging
    import logging.handlers
    logging.basicConfig(filename='shorah.log', level=logging.DEBUG,
                        format='%(levelname)s %(asctime)s %(filename)s: %(funcName)s() %(lineno)d: \t%(message)s',
                        datefmt='%Y/%m/%d %H:%M:%S')

    logging.info(' '.join(sys.argv))
    logging.info('shorah version:%s', __version__)
    # parse the args
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":  # and __package__ is None:
    main()
