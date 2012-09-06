#!/usr/bin/env python

# Copyright 2007-2012
# Niko Beerenwinkel,
# Nicholas Eriksson,
# Moritz Gerstung,
# Lukas Geyrhofer,
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

import subprocess
import os
import os.path
import sys
import dec

import logging
import logging.handlers

# Make a global logging object.
sholog = logging.getLogger("ShorahLog")

# set logging level
sholog.setLevel(logging.DEBUG)
# This handler writes everything to a file.
LOG_FILENAME = './shorah.log'
h = logging.handlers.RotatingFileHandler(LOG_FILENAME, 'w', maxBytes=100000, backupCount=5)
f = logging.Formatter("%(levelname)s %(asctime)s %(funcName)s %(lineno)d %(message)s")
h.setFormatter(f)
sholog.addHandler(h)

dn = os.path.dirname(__file__)

def run_f2r(filein):
    """
        2. translate to read format
        output: filein.read
        """

    my_prog = 'perl -I %s %s' % (dn+'/perllib', os.path.join(dn, 'fas2read.pl'))
    my_arg  = " -f %s" % filein
    assert os.path.isfile(filein), 'File %s not found' % filein
    
    # runs the fasta to read conversion
    try:
        retcode = subprocess.call(my_prog + my_arg, shell=True)
        if retcode < 0:
            sholog.error(my_prog + my_arg)
            sholog.error("Child %s was terminated by signal" % my_prog, -retcode)
        else:
            sholog.debug("Child %s returned %i" % (my_prog, retcode))
    except OSError, ee:
        sholog.error("Execution of %s failed:" % my_prog, ee)
    
    return retcode


def run_contain(filein):
    """
        3. eliminate redundant reads
        output: filein.rest
        """
    
    my_prog = os.path.join(dn, 'contain')
    my_arg  = " -f %s" % filein
    assert os.path.isfile('%s.fas' % filein), 'File %s not found' % filein
    
    #eliminates redundant reads
    try:
        retcode = subprocess.call(my_prog + my_arg, shell=True)
        if retcode < 0:
            sholog.error(my_prog + my_arg)
            sholog.error("Child %s was terminated by signal" % my_prog, -retcode)
        else:
            sholog.debug("Child %s returned %i" % (my_prog, retcode))
    except OSError, ee:
        sholog.error("Execution of %s failed:" % my_prog, ee)
    
    return retcode

def run_snv(reference,bam,sigma,incr):
    """
        x. Parse SNVs from local support, run strand bias
        output: raw_snv, snvs in 2+ windows, snvs + strand test
        """
    my_prog = os.path.join(dn, 'snv.py')
    my_arg = " -r %s -b %s -s %f -i %d" % (reference, bam, sigma, incr)
    try:
        retcode = subprocess.call(my_prog + my_arg, shell=True)
        if retcode < 0:
            sholog.error(my_prog + my_arg)
            sholog.error("Child %s was terminated by signal" % my_prog, -retcode)
        else:
            sholog.debug("Child %s returned %i" % (my_prog, retcode))
    except OSError, ee:
        sholog.error("Execution of %s failed:" % my_prog, ee)
    
    return retcode

def run_mm(filein, max_hap=200):
    """
        4. run maximum matching, output up to 200 haplotypes
        output: sample.$nr.geno
        """
    my_prog = os.path.join(dn, 'mm.py')
    my_arg  = " %s %i" % ('%s.rest' % filein, max_hap)
    
    assert os.path.isfile('%s.rest' % filein), 'File %s not found' % filein
    
    #minimal coverage of the read graph
    try:
        retcode = subprocess.call(my_prog + my_arg, shell=True)
        if retcode < 0:
            sholog.error(my_prog + my_arg)
            sholog.error("Child %s was terminated by signal" % my_prog, -retcode)
        else:
            sholog.debug("Child %s returned %i" % (my_prog, retcode))
    except OSError, ee:
        sholog.error("Execution of %s failed:" % my_prog, ee)
    
    return retcode


def run_freqEst(filein):
    """
        5. run EM
        output: sample.$nr.popl
        """
    my_prog = os.path.join(dn, 'freqEst')
    my_arg  = " -f %s" % filein
    assert os.path.isfile('%s.rest' % filein), 'File %s not found' % filein
    
    #estimate frequencies
    try:
        retcode = subprocess.call(my_prog + my_arg, shell=True)
        if retcode < 0:
            sholog.error(my_prog + my_arg)
            sholog.error("Child %s was terminated by signal" % my_prog, -retcode)
        else:
            sholog.debug("Child %s returned %i" % (my_prog, retcode))
    except OSError, ee:
        sholog.error("Execution of %s failed:" % my_prog, ee)
    
    return retcode


def main():
    """ Only called if run not interactively
        """
    import optparse
    import shutil
    import glob
    import dec
    
    # parse command line
    optparser = optparse.OptionParser()
    
    optparser.add_option("-b", "--bam", help="sorted bam format alignment file.",
                         default="", type="string", dest="b")
    optparser.add_option("-f", "--fasta", help="reference genome in fasta format.",
                         default="", type="string", dest="f")
    optparser.add_option("-a", "--alpha", help="alpha in dpm sampling <0.01>", default=0.01,
                         type="float", dest="a")
    optparser.add_option("-w", "--windowsize", help="window size <201>", default=201,
                         type="int", dest="w")
    optparser.add_option("-s", "--winshifts", help="number of window shifts <3>", default=3,
                         type="int", dest="s") # window shiftings, such that each base is covered up to win_shifts times
    optparser.add_option("-i", "--sigma", help="value of sigma to use when calling SNVs", default = 0.01,
                         type="float", dest="i")
    optparser.add_option("-x", "--maxcov", help="approximate maximum coverage allowed", default=10000,
                         type="int", dest="x")
    optparser.add_option("-r", "--region", help="region in format 'chr:start-stop', eg 'ch3:1000-3000'", default='',
                         type="string", dest="r")
    optparser.add_option("-k","--keep_files",help="keep intermediate files (Gibbs sampling)", default=False,
                         action="store_true", dest="k")
    (options, args) = optparser.parse_args()
    
    sholog.info(' '.join(sys.argv))
    in_bam = options.b
    in_fasta = options.f
    keep_all_files = options.k
    step = options.w
    win_shifts = options.s
    max_coverage = options.x
    sigma = options.i
    alpha = options.a
    region = options.r

    in_stem=os.path.split(in_bam)[1].split('.')[0]

    if not os.path.exists('snv/SNV.txt'):
        # run dec.py
        sholog.debug('running dec.py')
        dec.main(in_bam,in_fasta, step, win_shifts, max_coverage, sigma, region, keep_all_files, alpha)

        # copy file and run fas2reads
        shutil.move('%s.cor.fas' % in_stem, '%s_cor.fas' % in_stem)
        sholog.debug('running fas2reads')

        retcode = run_f2r('%s_cor.fas' % in_stem)
        if retcode is not 0:
            sholog.error('fas2reads did not return 0')
            sys.exit()

        # run contain
        retcode = run_contain('%s_cor' % in_stem)
        if retcode is not 0:
            sys.exit()

        # run mm.py (this might become a module)
        retcode = run_mm('%s_cor' % in_stem)
        if retcode is not 0:
            sys.exit()
        
        # run EM freqEst
        retcode = run_freqEst('%s_cor' % in_stem)
        if retcode is not 0:
            sys.exit()



    # run snv.py
    sholog.debug('running snv.py')
    retcode = run_snv(in_fasta, in_bam, sigma, step/win_shifts)
    if retcode is not 0:
        sholog.error('snv calling failed')

    # tidy snvs
    try:
        os.mkdir('snv')
    except:
        pass
    for snv_file in glob.glob('./SNV*'):
        shutil.move(snv_file, 'snv/')
# goodbye

if __name__ == "__main__":
    main()
