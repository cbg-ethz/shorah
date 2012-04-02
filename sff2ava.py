#!/usr/bin/env python

# Copyright 2007, 2008, 2009, 2012
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


import sys
import logging, logging.handlers

MID1='ACGAGTGCGT'
MID2='ACGCTCGACA'
MID3='AGACGCACTC'
MID4='AGCACTGTAG'
MID5='ATCAGACACG'
MID6='ATATCGCGAG'
MID7='CGTGTCTCTA'
MID8='CTCGCGTGTC'
MID9='TAGTATCAGC'
MID10='TCTCTATGCG'
MID11='TGATACGTCT'
MID12='TACTGAGCTA'
MID13='CATAGTAGTG'
MID14='CGAGAGATAC'

key = 'TCAG'
key_reverse='CTGA'
mid = 'ACGCTCGACA'
mid_reverse = 'TGTCGAGCGT'
f_primer = 'TTNTGGGAAGTTCAATTAGGAATACC'
r_primer = 'CATTCCTTTGGATGGGTTATGAAC'


LOG_FILENAME = './sff2ava.log'
# Make a global logging object.
x = logging.getLogger('logfun')
x.setLevel(logging.DEBUG)
# This handler writes everything to a file.
# h = logging.FileHandler('./log', 'w')
h = logging.handlers.RotatingFileHandler(LOG_FILENAME, 'w', maxBytes=100000, backupCount=5)
f = logging.Formatter('%(levelname)s %(asctime)s %(funcName)s line=%(lineno)d %(message)s')
h.setFormatter(f)
x.addHandler(h)
logfun = logging.getLogger('logfun')

def parse_com_line():
    
    options, args = None, None
    
    try:
        import argparse

        parser = argparse.ArgumentParser(description='Writes a sff file with the corrected reads',
                                        epilog='Input are mandatory')
        parser.add_argument('-c', '--corrected', dest='corrected', # default action is 'store',
                            help='file with the corrected reads (output of diri_sampler)')
        parser.add_argument('-s', '--sff', dest='sff', # default action is 'store',
                        help='original sff file')
        parser.add_argument('-o', '--output', dest='output', # default action is 'store',
                        help='output sff file')
        args = parser.parse_args()
    
    except ImportError:
        import optparse
        
        optparser = optparse.OptionParser()
        optparser.add_option("-c","--corrected", type="string", default="",
                             dest="corrected", help="file with the corrected reads (output of diri_sampler)")
        optparser.add_option("-s","--sff", type="string", default="",
                             dest="sff", help="original sff file")
        optparser.add_option("-o","--output", help="output sff file",
                             type="string", dest="output")
        
        (options, args) = optparser.parse_args()
    
    return options, args


if __name__ == '__main__':
    
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    
    logfun.info(' '.join(sys.argv))
    options, args = parse_com_line()
    if args:
        args = vars(args)
    else:
        args = vars(options)

    sffin = SeqIO.parse(args['sff'], 'sff')
    corr = SeqIO.parse(args['corrected'], 'fasta')
    output = open(args['output'], 'wb')
    
    wanted = {}
    for s in corr:
        if float(s.description.split('=')[1])> 0.9:
            wanted[s.id.split('#',1)[0]] = s.seq.tostring()
    records = (r for r in sffin if r.id.split('#',1)[0] in wanted)
    print "Found %i unique reads in %s" % (len(wanted), args['corrected'])
    cor_len = max([len(s) for s in wanted.values()])
    ava_recs = []
    for s in records:
        cql =  s.annotations['clip_qual_left']
        cqr = s.annotations['clip_qual_right']
        anns = s.annotations
        # seq_corr = s.seq[:cql] + wanted[s.id] + s.seq[cqr:]
        #seq_corr = Seq(s.seq[:cql].tostring() + mid + f_primer + wanted[s.id].replace('-', '') + r_primer + mid_reverse + key_reverse)
        seq_corr = Seq(s.seq[:cql].tostring() + f_primer + wanted[s.id].replace('-', '') + r_primer + key_reverse)
        seq_len = len(seq_corr)
        # print seq_len
        anns['flow_index'] = anns['flow_index'][:seq_len]
        anns['clip_qual_right'] = cql + cor_len
        # print anns
        s = SeqRecord(seq_corr, id=s.id, annotations=anns)
        s.letter_annotations = {"phred_quality": [40]*seq_len}
        ava_recs.append(s)

    SeqIO.write(ava_recs, output, 'sff')
'''
count = SeqIO.write(records, output_file, "sff")
print "Saved %i records from %s to %s" % (count, input_file, output_file)
if count < len(wanted):
    print "Warning %i IDs not found in %s" % (len(wanted)-count, input_file)

'''