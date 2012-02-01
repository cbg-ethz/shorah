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

#import os
import sys
#import pythonlib
import logging, logging.handlers

LOG_FILENAME = './s2f.log'
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

ins_check = True
acclength = 2.
amb_thresh = 2
min_length = 160
dna_code = ['A','C','T','G']
alphabet = ['A', 'C', 'G', 'T', 'N']
gen_ref_start = ''

f_code = {}
f_code['R'] = ['G', 'A']
f_code['Y'] = ['T', 'C']
f_code['K'] = ['G', 'T']
f_code['M'] = ['A', 'C']
f_code['S'] = ['G', 'C']
f_code['W'] = ['A', 'T']
f_code['B'] = ['G', 'T', 'C']
f_code['D'] = ['G', 'A', 'T']
f_code['H'] = ['A', 'C', 'T']
f_code['V'] = ['G', 'C', 'A']
f_code['N'] = ['A', 'G', 'C', 'T']


def parse_com_line():
    import optparse
    
    optparser = optparse.OptionParser()
    
    optparser.add_option("-f","--readfile", type="string", default="",
                         dest="input")
    optparser.add_option("-r","--ref", type="string", default="",
                         dest="ref")
    optparser.add_option("-o","--output", help="output suffix must be '.far' or none for stdout",
                         type="string", dest="o")
    optparser.add_option("-t","--threshold", help="if similarity is less, throw reads away... <default=0.7>",
                         type="float", dest="threshold", default=0.7)
    optparser.add_option("-a","--amplicon", help="include only reads whose overlap with the reference is larger than the threshold\
                         <default=0.0, i.e. deactivated>",
                         type="float", dest="amplicon", default=0.0)
    optparser.add_option("-d","--delete_files",help="delete temporary files <default=False>",
                         action="store_true", default=False, dest="d")
    optparser.add_option("-p","--pad_insert",help="insert padding gaps <default=not insert>",
                         action="store_true", default=False, dest="pad")
                         
    (options, args) = optparser.parse_args()
    
    return options, args


def check_indels(read, refseq):
    ''' Check indels on homopolymeric regions
    Randomize the gap position in case of homopolymer deletion, i.e.
    AAAC
    -AAC
    becomes (with probability 1/3)
    AAAC
    A-AC
    !!!!!!!!!!!!!!!!!    
    !!! WATCH OUT !!!
    !!!!!!!!!!!!!!!!!
    Not suffling at the moment
    '''
    
    import re
    #from random import shuffle
    
    max_homo = 10
    #ins_thresh = 0.0
    
    read_l = list(read)
    ref_l = list(refseq)
    ins = {}
    dels = {}
    start = len(read) - len(read.lstrip('-'))
        
    for b in dna_code:
        for i in range(max_homo+1, 2, -1):
            for m in re.finditer(b*i, read):
                a, e = m.start(), m.end()
                if refseq[a] == '-' and a not in ins and a > start and refseq[a+1:a+i] == ''.join(b*(i-1)): # read has an insertion
                    ins[a] = i
                if read[a-1] == '-' and a not in dels and a > start and read[a] == refseq[a-1]: # read has a deletion
                    dels[a] = i
                    
    
    for k, v in ins.iteritems(): # shuffle the reference, to gap in the reference
        tmp = ref_l[k:k+v]
        # shuffle(tmp)
        for p in enumerate(range(k, k+v)):
            ref_l[p[1]] = tmp[p[0]]
            
    for k, v in dels.iteritems():
        tmp = read_l[k-1:k+v]
        # shuffle(tmp)
        for p in enumerate(range(k-1, k+v)):
            read_l[p[1]] = tmp[p[0]]

    assert len(read_l) == len(ref_l), '%s %s' % (len(read_l), len(ref_l))
    
    newread, newref = ''.join(read_l), ''.join(ref_l)
    return newread, newref


def check_reference(fas_reference):
    ''' check the reference genome
    '''
    from Bio import SeqIO
    
    rh = open(fas_reference, 'rU')
    ref = list(SeqIO.parse(rh, 'fasta') )
    assert len(ref) == 1, 'One and only one sequence must be in the reference file'
    gen_length = len(ref[0].seq)
    if gen_length > 20:
        gen_ref_start = ref[0].seq.tostring()[:20]
    else:
        gen_ref_start = ref[0].seq.tostring()

    assert 'N' not in ref[0].seq, "Found an ambiguous position in the reference '%s'" % ref[0].id
    logfun.info('The reference genome length is %d' % gen_length )
    
    return gen_ref_start


def prepare_reads(f_fasta_filename):
    '''
    '''
    from Bio import SeqIO
    import math
    f_fasta = open(f_fasta_filename)
    seqlist_raw = list(SeqIO.parse(f_fasta, 'fasta'))
    
    seqlist = [ s for s in seqlist_raw if len(s) >= min_length ]
    countreads = len(seqlist)
    logfun.info('removed %d reads because shorter than %d' % (len(seqlist_raw) - countreads, min_length))
    
    # forward...
    f_fasta_forward_filename = 'tmp_reads_f.fas'
    f_fasta_forward = open(f_fasta_forward_filename, 'w')
    SeqIO.write(seqlist, f_fasta_forward, 'fasta')
    f_fasta_forward.close()
    
    # ...and reverse
    for seq in seqlist:
        seq.seq = seq.seq.reverse_complement()
    f_fasta.close()
    f_fasta_reverse_filename = 'tmp_reads_r.fas'
    f_fasta_reverse = open(f_fasta_reverse_filename,'w')
    SeqIO.write(seqlist,f_fasta_reverse,'fasta')
    f_fasta_reverse.close()
    
    logfun.info('Remaining %d reads', countreads)

    length = 0.
    length2 = 0.
    n = 0.

    for read in seqlist:
        if read.seq.count('N') < amb_thresh:
            len_seq = len(read.seq)
            length += float(len_seq)
            length2 += float(len_seq*len_seq)
            # readdict[read.] = [seq,len_seq]
            n += 1.

    meanlr = length/n
    stdlr = math.sqrt((n*length2 - length*length)/(n*n - n))
    #allowed_length = [meanlr - acclength * stdlr, meanlr + (1+acclength) * stdlr]
    logfun.info('mean length=%d, stderr=%d' % (meanlr, stdlr))


def align_strand(al_info):
    '''
    '''
    from Bio.Emboss.Applications import NeedleCommandline
    import subprocess
    
    ref_file = al_info['ref_file']
    in_file = al_info['in_file']
    out_file = al_info['out_file']
    
    cline = NeedleCommandline(gapopen=6.0, gapextend=3.0)
    cline.asequence = ref_file
    cline.bsequence = in_file
    cline.outfile = out_file
    cline.aformat = 'markx10'
    cml = str(cline) + ' -adesshow3 -auto'
    logfun.info(cml)
    try:
        retcode = subprocess.call(cml, shell=True)
        if retcode < 0:
            logfun.info('Child diri_sampler was terminated by signal %d' -retcode)
        else:
            logfun.info('Child diri_sampler returned %d' % retcode)
    except OSError, ee:
        logfun.exception('Execution of diri_sampler failed:' + ee)
    
    return


def parse_alignments(threshold):
    '''
    '''
    
    from pythonlib.MarkxIO import Markx10Iterator
    counted_ident = []
    
    countreads_afterreverse = 0
    count_forward = 0
    count_reverse = 0
    
    reads = {}
    ref = {}
    descriptions = {}
    global_ins = {}
    insertions = {}
    ID = 1
    disc_seq = []
    
    f_forward = open('tmp_align_f.needle')
    f_reverse = open('tmp_align_r.needle')
    
    forwardaligniter = Markx10Iterator(f_forward)
    reversealigniter = Markx10Iterator(f_reverse)
    ref_start = 100000
    ref_end = 0
    pos = 0

    logfun.info('parsing the alignments')
    
    # iterates through the alignments
    while True:
        
        try:
            f_align = forwardaligniter.next()
            r_align = reversealigniter.next()
        except:
            break
        if f_align is None or r_align is None:
            break
        
        #assert f_align.get_all_seqs()[1].id == r_align.get_all_seqs()[1].id, 'same seq back and forward'
        #descr = f_align.get_all_seqs()[1].id
        descr = list(f_align)[0].id
        assert descr == list(r_align)[0].id, 'same seq back and forth'

        basecount = 0.
        idencount = 0.
        aligned_gaps = 0.
        
        accept_read = False
        
        if float(f_align._annotations['sw_score']) > float(r_align._annotations['sw_score']):
            #tmp = f_align.get_seq_by_num(1).tostring().upper()
            #refseq = f_align.get_seq_by_num(0).tostring().upper()
            pair_here = list(f_align)
            tmp = pair_here[1].seq.tostring().upper()
            refseq = pair_here[0].seq.tostring().upper()
            count_forward += 1
        else:
            #tmp = r_align.get_seq_by_num(1).tostring().upper()
            #refseq = r_align.get_seq_by_num(0).tostring().upper()
            pair_here = list(r_align)
            tmp = pair_here[1].seq.tostring().upper()
            refseq = pair_here[0].seq.tostring().upper()
            count_reverse += 1
    
        assert len(tmp) == len(refseq)
        
        # if ins_check:
        tmp, refseq = check_indels(tmp, refseq)
        
        # get start and stop values for the alignment,
        # only where ref and read overlap
        start = False
        for c in enumerate(zip(refseq, tmp)):
            if not start and c[1].count('-') == 0:
                align_start = c[0]
                start = True
            if start == True and c[1].count('-') == 0:
                align_end = c[0]
        
        if align_start < ref_start:
            ref_start = align_start
        if align_end > ref_end:
            ref_end = align_end
        
        tmpl = list(tmp)
    
        # take only the overlap with the reference
        # if the read begins before the ref_genome...
        for i in range(align_start):
            tmpl[i] = '-'
    
        # ...or if it continues also after it
        tmpl = tmpl[:align_end+1]
        refseq = refseq[:align_end+1]
        
        tmp = ''.join(tmpl)
        
        if len(tmp) != len(refseq):
            logfun.error('Read and match lengths are %d %d' % (len(tmp), len(refseq)))
            sys.exit()
        # alignment must always start from the beginning of the reference
        if not refseq.replace('-','').startswith(gen_ref_start):
            logfun.error('should consider from the beginning of the reference only')
            sys.exit()
        # last character of the matched reference cannot be a gap
        if refseq.endswith('-'):
            logfun.error('last character of the matched reference cannot be a gap\nrefseq:\n%s' % refseq)
            sys.exit()
            
        # detect positions of insertions
        pos_ins = {}
        i = 0
        count = 0
        while i <= align_end:
            if refseq[i] == '-':
                # there's an insertion, how long?
                ext = len(refseq[i:]) - len(refseq[i:].lstrip('-'))
                pos_ins[count] = ext
                i += ext-1
            elif refseq[i] != '-':
                count += 1
            i += 1
        
        for i in range(align_start, align_end+1):
            if tmp[i] == refseq[i]:
                idencount += 1.
            if tmp[i] in dna_code:
                basecount += 1.
        try:
            counted_ident.append(idencount/basecount)
        except:
            logfun.error('Read %s is empty' %  list(f_align)[1].id)
            sys.exit()
            
        read_id = descr + '_' + str(ID)
        len_no_flank = len(tmp.strip('-'))
        if threshold < idencount/basecount:
            reads[read_id] = tmp
            ref[read_id] = refseq
            insertions[read_id] = pos_ins
            descriptions[read_id] = descr
            
            for pi in pos_ins:
                try:
                    if global_ins[pi] < pos_ins[pi]:
                        global_ins[pi] = pos_ins[pi]
                except:
                    global_ins[pi] = pos_ins[pi]
                    
            ID += 1
            countreads_afterreverse += 1
        else:
            pass
        
    logfun.info('%d reads were above the threshold (included)' % countreads_afterreverse)
    logfun.info('forward: %d, reverse: %d' % (count_forward, count_reverse))
    
    if len(reads) != len(ref):
        logfun.exception('%d reads, %d references' % (len(reads), len(ref)))
        sys.exit()
    
    return ref, reads, insertions, descriptions, global_ins


def pad_alignment(ref, reads, insertions, descriptions, global_ins, amp_thresh, f_output):
    '''
    '''
    # ##############################
    # propagating the gaps         #
    # from the pairwise alignments #
    # ##############################
    import operator
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    
    ids = ref.keys()

    new_reads = {}
    begin = {}
    end = []
    pad = '-'
    # fasta_length = 80

    for ID in ids:
        to_insert = {}
    
        for ti in global_ins:
            try:
                already = insertions[ID][ti]
            except KeyError:
                already = 0
            to_insert[ti] = global_ins[ti] - already
            
        i = 0
        this_read = []
        start = False
        padded = 0
        
        for c in enumerate(zip(reads[ID], ref[ID])):
            if i in to_insert and c[1][1] != '-':
                n = to_insert[i]
                this_read.extend(pad*n)
                padded += n
                
            this_read.append(c[1][0])
            
            if c[1][1] != '-':
                i += 1
                
            if start == False and c[1].count('-') == 0:
                start = True
                begin[ID] = c[0] + len(this_read) + 1
            
        new_reads[ID] = ''.join(this_read)
        end.append(max(new_reads[ID].rfind(b) for b in alphabet))
        
    logfun.info('sorting')
    items = sorted(begin.iteritems(), key=operator.itemgetter(1), reverse=False)
    logfun.info('writing reads to file')
    min_begin = min(begin.values())
    max_end = max(end)
    records = []
    too_many_flank = 0
    for i in items:
        ID = i[0]
        if len(new_reads[ID]) < max_end:
            new_reads[ID] += '-'*(max_end - len(new_reads[ID]))        
        len_no_flank = float(len(new_reads[ID].strip('-')))
        
        if len_no_flank/len(new_reads[ID]) > amp_thresh:
            seq_here = Seq(new_reads[ID])
            records.append(SeqRecord(seq_here, id=ID, description=descriptions[ID]))
        else:
            too_many_flank += 1
    logfun.info('%d reads with too many flanking gaps (only useful in amplicon analysis)' % too_many_flank)
    written = SeqIO.write(records, f_output, 'fasta')
    logfun.info('%d reads written to output' % written)

def merge_alignment(ref, reads, descriptions, amp_thresh, f_output):
    '''
    '''
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    ids = ref.keys()

    new_reads = {}
    begin = {}
    end = []
    pad  = '-'
    fasta_length = 80
    for ID in ids:
        this_read = []
        start = False
        i = 0
        for c in zip(reads[ID], ref[ID]):
            if c[1] != '-':
                this_read.append(c[0])
                i += 1
            
            if start == False and c[0] != '-':
                start = True
                begin[ID] = i
        new_reads[ID] = ''.join(this_read)
        end.append(len(this_read))
    logfun.info('sorting')
    items = [ (v, k) for k, v in begin.items() ]
    
    items.sort()
    
    logfun.info('writing reads to file')
    max_end = max(end)
    records = []
    too_many_flank = 0
    for i in items:
        ID = i[1]
        if len(new_reads[ID]) < max_end:
            new_reads[ID] += '-'*(max_end - len(new_reads[ID]))
        
        len_no_flank = float(len(new_reads[ID].strip('-')))
        
        if len_no_flank/len(new_reads[ID]) > amp_thresh:
            seq_here = Seq(new_reads[ID])
            records.append(SeqRecord(seq_here, id=ID, description=descriptions[ID]))
        else:
            too_many_flank += 1
    logfun.info('%d reads with too many flanking gaps (only useful in amplicon analysis)' % too_many_flank)
    written = SeqIO.write(records, f_output, 'fasta')
    logfun.info('%d reads written to output' % written)
    
    return


def main(ref_file, in_file, out_file, thresh, pad_insert, keep_files):
    '''
    '''
    
    from multiprocessing import Pool
    import os.path
    
    logfun.info('\n\n')
    logfun.info('running with:')
    logfun.info(' '.join(sys.argv) + '\n')
    
    if in_file.split('.')[-1] in ['fas', 'fasta', 'fa', 'fna']:
        if os.path.isfile(in_file):
            f_fasta_filename = in_file
        else:
            logfun.error('fasta file "%s" not found...' % in_file)
            sys.exit()
    else:
        logfun.error('format of input-file not supported...')
        sys.exit()
        
    if out_file == None:
        f_output = sys.stdout
    else:
        if(out_file.split('.')[-1] != 'far'):
            logfun.error('The suffix of output file must be .far')
            sys.exit()
        f_output = open(out_file, 'w')

    gen_ref_start = check_reference(ref_file)
    
    prepare_reads(f_fasta_filename)
    
    af = {
        'in_file': 'tmp_reads_f.fas',
        'ref_file': ref_file,
        'out_file': 'tmp_align_f.needle'
        }

    ar = {
        'in_file': 'tmp_reads_r.fas',
        'ref_file': ref_file,
        'out_file': 'tmp_align_r.needle'
        }
    
    if not( os.path.exists('tmp_align_f.needle') and os.path.exists('tmp_align_r.needle') ):
        logfun.info('running the alignments')
        strands = [af, ar]
        pool = Pool(processes=2)
        pool.map(align_strand, strands)
    else:
        logfun.warning('reading alignments from files tmp_align_[fr].needle')
    
    ref, reads, insertions, descriptions, global_ins = parse_alignments(thresh)
    
    if pad_insert:
        logfun.info('pad alignment, insertions are propagated')
        pad_alignment(ref, reads, insertions, descriptions, global_ins, amp_thresh, f_output)
    else:
        logfun.info('merge alignment, insertions are not propagated')
        merge_alignment(ref, reads, descriptions, amp_thresh, f_output)
        
    try:
        os.remove('tmp_reads_f.fas')
    except:
        pass
    try:
        os.remove('tmp_reads_r.fas')
    except:
        pass
    
    if not keep_files:
        logfun.info('removing alignment files')
        try:
            os.remove('tmp_align_f.needle')
        except:
            pass
        try:
            os.remove('tmp_align_r.needle')
        except:
            pass

if __name__ == "__main__":
    
    options, args = parse_com_line()
    
    ref_file = options.ref
    in_file = options.input
    out_file = options.o
    thresh = options.threshold
    pad_insert = options.pad
    amp_thresh = options.amplicon
    try:
        del_files = options.d
        keep_files = False
    except:
        keep_files = True
    keep_files = True
    main(ref_file, in_file, out_file, thresh, pad_insert, keep_files)
