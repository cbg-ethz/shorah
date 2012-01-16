#!/usr/bin/env python
# Thanks to Moritz Gerstung
# and AB

import pysam
import sys
import Bio.Seq
import Bio.SeqRecord
import Bio.SeqIO
import Bio.Align.Generic
from Bio.Alphabet import IUPAC, Gapped
alphabet = Gapped(IUPAC.ambiguous_dna)
SANGER_SCORE_OFFSET = ord("!")
q_mapping = dict()
for letter in range(0,255):
    q_mapping[chr(letter)] = letter-SANGER_SCORE_OFFSET
QC = 6



def getseq(AlignedRead, start=0, stop=0, qcutoff=QC):
    """Retrieve the sequence of an AlignedRead object between start and stop of the reference position.
    The output will be padded by N's if the region exceeds the read length.
    """
    rpos = 0 # position in the read
    gaps = 0 # number of gaps added
    fasta = str() # will hold the alignment
    
    # Aligned?
    if AlignedRead.is_unmapped:
        return

    rseq = AlignedRead.seq
    pos = AlignedRead.pos
    seq = ""

    for i,s in enumerate(rseq):
        if q_mapping[AlignedRead.qual[i]] >= qcutoff:
            seq += s
        else:
            seq += 'N'
    
    if AlignedRead.cigar != None:
        #print "Crap in:", AlignedRead.qname 
        for op in AlignedRead.cigar:
            if op[0] == 0:
                fasta = fasta + seq[rpos:(rpos+op[1])]
                rpos = rpos + op[1]
            else: 
                if op[0] == 2: # Deletion
                    for i in range(op[1]):
                        fasta = fasta + '-' #Insert gaps
                    gaps = gaps + op[1]
                else:
                    if op[0] == 1 or op[0] == 4: # HACK! Insertion, may cause trouble later..
                        #print "Insertions in:", AlignedRead.qname
                        rpos = rpos + op[1] # Skip
                    else:
                        pass
                        #print "Crap in: ", op[0], AlignedRead.qname
    else:
        fasta = fasta + seq
    
    # Pad the ends
    for i in range(pos - start):
        fasta = 'n' + fasta
    for i in range(stop - pos - AlignedRead.rlen - gaps):
        fasta = fasta + 'n'
    
    # Compute range to output
    begin = max(0, start - AlignedRead.pos)
    end = stop-start +1 + begin
    if not len(fasta[begin:end]) == end - begin:
        for i in range( end - begin - len(fasta[begin:end])):
            fasta += 'n' # Pad end because we omitted insertions
        #print fasta[begin:end]
    return fasta[begin:end]

def printfasta(AlignedRead, start=0, stop=0, chrom=None, out = sys.stdout, Nmax = 0.1, qcutoff = QC):
    """Retrieve the sequence of an aligned read and print in fasta format.
    """
    fasta = getseq(AlignedRead, start=start, stop=stop, qcutoff=qcutoff)
    STRND = ['+', '-']
    if fasta.count('N') < Nmax * len(fasta):
        print >>out, ">%s|%s:%i%s" %(AlignedRead.qname, chrom, AlignedRead.pos, STRND[AlignedRead.is_reverse])
        #print >>out, ">%s%s" %(AlignedRead.cigar, STRND[AlignedRead.is_reverse])
        print >>out, fasta
        return 1
    return 0

def getSeqRecord(AlignedRead, start=0, stop=0):
    """Get a SeqRecord object from AlignedRead between positions start and stop.
    """
    seq = getseq(AlignedRead, start=start, stop=stop)
    return Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq, alphabet), id = AlignedRead.qname, name = "%s:%i-%i" % (AlignedRead.rname, start, stop), description="", annotations = {"strand": AlignedRead.is_reverse, "CIGAR": AlignedRead.cigar})
    
def bam2fasta(samfile, chrom = None, start = None, stop = None, minlen = 1, out = sys.stdout, qcutoff = QC, strand = 2):
    """Extract the reads from an samfile and print as fasta.
    """
    iter = samfile.fetch(chrom, start, stop)
    i = 0
    for read in iter:
        soft_clipped = sum([op[1] for op in read.cigar if op[0] in (4,1)])
        if read.rlen - start + read.pos + 1 > minlen +soft_clipped and stop - read.pos +1 >= minlen +soft_clipped and (strand == 2 or read.is_reverse == strand):
            i += printfasta(read, chrom=chrom, start=start, stop=stop, out = out, qcutoff=qcutoff)
    print >>sys.stderr, "[samtools] Fetched %i reads from %s:%i-%i." %(i, chrom, start, stop)

#            i = i + 1
#        if i > 50:
#            break

def bam2Alignment(samfile, chrom = None, start = None, stop = None, minlen = 1, out = sys.stdout):
    """Read alignment from samfile and return Alignment object.
    """
    iter = samfile.fetch(chrom, start, stop)
    from Bio.Align.Generic import Alignment 
    from Bio.Alphabet import IUPAC, Gapped
    alphabet = Gapped(IUPAC.ambiguous_dna)
    aln = Alignment(alphabet)
    for read in iter:
        soft_clipped = sum([op[1] for op in read.cigar if op[0] in (4,1)])
        #print soft_clipped, read_cigar
        if read.rlen - start + read.pos + 1 > minlen + soft_clipped and stop - read.pos +1 >= minlen + soft_clipped:
            aln._records.append(getSeqRecord(read, start=start, stop=stop))
    return aln

def getHaplotypes(aln, n = 10, fmin = 0.0):
    """Get the haplotypes of the aligment aln.
    """
    count = {}
    from Bio.Align.Generic import Alignment 
    haplotypes = Alignment(alphabet)
    for record in aln:
        count[record.seq.tostring()] = count.get(record.seq.tostring(),0) + 1
    for i,seq in enumerate(sorted( count.keys(), key=lambda x: count[x], reverse=True)[:n]):
        f = count[seq]/float(len(aln))
        if f > fmin:
            haplotypes._records.append(Bio.SeqIO.SeqRecord(Bio.Seq.Seq(seq, alphabet), id = "Hap%04i" % (i +1), name =  "%f" % (f)))
    return haplotypes
    
if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage = "%prog [options] <bamfile>")
    parser.add_option("-c", '--chromosome', type = 'string', dest = 'chrom', default = None, help = "Chromosome [None]")
    parser.add_option("-l", '--left-end', type = 'int', dest = 'start', default = 0, help = "Left end of region [0]")
    parser.add_option("-r", '--right-end', type = 'int', dest = 'stop', default = 100, help = "Right end of region [100]")
    parser.add_option("-m", '--min-length', type = 'int', dest = 'minlen', default = 50, help = "Minimal length of unpadded reads [50].")
    parser.add_option("-q", '--quality-cutoff', type = 'int', dest = 'qcutoff', default = QC, help = "Mask bases with quality lower than QCUTOFF by 'N'.")
    parser.add_option("-s", '--strand', type = 'int', dest = 'strand', default = 2, help = "Only report alignment from strand (0: forward, 1: reverse, 2: both) [2].")
    parser.add_option("-o", '--outfile', type = 'string', dest = 'outfile', default = None, help = "Print output [stdout].")
    (options, args) = parser.parse_args(sys.argv[1:])
    if len(args) != 1:
        parser.error("Incorrect number of arguments.")
    samfile = pysam.Samfile(args[0], 'rb')
    if options.outfile:
        outfile = open(options.outfile, 'w')
    else:
        outfile = sys.stdout
    bam2fasta(samfile, chrom = options.chrom, start = options.start, stop = options.stop, minlen = options.minlen, out = outfile, qcutoff=options.qcutoff, strand = options.strand)
    samfile.close()
    outfile.close()
    sys.exit()
