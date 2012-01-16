# Copyright Osvaldo Zagordi 2008, 2012
"""
This module contains a parser for the pairwise alignments produced by
EMBOSS software with format markx10
"""
# Adapted from FastaIO.py by Peter Cook
import os
import sys
from StringIO import StringIO
from Bio.Alphabet import generic_alphabet, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Generic import Alignment
from Bio.Align import MultipleSeqAlignment
from Bio.AlignIO.Interfaces import AlignmentIterator

class Markx10Iterator(AlignmentIterator) :
    """Alignment iterator for the EMBOSS tool's pairwise alignment output.

    This is for reading the pairwise alignments output by EMBOSS
    program when called with the -aformat markx10 command line option for machine
    readable output.  For more details about the FASTA tools, see the website:
    http://emboss.sourceforge.net/docs/themes/AlignFormats.html
    http://emboss.sourceforge.net/docs/themes/alnformats/align.markx10

    This class is intended to be used via the Bio.AlignIO.parse() function
    by specifying the format as "markx10" as shown in the following code:

    for a in MyAlignIO.parse(handle, "markx10") :
            assert len(a.get_all_seqs()) == 2, "Should be pairwise!"
            print "Alignment length %i" % a.get_alignment_length()

    Note that this is not a full blown parser for all the information
    in the output - for example, most of the header and all of the
    footer is ignored.  Also, the alignments are not batched according to
    the input queries.
    """
    
    def next(self) :
        """Reads from the handle to construct and return the next alignment.

        This returns the pairwise alignment of query and match/library
        sequences as an Alignment object containing two rows."""

        handle = self.handle
        try :
            #Header we saved from when we were parsing
            #the previous alignment.
            line = self._header
            print self._header.strip(), '--> self_header'
            del self._header
        except AttributeError:      
            line = handle.readline()
        if not line:
            return None

        if line.startswith('#-') :
            #Reached the end of the alignments, no need to read the footer...
            return None
        if line.startswith("##") :
            #Skip the file header before the alignments.  e.g.
#            print line.strip()
            line = self._skip_file_header(line)
#        print 'Back from file header skip'
        assert line.startswith('#'), line

        while not line.startswith('#=') :
            line = self.handle.readline()

        if line.startswith('#='):
            #Moved onto the next query sequence!
            self._query_descr = ""
            self._query_header_annotation = {}
            #Read in the query header
            line = self._parse_query_header(line)
        if not line :
            #End of file
            return None

        
        assert line.startswith(">>") and not line.startswith(">>>"), line

        query_seq_parts, match_seq_parts = [], []
        query_annotation, match_annotation = {}, {}
        match_descr = ""
        alignment_annotation = {}

        #This should be followed by the target match numbering line, then more tags.
        #e.g.
        """
        >>#2
        ; sw_score: 41.0
        ; sw_ident: 0.846
        ; sw_overlap: 13
        """
        
        if not line.startswith(">>") and not line.startswith(">>>") :
            raise ValueError("Expected target line starting '>>'")
        match_descr = line[2:].strip()
        #print match_descr, 'match'
        #Handle the following "alignment hit" tagged data, e.g.
        line = handle.readline()
        line = self._parse_tag_section(line, alignment_annotation)
        assert not line.startswith("; ")

        #Then we have the alignment numbers and sequence for the query
        """
        >gi|10955265| ..
        ; sq_len: 346
        ; sq_offset: 1
        ; sq_type: p
        ; al_start: 197
        ; al_stop: 238
        ; al_display_start: 167
        DFMCSILNMKEIVEQKNKEFNVDIKKETIESELHSKLPKSIDKIHEDIKK
        QLSC-SLIMKKIDVEMEDYSTYCFSALRAIEGFIYQILNDVCNPSSSKNL
        GEYFTENKPKYIIREIHQET
        """
        if not (line.startswith(">") and line.strip().endswith("..")):
            raise ValueError("Expected line starting '>' and ending '..'")
        assert self._query_descr.startswith(line[1:].split()[0])
        
        #Handle the following "query alignment" tagged data
        line = handle.readline()
        line = self._parse_tag_section(line, query_annotation)
        assert not line.startswith("; ")

        #Now should have the aligned query sequence (with leading flanking region)
        while not line.startswith(">") :
            query_seq_parts.append(line.strip())
            line = handle.readline()
#            print 'queryseq', line.strip()
        #Handle the following "match alignment" data
        """
        >gi|152973545|ref|YP_001338596.1| ..
        ; sq_len: 242
        ; sq_type: p
        ; al_start: 52
        ; al_stop: 94
        ; al_display_start: 22
        IMTVEEARQRGARLPSMPHVRTFLRLLTGCSRINSDVARRIPGIHRDPKD
        RLSSLKQVEEALDMLISSHGEYCPLPLTMDVQAENFPEVLHTRTVRRLKR
        QDFAFTRKMRREARQVEQSW
        """
        #Match identifier
        if not (line.startswith(">") and line.strip().endswith("..")):
            raise ValueError("Expected line starting '>' and ending '..', got '%s'" % repr(line))
        #print '----->', line.strip(), match_descr
        match_descr = line[1:].split()[0] + match_descr
        
        #assert match_descr.startswith(line[1:].split()[0])
#        assert self._match_descr.startswith(line[1:].split()[0])

        #Tagged data,
        line = handle.readline()
        line = self._parse_tag_section(line, match_annotation)
        assert not line.startswith("; ")
        
        #Now should have the aligned query sequence with flanking region...
        while not (line.startswith(">") or ">>>" in line) and not line.startswith('#'):
            match_seq_parts.append(line.strip())
            line = handle.readline()
            if not line:
                #End of file
                return None
        if line.startswith('>') or '>>>' in line:
            self._header = line

        #We built a list of strings and then joined them because
        #its faster than appending to a string.
        query_seq = "".join(query_seq_parts)
        match_seq = "".join(match_seq_parts)
        del query_seq_parts, match_seq_parts
        #Note, query_seq and match_seq will usually be of different lengths, apparently
        #because in the m10 format leading gaps are added but not trailing gaps!

        #Remove the flanking regions,
        query_align_seq = self._extract_alignment_region(query_seq, query_annotation)
        match_align_seq = self._extract_alignment_region(match_seq, match_annotation)

        #The "sq_offset" values can be specified with the -X command line option.
        #The appear to just shift the origin used in the calculation of the coordinates.
        
        if ("sq_offset" in query_annotation and query_annotation["sq_offset"] != "1") \
        or ("sq_offset" in match_annotation and match_annotation["sq_offset"] != "1") :
            #Note that until some point in the v35 series, FASTA always recorded one
            #for the query offset, and ommitted the match offset (even when these were
            #query_seq the -X command line option).
            #TODO - Work out how exactly the use of -X offsets changes things.
            #raise ValueError("Offsets from the -X command line option are not (yet) supported")
            pass

# this is not useful when using stretcher
#        if len(query_align_seq) != len(match_align_seq) :
#            raise ValueError("Problem parsing the alignment sequence coordinates")
        if "sw_overlap" in alignment_annotation :
            if int(alignment_annotation["sw_overlap"]) != len(query_align_seq) :
                raise ValueError("Specified sw_overlap = %s does not match expected value %i" \
                                 % (alignment_annotation["sw_overlap"],
                                    len(query_align_seq)))

        #TODO - Look at the "sq_type" to assign a sensible alphabet?
        q_rec = SeqRecord(Seq(query_align_seq), id=self._query_descr.split()[0].strip(","))
        m_rec = SeqRecord(Seq(match_align_seq), id=match_descr.split()[0].strip(","))
        
        alignment = MultipleSeqAlignment([q_rec, m_rec], self.alphabet)
        #TODO - Introduce an annotated alignment class?
        #For now, store the annotation a new private property:
        alignment._annotations = {}
        
        #Want to record both the query header tags, and the alignment tags.
        for key, value in self._query_header_annotation.iteritems() :
            alignment._annotations[key] = value
        for key, value in alignment_annotation.iteritems() :
            alignment._annotations[key] = value
                #        for k, v in query_annotation.iteritems():
                #            print k, v
        
        record = list(alignment)[0]
                #record.name = "query"
        record.annotations["original_length"] = int(query_annotation["sq_len"])
        # Roba mia
        for k in query_annotation.keys():
            record.annotations[k] = query_annotation[k]

        record = list(alignment)[1]
        assert record.id == match_descr.split()[0].strip(",") or record.description == match_descr
        # assert record.seq.tostring() == match_align_seq
        # record.id = match_descr.split()[0].strip(",")
        #record.name = "match"
        record.annotations["original_length"] = int(match_annotation["sq_len"])
        # Roba mia
        for k in query_annotation.keys():
            record.annotations[k] = match_annotation[k]
        return alignment

    def _skip_file_header(self, line) :
        """Helper function for the main parsing code.

        Skips over the file header region.
        """
        #e.g. This region:
        """
        ########################################
        # Program: water
        # Rundate: Thu  5 Jun 2008 16:56:06
        # Commandline: water
        #    -asequence short1.fas
        #    -bsequence short2.fas
        #    -gapopen 10.0
        #    -gapextend 0.5
        #    -outfile short1.faswater
        #    -aformat3 markx10
        #    -outfile out
        # Align_format: markx10
        # Report_file: out
        ########################################
        """
        #Note that it is not recording the command line here
        #from the # line
        line = self.handle.readline()
        while not line.startswith('##') :
            line = self.handle.readline()

        assert line.startswith('##'), 'The header file should finish here'
        return line
    
    def _parse_query_header(self, line) :
        """Helper function for the main parsing code.

        Skips over the free format query header, extracting the tagged parameters.

        If there are no hits for the current query, it is skipped entirely."""
#        print '------>>>>>NOW IN parse_quesry_header<<<<<<<<<----------'
        
        #e.g. this region (where there is often a histogram too):
        """
        #=======================================
        #
        # Aligned_sequences: 2
        # 1: short1
        # 2: short2
        # Matrix: EDNAFULL
        # Gap_penalty: 10.0
        # Extend_penalty: 0.5
        #
        # Length: 13
        # Identity:      12/13 (92.3%)
        # Similarity:    12/13 (92.3%)
        # Gaps:           1/13 ( 7.7%)
        # Score: 50.0
        # 
        #
        #=======================================
        """
        #Sometime have queries with no matches, in which case we continue to the next query block:
        
        self._query_header_annotation = {}
        self._query_descr = ""
        
        while not line.startswith('#=') :
            line = self.handle.readline()
        
        line = self.handle.readline()
        #We ignore the free form text...
        while not line.startswith(">>>") :
            #print "Ignoring %s" % line.strip()
            line = self.handle.readline()
            if not line :
                raise ValueError("Premature end of file!")
            if line.startswith('#-') :
                #End of alignments, looks like the last query
                #or queries had no hits.
                return line

        #Now want to parse this section:
        """
        >>>short1, 14 nt vs short3, 13 nt
        ; mp_name: EMBOSS
        ; mp_ver: 5.0.0
        ; pg_name: water
        ; pg_ver: 5.0.0
        ; pg_matrix: EDNAFULL
        ; pg_gap-pen: -10.0 -0.5
        """
        assert line.startswith(">>>"), line
        self._query_descr = line[3:].strip()

        #Handle the following "program" tagged data,
        line = self.handle.readline()
        line = self._parse_tag_section(line, self._query_header_annotation)

        assert not line.startswith("; "), line
        assert line.startswith(">>"), line
        return line


    def _extract_alignment_region(self, alignment_seq_with_flanking, annotation) :
        """Helper function for the main parsing code.

        To get the actual pairwise alignment sequences, we must first
        translate the un-gapped sequence based coordinates into positions
        in the gapped sequence (which may have a flanking region shown
        using leading - characters).  To date, I have never seen any
        trailing flanking region shown in the m10 file, but the
        following code should also cope with that.

        Note that this code seems to work fine even when the "sq_offset"
        entries are prsent as a result of using the -X command line option.
        """
        align_stripped = alignment_seq_with_flanking#.strip("-")
        start = int(annotation['al_start']) \
              - int(annotation['al_display_start'])
        end   = int(annotation['al_stop']) \
              - int(annotation['al_display_start']) \
              + align_stripped.count("-") + 1
        return align_stripped[start:end]


    def _parse_tag_section(self, line, dictionary) :
        """Helper function for the main parsing code.

        line - supply line just read from the handle containing the start of
               the tagged section.
        dictionary - where to record the tagged values

        Returns a string, the first line following the tagged section."""
#        print '------>>>>>NOW IN parse_tag_section<<<<<<<<<----------', line


        if not line.startswith("; ") :
            raise ValueError("Expected line starting '; '")
        while line.startswith("; ") :
            tag, value = line[2:].strip().split(": ",1)
            #fasta34 and fasta35 will reuse the pg_name and pg_ver tags
            #for the program executable and name, and the program version
            #and the algorithm version, respectively.  This may be a bug.
            #if tag in dictionary :
            #    raise ValueError("Repeated tag '%s' in section" % tag)
            dictionary[tag] = value
#            print line.strip(), tag, value
            line = self.handle.readline()
        return line

def main():
    import sys
    f = sys.argv[1]
    forwardaligniter = Markx10Iterator(open(f))
    print >> sys.stderr, 'here'
    print forwardaligniter
    f_align = forwardaligniter.next()
        #assert f_align.get_all_seqs()[1].id == r_align.get_all_seqs()[1].id, 'same seq back and forward'
        #descr = f_align.get_all_seqs()[1].id
    descr = list(f_align)[0].id
    print f_align[0]
    print descr
    

if __name__ == '__main__':
    main()