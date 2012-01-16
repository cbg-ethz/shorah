#!/usr/bin/env python
# Python implementation of an ABIF file reader according to Applied Biosystems' specificatons,
# see http://www.appliedbiosystems.com/support/software_community/ABIF_File_Format.pdf
#
# This code is published by Interactive Biosoftware, France,
# see http://www.interactive-biosoftware.com/
# under GPL license,
# see http://www.gnu.org/licenses/gpl.html
#
# Author: Francis Wolinski
# Version: 1.0, March 2007
# Copyright (c) Francis Wolinski 2007
#
# User Manual
#
# Conversion of ABIF data types to Python types (see struct.unpack method):
# type 1 = byte -> integer
# type 2 = char -> string
# type 3 = word -> long
# type 4 = short -> integer
# type 5 = long -> integer
# type 7 = float -> float
# type 8 = double -> float
# type 10 = date -> datetime.date instance
# type 11 = time -> datetime.time instance
# type 12 = thumb -> tuple
# type 13 = bool -> True or False
# type 18 = pString -> string
# type 19 = cString -> string
# type = 1024+ = user -> NotImplemented: to be overwritten in user's code in ABIFReader.readNextUserData method
# type = other -> NotImplemented
#
# from ABIFReader import *
# reader = ABIFReader(<filename>) # creates an instance of ABIFReader
# reader.version # version of ABIF file
# reader.showEntries() # print all entries of ABIF file "<name> (<num>) / <type> (<size>)"
# data = reader.getData(<name>[, <num>]) # read data for entry named <name> with number <num>, by default <num> is 1
# reader.close() # close the file, since it is kept open
#

#from numpy import *
import struct
import datetime
import sys
#import Gnuplot, Gnuplot.funcutils

SFF_TYPES = {1: 'byte', 2: 'char', 3: 'word', 4: 'short', 5: 'long', 7: 'float', 8: 'double',\
        10: 'date', 11: 'time', 12: 'thumb', 13: 'bool', 18: 'pString', 19: 'cString'}
verbose = False

class SFFReader:
    def __init__(self, fn):
        self.filename = fn
        self.file = open(fn, 'rb')
        self.type = self.readNextString(4)
        #self.type = self.readNextUnsignedInt_32()
        #self.type = self.readNextDouble()

        if self.type != '.sff':
            self.close()
            raise SystemExit("error: No SFF file '%s'" % fn)
        self.version = self.readNextString(4)
	self.index_offset = self.readNextUnsignedInt_64()
	self.index_length = self.readNextUnsignedInt_32()
	self.number_of_reads = self.readNextUnsignedInt_32()
	self.header_length = self.readNextUnsignedInt_16()
	self.key_length = self.readNextUnsignedInt_16() 
	self.number_of_flows_per_read = self.readNextUnsignedInt_16()
	self.flowgram_format_code = self.readNextUnsignedInt_8()
	self.flow_chars = self.readNextString(self.number_of_flows_per_read)
	self.key_sequence = self.readNextString(self.key_length)
	self.eight_byte_padding = self.readNextByte()
	while self.file.tell() % 8 != 0:
		self.eight_byte_padding = self.readNextByte()

        if verbose:
            print "version:\t" + str(self.version)
            print "type:\t"+ str(self.type)
            print "index_offset:\t" + str(self.index_offset)
            print "index_length:\t" + str(self.index_length)
            print "number_of_reads:\t" + str(self.number_of_reads)
            print "header_length:\t" + str(self.header_length )
            print "key_length:\t" + str(self.key_length )
            print "number_of_flows_per_read:\t" + str(self.number_of_flows_per_read)
            print "flowgram_format_code:\t" + str(self.flowgram_format_code)
            print "flow_chars:\t" + str(self.flow_chars)
            print "key_sequence:\t" + str(self.key_sequence)
            print "eight_byte_padding:\t" + str(self.eight_byte_padding)


	self.curPos = self.file.tell()
	if self.index_offset != 0 and self.index_length != 0:
		self.seek(self.index_offset)
		self.index_magic_number = self.readNextUnsignedInt_32()
        
                if verbose:
                    print "index_magic_number:\t" + str(self.index_magic_number)
        
		self.index_version = self.readNextString(4)
                
                if verbose:
                    print "index_version:\t" + str(self.index_version)
                
		while self.file.tell() % 8 != 0:
			self.eight_byte_padding = self.readNextByte()
                        
                        if verbose:
                            print "eight_byte_padding:\t" + str(self.eight_byte_padding)

	self.file.seek(self.curPos)



       	self.reads = self.getReads()

    def getReads(self):
    	reads = {}
	for n in range(self.number_of_reads):

		reads[n] = readEntry(self)
	
	return reads.values()
    
    
    def getData(self, name, num = 1):
        entry = self.getEntry(name, num)
        if not entry:
            raise SystemExit("error: Entry '%s (%i)' not found in '%s'" % (name, num, self.filename))
        self.seek(entry.mydataoffset())
        data = self.readData(entry.elementtype, entry.numelements)
        if data != NotImplemented and len(data) == 1:
            return data[0]
        else:
            return data

    def showEntries(self):
        for e in self.entries:
            print e

    def getEntry(self, name, num):
        for e in self.entries:
            if e.name == name and e.number == num:
                return e
        return None

    def readData(self, type, num):
        if type == 1:
            return [self.readNextByte() for i in range(num)]
        elif type == 2:
            return self.readNextString(num)
        elif type == 3:
            return [self.readNextUnsignedInt() for i in range(num)]
        elif type == 4:
            return [self.readNextShort() for i in range(num)]
        elif type == 5:
            return [self.readNextLong() for i in range(num)]
        elif type == 7:
            return [self.readNextFloat() for i in range(num)]
        elif type == 8:
            return [self.readNextDouble() for i in range(num)]
        elif type == 10:
            return [self.readNextDate() for i in range(num)]
        elif type == 11:
            return [self.readNextTime() for i in range(num)]
        elif type == 12:
            return [self.readNextThumb() for i in range(num)]
        elif type == 13:
            return [self.readNextBool() for i in range(num)]
        elif type == 18:
            return self.readNextpString()
        elif type == 19:
            return self.readNextcString()
        elif type >= 1024:
            return self.readNextUserData(type, num)
        else:
            return NotImplemented

    def readNextBool(self):
        return readNextByte(self) == 1

    def readNextByte(self):
        return self.primUnpack('B', 1)

    def readNextChar(self):
        return self.primUnpack('c', 1)

    def readNextcString(self):
        chars = []
        while True:
            c = self.readNextChar()
            if ord(c) == 0:
                return ''.join(chars)
            else:
                chars.append(c)

    def readNextDate(self):
        return datetime.date(self.readNextShort(), self.readNextByte(), self.readNextByte())

    def readNextDouble(self):
        return self.primUnpack('>d', 8)

    def readNextInt(self):
        return self.primUnpack('>i', 4)

    def readNextFloat(self):
        return self.primUnpack('>f', 4)

    def readNextLong(self):
        return self.primUnpack('>l', 4)

    def readNextpString(self):
        nb = self.readNextByte()
        chars = [self.readNextChar() for i in range(nb)]
        return ''.join(chars)

    def readNextShort(self):
        return self.primUnpack('>h', 2)

    def readNextString(self, size):
        chars = [self.readNextChar() for i in range(size)]
        return ''.join(chars)
    
    def readNextThumb(self):
        return (self.readNextLong(), self.readNextLong(), self.readNextByte(), self.readNextByte())

    def readNextTime(self):
        return datetime.time(self.readNextByte(), self.readNextByte(), self.readNextByte(), self.readNextByte())

    def readNextUnsignedInt(self):
        return self.primUnpack('>I', 4)

    def readNextUnsignedInt_64(self):
        return self.primUnpack('>Q', 8)

    def readNextUnsignedInt_32(self):
        return self.primUnpack('>I', 4)

    def readNextUnsignedInt_16(self):
        return self.primUnpack('>H', 2)

    def readNextUnsignedInt_8(self):
        return self.primUnpack('>B', 1)
    
    def readNextUserData(self, type, num):
        # to be overwritten in user's code
        return NotImplemented

    def primUnpack(self, format, nb):
        x = struct.unpack(format, self.file.read(nb))
        return x[0]
    
    def close(self):
        self.file.close()

    def seek(self, pos):
        self.file.seek(pos)

    def tell(self):
        return self.file.tell()

class readEntry:
    def __init__(self, reader):
        self.read_header_length = reader.readNextUnsignedInt_16()
	#print "read_header_length:\t" + str(self.read_header_length)
	self.name_length = reader.readNextUnsignedInt_16()
	#print "name_length:\t" + str(self.name_length )
	self.number_of_bases = reader.readNextUnsignedInt_32()
	#print "number_of_bases:\t" + str(self.number_of_bases )
	self.clip_qual_left = reader.readNextUnsignedInt_16()
	#print "clip_qual_left:\t" + str(self.clip_qual_left)
	self.clip_qual_right = reader.readNextUnsignedInt_16()
	#print "clip_qual_right:\t" + str( self.clip_qual_right )
	self.clip_adapter_left = reader.readNextUnsignedInt_16()
	#print "clip_adapter_left:\t" + str(self.clip_adapter_left)
	self.clip_adapter_right = reader.readNextUnsignedInt_16()
	#print "clip_adapter_right:\t" + str(self.clip_adapter_right)
	self.name = reader.readNextString(self.name_length)
	#print "name:\t" + str(self.name)
	while reader.file.tell() % 8 != 0:
		self.eight_byte_padding = reader.readNextByte()
		#print "eight_byte_padding:\t" + str(self.eight_byte_padding)


       

	self.flowgram_values = [0] * reader.number_of_flows_per_read
	for i in range(reader.number_of_flows_per_read):
	    	if reader.flowgram_format_code == 1:
		    	self.flowgram_values[i] = reader.readNextUnsignedInt_16()
		        #print "flowgram_values:\t" + str(self.flowgram_values)
	
        
        
	self.flow_index_per_base = [0] * self.number_of_bases
	for i in range(self.number_of_bases):
		self.flow_index_per_base[i] = reader.readNextUnsignedInt_8()
		#print "flow_index_per_base:\t" + str(self.flow_index_per_base[i] )
        
        
	self.bases = reader.readNextString(self.number_of_bases)
	#print "bases:\t" + str(self.bases )

        
	self.quality_scores = [0] * self.number_of_bases
	for i in range(self.number_of_bases):
	        self.quality_scores[i] = reader.readNextUnsignedInt_8()
		#print "quality_scores:\t" + str(self.number_of_bases )

	while reader.file.tell() % 8 != 0:
	        	self.eight_byte_padding = reader.readNextByte()
	
	
		

    def __str__(self):
        return "%s (%i) / %s (%i)" % (self.name, self.number, self.mytype(), self.numelements)

    def mydataoffset(self):
        if self.datasize <= 4:
            return self.dataoffsetpos
        else:
            return self.dataoffset

    def mytype(self):
        if self.elementtype < 1024:
            return ABIF_TYPES.get(self.elementtype, 'unknown')
        else:
            return 'user'

def histogram(data, start, stop, mesh):
    """Returns the histogram vector for data"""

    histo = [0]*mesh
    step = float(stop - start)/mesh
    for d in data:
        for i in range(mesh):
            if start + i*step < d and d <= start + (i+1)*step:
                histo[i] = histo[i] + 1
    # for i in range(mesh):
    #     print start + i*step, histo[i]
    hdata = [ [ start + i*step, histo[i] ] for i in range(mesh) ]
    return hdata

###############################################################
def main():
    filename = "E9K9WQI16.sff" 
    reader = SFFReader(filename) 
    reads = reader.reads
    bb = [i.number_of_bases for i in reads]
    th = histogram(bb, 200, 350, 50)
    print th[:10]
    g = Gnuplot.Gnuplot(debug=1)
    g.title('Read length from experiment %s, total=%i reads' % (filename, reader.number_of_reads) )
    g('set data style boxes')
    g.xlabel('read length')
    g.ylabel('counts')
    g.plot(th)
    
# when executed, just run demo():
if __name__ == '__main__':
    main()
