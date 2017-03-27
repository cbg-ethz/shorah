# utilities for reading files of reads
#
# Copyright 2007, 2008, 2009
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

package ReadUtil;

use strict;
use NikoUtil;
use DNA;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;

$VERSION = 1.00;
@ISA = qw(Exporter);

# autoexport:
@EXPORT = qw(read_ReadSim
readReadFile readFastaFile
compareToCons compareReadToCons
);
# export on request:
@EXPORT_OK = qw();



sub read_ReadSim {
	my $ReadSimOut = shift;
	my @readList;

	my ($h, $s) = read_fastas($ReadSimOut);
	foreach my $i (0 .. @$s-1) {
		# process $h->[$i] =  read_0001|beg|2|length|92|forward|HIV
		# and     $s->[$i] = TCAAATCACTCTTTGGCAACGACCCCTCGTCGACAATAA

		#split on '|' doesn't work well, so first substitute spaces
		$h->[$i] =~ s/\|/ /g;
		my @head = split ' ', $h->[$i];

		die "missing header $h->[$i]\n" if (@head <= 1);
		my $seqID = $head[0];
		my $startpos = $head[2];
		my $dir = $head[5];
		my $source = $head[6];

		my %read;
		$read{id} = $seqID;
		$read{st} = $startpos;
		if ($dir eq 'forward') {
			$read{seq} = $s->[$i];
		} else {
			$read{seq} = revComp($s->[$i]);
		}
		$read{len} = length($s->[$i]);
		$read{end} = $startpos + $read{len} - 1;
		$read{dir} = $dir;
		$read{source} = $source;
		push(@readList, \%read);
	}

	return \@readList;
}

sub readReadFile {
	# read a file of reads given by 
	# startpos sequence
	# startpos sequence
	# startpos sequence
	# ...
	my $file = $_[0];
	my @readList;

	if (! defined $file) {
		die "Error: No fasta file specified\n";
	} 
	
	if (! -r $file) {
		die "Fasta file \"$file\" not readable\n";
	}
    open(IN, "<$file") or die "Can't open $file: $!";
    while (<IN>) {
		chomp;
		next if (/^#/);
		my ($start, $seq) = split;
		my %read;
		$read{st} = $start;
		$read{seq} = $seq;
		$read{len} = length($seq);
		$read{end} = $start + $read{len} - 1;
		push(@readList, \%read);
	}
	return \@readList;
}




sub readFastaFile {
    my $file = $_[0];
	my $seqID;
	my $startpos;
	my @readList;
	my $seq = '';
	my @head;

	if (! defined $file) {
		die "Error: No fasta file specified\n";
	} 
	
	if (! -r $file) {
		die "Fasta file \"$file\" not readable\n";
	}
    open(IN, "<$file") or die "Can't open $file: $!";
    while (<IN>) {
		chomp;
		if (/^\#/) {
			next;
		} elsif ($_ eq '') {
			next;
		} elsif (/^>/) {
			# header looks like ">seqID startpos ..."
			if (length($seq) > 0) {
				#we've found a sequence, so store it

				my %read;
				$read{id} = $seqID;
				$read{st} = $startpos;
				$read{seq} = $seq;
				$read{len} = length($seq);
				$read{end} = $startpos + $read{len} - 1;
				#my @tmp = @head;
				#$read{qs} = \@tmp;
				push(@readList, \%read);
				$seq = '';
			}
			@head = split;
			if (@head > 1) {
				$seqID = shift @head; #$head[0];
				$seqID =~ s/>//;
				$startpos = shift @head; #$head[1];
			} else {
				print "Error - @head missing characters\n";
			}
		} else {
			#a sequence line
			#remove spaces at beginning and end
			s/^\s+//;
			s/\s+$//;
			$seq .= $_;
		}

	}
	#still have sequence in the queue, need to put the last one in
	#this should be done in a better fashion...
	my %read;
	$read{id} = $seqID;
	$read{st} = $startpos;
	$read{seq} = $seq;
	$read{len} = length($seq);
	$read{end} = $startpos + $read{len} - 1;
#	$read{qs} = \@head;
	push(@readList, \%read);

	close IN;
	return \@readList;
}

sub compareToCons {
	my $seq = shift;
	my $cons = shift;
	my @ans;
	#die "$seq\n$cons\n" if (length($seq) != length($cons));
	my $minl = Nmin(length($seq), length($cons));
	#input: two strings of same length
	#output: list of differences in format A90I
	foreach my $i (1..$minl) {
		my $c1 = substr($seq,$i-1,1);
		my $c0 = substr($cons,$i-1,1);
		if (!($c1 eq $c0)) {
			push(@ans, $c0 . $i . $c1);
		}
	}
	return \@ans;
}

sub compareReadToCons {
	# not really right, since  compareToCons
	# prints in local coordinates...
	my $read = shift;
	my $cons = shift;
	#input a read %r = {id:AHAHA seq : ACGTAT, len : 6, start : 80}
	#and a string
	#output: %r = {id: AHAHA, ..., muts = join(',', compareToCons(read, cons) }
	#my %r;
	#$r{id} = $read->{id};
	#$r{len} = $read->{len};
	$read->{muts} = join(',', compareToCons($read->{seq}, substr($cons, $read->{st}, $read->{len})));
}

1;
