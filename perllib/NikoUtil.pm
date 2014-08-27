#
# NikoUtil.pm - Utilities
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

package NikoUtil;

use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;

$VERSION = 1.00;
@ISA = qw(Exporter);

# autoexport:
@EXPORT = qw(
		trim
		is_numeric
		round
		intrand  discr_rand
		concat_list
		read_fasta  read_fastas  print_fasta
		print_matrix  transpose  exp_matrix
		where_clause_list
		hamming
		hamming_distance  alphabet
		makeDistMatrix
		random
		Nmin Nmax
);

# export on request:
@EXPORT_OK = qw();

sub concat_list {  # concatenate lists

	my @lists = @_;  # list of refs. to lists to be concatenated

# concatenate all lists
	my @all; 
	for (@lists) {
		push @all, @{$_};
	}

# count occurences in concatenation
	my (@union, %count);
	for (@all) {
		$count{"@{$_}"}++;
		push @union, $_ unless ($count{"@{$_}"} > 1);
	}

	return \@union;  # multiple occurences appear only once
}



sub is_numeric {  # Cookbook p. 44
	my $val = shift || '';  # return FALSE if $val is FALSE

	return ($val =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)$/);
}



sub round {  # round to nearest integer

	my @input = @_;

	for (@input) {
		$_ = int($_ + 0.5);
	}

	return wantarray ? @input : $input[0];
}



sub intrand {

# random integer number between and including $a and $b

	my $a = shift;
	my $b = shift;

	return int($a) + int(rand($b-$a+1));
}



sub discr_rand {

# random integers in given proportion

	my @proportion = @_;
	my $n = scalar(@proportion);

	my @cum;
	my $i;

	my $sum = 0.0;
	for $i (0 .. $n - 1) {
		$sum += $proportion[$i];
	}
	$cum[0] = ($proportion[0] / $sum);
	for $i (1 .. $n - 1) {
		$cum[$i] = $cum[$i-1] + ($proportion[$i] / $sum);
	}

	my $r = rand();
#print "@cum  --->  $r\n";

	$i = 0;
	while (($i < $n-1) && ($cum[$i] < $r)) {
		$i++;
	}

	return $i;
}



sub trim {  # remove white spaces to the left and right

	my @input = @_;

	for (@input) {
		s/^\s+//;
		s/\s+$//;
	}

	return wantarray ? @input : $input[0];
}



sub read_fasta {

	my $filename = shift;

	my $sep = chr(13);  # CR

	open FASTA, $filename 
		or die "Can't open fasta file $filename!\n";
	my @lines = <FASTA>;
	close FASTA;

	for my $i (0 .. $#lines) {
		splice (@lines, $i, 1, (split /$sep/, $lines[$i]));
	}

	return (substr($lines[0], 1), join('', @lines[1..$#lines]));
}



sub read_fastas {

	my $filename = shift;

	my (@header, @seq);
	my $i = -1;

	open FASTA, $filename 
		or die "Can't open fasta file -- $filename!\n";

	while (my $line = <FASTA>) {
		next if ($line =~ /^#/);
		if (substr($line, 0, 1) eq ">") {  # header
			$i++;
			$header[$i] = trim(substr($line, 1));
		}
		else {  # sequence
			die "Error - fasta file without header lines @ $line\n" if ($i == -1);
			$seq[$i] .= trim($line);
		}
	}

	close FASTA;

	return (\@header, \@seq);

}



sub print_fasta {

	my $head = shift;
	my $seq  = shift;

	my $width = 60;

	my @seq = split //, $seq;

	my $fasta = ">$head\n";
	while (scalar(@seq) > 0) {
		$fasta .= join('', splice(@seq, 0, $width)) . "\n";
	}

#return $fasta . "\n";
	return $fasta;
}



sub print_matrix {

# print matrix (tab- and newline-separated)

	my @mat       = @{shift @_};    # data im $mat[][]
	my $row_label = shift;          # row labels == first column to print
	my $col_label = shift;          # column labels == first row to print
	my $mat_label = shift || "\\";  # upper leftmost cell

	print $mat_label, "\t" if ($row_label);
	if ($col_label) {
		print join("\t", @$col_label), "\n";
	}

	for my $i (0 .. $#mat) {
		print $row_label->[$i], "\t" if ($row_label);
		for my $j (0 .. $#{$mat[$i]}) {
			if (defined($mat[$i][$j])) {
				print $mat[$i][$j], "\t";
			}
			else {
				print "undef\t";
			}
		}
		print "\n";
	}
}



sub transpose {

# transpose matrix

	my @M = @{shift @_};

	my @Mt;
	for my $i (0 .. $#M) {
		for my $j (0 .. $#{$M[0]}) {
			$Mt[$j][$i] = $M[$i][$j];
		}
	}

	return \@Mt;
}



sub exp_matrix {

	my $mat  = shift;
	my $base = shift || 10;

	my $emat;

	for my $i (0 .. scalar(@$mat)-1) {
		for my $j (0 .. scalar(@{$mat->[0]})-1) {
			$emat->[$i][$j] = $base ** $mat->[$i][$j];
		}
	}

	return $emat;
}



sub where_clause_list {

	my $field  = shift;
	my @values = @_;

	my @cond = ();

	for (@values) {
		if (is_numeric($_)) {
			push @cond, sprintf "%s=%d", $field, $_;
		}
		else {
			push @cond, sprintf "%s='%s'", $field, $_;
		}
	}

	return join(' or ', @cond);
}




sub Nmax {

# Maximum of a vector containing undef's

	my @input = @_;

	my ($max, $argmax);

	for my $i (0 .. $#input) {
		if (defined($input[$i])) {
			if (! defined($argmax)) {
				$max = $input[$i];
				$argmax = $i;
			}
			elsif ($input[$i] > $max) {  # same action, but avoids warning
				$max = $input[$i];
				$argmax = $i;
			}
		}
	}

	return wantarray ? ($max, $argmax) : $max;
}


sub Nmin {

# Minimum of a vector containing undef's

	my @input = @_;

	my ($min, $argmin);

	for my $i (0 .. $#input) {
		if (defined($input[$i])) {
			if (! defined($argmin)) {
				$min = $input[$i];
				$argmin = $i;
			}
			elsif ($input[$i] < $min) {  # same action, but avoids warning
				$min = $input[$i];
				$argmin = $i;
			}
		}
	}
	return wantarray ? ($min, $argmin) : $min;
}

sub hamming_distance {

	my $X = uc(shift);
	my $Y = uc(shift);

	my $min;
	if (length($X) > length($Y)) {
		$min = length($Y);
	} else {
		$min = length($X);
	}

	my $comparisons = 0;
	my $dist = 0;
	my ($x, $y);
	for my $pos (0 .. $min - 1) {
		$x = substr($X,$pos,1);
		$y = substr($Y,$pos,1);
		if (($x ne '.') && ($y ne '.')) {
			$comparisons++;
			if ($x ne $y) {
				$dist++;
			}
		}
	}

	#die "Incomparable strings:\n$X\n$Y\n" if ($comparisons == 0);
	if ($comparisons == 0) {
		#print "Error: incomparable strings:\n$X\n$Y\n";
		return 100;
	}

	return ($dist / $comparisons);
}

sub makeDistMatrix  {
	my $seqs = shift;
	my @ans;
	my $N = @$seqs;
	foreach my $i (0..$N-1) {
		foreach my $j (0..$N-1) {
			$ans[$i][$j] = hamming($$seqs[$i],$$seqs[$j]);
		}
		#print "$i\n" if ($i % 10 == 0);
	}
	#print "distance matrix:\n";
	#print_matrix(\@ans);
	#print "\n\n";
	return \@ans;
}

sub hamming {
	#hamming distance, absolute, not caring about '.'
	my $s = shift;
	my $t = shift;
	my $diff = 0;
	# compute numer of differences between $s and $t
	#die "ERROR in hamming $s\n$t\n" if (length($s) != length($t));
	my $len = Nmin(length($s), length($t));
	foreach my $i (0.. $len-1) {
		if (!(substr($s,$i,1) eq substr($t,$i,1))) {
			$diff++;
		}
	}
	return $diff;
}



sub alphabet {

	my $type = uc(shift);

	my (@Alphabet, @Special);

	if ($type eq 'AA') {
		@Alphabet = qw(A C D E F G H I K L M N P Q R S T V W Y);
		@Special = ("-", ".", "*");
	}
	elsif ($type eq 'NT') {
		@Alphabet = qw(A C G T);
		@Special = ("-", ".", "N");
	}
	elsif ($type eq '01') {
		@Alphabet = qw(0 1);
		@Special = ("-", ".");
	}
	else {
		die "Unknown alphabet -- $type!\n";
	}

	# indices:
	my %LtrIdx;
	for my $k (0 .. $#Special) {
		$LtrIdx{$Special[$k]} = -$k - 1;
	}
	for my $k (0 .. $#Alphabet) {
		$LtrIdx{$Alphabet[$k]} = $k;
	}

	return (\@Alphabet, \@Special, \%LtrIdx);
}

sub random {  # generate an integer random number                                                    

	my $x = shift;  # between $x
	my $y = shift;  # and $y (including both)

	return int(rand($y - $x + 1)) + $x;
}


1;
