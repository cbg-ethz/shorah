package Util;
# various statistical utilities

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

use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;
use NikoUtil;

$VERSION = 1.00;
@ISA = qw(Exporter);

# autoexport:
@EXPORT = qw(dec2bin findProg Psystem average median countChar
normalize maxind makeRandomDist makeRandomDistMin binomial fact rand_perm
computeInDelMut);
# export on request:
@EXPORT_OK = qw();

sub dec2bin {
	my $arg = shift;
	my $N = shift;
	substr(unpack("B32", pack('C', $arg)),-$N);
}

sub findProg {
	my $name = shift;
	my $ans;

	my @pathList = split(':',$ENV{'PATH'});
	push(@pathList,"$ENV{'HOME'}/bin");
	push(@pathList,'.');

	my $found = 0;
	foreach my $dir (@pathList) {
		if (-r "$dir/$name") {
			$ans = "$dir/$name";
			$found = 1;
			print "found $name as $ans\n";
			last;
		}
	}
	if ($found == 0) {
		print "error, can't find $name in PATH = @pathList\n";
		exit 1;
	}
	return $ans;
}

sub Psystem {
        my $s = shift;
        print "exec : $s\n";
        system ("$s");
}

sub Fsystem {
        my $s = shift;
        print "$s\n";
}

sub average {
	my $l = shift;
	if (@$l == 0) {
		return undef;
	}
	my $ave =  0;
	foreach (@$l) {
		$ave += $_;
	} 
	$ave /= @$l;
	return $ave;
}

sub median {
	#actually computes the mode
	my $l = shift;
	my %c;
	foreach (@$l) {
		$c{$_}++;
	}
	my $maxcnt = 0;
	my $med;
	foreach (keys(%c)) {
		if ($c{$_} > $maxcnt) {
			$maxcnt = $c{$_};
			$med = $_;
		}
	}
	return $med;
}

sub countChar {
	#count occurances of $char in $s
	#
	my $s = shift;
	my $char = shift;
	#print "counting $char in $s\n";
	# this little trick doesn't work if char = '.'
	#$_ = $s;
	#return tr/$char/$char/;

	my $count = 0;
	foreach (0.. length($s)) {
		if (substr($s,$_,1) eq $char) {
			$count++;
		}
	}
	#print "get $count\n";
	return $count;
}

sub normalize {
	#input: list of positive integers
	#output: list / sum
	my $l = shift;
	#print "Calling norm with @$l\n";
	my $sum = 0;
	foreach my $a (@$l) {
		if (defined $a) {
			$sum += $a;
		}
	}
	if ($sum > 0) {
		foreach my $a (@$l) {
			if (defined $a) {
				$a = $a / $sum;
			} else {
				$a = 0;
			}
		}
	}
}

sub maxind {
	my $l = shift;
	my $max = -1e100;
	my $maxind = 0;
	foreach my $i (0..@$l-1) {
		if ($l->[$i] > $max) {
			$max = $l->[$i];
			$maxind = $i;
		}
	}
	return $maxind;
}

sub makeRandomDist {
	my $n = int(shift);
	return 0 if ($n < 1);
	my @ans;

	foreach my $i (1 .. $n) {
		push(@ans, int(rand(1000)));
	}
	normalize(\@ans);
	return \@ans;
}

sub makeRandomDistMin {
	my $n = int(shift);
	my $minFreq = shift;
	return 0 if ($n < 1);
	die "$n * $minFreq > 1\n" if ($n * $minFreq > 1);
	my @ans;

	push @ans, 0;
	foreach my $i (1 .. $n-1) {
		push(@ans, int(rand(1000)));
	}
	# normalize so that sum @ans = 1 - $minFreq*$n
	my $num = 1 - $minFreq*$n;
	my $sum = 0;
	$sum += $_ foreach (@ans);
	$_ *= $num/$sum foreach (@ans);
	$_ += $minFreq foreach (@ans);
	return \@ans;
}

sub binomial {
	my ($n, $k) = round(@_);
	if ($k > $n) {
		return 0;
	}
	return (fact($n) / (fact($k) * fact($n - $k)));
}

sub fact {
	my $n = shift;
	if ($n == 0) {
		return 1;
	} else {
		my $ans = 1;
		foreach (1..$n) {
			$ans *= $_;
		}
		return $ans;
	}
}


sub rand_perm {

	# returns a random permutation of the $N integers 0, ..., $N-1

	my $N = shift;

	my @list = (0 .. $N-1);
	my @perm;

	while (@list) {
		my $r = random(0, $#list);
		push @perm, splice(@list, $r, 1);
	}

	return @perm;
}

sub computeInDelMut {
	my $a = shift;
	my $b = shift;
	my $del = countChar($a, '-');
	my $ins = countChar($b, '-');
	my $Mut = 0;
	foreach my $i (0 .. length($a) - 1) {
		my $x = substr($a,$i,1);
		my $y = substr($b,$i,1);
		if (($x ne $y) and ($x ne '-') and ($y ne '-')) {
			$Mut++;
		}
	}
	return ($ins, $del, $Mut);
}

1;
