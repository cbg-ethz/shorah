# DNA.pm - Utilities

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

package DNA;

use strict;
use NikoUtil;
use Util;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;

$VERSION = 1.00;
@ISA = qw(Exporter);

# autoexport:
@EXPORT = qw(
	translate
	makeCodons
	reverseCodons
	revComp
	NW SW
	generalizedSW
	blastAlign
	correctGaps
	correctGaps2
	ambigString
	);
# export on request:
@EXPORT_OK = qw();

sub translate {
	# input:
	# dna string
	# codon hash
	# optional: if ambiguous, should we output just one aa or
	# a whole block [CKEIF] ??

	my $dna = shift;
	my $codon = shift;
	my $JUSTONEAA = 0;
	if (@_ > 0) {
		$JUSTONEAA = shift;
	}
	my $pro = '';
	my $cod;
	my $choices = 1;
	my $ambig = 0;

	for(my $i=0; $i < length($dna); $i+=3) {
		$cod = substr($dna,$i,3);
		if (length($cod) == 3) {
			if (defined $codon->{$cod}) {
				$pro .= $codon->{$cod};
			} else {
				my $l = resolveCodon($cod);
				my %poss;
				foreach (@$l) {
					if (defined ($codon->{$_})) {
						$poss{$codon->{$_}} = 1;
					} else {
						$poss{' '} = 1;
					}
				}
				my $tmp = join('', sort(keys(%poss)));
				if (keys(%poss) > 1) {
						if ($JUSTONEAA) {
								$pro .= substr($tmp,0,1);
						} else {
								$pro .= '[' . $tmp . ']';
						}
						$choices *= keys(%poss);
						$ambig += 1;
				} else {
						$pro .= $tmp;
				}
#print "#error : don't know codon ", substr($dna,$i,3),"\n";
			}
		}
	}
	#print "$choices choices for protein\n";
	#print "$ambig ambiguous codons\n";
	#print "$pro\n";
	return $pro;
}

sub resolveCodon {
	#take codon with ambiguity codes and output list of all codons that could
	#give it
	my $c = shift;
	my @l = split '', $c;
	my @ans = ('', '', '');
	my @codons;
	foreach my $p (0..2) {
		if ($l[$p] eq 'A') {
			$ans[$p] .= 'A';
		} elsif ($l[$p] eq 'C') {
			$ans[$p] .= 'C';
		} elsif ($l[$p] eq 'G') {
			$ans[$p] .= 'G';
		} elsif ($l[$p] eq 'T') {
			$ans[$p] .= 'T';
		} elsif ($l[$p] eq 'R') {
			$ans[$p] .= 'AG';
		} elsif ($l[$p] eq 'Y') {
			$ans[$p] .= 'CT';
		} elsif ($l[$p] eq 'K') {
			$ans[$p] .= 'GT';
		} elsif ($l[$p] eq 'M') {
			$ans[$p] .= 'AC';
		} elsif ($l[$p] eq 'S') {
			$ans[$p] .= 'CG';
		} elsif ($l[$p] eq 'W') {
			$ans[$p] .= 'AT';
		} elsif ($l[$p] eq 'D') {
			$ans[$p] .= 'AGT';
		} elsif ($l[$p] eq 'H') {
			$ans[$p] .= 'ACT';
		} elsif ($l[$p] eq 'B') {
			$ans[$p] .= 'CGT';
		} elsif ($l[$p] eq 'V') {
			$ans[$p] .= 'ACG';
		} elsif ($l[$p] eq 'X') {
			$ans[$p] .= 'ACGT';
		} elsif ($l[$p] eq 'N') {
			$ans[$p] .= 'ACGT';
		}
	}
	#print "$c -> @ans\n";
	foreach my $i (0 .. length($ans[0])-1) {
		foreach my $j (0 .. length($ans[1])-1) {
			foreach my $k (0 .. length($ans[2])-1) {
				push(@codons, substr($ans[0],$i,1) . substr($ans[1],$j,1) . substr($ans[2],$k,1));
			}
		}
	}
	#print "Codon $c -> @codons\n";

	return \@codons;
}

sub makeCodons {
	my %codon;
	$codon{GCT} = $codon{GCC} = $codon{GCA} = $codon{GCG} = 'A';
	$codon{TTA} = $codon{TTG} = $codon{CTT} = $codon{CTC} = $codon{CTA} = $codon{CTG} = 'L';
	$codon{CGT} = $codon{CGC} = $codon{CGA} = $codon{CGG} = $codon{AGA} = $codon{AGG} = 'R';
	$codon{AAA} = $codon{AAG} = 'K';
	$codon{AAT} = $codon{AAC} = 'N';
	$codon{ATG} = 'M';
	$codon{GAT} = $codon{GAC} = 'D';
	$codon{TTT} = $codon{TTC} = 'F';
	$codon{TGT} = $codon{TGC} = 'C';
	$codon{CCT} = $codon{CCC} = $codon{CCA} = $codon{CCG} = 'P';
	$codon{CAA} = $codon{CAG} = 'Q';
	$codon{TCT} = $codon{TCC} = $codon{TCA} = $codon{TCG} = $codon{AGT} = $codon{AGC} = 'S';
	$codon{GAA} = $codon{GAG} = 'E';
	$codon{ACT} = $codon{ACC} = $codon{ACA} = $codon{ACG} = 'T';
	$codon{GGT} = $codon{GGC} = $codon{GGA} = $codon{GGG} = 'G';
	$codon{TGG} = 'W';
	$codon{CAT} = $codon{CAC} = 'H';
	$codon{TAT} = $codon{TAC} = 'Y';
	$codon{ATT} = $codon{ATC} = $codon{ATA} = 'I';
	$codon{GTT} = $codon{GTC} = $codon{GTA} = $codon{GTG} = 'V';
	#$codon{ATG} = '+';
	$codon{TAG} = $codon{TGA} = $codon{TAA} = '*';
	return \%codon;
}
sub reverseCodons {
	my $codon = shift;
	my %a2c;
	#$a2c{A} = [GCT, GCC, GCA, GCG]
	foreach my $cod (keys(%$codon)) {
		push(@{$a2c{$$codon{$cod}}},$cod);
	}
	return \%a2c;
}

sub revComp {
	my $seq = $_[0];
	my $seqRC = "";
	my $I;

	for($I = length($seq)-1; $I >= 0; $I--) {
		my $base = substr($seq, $I, 1);
		if($base eq "A") { 
			$seqRC .= "T"; 
		} elsif($base eq "C") { 
			$seqRC .= "G"; 
		} elsif($base eq "G") { 
			$seqRC .= "C"; 
		} elsif($base eq "T") { 
			$seqRC .= "A"; 
		} else {
			$seqRC .= $base;
		}
	}
   return($seqRC);
}


sub trans {
	my $x = shift;
	if ($x eq 'A') {
		return 0;
	} elsif ($x eq 'C') {
		return 1;
	} elsif ($x eq 'G') {
		return 2;
	} else {
		return 3;
	}
}

sub NW {
	# be careful if using... .fixes went into SW
	my $a = shift;
	my $b = shift;
	my $d;
	my $S;
	if (@_ > 0) {
		$d = shift; #linear gap
		$S = shift; #[][] = costs
	} else {
		$d = -400;
		$S = [[91,-114,-31,-123],[-114,100,-125,-31],[-31,-125,100,-114],[-123,-31,-114,91]];
	}

	my @A = split('', $a);
	my @B = split('', $b);
	my ($i, $j, $choice1, $choice2, $choice3);
	my $AlignmentA;
	my $AlignmentB;
	my @F;

	foreach $i (0..@A) {
		$F[$i][0] = $d * $i;
	}
	foreach $j (0..@B) {
		$F[0][$j] = $d * $j;
	}

	foreach $i (1..@A) {
		foreach $j (1..@B) {
			$choice1 = $F[$i-1][$j-1] + $S->[trans($A[$i-1])][trans($B[$j-1])];
			$choice2 = $F[$i-1][$j] + $d;
			$choice3 = $F[$i][$j-1] + $d;
			$F[$i][$j] = max($choice1,$choice2,$choice3);
		}
	}

	#foreach $i (0..@A) {
	#	foreach $j (0..@B) {
	#		print "$F[$i][$j] ";
	#	}
	#	print "\n";
	#}

	$AlignmentA = '';
	$AlignmentB = '';
	$i = @A;
	$j = @B;

	while ($i > 0 and  $j > 0) {
		my $Score = $F[$i][$j];
		my $ScoreDiag = $F[$i - 1][$j - 1];
		my $ScoreUp = $F[$i][$j - 1];
		my $ScoreLeft = $F[$i - 1][$j];
		my $Ai = trans($A[$i-1]);
		my $Bj = trans($B[$j-1]);

		if ($Score == $ScoreDiag + $S->[$Ai][$Bj]) {
			$AlignmentA = $A[$i-1] . $AlignmentA;
			$AlignmentB = $B[$j-1] . $AlignmentB;
			$i--;
			$j--;
		} elsif ($Score == $ScoreLeft + $d) {
			$AlignmentA = $A[$i-1] . $AlignmentA;
			$AlignmentB = '-' . $AlignmentB;
			$i--;
		} else { # (Score == ScoreUp + d)
			$AlignmentA = '-' . $AlignmentA;
			$AlignmentB = $B[$j-1] . $AlignmentB;
			$j--;
		}
	}
	while ($i > 0) {
		$AlignmentA = $A[$i-1] . $AlignmentA;
		$AlignmentB = '-' . $AlignmentB;
		$i--;
	}
	while ($j > 0) {
		$AlignmentA = '-' . $AlignmentA;
		$AlignmentB = $B[$j-1] . $AlignmentB;
		$j--;
	}
	return ($AlignmentA,$AlignmentB);
}


sub SW {

	#Smith Waterman local alignment
	# assume that first sequence is the short one
	my $a = shift;
	my $b = shift;
	my $d;
	my $S;
	if (@_ > 0) {
		$d = shift; #linear gap
		$S = shift; #[][] = costs
	} else {
		#$d = -100;
		#$S = [[91,-114,-31,-123],[-114,100,-125,-31],[-31,-125,100,-114],[-123,-31,-114,91]];
		$d = -4;
		$S = [[4,-4,-4,-4],[-4,4,-4,-4],[-4,-4,4,-4],[-4,-4,-4,4]];
	}
	my $DEBUG = 0;
	if (@_ > 0) {
		$DEBUG = shift;
	}


	my @A = split('', $a);
	my @B = split('', $b);
	my @tA = map{trans($_)} @A;
	my @tB = map{trans($_)} @B;
	#my @tA = map {trans($_)} split ('', $a);
	#my @tB = map {trans($_)} split ('', $b);
	my ($i, $j);
	my $AlignmentA;
	my $AlignmentB;
	my @F;
	my ($largestF, $largesti,$largestj);
	$largestF = 1;

	foreach $i (0..@A) {
		$F[$i][0] = 0;
	}
	foreach $j (0..@B) {
		$F[0][$j] = 0;
	}


	my ($c1, $c2, $c3, $m);

	foreach $i (1..@A) {
		foreach $j (1..@B) {
			$c1 = $F[$i-1][$j-1] + $S->[$tA[$i-1]][$tB[$j-1]];
			$c2 = $F[$i-1][$j] + $d;
			$c3 = $F[$i][$j-1] + $d;
			$m = $c1 > $c2 ? ($c1 > $c3 ? $c1 : $c3) : ($c2 > $c3 ? $c2 : $c3);
			$F[$i][$j] = $m > 0 ? $m : 0;
			#$F[$i][$j] = max($choice1,$choice2,$choice3,0);

			#$F[$i][$j] = max($F[$i-1][$j-1] + $S->[$tA[$i-1]][$tB[$j-1]],
			#	$F[$i-1][$j] + $d, $F[$i][$j-1] + $d, 0);

			if ($F[$i][$j] >= $largestF) {
				#print "max at $i $j -- $F[$i][$j]\n";
				$largestF = $F[$i][$j];
				$largesti = $i;
				$largestj = $j;
			}
		}
	}
	if ($DEBUG) {
		foreach my $row (@F) {
			print join("\t", @$row), "\n";
		}
	}


	#now trace back from the largest value in the matrix and build up the alignment

	$AlignmentA = '';
	$AlignmentB = '';
	$i = $largesti;
	$j = $largestj;
	#print "Largest F at $i, $j\n";
	#
	#	my $gapend = (length($a) - $largesti);
	#	if ($gapend > 0) {
	#		print "adding $gapend gaps at end\n";
	#		$AlignmentA = 'x' x $gapend;
	#		#foreach my $x ($j .. ($j + $gapend - 1)) {
	#		#	print "$x $B[$x]\n";
	#		#}
	#		#$AlignmentB = lc(join('',@B[$j .. ($j + $gapend-1)]));
	#		$AlignmentB = 'x' x $gapend;
	#	}
	#
	while ($F[$i][$j] > 0) {
		if ($DEBUG) {
			print "$F[$i-1][$j-1] $F[$i-1][$j]\n$F[$i][$j-1] $F[$i][$j]\n";
		}
		if ($F[$i][$j] == $F[$i-1][$j-1] + $S->[$tA[$i-1]][$tB[$j-1]]) {
			$AlignmentA = $A[$i-1] . $AlignmentA;
			$AlignmentB = $B[$j-1] . $AlignmentB;
			$i--; $j--;
		} elsif ($F[$i][$j] == $F[$i-1][$j] + $d) {
			$AlignmentA = $A[$i-1] . $AlignmentA;
			$AlignmentB = '-' . $AlignmentB;
			$i--;
		} else {
			$AlignmentA = '-' . $AlignmentA;
			$AlignmentB = $B[$j-1] . $AlignmentB;
			$j--;
		}
		if ($DEBUG) {
			print "$AlignmentA\n$AlignmentB\n";
			print '-' x 20, "\n";
		}
	}

	if (0) {

	while ($F[$i][$j] > 0) {
		if ($F[$i][$j] == $F[$i-1][$j] + $d) {
			$AlignmentA = $A[$i-1] . $AlignmentA;
			$AlignmentB = '-' . $AlignmentB;
			$i--;
		} elsif ($F[$i][$j] == $F[$i][$j-1] + $d) {
			$AlignmentA = '-' . $AlignmentA;
			$AlignmentB = $B[$j-1] . $AlignmentB;
			$j--;
		} elsif ($F[$i][$j] == $F[$i-1][$j-1] + $S->[$tA[$i-1]][$tB[$j-1]]) {
			$AlignmentA = $A[$i-1] . $AlignmentA;
			$AlignmentB = $B[$j-1] . $AlignmentB;
			$i--; $j--;
		} 
		print "$AlignmentA\n$AlignmentB\n";
		print '-' x 20, "\n";
	}
}
	#while ($i > 0) {
	#	$AlignmentA = '-' . $AlignmentA;
	#	$AlignmentB = $B[$i] . $AlignmentB;
	#	#$AlignmentB = 'x' . $AlignmentB;
	#	$i--;
	#}

	return ($AlignmentA,$AlignmentB, $j);
}

sub generalizedSW {

	#Smith Waterman local alignment
	#align a sequence to a quasi-sequence:
	#
	# so $a = 'AACCTT'
	# and $b = [ {A => .5, C => .5}, ... ]
	# assume that first sequence is the short one
	# b is a list of lists, actually, in order ACGT
	# b = [ [.5, .5, 0, 0], ... ]
	my $a = shift;
	my $b = shift;
	# ltrIdx is hash A => 0, ... T => 3
	# alphabet is list [A C G T]
	my $LtrIdx = shift;
	my $alphabet = shift;
	my $d;
	my $S;
	if (@_ > 0) {
		$d = shift; #linear gap
		$S = shift; #[][] = costs
	} else {
		#$d = -100;
		#$S = [[91,-114,-31,-123],[-114,100,-125,-31],[-31,-125,100,-114],[-123,-31,-114,91]];
		$d = -4;
		#$S = [[4,-8,-8,-8],[-8,4,-8,-8],[-8,-8,4,-8],[-8,-8,-8,4]];
		$S = [[4,-4,-4,-4],[-4,4,-4,-4],[-4,-4,4,-4],[-4,-4,-4,4]];
	}
	my $DEBUG = 0;
	if (@_ > 0) {
		$DEBUG = shift;
	}


	if ($DEBUG) {
		print "Aligning\n$a\n";
		print "[@$_], " foreach (@$b);
		print "\n";
	}
	my @A = split('', $a);
	my @tA = map{$LtrIdx->{$_}} @A;

	my $numB = @$b;
	#my @B = split('', $b);
	#my @tB = map{trans($_)} @B;

	my @B = map {$alphabet->[maxind($_)]} @$b;

	my ($i, $j);
	my $AlignmentA;
	my $AlignmentB;
	my @F;
	my ($largestF, $largesti,$largestj);
	$largestF = 1;

	foreach $i (0..@A) {
		$F[$i][0] = 0;
	}
	foreach $j (0..$numB) {
		$F[0][$j] = 0;
	}


	my ($c1, $c2, $c3, $m);

	foreach $i (1..@A) {
		foreach $j (1..$numB) {
			#$c1 = $F[$i-1][$j-1] + $S->[$tA[$i-1]][$tB[$j-1]];
			#$c1 = F[i-1][j-1] + sum_alphabet pr(b_j-1 = K)*S->[a_i-1][K]
			$c1 = $F[$i-1][$j-1] + score($tA[$i-1], $b->[$j-1], $S);

			$c2 = $F[$i-1][$j] + $d;
			$c3 = $F[$i][$j-1] + $d;
			#print "$i $j = max($c1 $c2 $c3)\n";
			$m = $c1 > $c2 ? ($c1 > $c3 ? $c1 : $c3) : ($c2 > $c3 ? $c2 : $c3);
			$F[$i][$j] = $m > 0 ? $m : 0;
			#$F[$i][$j] = max($choice1,$choice2,$choice3,0);

			#$F[$i][$j] = max($F[$i-1][$j-1] + $S->[$tA[$i-1]][$tB[$j-1]],
			#	$F[$i-1][$j] + $d, $F[$i][$j-1] + $d, 0);

			if ($F[$i][$j] >= $largestF) {
				#print "max at $i $j -- $F[$i][$j]\n";
				$largestF = $F[$i][$j];
				$largesti = $i;
				$largestj = $j;
			}
		}
	}
	if ($DEBUG) {
		foreach my $row (@F) {
			print join("\t", @$row), "\n";
		}
	}


	#now trace back from the largest value in the matrix and build up the alignment

	$AlignmentA = '';
	$AlignmentB = '';
	$i = $largesti;
	$j = $largestj;
	$i = @A;
	$j = $numB;
	#print "Largest F at $i, $j\n";
	#
	#	my $gapend = (length($a) - $largesti);
	#	if ($gapend > 0) {
	#		print "adding $gapend gaps at end\n";
	#		$AlignmentA = 'x' x $gapend;
	#		#foreach my $x ($j .. ($j + $gapend - 1)) {
	#		#	print "$x $B[$x]\n";
	#		#}
	#		#$AlignmentB = lc(join('',@B[$j .. ($j + $gapend-1)]));
	#		$AlignmentB = 'x' x $gapend;
	#	}
	#
	while ($F[$i][$j] > 0) {
		if ($DEBUG) {
			print "$F[$i-1][$j-1] $F[$i-1][$j]\n$F[$i][$j-1] $F[$i][$j]\n";
		}
		if ($F[$i][$j] == $F[$i-1][$j-1] + score($tA[$i-1], $b->[$j-1], $S)) {
			$AlignmentA = $A[$i-1] . $AlignmentA;
			$AlignmentB = $B[$j-1] . $AlignmentB;
			$i--; $j--;
		} elsif ($F[$i][$j] == $F[$i-1][$j] + $d) {
			$AlignmentA = $A[$i-1] . $AlignmentA;
			$AlignmentB = '-' . $AlignmentB;
			$i--;
		} else {
			$AlignmentA = '-' . $AlignmentA;
			$AlignmentB = $B[$j-1] . $AlignmentB;
			$j--;
		}
		if ($DEBUG) {
			print "$AlignmentA\n$AlignmentB\n";
			print '-' x 20, "\n";
		}
	}

	return ($AlignmentA,$AlignmentB, $j);
}

sub blastAlign {
#s2 should be a subsequence of s1
	my $s1 = shift;
	my $s2 = shift;
	system("echo $s1 > /tmp/seq1");
	system("echo $s2 > /tmp/seq2");
	system("bl2seq -i /tmp/seq1 -j /tmp/seq2 -p blastn |grep -A 2 Query:");
}

sub score {
	my $a = shift;
	my $b = shift;
	my $S = shift;
	my $ans = 0;
	foreach my $baseInd (0 .. @$b-1) {
		$ans += $S->[$a][$baseInd] * $b->[$baseInd];
	}
	#print "Score $a, [@$b] = $ans\n";
	return $ans;
}

sub correctGaps {
	my $reads = shift;

	my @AAlist;
	foreach my $rd (@$reads) {
		$rd->{seq} =~ s/N/-/g;
		foreach my $pos ($$rd{st} .. $$rd{end}-1) {
			${$AAlist[$pos]}{substr($$rd{seq},$pos - $$rd{st}, 1)}++;
		}
	}
	my @corInd;
	my @cons;

	foreach my $i (0..@AAlist-1) {
		#my @chars = keys(%{$AAlist[$i]});
		if ((defined ${$AAlist[$i]}{'-'}) and (keys(%{$AAlist[$i]}) == 2)) {
			#foreach (sort{$AAlist[$i]{$b} <=> $AAlist[$i]{$a}}keys(%{$AAlist[$i]})) {
			delete ${$AAlist[$i]}{'-'};
			push(@corInd, $i);
			#my @c = keys(%{$AAlist[$i]});
			my @c = sort {$AAlist[$i]{$b} <=> $AAlist[$i]{$a}} keys(%{$AAlist[$i]});
			push(@cons, $c[0]);
		}
	}

	foreach my $r (@$reads) {
		foreach my $i (0 .. @corInd-1) {	
			#if r->st <= corInd[i]  <= r->end
			#and r->seq . substr(i-start) == '-'
			#correct to cons
			if (($r->{st} <= $corInd[$i]) and ($r->{end} >= $corInd[$i])) {
				if (substr($r->{seq},$corInd[$i] - $r->{st},1) eq '-') {
					substr($r->{seq},$corInd[$i] - $r->{st},1)  = $cons[$i];
				}
			}
		}
	}
}



sub correctGaps2 {
	#take a set of reads (given as two lists, @head and @seq)
	#and correct all gaps for which there is only one choice
	my $head = shift;
	my $seq = shift;
	my $starts;
	my $ends;
	#fix starts and ends
	my $n = @$head;
	my $num = 0;
	foreach my $i (0 .. $n-1) {
		my ($left, $right) = split /\$/, $head->[$i];
		my ($id, $start, undef) = split / /, trim($left);
		my $len = length($seq->[$i]);
		push(@$starts, $start);
		push(@$ends, $start + $len);
	}

	my @AAlist;
	foreach my $rd (0 .. $n-1) {
		$seq->[$rd] =~ s/N/-/g;
		#if ($seq->[$rd] =~ / /) {
		#	print "SPACE: read $rd == .$seq->[$rd].\n";
		#}
		#print "read $rd starts at $starts->[$rd] ends $ends->[$rd]\n";
		#print "length ", length($seq->[$rd]), "\n";
		foreach my $pos ($starts->[$rd] .. $ends->[$rd] - 1) {
			${$AAlist[$pos]}{substr($seq->[$rd],$pos - $starts->[$rd], 1)}++;
		}
	}
	my @corInd;
	my @cons;

	foreach my $i (0..@AAlist-1) {
		#my @chars = keys(%{$AAlist[$i]});
		#print "pos $i, keys = [", join(',', keys(%{$AAlist[$i]})), "]\n";
		if ((defined ${$AAlist[$i]}{'-'}) and (keys(%{$AAlist[$i]}) == 2)) {
			#foreach (sort{$AAlist[$i]{$b} <=> $AAlist[$i]{$a}}keys(%{$AAlist[$i]})) 
			delete ${$AAlist[$i]}{'-'};
			push(@corInd, $i);
			my @c = sort {$AAlist[$i]{$b} <=> $AAlist[$i]{$a}} keys(%{$AAlist[$i]});
			push(@cons, $c[0]);
		}
	}
	#print "correcting at @corInd\n";

	foreach my $r (0 .. $n-1) {
		foreach my $i (0 .. @corInd-1) {	
			#if r->st <= corInd[i]  <= r->end
			#and r->seq . substr(i-start) == '-'
			#correct to cons
			if (($starts->[$r] <= $corInd[$i]) and ($ends->[$r] >= $corInd[$i])) {
				if (substr($seq->[$r],$corInd[$i] - $starts->[$r],1) eq '-') {
					substr($seq->[$r],$corInd[$i] - $starts->[$r],1)  = $cons[$i];
					#print "Correct at $corInd[$i]\n";
					$num++;
				}
			}
		}
	}
	return $num;
}

sub ambigString {
	#input: $a = asdf[123]
	#output  @b = ({a=>1} {s=>1} {d=>1} {f=>1} {1=>1, 2=>1, 3=>1})
	#
	my $a = shift;
	my @b;

	my $i = 0;

	while ($i < length($a)) {
		my $chunk = '';
		if (substr($a,$i,1) ne '[') {
			$chunk = substr($a,$i,1);
		} else {
			# find end of [...] block
			while (substr($a,$i,1) ne ']') {
				$chunk .= substr($a,$i,1);
				$i++;
			}
		}
		$i++;
		my %tmp = map {$_ => 1.0/length($chunk)} split '', $chunk;
		push(@b, \%tmp);
	}
	return \@b;
}


1;
