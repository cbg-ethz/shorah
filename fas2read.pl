#!/usr/bin/perl -w
#translate a fasta file to read format

use lib "perllib";
use strict;
use Getopt::Long;
use DNA;
use ReadUtil;

my ($HELP, $VERBOSE, $FASTAFILE);

GetOptions(
	   "help" => \$HELP,
	   "fasta=s" => \$FASTAFILE,
	  );
help() if (($HELP) or (! defined $FASTAFILE));
my $outfile = $FASTAFILE;
$outfile =~ s/fas$/read/;
print "Reading $FASTAFILE\nOutput to $outfile\n";
open OUT, ">$outfile";
my $reads = readFastaFile($FASTAFILE);

foreach my $rd (@$reads) {
  print OUT "$rd->{st} $rd->{seq}\n";
}
close OUT;

sub help {
  print <<END;
usage: $0 -f file.fas [ -h -v -o]
	reads from file.fas, outputs to file.read
END
  exit;
}
