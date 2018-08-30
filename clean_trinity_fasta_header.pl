#!/usr/bin/perl

use warnings;
use strict;

# script usage
my $usage = "\nUsage: $0 <trinity_assembly.fa>\n\n";

# make sure a FASTA file has been specified
$ARGV[0] or die $usage;

my $assembly = shift @ARGV;

open( FAS, $assembly ) || die "Couldn't open $assembly: $!\n";
while( my $line = <FAS> ){
  if( $line =~ /^>/ ){
    my @fasta_header = split /\s+/, $line;
    print $fasta_header[0], "\n";
  }else{
    print $line;
  }
}
close FAS;

exit;
