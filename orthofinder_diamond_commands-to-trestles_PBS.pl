#!/usr/bin/perl

use warnings;
use strict;

my $usage = "$0 diamond-commands.list outdir\n";
$ARGV[1] or die $usage;

# set output directory
my $outdir = chomp($ARGV[1]);

# wall time
my $wall = 6;

# print top of PBS script
print '#PBS -N ' . $ARGV[0] . "\n";
print '#PBS -q q06h32c', "\n";
print '#PBS -j oe', "\n";
print '#PBS -o ' . $ARGV[0] . '.$PBS_JOBID' . "\n";
print '#PBS -l nodes=1:ppn=32', "\n";
print '#PBS -l walltime=' . $wall . ':00:00', "\n\n";

print 'cd /scratch/$PBS_JOBID', "\n\n";

print "# run the Diamond searches\n\n";
  
# read and print Diamond commands
open( IN, $ARGV[0] ) || die "Couldn't open $ARGV[0]: $!\n";
while(<IN>){
  print "/share/apps/bioinformatics/diamond/0.9.21.122/";
  print;
}
close IN;

print "\n";

print "# move the output files\n";
print "mv * $outdir\n\n";

exit;
