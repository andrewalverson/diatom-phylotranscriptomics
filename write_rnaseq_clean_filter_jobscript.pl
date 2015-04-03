#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

our $SCRIPT_WALL   = 18;    # walltime for running 'rnaseq_clean_filter.pl'
our $SCRIPT_QUEUE  = 'aja'; # queue for Trinity PBS script
our $TRINITY_WALL  = 12;    # walltime for Trinity, which will be written to the Trinity PBS job script
our $TRINITY_QUEUE = 'mem512GB64core'; # queue for Trinity PBS script
my $MIN_KMER_RNA   = 1;     # for Trinity, min_kmer_cov parameter (rRNA assembly)
my $MIN_KMER_ORG   = 1;     # for Trinity, min_kmer_cov parameter (organelle assembly)
my $MIN_KMER_NUC   = 2;     # for Trinity, min_kmer_cov parameter (nuclear assembly)
our $NORMALIZE;

parseArgs();

# of processors in queue
my $ppn;

# get queue information
if( $SCRIPT_QUEUE eq 'aja' ){
  $ppn = 16;
}elsif( $SCRIPT_QUEUE eq 'med12core' ){
  $ppn = 12;
}elsif( $SCRIPT_QUEUE eq 'med16core' ){
  $ppn = 16;
}else{
  die "Queue \'$SCRIPT_QUEUE\' not recognized\n";
}

# normalize or not
if( $NORMALIZE ){
  $NORMALIZE = "--normalize";
}else{
  $NORMALIZE = "";
}

# get RNA ID
my $rnaID = shift @ARGV;

# if this is a re-sequenced RNA, separate numbers and letters
my $suffix = "";
if( $rnaID =~ /(\d+)([a-zA-Z])/ ){
  $rnaID  = $1;
  $suffix = $2;
}

# print "rnaID:  $rnaID\n";
# print "suffix: $suffix\n";

my $outfile = "run_rnaseq2clean_filter_R$rnaID$suffix.pbs";

open( PBS, '>', $outfile ) || die "Can't open $outfile: $!\n";

# print PBS lines
print PBS "#PBS -N rnaseq2clean_filter_R$rnaID$suffix\n";
print PBS "#PBS -q $SCRIPT_QUEUE\n";
print PBS '#PBS -j oe', "\n";
print PBS '#PBS -m abe', "\n";
print PBS '#PBS -M aja@uark.edu ', "\n";
print PBS "#PBS -o rnaseq2clean_filter_R$rnaID$suffix\." . '$PBS_JOBID' . "\n";
print PBS "#PBS -l nodes=1:ppn=$ppn\n";
print PBS '#PBS -l walltime=' . $SCRIPT_WALL . ':00:00', "\n\n";
print PBS 'cd $PBS_O_WORKDIR', "\n\n";

print PBS 'module purge', "\n";
print PBS 'module load bowtie', "\n";
print PBS 'module load samtools/0.1.19', "\n";
print PBS "module load bbmap\n\n";

print PBS "/home/aja/local/src/scripts/rnaseq_clean_filter.pl $rnaID --wall=$TRINITY_WALL --queue=$TRINITY_QUEUE $NORMALIZE --min_kmer_rna=$MIN_KMER_RNA --min_kmer_org=$MIN_KMER_ORG --min_kmer_nuc=$MIN_KMER_NUC\n";

close PBS;

exit;

############################################SUBROUTINES#############################################
sub parseArgs{

  my $usage = "\nUsage: $0 RNA_ID [options]

   options for running this script
          --script_wall   - walltime for running 'rnaseq_clean_filter.pl' (integer, default = 12 hrs)
          --script_queue  - queue for running 'rnaseq_clean_filter.pl' ('aja' [default], 'med12core', 'med16core')
          --normalize     - perform kmer normalization with BBNorm (default: false)

   options for Trinity run
          --trinity_wall  - walltime for Trinity PBS script (integer, default = 12 hrs)
          --trinity_queue - queue for Trinity PBS script ('mem256GB48core', 'mem512GB64core' [default], 'mem768GB32core', 'random')
          --min_kmer_rna   - Trinity min_kmer_cov parameter for rRNA assembly (default: 1)
          --min_kmer_org   - Trinity min_kmer_cov parameter for organelle assembly (default: 1)
          --min_kmer_nuc   - Trinity min_kmer_cov parameter for nuclear assembly (default: 2)

\n";

  my $result = GetOptions
	(
	 'script_wall=i'   => \$SCRIPT_WALL,
	 'script_queue=s'  => \$SCRIPT_QUEUE,
	 'trinity_wall=i'  => \$TRINITY_WALL,
	 'trinity_queue=s' => \$TRINITY_QUEUE,
	 'normalize!'      => \$NORMALIZE,
	 'min_kmer_rna=i'  => \$MIN_KMER_RNA,
	 'min_kmer_org=i'  => \$MIN_KMER_ORG,
	 'min_kmer_nuc=i'  => \$MIN_KMER_NUC,
	);
  
  $ARGV[0] or die $usage;
}
####################################################################################################
