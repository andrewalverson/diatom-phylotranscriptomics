#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

# settings for this script
my $SCRIPT_WALL   = 18;    # walltime for running 'rnaseq_clean_filter.pl'
my $SCRIPT_QUEUE  = 'aja'; # queue for Trinity PBS script

# settings for writing Trinity job script within 'rnaseq_clean_filter.pl'
my $TRINITY_QUEUE  = 'mem512GB64core'; # queue for Trinity PBS script
my $TRINITY_WALL   = 12; # walltime for Trinity, which will be written to the Trinity PBS job script
my $MIN_KMER_RNA   = 1;  # for Trinity, min_kmer_cov parameter (rRNA assembly)
my $MIN_KMER_ORG   = 1;  # for Trinity, min_kmer_cov parameter (organelle assembly)
my $MIN_KMER_NUC   = 1;  # for Trinity, min_kmer_cov parameter (nuclear assembly)

# whether to normalize with BBNORM while running 'rnaseq_clean_filter.pl'
my $NORMALIZE;

# settings for running Trimmomatic within 'rnaseq_clean_filter.pl'
my $SLIDINGWINDOW = '4:5'; # Trimmomatic SLIDINGWINDOW paramter - number of 5' bp to trim
my $TRIM_PHRED    =  2;    # Trimmomatic LEADING and TRAILING partameters - bases below this Phred score will be trimmed from 5' and 3' ends of reads
my $MIN_LEN       = 25;    # Trimmomatic MINLEN parameter - remove reads shorter than this


parseArgs();

# get RNA ID
my $rnaID = shift @ARGV;
my $suffix = "";

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

# determine whether this an Alverson lab RNA ID
if( $rnaID =~ /(\d+)([a-zA-Z])?/ ){
  $rnaID  = "R$1";
  $2 and $suffix = $2;
}

# print "rnaID:  $rnaID\n";
# print "suffix: $suffix\n";

my $outfile = "run_rnaseq2clean_filter_$rnaID$suffix.pbs";

open( PBS, '>', $outfile ) || die "Can't open $outfile: $!\n";

# print PBS lines
print PBS "#PBS -N rnaseq2clean_filter_$rnaID$suffix\n";
print PBS "#PBS -q $SCRIPT_QUEUE\n";
print PBS '#PBS -j oe', "\n";
print PBS '#PBS -m abe', "\n";
print PBS '#PBS -M aja@uark.edu ', "\n";
print PBS "#PBS -o rnaseq2clean_filter_$rnaID$suffix\." . '$PBS_JOBID' . "\n";
print PBS "#PBS -l nodes=1:ppn=$ppn\n";
print PBS '#PBS -l walltime=' . $SCRIPT_WALL . ':00:00', "\n\n";
print PBS 'cd $PBS_O_WORKDIR', "\n\n";

print PBS 'module purge', "\n";
print PBS 'module load bowtie', "\n";
print PBS 'module load samtools/0.1.19', "\n";
print PBS "module load bbmap\n\n";

print PBS "/home/aja/local/src/scripts/rnaseq_clean_filter.pl $rnaID --wall=$TRINITY_WALL --queue=$TRINITY_QUEUE $NORMALIZE \
          --min_kmer_rna=$MIN_KMER_RNA --min_kmer_org=$MIN_KMER_ORG --min_kmer_nuc=$MIN_KMER_NUC
          --window=$SLIDINGWINDOW --trim_phred=$TRIM_PHRED --min_len=$MIN_LEN\n";

close PBS;

exit;

############################################SUBROUTINES#############################################
sub parseArgs{

  my $usage = "\nUsage: $0 RNA_ID [options]


   options for running the 'rnaseq_clean_filter.pl' pipeline
          --script_wall   - walltime for running 'rnaseq_clean_filter.pl' (integer, default = 12 hrs)
          --script_queue  - queue for running 'rnaseq_clean_filter.pl' ('aja' [default], 'med12core', 'med16core')


   option to normalize within 'rnaseq_clean_filter.pl' pipeline
          --normalize     - perform kmer normalization with BBNorm (default: false)


   options for running Trimmomatic within 'rnaseq_clean_filter.pl' pipeline
          --window     - Trimmomatic SLIDINGWINDOW parameter specified as <windowSize><requiredQuality> (default: '4:5')
          --trim_phred - Trimmomatic LEADING and TRAILING partameters; bases below this minimum Phred quality
                             will be trimmed from 5' and 3' ends (default: 2)
          --min_len    - filter reads shorter than this length threshold (default: 25)


   options for writing Trinity PBS job script within 'rnaseq_clean_filter.pl' pipeline
          --trinity_wall  - walltime for Trinity PBS script (integer, default = 12 hrs)
          --trinity_queue - queue for Trinity PBS script ('mem256GB48core', 'mem512GB64core' [default], 'mem768GB32core', 'random')
          --min_kmer_rna   - Trinity min_kmer_cov parameter for rRNA assembly (default: 1)
          --min_kmer_org   - Trinity min_kmer_cov parameter for organelle assembly (default: 1)
          --min_kmer_nuc   - Trinity min_kmer_cov parameter for nuclear assembly (default: 1)

\n";

  my $result = GetOptions
	(
	 'script_wall=i'   => \$SCRIPT_WALL,
	 'script_queue=s'  => \$SCRIPT_QUEUE,
	 'normalize!'      => \$NORMALIZE,
	 'window=s'        => \$SLIDINGWINDOW,
	 'trim_phred=s'    => \$TRIM_PHRED,
	 'min_len=s'       => \$MIN_LEN,
	 'trinity_wall=i'  => \$TRINITY_WALL,
	 'trinity_queue=s' => \$TRINITY_QUEUE,
	 'min_kmer_rna=i'  => \$MIN_KMER_RNA,
	 'min_kmer_org=i'  => \$MIN_KMER_ORG,
	 'min_kmer_nuc=i'  => \$MIN_KMER_NUC,
	);
  
  $ARGV[0] or die $usage;
}
####################################################################################################
