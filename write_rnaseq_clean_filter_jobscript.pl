#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

# settings for this script
my $SCRIPT_WALL   = 18;    # walltime for running 'rnaseq_clean_filter.pl'
my $SCRIPT_QUEUE  = 'aja'; # queue for Trinity PBS script

# settings for writing Trinity job script within 'rnaseq_clean_filter.pl'
my $TRINITY_QUEUE  = 'mem512GB64core'; # queue for Trinity PBS script
my $TRINITY_WALL   = 120; # walltime for Trinity, which will be written to the Trinity PBS job script
my $MIN_KMER_RNA   = 1;   # for Trinity, --min_kmer_cov parameter (rRNA assembly)
my $MIN_KMER_ORG   = 1;   # for Trinity, --min_kmer_cov parameter (organelle assembly)
my $MIN_KMER_NUC   = 1;   # for Trinity, --min_kmer_cov parameter (nuclear assembly)
my $STRANDED;             # Trinity, --SS_lib_type parameter to specify a stranded library
my $NORMALIZE;            # Trinity, --no_normalize_reads parameter (Trinity normalizes reads by default after 21 Sept 2016)
my $REFERENCE      = 'phaeo';  # Transrate 'reference' parameter, whether to use Phaeo or Thaps as the reference proteome

# settings for running Trimmomatic within 'rnaseq_clean_filter.pl'
my $SLIDINGWINDOW = '4:2'; # Trimmomatic SLIDINGWINDOW paramter - number of 5' bp to trim
my $MIN_LEN       = 30;    # Trimmomatic MINLEN parameter - remove reads shorter than this
my $PHRED_OUT     = 64;    # Trimmomatic TOPHRED parameter

parseArgs();

# get RNA ID
my $ID     = shift @ARGV;
my $type   = "";
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
}elsif( $SCRIPT_QUEUE eq 'mem96GB12core' ){
  $ppn = 12;
}elsif( $SCRIPT_QUEUE eq 'mem512GB64core' ){
  $ppn = 32;
}elsif( $SCRIPT_QUEUE eq 'mem768GB32core' ){
  $ppn = 32;
}else{
  die "Queue \'$SCRIPT_QUEUE\' not recognized\n";
}

# determine whether this an Alverson lab RNA or Library ID
if( $ID =~ /\A([RL])(\d+)([a-zA-Z])?/i ){
  $type = $1;
  $ID   = $2;
  $3 and $suffix = $3;
}

my $outfile = "$type$ID$suffix\_rnaseq2clean_filter.pbs";

# parse stranded option
if( $STRANDED ){
  $STRANDED = '--stranded';
}else{
  $STRANDED = '';
}

# parse normalize option
if( $NORMALIZE ){
  $NORMALIZE = '--normalize ';
}else{
  $NORMALIZE = '';
}

open( PBS, '>', $outfile ) || die "Can't open $outfile: $!\n";

# print PBS lines
print PBS "#PBS -N rnaseq2clean_filter_$type$ID$suffix\n";
print PBS "#PBS -q $SCRIPT_QUEUE\n";
print PBS '#PBS -j oe', "\n";
print PBS '#PBS -m abe', "\n";
print PBS '#PBS -M aja@uark.edu ', "\n";
print PBS "#PBS -o rnaseq2clean_filter_$type$ID$suffix\." . '$PBS_JOBID' . "\n";
print PBS "#PBS -l nodes=1:ppn=$ppn\n";
print PBS '#PBS -l walltime=' . $SCRIPT_WALL . ':00:00', "\n\n";
print PBS 'cd $PBS_O_WORKDIR', "\n\n";

print PBS 'module purge', "\n";
print PBS 'module load bowtie', "\n";
print PBS 'module load samtools/0.1.19', "\n";
print PBS "module load bbmap\n\n";

print PBS "/home/aja/local/src/scripts/rnaseq/rnaseq_clean_filter.pl $type$ID$suffix --wall=$TRINITY_WALL --trinity_queue=$TRINITY_QUEUE ", '\\', "\n",
          "          --min_kmer_rna=$MIN_KMER_RNA --min_kmer_org=$MIN_KMER_ORG --min_kmer_nuc=$MIN_KMER_NUC $STRANDED $NORMALIZE", ' \\', "\n",
          "          --window=$SLIDINGWINDOW --min_len=$MIN_LEN --reference=$REFERENCE --phred_out=$PHRED_OUT\n";

close PBS;

exit;

############################################SUBROUTINES#############################################
sub parseArgs{

  my $usage = "\nUsage: $0 RNA/Library_ID [options]


   options for running the 'rnaseq_clean_filter.pl' pipeline
          --script_wall   - walltime for running 'rnaseq_clean_filter.pl' (integer, default = 12 hrs)
          --script_queue  - queue for running 'rnaseq_clean_filter.pl' ('aja' [default], 'med12core', 'med16core', 'mem96GB12core', 'mem512GB64core', 'mem768GB32core')


   options for running Trimmomatic within 'rnaseq_clean_filter.pl' pipeline
          --window     - Trimmomatic SLIDINGWINDOW parameter specified as <windowSize><requiredQuality> (default: '4:2')
          --min_len    - filter reads shorter than this length threshold (default: 30)
          --phred_out  - output Phred scale (default: 64)

   options for writing Trinity PBS job script within 'rnaseq_clean_filter.pl' pipeline
          --trinity_wall  - walltime for Trinity PBS script (integer, default = 120 hrs)
          --trinity_queue - queue for Trinity PBS script ('mem256GB48core', 'mem512GB64core' [default], 'mem768GB32core', 'random')
          --min_kmer_rna   - Trinity --min_kmer_cov parameter for rRNA assembly (default: 1)
          --min_kmer_org   - Trinity --min_kmer_cov parameter for organelle assembly (default: 1)
          --min_kmer_nuc   - Trinity --min_kmer_cov parameter for nuclear assembly (default: 1)
          --stranded       - Trinity --SS_lib_type parameter, set to RF when true (default: false [not stranded])
          --normalize      - Trinity --no_normalize_reads parameter, Trinity normalizes by default (boolean, default: do not normalize)


   Transrate options
          --reference  - specify whether to use Phaeodactylum (phaeo) or T. pseudonana (thaps)
                             as the reference proteome (default: phaeo)


\n";

  my $result = GetOptions
	(
	 'script_wall=i'   => \$SCRIPT_WALL,
	 'script_queue=s'  => \$SCRIPT_QUEUE,

	 'window=s'        => \$SLIDINGWINDOW,
	 'min_len=s'       => \$MIN_LEN,
	 'phred_out=s'     => \$PHRED_OUT,
	 
	 'trinity_wall=i'  => \$TRINITY_WALL,
	 'trinity_queue=s' => \$TRINITY_QUEUE,
	 'min_kmer_rna=i'  => \$MIN_KMER_RNA,
	 'min_kmer_org=i'  => \$MIN_KMER_ORG,
	 'min_kmer_nuc=i'  => \$MIN_KMER_NUC,
	 'stranded!'       => \$STRANDED,
	 'normalize!'      => \$NORMALIZE,

	 'reference=s'     => \$REFERENCE,

	);

  $ARGV[0] or die $usage;
}
####################################################################################################
