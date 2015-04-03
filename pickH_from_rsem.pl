#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

my $RSEM       = 'RSEM.isoforms.results';
my $MIN_LENGTH = 300;
my $CLEAN_NAME;
my $MOD_NAME;

parseArgs();
my $assembly = shift @ARGV;

my %grab;            # key = isoform ID with highest coverage, value = 1 (placeholder)
my %pcnt_counter;    # key = gene ID ($gene_ID), value = sum of IsoPct values for that gene (will sum to 100 when all isoforms have been found)
my %genes;           # key = gene ID ($gene_ID), value = array of IsoPct values for that gene where the index corresponds to the isoform number (zero array position is null)
my %isoform_lengths; # key = gene ID ($gene_ID), value = array of IsoPct values for that gene, where the array index corresponds to the isoform number (zero array position is null)

# read in RSEM isoform results
open( RSEM_ISO, $RSEM ) || die "Couldn't open $RSEM: $!\n";
while ( my $line = <RSEM_ISO> ){
  chomp $line;
  $line =~ /^transcript_id/ and next;
  my ( $transcript_id, $gene_id, $length, $effective_length, $expected_count, $tpm, $fpkm, $IsoPct ) = split( /\t/, $line );
  
  my $isoform_num = chop $transcript_id;
  $pcnt_counter{$gene_id} += $IsoPct;
  $genes{$gene_id}->[$isoform_num] = $IsoPct;
  $isoform_lengths{$gene_id}->[$isoform_num] = $length;
}
close RSEM_ISO;

# number of genes with IsoPct = 0
my $zero_coverage = 0;

# loop over genes, determine isoform with highest $IsoPcnt and grab its ID
foreach my $gene_id ( keys %pcnt_counter ){

  # get array index of isoform with the largest IsoPct
  my( $maxIndex, $maxValue ) = getIndexOfMaxValue( @{$genes{$gene_id}} );

  # skip genes (1 isoform/gene at this point) with zero total coverage
  if( $maxValue == 0 ){
    $zero_coverage++;
    # print $gene_id, "\n";
  }else{
    # make sure I've captured all the IsoPct values correctly
    round( $pcnt_counter{$gene_id} ) == 100 or die "IsoPct values sum to $pcnt_counter{$gene_id} for $gene_id: problem with script.\n";
    
    # proceed if isoform meets the minimum length requirement
    if( $isoform_lengths{$gene_id}->[$maxIndex] >= $MIN_LENGTH ){

      # assemble name of the isoform as it appears in the FASTA header
      my $maxIsoformID = $gene_id . "_i" . $maxIndex;

      # store isoform with the highest coverage
      $grab{$maxIsoformID} = 1;

      # print $maxIsoformID, "\t", $maxIndex, "\t", $maxValue, "\t", $isoform_lengths{$gene_id}->[$maxIndex], "\n";
    }
  }
}

# for troubleshooting and checking results
# my $c = 0;
# foreach( keys %grab ){
#   $c++;
# }
# print $c, "\n";
# print $zero_coverage, "\n";

fetchSeqs( $assembly, \%grab );

exit;


############################################SUBROUTINES#############################################
sub getIndexOfMaxValue{
  my @array = @_;

  my $max_value = 0;
  my $max_index = 0;

  # start at $i=1 (isoform 1, 0 position is null)
  for( my $i = 1; $i < @array; $i++ ){
    if($array[$i] > $max_value){
      $max_index = $i;
      $max_value = $array[$i];
    }
  }
  
return( $max_index, $max_value );

}
####################################################################################################
sub fetchSeqs{
  my( $assembly, $grab_ref ) = @_;
  my @fasta_header;
  
  open( FAS, $assembly ) || die "Couldn't open $assembly: $!\n";
  while( my $line = <FAS> ){
    if( $line =~ /^>/ ){
      @fasta_header = split /\s+/, $line;
      $fasta_header[0] =~ s/^>//;
      if( $grab_ref->{$fasta_header[0]} ){
	if( $CLEAN_NAME ){
	  print ">$fasta_header[0]\n";
	}else{
	  print $line;
	}
      }
    }else{
      $grab_ref->{$fasta_header[0]} and print $line;
    }
  }
  close FAS;
}
####################################################################################################
sub round {
  my $number = shift;
  return int($number + .5);
}
####################################################################################################
sub parseArgs{

  my $usage = "\n\nUsage: $0 trinity_assembly.fa [options]

   options
          --rsem        - isoforms output from RSEM (default = 'RSEM.isoforms.results')
          --min_length  - only keep contigs >= than <integer> bp in length (default = 300)
          --clean_name  - clean up the long FASTA header output by Trinity (default = false)
          --modify_name - add suffix to the FASTA header (e.g., species name) to make it more informative (string)
\n\n";

  
  my $result = GetOptions
	(
	 'rsem=s'        => \$RSEM,
	 'min_length=i'  => \$MIN_LENGTH,
	 'clean_name!'   => \$CLEAN_NAME,
	 'modify_name=s' => \$MOD_NAME,
	);
  
  $ARGV[0] or die $usage;
}
####################################################################################################
