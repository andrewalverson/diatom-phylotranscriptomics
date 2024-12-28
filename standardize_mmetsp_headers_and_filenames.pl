#!/usr/bin/perl

# standardize FASTA headers and filenames for MMETSP Transdecoder files

use strict;
use warnings;
use Getopt::Long;

# read command line options
parseArgs();
my $infile = $ARGV[0];
my $source = $ARGV[1];
my($taxon, $strain, $outfile, $header, @infile_name_parts);

if($source eq 'mmetsp'){
  # assemble name of output file (MMETSP)
  @infile_name_parts = split( /_nuclear\.good_transcripts\.short_name\.filtered_transcripts\.largest_cluster_transcripts\.fa\.transdecoder\./, $infile );
}else{
  # assemble name of output file (NIES-3581)
  @infile_name_parts = split( /_nuclear\.largest_cluster_transcripts\.fa\.transdecoder\./, $infile );
}

my $file_suffix = pop @infile_name_parts;
if($infile_name_parts[0] =~ /\./){
  ($taxon, $strain) = split /\./, $infile_name_parts[0];
  $taxon =~ s/Thalassiosira_nitzschioides/Thalassionema_nitzschioides/;
  $taxon =~ s/Skeletonema_grathae/Skeletonema_grethae/;
  $outfile = $taxon . "_" . $strain . "." . $file_suffix;
  $header  = $taxon . "_" . $strain;
}else{
  $taxon = shift @infile_name_parts;
  $taxon =~ s/Thalassiosira_nitzschioides/Thalassionema_nitzschioides/;
  $taxon =~ s/Skeletonema_grathae/Skeletonema_grethae/;
  $outfile = $taxon . "." . $file_suffix;
  $header  = $taxon;
}

open( OUT, ">$outfile") || die "Couldn't open $outfile: $!\n";

# open and parse FASTA file
open( FASTA, "$infile" ) || die "Couldn't open $infile: $!\n";
while( my $line = <FASTA>){
  my $trinity_contig;
  if( $line =~ /^>/ ){ # FASTA header
    # transdecoder output - need to include the "m.\d+" translation ID or there will be redundant FASTA headers in the output
    if( $line =~ /(TRINITY_DN\d+_c\d+_g\d+_i\d+)\|(m.\d+)/ ){
      $trinity_contig = "$1\_$2";
    }elsif( $line =~ /(TRINITY_DN\d+_c\d+_g\d+_i\d+(\.p\d+)*)/ ){
      $trinity_contig = "$1";
    }else{
      die "Did not parse Trinity contig id\n$line\n";
    }

    # print FASTA header
    print OUT ">", $header, "_", $trinity_contig, "\n";
  }else{ # FASTA sequence
    print OUT $line;
  }
}
close FASTA;
close OUT;

############################################SUBROUTINES#############################################
sub parseArgs{

  my $usage = "\nUsage: $0 file.fasta source

   Required: name of FASTA file <file.fasta>
   Required: source ('mmetsp' or 'nies')

   Options: none

";

  my $result = GetOptions(
			 );
  
  $ARGV[1] or die $usage;

}
####################################################################################################
