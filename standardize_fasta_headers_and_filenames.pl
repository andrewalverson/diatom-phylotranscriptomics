#!/usr/bin/perl

# standardize FASTA headers among CroCo output and Transdecoder cds, mRNA, and pep files

use strict;
use warnings;
use Getopt::Long;
use DBI;

# read command line options
parseArgs();
my $infile = shift @ARGV;
my ( $sample_type, $sample_num, $taxon_id );

# parse file name to get RNA or Library ID
if($infile =~ /^nitz2144/){
  $taxon_id = 'Nitzschia_sp_CCMP2144';
}elsif($infile =~ /^lepto/){
  $taxon_id = 'Leptocylindrus_danicus_ECT3929';
}elsif($infile =~ /^bolido/){
  $taxon_id = 'Bolidomonas_pacifica_CCMP1866';
}elsif($infile =~ /^hemi/){
  $taxon_id = 'Hemiaulus_sinensis_24i10-1A';
}elsif($infile =~  /^([LR])(\d+)/){
  $sample_type = $1;
  $sample_num  = $2;
  
  # assemble taxon ID for FASTA header
  $taxon_id = sql_fetch( $sample_type, $sample_num );
}

# assemble name of output file
my @infile_name_parts = split( /\./, $infile );
my $file_suffix = pop @infile_name_parts;
my $outfile = $taxon_id . "." . $file_suffix;

my $trinity_contig;

open( OUT, ">$outfile") || die "Couldn't open $outfile: $!\n";

# open and parse FASTA file
open( FASTA, "$infile" ) || die "Couldn't open $infile: $!\n";
while( my $line = <FASTA>){
  if( $line =~ /^>/ ){ # FASTA header
    # transdecoder output - need to include the "m.\d+" translation ID or there will be redundant FASTA headers in the output
    if( $line =~ /(TRINITY_DN\d+_c\d+_g\d+_i\d+)\|(m.\d+)/ ){
      $trinity_contig = "$1\_$2";
    # CroCo output
    }elsif( $line =~ /(TRINITY_DN\d+_c\d+_g\d+_i\d+)/ ){
      $trinity_contig = $1;
    }else{
      die "Did not parse Trinity contig id\n";
    }

    # print FASTA header
    print OUT ">", "$taxon_id", "_", $trinity_contig, "\n";
  
  }else{ # FASTA sequence
    print OUT $line;
  }
}
close FASTA;
close OUT;

############################################SUBROUTINES#############################################
sub sql_fetch{

  my ( $sample_type, $sample_num ) = @_;
  my $HOST      = 'hpceq';          # the MySQL server on UARK RAZOR cluster
  my $DB        = 'transcriptome';  # MySQL datbase
  my $table;

  # for connecting to the MySQL server
  my $ds     = "DBI:mysql:$DB:$HOST";
  my $user   = 'aja';
  my $passwd = '123aja321';

  my $dbh = DBI->connect("DBI:mysql:database=transcriptome;host=hpceq",
			 "aja", "$passwd",
			{ 'AutoCommit' => 0, 'RaiseError' => 1 }) or die $DBI::errstr;

  my $sth;  # holds the SQL prepare statement

  my( $query_columun, $type );
  if( $sample_type eq 'R' ){
    $table = 'RNA_Book';
    $query_columun = 'RNA_ID';
  }elsif( $sample_type eq 'L' ){
    $table = 'Library_Book';
    $query_columun = 'Library_ID';
  }
  
  # prepare and execute the SQL SELECT statement for this RNA_ID
  $sth = $dbh->prepare( "SELECT * FROM $table WHERE $query_columun = $sample_num" );
  $sth->execute();
  
  # fetch the SQL entry for this RNA_ID into a hash
  my $hash = $sth->fetchrow_hashref();
  
  # finish the query
  $sth->finish;

  # build sequence identifier (e.g., 'Rhizosolenia_setigera_ECT2AJA-026_R20b' )
  my $taxon_id;

  if( $hash->{Genus} ){
    $taxon_id = ucfirst $hash->{Genus};
    if( $hash->{Species} ){
      $hash->{Species} =~ s/\s+/_/g;
      $taxon_id .= "_" . lc $hash->{Species};
    }
    $taxon_id .= "_" .  $hash->{Culture_ID};
    $taxon_id .= "_" .  $sample_type . $sample_num;
  }else{
    $taxon_id = "Unidentified_" . $hash->{Culture_ID} . "_" . $sample_type . $sample_num;
  }
  
  # remove '?' and '.' charcters from taxon ID
  $taxon_id =~ s/\?//g;
  $taxon_id =~ s/\.//g;
  
  # disconnect from the MySQL database
  $dbh->disconnect;

  return $taxon_id;
  
}
####################################################################################################
sub parseArgs{

  my $usage = "\nUsage: $0 file.fasta [options]

   Required: name of FASTA file <file.fasta>
   Options: none

";

  my $result = GetOptions(
			 );
  
  $ARGV[0] or die $usage;

}
####################################################################################################
