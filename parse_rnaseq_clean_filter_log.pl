#!/usr/bin/perl

# parse log file created by rnaseq_clean_filter.pl

use strict;
use warnings;
use DBI;

$ARGV[0] or die "usage: $0 assembly-directory\n";

my $dir = shift;
my $log_file = `find $dir -name \"*_clean_filter.log\" -print`;

open( LOG, "$log_file" ) || die "Couldn't open $log_file: $!\n";

my ( $sample_ID, $type, $table );
my ( $genus, $species );
my $num_trimmomatic_dropped = 0;
my $num_raw_reads = 0;
my $num_corrected_bases = 0;
my $num_nuclear_pairs = 0;
my $num_nuclear_merged = 0;
my $bbmerge_setting;
my $percent_diatom_LSU = 0;
my $percent_diatom_SSU = 0;
my $percent_tol_SSU = 0;

################ PARSE LOG FILE ################
while( <LOG>) {
  chomp;
  
  if( /Sample ID: \A([LR])(\d+)/ ){
    $type      = $1;
    $sample_ID = $2;

    $type eq 'L' and $table = 'Library_Book';
    $type eq 'R' and $table = 'RNA_Book';
  }

  ( $genus, $species) = sql_fetch( $sample_ID, $table );

  if( /BEGIN NUCLEAR BBMERGE/ ){
    while( <LOG> ){
      /^END/ and last;
      /(v*strict)=t/    and $bbmerge_setting    = $1;
      /Pairs:\s+(\d+)/  and $num_nuclear_pairs  = $1;
      /Joined:\s+(\d+)/ and $num_nuclear_merged = $1;
    }
  }elsif( /BEGIN BOWTIE2 LSU RRNA/ ){
    while( <LOG> ){
      /^END/ and last;
      /([\d\.]+)% overall alignment rate/ and $percent_diatom_LSU = $1
    }
  }elsif( /BEGIN BOWTIE2 SSU_TOL RRNA/ ){
    while( <LOG> ){
      /^END/ and last;
      /([\d\.]+)% overall alignment rate/ and $percent_tol_SSU = $1
    }
  }elsif( /BEGIN BOWTIE2 SSU_DIATOM RRNA/ ){
    while( <LOG> ){
      /^END/ and last;
      /([\d\.]+)% overall alignment rate/ and $percent_diatom_SSU = $1
    }
  }elsif( /BEGIN TRIMMOMATIC/ ){
    while( <LOG> ){
      /^END/ and last;
      /Dropped: (\d+)/ and $num_trimmomatic_dropped = $1;
    }
  }elsif( /BEGIN RCORRECTOR/ ){
    while( <LOG> ){
      /^END/ and last;
      /^Processed (\d+) reads/ and $num_raw_reads = $1;
      /Corrected (\d+) bases./ and $num_corrected_bases = $1;
    }
  }
}

close LOG;


################ PARSE BUSCO SHORT SUMMARY ################

my $complete_buscos = 0;
my $fragmented_buscos = 0;
my $total_searched_buscos = 0;

open( BUSCO, "$dir/run_$sample_ID\_BUSCO/short_summary_$sample_ID\_BUSCO" ) || die "Couldn't open $dir/run_$sample_ID\_BUSCO/short_summary_$sample_ID\_BUSCO: $!\n";
while( <BUSCO> ){

  /(\d+)\s+Complete BUSCOs/             and $complete_buscos = $1;
  /(\d+)\s+Fragmented BUSCOs/           and $fragmented_buscos = $1;
  /(\d+)\s+Total BUSCO groups searched/ and $total_searched_buscos = $1;

}
close BUSCO;

################ PARSE TRANSRATE SUMMARY ################

my $num_contigs      = 0;
my $num_bases        = 0;
my $contig_N50       = 0;
my $assembly_gc      = 0;
my $p_good_mapping   = 0;
my $p_refs_with_crbb = 0;
my $transrate_score  = 0;

open( TRANSRATE, "$dir/transrate.csv" ) || die "Couldn't open $dir/transrate.csv: $!\n";
while( my $line = <TRANSRATE> ){

  $line =~ /^assembly/ and next;
  chomp $line;
  
  my @cols = split ',', $line;
  
  $num_contigs      = $cols[1];
  $num_bases        = $cols[4];
  $contig_N50       = $cols[13];
  $assembly_gc      = $cols[16];
  $p_good_mapping   = $cols[27];
  $p_refs_with_crbb = $cols[45];
  $transrate_score  = $cols[57];
  
}

close TRANSRATE;


##################### PRINT RESULTS #####################

#print join "\t", "sample_ID", "num_reads (F+R)", "num_corrected_bases", "num_dropped_reads", "percent_diatom_LSU", "percent_diatom_SSU", "percent_TOL_SSU", "num_nuclear_read_pairs", "bbmerge_setting", "percent_merged", "Assembled contigs", "Assembled bases", "Assembly GC", "Contig N50", "Percent good mappings", "Percent Phaeo CRBB", "Transrate score", "BUSCOs (complete+fragmented)", "BUSCOs (percent)", "\n";

print  "$sample_ID\t$genus\t$species\t$num_raw_reads\t$num_corrected_bases\t$num_trimmomatic_dropped\t";
print  $percent_diatom_LSU, "\t", $percent_diatom_SSU, "\t", $percent_tol_SSU, "\t";
print  $num_nuclear_pairs, "\t", $bbmerge_setting, "\t";
printf "%.2f\t", ($num_nuclear_merged/$num_nuclear_pairs)*100;
printf "$num_contigs\t$num_bases\t%.2f\t$contig_N50\t", $assembly_gc;
printf "%.2f\t%.2f\t%.2f\t", $p_good_mapping*100, $p_refs_with_crbb*100, $transrate_score;
print  $complete_buscos + $fragmented_buscos, "\t";
printf "%.2f", (($complete_buscos + $fragmented_buscos)/$total_searched_buscos)*100;
print  "\n";

exit;


############################################SUBROUTINES#############################################
sub sql_fetch{

  my ($id_num, $table) = @_;
  my $HOST      = 'hpceq';          # the MySQL server on UARK RAZOR cluster
  my $DB        = 'transcriptome';  # MySQL datbase
  my $TABLE     = 'RNA_Book';       # table within MySQL database

  # for connecting to the MySQL server
  my $ds     = "DBI:mysql:$DB:$HOST";
  my $user   = 'aja';
  my $passwd = '123aja321';

# my $HOST      = 'hpceq';          # the MySQL server on UARK RAZOR cluster
# my $DB        = 'transcriptome';  # MySQL datbase
# my $TABLE     = 'RNA_Book';       # table within MySQL database
  
  my $dbh = DBI->connect("DBI:mysql:database=transcriptome;host=hpceq",
			 "aja", "$passwd",
			{ 'AutoCommit' => 0, 'RaiseError' => 1 }) or die $DBI::errstr;

  my $sth;  # holds the SQL prepare statement

  my( $query_columun, $type );
  if( $table eq 'RNA_Book' ){
    $query_columun = 'RNA_ID';
    $type = 'R';
  }elsif( $table eq 'Library_Book' ){
    $query_columun = 'Library_ID';
    $type = 'L';
  }
  
  # prepare and execute the SQL SELECT statement for this RNA_ID
  $sth = $dbh->prepare( "SELECT * FROM $table WHERE $query_columun = $id_num" );
  $sth->execute();
  
  # fetch the SQL entry for this RNA_ID into a hash
  my $hash = $sth->fetchrow_hashref();
  
  # finish the query
  $sth->finish;

  # build directory name (e.g., 'R20b_rhizosolenia_setigera_ECT2AJA-026' )
  my $dir_id;

  if( $hash->{Genus} ){
    $dir_id = $type . $id_num . "_" . lcfirst $hash->{Genus};
    if( $hash->{Species} ){
      $hash->{Species} =~ s/\s+/_/g;
      $dir_id .= "_" . lc $hash->{Species};
    }
    $dir_id .= "_" .  $hash->{Culture_ID};
  }else{
    $dir_id = $type . $id_num . "_" . $hash->{Culture_ID};
  }
  
  # remove '?' charcters from directory ID
  $dir_id =~ s/\?//g;
  $dir_id =~ s/\.//g;
  
  # build name for log file name
  my $logname = "$type$id_num\_clean_filter.log";

  # disconnect from the MySQL database
  $dbh->disconnect;

  # print $dir_id, "\n" and exit;
  
  return( $hash->{Genus}, $hash->{Species} );
  
}
####################################################################################################
