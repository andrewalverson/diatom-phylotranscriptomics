#!/usr/bin/perl

# parse log file created by rnaseq_clean_filter.pl

use strict;
use warnings;

die "usage: $0 assembly-directory\n" unless $#ARGV == 0;

my $dir = shift;
my $log_file = `find $dir -name \"*_clean_filter.log\" -print`;

open( LOG, "$log_file" ) || die "Couldn't open $log_file: $!\n";

my $sample_ID;
my $num_trimmomatic_dropped;
my $num_raw_reads;
my $num_corrected_bases;
my $num_organelle_pairs;
my $num_rrna_pairs;
my $num_nuclear_pairs;
my $num_nuclear_merged;
my $bbmerge_setting;


################ PARSE LOG FILE ################
while( <LOG>) {
  # print;
  chomp;
  /Sample ID: (\S+)/ and $sample_ID = $1;

  if( /BEGIN NUCLEAR BBMERGE/ ){
    while( <LOG> ){
      /^END/ and last;
      /(v*strict)=t/    and $bbmerge_setting    = $1;
      /Pairs:\s+(\d+)/  and $num_nuclear_pairs  = $1;
      /Joined:\s+(\d+)/ and $num_nuclear_merged = $1;
    }
  }elsif( /BEGIN ORGANELLE BBMERGE/ ){
    while( <LOG> ){
      /^END/ and last;
      /Pairs:\s+(\d+)/ and $num_organelle_pairs = $1;
    }
  }elsif( /BEGIN RRNA BBMERGE/ ){
    while( <LOG> ){
      /^END/ and last;
      /Pairs:\s+(\d+)/ and $num_rrna_pairs = $1;
    }
  }elsif( /BEGIN BOWTIE2 ORGANELLE/ ){
    while( <LOG> ){
      /^END/ and last;
    }
  }elsif( /BEGIN BOWTIE2 RRNA/ ){
    while( <LOG> ){
      /^END/ and last;
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

my $complete_buscos;
my $fragmented_buscos;
my $total_searched_buscos;

open( BUSCO, "$dir/run_$sample_ID\_BUSCO/short_summary_$sample_ID\_BUSCO" ) || die "Couldn't open $dir/run_$sample_ID\_BUSCO/short_summary_$sample_ID\_BUSCO: $!\n";
while( <BUSCO> ){

  /(\d+)\s+Complete BUSCOs/             and $complete_buscos = $1;
  /(\d+)\s+Fragmented BUSCOs/           and $fragmented_buscos = $1;
  /(\d+)\s+Total BUSCO groups searched/ and $total_searched_buscos = $1;

}
close BUSCO;

################ PARSE TRANSRATE SUMMARY ################

my $num_contigs;
my $num_bases;
my $contig_N50;
my $assembly_gc;
my $p_good_mapping;
my $p_refs_with_crbb;
my $transrate_score;

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

print join "\t", "sample_ID", "num_reads (F+R)", "num_corrected_bases", "num_dropped_reads", "percent_organelle_reads", "percent_rrna_reads", "num_nuclear_reads", "bbmerge_setting", "percent_merged", "Assembled contigs", "Assembled bases", "Assembly GC", "Contig N50", "Percent good mappings", "Percent Phaeo CRBB", "Transrate score", "BUSCOs (complete+fragmented)", "BUSCOs (percent)", "\n";

my $num_raw_reads_commify     = commify($num_raw_reads);
my $num_nuclear_reads_commify = commify($num_nuclear_pairs*2);
$num_corrected_bases          = commify($num_corrected_bases);
$num_bases                    = commify($num_bases);
$num_contigs                  = commify($num_contigs);

printf "$sample_ID\t$num_raw_reads_commify\t$num_corrected_bases\t$num_trimmomatic_dropped\t%.2f\t%.2f\t", ($num_organelle_pairs*200)/$num_raw_reads, ($num_rrna_pairs*200)/$num_raw_reads;
print  $num_nuclear_reads_commify, "\t", $bbmerge_setting, "\t";
printf "%.2f\t", ($num_nuclear_merged/$num_nuclear_pairs)*100;
printf "$num_contigs\t$num_bases\t%.2f\t$contig_N50\t", $assembly_gc;
printf "%.2f\t%.2f\t%.2f\t", $p_good_mapping*100, $p_refs_with_crbb*100, $transrate_score;
print  $complete_buscos + $fragmented_buscos, "\t";
printf "%.2f", (($complete_buscos + $fragmented_buscos)/$total_searched_buscos)*100;
print "\n";

exit;


sub commify {
  my $text = reverse $_[0];
  $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return scalar reverse $text
}
