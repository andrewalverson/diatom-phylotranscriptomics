#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Benchmark;
use DBI;

#----------------------------------------------------------------------
# set global variables
#----------------------------------------------------------------------
my $QUEUE     = 'mem512GB64core'; # job queue for Trinity PBS script
my $USE_SQL;                      # boolean (defined and query the MySQL database for genus, species, etc.)
my $HOST      = 'hpceq';          # the MySQL server on UARK RAZOR cluster
my $DB        = 'transcriptome';  # MySQL datbase
my $TABLE     = 'RNA_Book';       # table within MySQL database
my $PBS_ONLY;                # skip analyses and just generate a Trinity job script
my $WALLTIME      = 48;      # wall time for Trinity run
my $SLIDINGWINDOW = '4:2';   # Trimmomatic SLIDINGWINDOW paramter (window-size:avg-quality)
my $TRIM_PHRED    =  2;      # Trimmomatic LEADING and TRAILING partameters - bases below this Phred score will be trimmed from 5' and 3' ends of reads
my $MIN_LEN       = 30;      # Trimmomatic MINLEN parameter - remove reads shorter than this
my $AVG_QUAL      = 28;      # Trimmomatic AVGQUAL parameter - drop the read if the average quality is below the specified level
my $PHRED_OUT     = 64;      # Trimmomatic TOPHRED33 or TOPHRED64 base quality coding of the output sequences
my $STRANDED      =  1;      # Trinity, SS_lib_type parameter to specify a stranded library
my $NORMALIZE     =  1;      # Trinity, --no_normalize_reads parameter (Trinity normalizes reads by default after 21 Sept 2016)
my $MIN_KMER_RNA  =  1;      # Trinity, min_kmer_cov parameter (rRNA assembly)
my $MIN_KMER_NUC  =  1;      # Trinity, min_kmer_cov parameter (nuclear assembly)
my $REFERENCE     = 'phaeo'; # Transrate 'reference' parameter, whether to use Phaeo or Thaps as the reference proteome


#----------------------------------------------------------------------
# set paths to directories, bowtie2 databases, and executables
#----------------------------------------------------------------------
my $parent                = '/scratch/aja/rnaseq';
my $rcorrector_executable = '/share/apps/bioinformatics/Rcorrector/run_rcorrector.pl';
my $trinity_plugins       = '/share/apps/trinity/trinityrnaseq-2.3.2/trinity-plugins';
my $trinity_utils         = '/share/apps/trinity/trinityrnaseq-2.3.2/util';
my $transrate_executable  = '/share/apps/transrate/transrate-1.0.1-linux-x86_64/transrate';
my $parallel_executable   = '/share/apps/parallel/20150822/bin/parallel';
my $bowtie_univec_db      = '/storage/aja/bowtie_univec_database/UniVec_Core';
my $bowtie_lsu_db         = '/storage/aja/bowtie_lsu_database/diatom_lsu';
my $bowtie_ssu_db         = '/storage/aja/bowtie_ssu_database/diatom_ssu';
my $bowtie_TOL_ssu_db     = '/storage/aja/bowtie_ssu_tree_of_life/ssu_tree_of_life';
my $blast_diatom_rrna_db  = '/storage/aja/blast_ssu_rrna_database/diatom_ssu.fa';
my $blast_tol_rrna_db     = '/storage/aja/bowtie_ssu_tree_of_life/ssu_tree_of_life';
my $adapterDB             = "/storage/aja/illumina_adapter_database/TruSeq_adapters.fa";
my $swissprot_db          = '/storage/aja/swissprot_db/uniprot_sprot.fasta';
my $pfam_db               = '/storage/aja/pfam_db/Pfam-A.hmm';
my $phaeoProteins         = "/storage/aja/diatom_protein_dbs/phaeoAA.fa";
my $thapsProteins         = "/storage/aja/diatom_protein_dbs/thapsAA.fa";
my $copy_output_to        = '/storage/aja/rnaseq_assemblies';
my $copy_output_to2       = 'storage07:/data/data/rnaseq_assemblies';
my $trinity_PBS_script;      # name of Trinity PBS script
my $transdecoder_PBS_script; # name of Transdecoder PBS script


# read command line options
parseArgs();

# set the reference proteome for Transrate
$REFERENCE eq 'phaeo' and $REFERENCE = $phaeoProteins;
$REFERENCE eq 'thaps' and $REFERENCE = $thapsProteins;

# get the RNA ID and declare related variables
my $ID  = shift @ARGV;
my $suffix = "";
my( $logfile, $directory_id, $genus, $species, $cultureID, $fwd_reads, $rev_reads );

#----------------------------------------------------------------------
# parse the RNA ID
#----------------------------------------------------------------------
# determine whether this has an Alverson lab RNA or Library ID
if( $ID =~ /\A([LR])(\d+)([a-zA-Z])?/ ){ # Alverson lab RNA or Library ID
  my $type   = $1;
  my $id_num = $2;
  $3 and $suffix = $3;

  $type eq 'L' and $TABLE = 'Library_Book';
  $type eq 'R' and $TABLE = 'RNA_Book';

  # extract metadata from MySQL server
  ( $genus, $species, $cultureID, $logfile, $directory_id ) = sql_fetch( $id_num, $suffix, $TABLE );

  # remove '?' characters from genus and species names
  $genus   =~ s/\?//g;
  $species =~ s/\?//g;
  
  # establish file names for forward and reverse raw reads
  $fwd_reads = "$ID\_1.fq.gz";
  $rev_reads = "$ID\_2.fq.gz";

}else{  # generic/other RNA ID
  # establish names for log file and output directory
  $logfile = "$ID\_clean_filter.log";
  $directory_id = $ID;

  # establish file names for forward and reverse raw reads
  $fwd_reads = "$ID\_1.fq.gz";
  $rev_reads = "$ID\_2.fq.gz";
}


#----------------------------------------------------------------------
# make output directory if it doesn't already exist
#----------------------------------------------------------------------
-d "$parent/$directory_id" or system( "mkdir $parent/$directory_id" );


#----------------------------------------------------------------------
# Trinity PBS file only? If so, make it and exit.
#----------------------------------------------------------------------
if( $PBS_ONLY ){

  # accessing these subroutines only to get file names for rRNA and nuclear reads
  my( $fwd_trimmed, $rev_trimmed ) = run_Trimmomatic( $ID );

  my( $fwd_ssu_diatom, $rev_ssu_diatom )                                 = filter_rRNA( $ID, "ssu_diatom" );
  my( $fwd_ssu_tol, $rev_ssu_tol )                                       = filter_rRNA( $ID, "ssu_tol" );

  my( $fwd_lsu, $rev_lsu, $fwd_trimmed_filtered, $rev_trimmed_filtered ) = filter_rRNA( $ID, "lsu" );
  my( $fwd_nuc_org_merged, $rev_nuc_org_merged ) = merge_reads_with_BBMerge( $ID, $fwd_trimmed_filtered, $rev_trimmed_filtered, $logfile, "nuclear" );

  # write Trinity PBS job script
  ( $trinity_PBS_script, $transdecoder_PBS_script ) = writeTrinity_Transdecoder_PBS( $logfile, $fwd_trimmed, $rev_trimmed, $fwd_trimmed_filtered, $rev_trimmed_filtered, $fwd_nuc_org_merged, $rev_nuc_org_merged, $fwd_ssu_diatom, $rev_ssu_diatom, $fwd_ssu_tol, $rev_ssu_tol, $fwd_lsu, $rev_lsu, $ID, $genus, $species, $cultureID, $parent, $directory_id, $trinity_plugins, $trinity_utils );
  
  # exit the program
  exit;

}else{
  
  # log the current time
  my $timestamp = localtime;
  
  
  #----------------------------------------------------------------------
  # open log file and record day, time, and RNA ID
  #----------------------------------------------------------------------
  open( LOGFILE, ">", $logfile  ) || die "Can't open $logfile: $!\n";
  print LOGFILE $timestamp, "\n\n";
  print LOGFILE $0, "\n\n";
  if( $genus ){
    print LOGFILE "$genus, ";
    $species =~ /\S+/ and print lc $species;
    print LOGFILE "Sample ID: $ID\n\n";
  }else{
    print LOGFILE "Sample ID: $ID\n\n";
  }
  close LOGFILE;


  #----------------------------------------------------------------------
  # use rcorrector for error correction
  #----------------------------------------------------------------------
  my( $fwd_corrected, $rev_corrected ) = run_rcorrector( $ID, $logfile, $fwd_reads, $rev_reads );


  #----------------------------------------------------------------------
  # use Trimmomatic to remove adapters, trim reads, and remove low-quality read pairs
  #----------------------------------------------------------------------
  my( $fwd_trimmed, $rev_trimmed ) = run_Trimmomatic( $ID, $logfile, $fwd_corrected, $rev_corrected );

  
  #----------------------------------------------------------------------
  # filter UniVec and rRNA sequences with Bowtie2
  #----------------------------------------------------------------------
  # Bowtie2 options
  # -1 -> fwd reads input
  # -2 -> rev reads input
  # --phred33    -> quality scores are PHRED+33
  # --local      -> perform local alignments
  # --un-conc-gz -> output unmapped reads (i.e. non-organelle, this is what I want)
  # --al-conc-gz -> output mapped reads (i.e. organelle reads, might want these later on)
  
  #----------------------------------------------------------------------
  # Get Bowtie2 version
  #----------------------------------------------------------------------
  open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
  print LOGFILE "# ----------------------------------------------------------------------\n";
  print LOGFILE "# Bowtie2 version\n";
  print LOGFILE "# ----------------------------------------------------------------------\n\n";
  close LOGFILE;
  system ( "bowtie2 --version >> $logfile" );
  system ( "echo >> $logfile" );


  #---------------------------------------------------------------------
  # run Bowtie2 to filter out vector contaminants
  #----------------------------------------------------------------------
  my( $fwd_noVec, $rev_noVec ) = filter_univec( $ID, $logfile, $fwd_trimmed, $rev_trimmed, $bowtie_univec_db );


  #----------------------------------------------------------------------
  # run Bowtie2 to filter out SSU and LSU rRNA sequences
  #----------------------------------------------------------------------
  # filter diatom SSU rRNA reads - no need to keep unmapped reads, will make an unmapped SSU reads file based on mapping to the TOL SSU db
  my( $fwd_ssu_diatom, $rev_ssu_diatom ) = filter_rRNA( $ID, "ssu_diatom", $logfile, $fwd_noVec, $rev_noVec, $bowtie_ssu_db );

  # filter all (TOL) SSU rRNA reads - unmapped read files should be devoid of *ALL* SSU reads (diatoms + everything else)
  my( $fwd_ssu_tol, $rev_ssu_tol, $fwd_noVec_noSSU, $rev_noVec_noSSU ) = filter_rRNA( $ID, "ssu_tol", $logfile, $fwd_noVec, $rev_noVec, $bowtie_TOL_ssu_db );

  # filter diatom LSU rRNA reads - unmapped read files should be devoid of all rRNA reads (SSU + most LSU)
  my( $fwd_lsu, $rev_lsu, $fwd_trimmed_filtered, $rev_trimmed_filtered ) = filter_rRNA( $ID, "lsu", $logfile, $fwd_noVec_noSSU, $rev_noVec_noSSU, $bowtie_lsu_db );

  
  #----------------------------------------------------------------------
  # use BBMerge to merge nuclear+organelle reads
  #----------------------------------------------------------------------
  # merge nuclear reads wth BBMerge
  my( $fwd_nuc_org_merged, $rev_nuc_org_merged ) = merge_reads_with_BBMerge( $ID, $fwd_trimmed_filtered, $rev_trimmed_filtered, $logfile, "nuclear" );

  
  #----------------------------------------------------------------------
  # generate Trinity PBS job script
  #----------------------------------------------------------------------
  ( $trinity_PBS_script, $transdecoder_PBS_script ) = writeTrinity_Transdecoder_PBS( $logfile, $fwd_trimmed, $rev_trimmed, $fwd_trimmed_filtered, $rev_trimmed_filtered, $fwd_nuc_org_merged, $rev_nuc_org_merged, $fwd_ssu_diatom, $rev_ssu_diatom, $fwd_ssu_tol, $rev_ssu_tol, $fwd_lsu, $rev_lsu, $ID, $genus, $species, $cultureID, $parent, $directory_id, $trinity_plugins, $trinity_utils );

  
  #----------------------------------------------------------------------
  # move or delete output files
  #----------------------------------------------------------------------
  open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
  print LOGFILE "\n# ----------------------------------------------------------------------\n";
  print LOGFILE "# Cleaning up and moving output files\n";
  print LOGFILE "# ----------------------------------------------------------------------\n\n";
  
  # remove Rcorrector, Trimmomatic, and Bowtie2 UniVec output files - don't need these
  print LOGFILE "Deleting intermediate files:\n";
  print LOGFILE "    $fwd_corrected\, $rev_corrected\, $fwd_noVec\, $rev_noVec\, $fwd_noVec_noSSU\, $rev_noVec_noSSU\n\n";
  system( "rm $fwd_corrected $rev_corrected $fwd_noVec $rev_noVec $fwd_noVec_noSSU $rev_noVec_noSSU" );

  # move files to output directory created for this species/RNA_ID
  print LOGFILE "Moving log file; Trinity PBS script; nuclear and rRNA reads to:\n";
  print LOGFILE "    $parent/$directory_id/\n\n";

  system( "mv $logfile $trinity_PBS_script $transdecoder_PBS_script $fwd_nuc_org_merged $rev_nuc_org_merged $fwd_trimmed $rev_trimmed $fwd_trimmed_filtered $rev_trimmed_filtered $fwd_ssu_diatom $rev_ssu_diatom $fwd_ssu_tol $rev_ssu_tol $fwd_lsu $rev_lsu $ID\_rnaseq_clean_filter.pbs $parent/$directory_id/" );
  
  # create shell script to:
  # (1) remove processed read files (can recreate these), (2) rsync output to storage07, (3) delete output directory from scratch
  # my $rm_outdir = "rm_$ID\_outdir.sh";
  # open( RMD, ">", "$parent/$rm_outdir"  ) || die "Can't open $rm_outdir: $!\n";
  # print RMD "rm $parent/$directory_id/*.fq.gz\n";
  # print RMD "rsync -a $parent/$directory_id storage07:/data/data/rnaseq/\n";
  # print RMD "rm -r $parent/$directory_id/\n";
  # close RMD;
  
  # make these rm shell scripts executable
  # system( "chmod u+x $parent/$rm_outdir" );

}

# exit the program
exit;


############################################SUBROUTINES#############################################
sub sql_fetch{

  my ($id_num, $suffix, $table) = @_;
  
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
    $dir_id = $type . $id_num . $suffix . "_" . lcfirst $hash->{Genus};
    if( $hash->{Species} ){
      $hash->{Species} =~ s/\s+/_/g;
      $dir_id .= "_" . lc $hash->{Species};
    }
    $dir_id .= "_" .  $hash->{Culture_ID};
  }else{
    $dir_id = $type . $id_num . $suffix . "_" . $hash->{Culture_ID};
  }
  
  # remove '?' charcters from directory ID
  $dir_id =~ s/\?//g;
  $dir_id =~ s/\.//g;
  
  # build name for log file name
  my $logname = "$type$id_num$suffix\_clean_filter.log";

  # disconnect from the MySQL database
  $dbh->disconnect;

  # print $dir_id, "\n" and exit;
  
  return( $hash->{Genus}, $hash->{Species}, $hash->{Culture_ID}, $logname, $dir_id );
  
}
####################################################################################################
sub run_rcorrector{

  # this subroutine uses rcorrector to correct errors in the raw reads

  my ( $ID, $logfile, $fwd_in, $rev_in ) = @_;

  my $fwd_out = "$ID\_1.cor.fq.gz"; # (fwd)
  my $rev_out = "$ID\_2.cor.fq.gz"; # (rev)

  # rcorrector settings based on: http://oyster-river-protocol.readthedocs.org/en/latest/
  my $rcorrector_command = "$rcorrector_executable -k 31 -t 15 -1 $fwd_in -2 $rev_in 2>>$logfile";
  
  open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
  print LOGFILE "# ----------------------------------------------------------------------\n";
  print LOGFILE "# Correcting errors in raw reads with rcorrector\n";
  print LOGFILE "# ----------------------------------------------------------------------\n\n";
  print LOGFILE "rcorrector run with: \'$rcorrector_command\'\n\n";
  print LOGFILE "BEGIN RCORRECTOR;\n";
  close LOGFILE;
  
  my $start_time = new Benchmark;

  # run rcorrector
  system( $rcorrector_command );
  
  my $end_time   = new Benchmark;
  my $difference = timediff( $end_time, $start_time );

  # re-open log file for printing
  open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
  print LOGFILE "\nrcorrector CPU: ", timestr( $difference ), "\n";
  print LOGFILE "Output files: $fwd_out, $rev_out\n";
  print LOGFILE "END RCORRECTOR;\n\n\n\n";
  close LOGFILE;
  
  return( $fwd_out, $rev_out );
  
}
####################################################################################################
sub run_Trimmomatic{

  # this subroutine uses Trimmomatic to trim reads and remove adapter sequences

  my ( $ID, $logfile, $fwd_in, $rev_in ) = @_;

  # 1 - remove Illumina adapters (Trimmomatic: ILLUMINACLIP using paramaters 2:40:15; chosen based on http://dx.doi.org/10.3389/fgene.2014.00013)
  # 2 - perform sliding window trimming (Trimmomatic: SLIDINGWINDOW; default parameters = 4:2; chosen based on http://dx.doi.org/10.3389/fgene.2014.00013)
  # 3 - remove reads <= $MIN_LEN bp in length  (Trimmomatic: MINLEN;  default = 25; chosen based on http://dx.doi.org/10.3389/fgene.2014.00013)

  my $fwd_out = "$ID\_trimmed_1.fq.gz"; # (fwd)
  my $rev_out = "$ID\_trimmed_2.fq.gz"; # (rev)

  if( $PBS_ONLY ){
    return( $fwd_out, $rev_out );
  }
  
  $AVG_QUAL =~ /none/i and $AVG_QUAL = '';
  $AVG_QUAL =~ /\d+/   and $AVG_QUAL = "AVGQUAL:$AVG_QUAL";
    
  my $trimmomatic_command = "java -jar /share/apps/trimmomatic/Trimmomatic-0.32/trimmomatic-0.32.jar PE $fwd_in $rev_in $fwd_out $ID\_junk\_1.fq.gz $rev_out $ID\_junk_2.fq.gz ILLUMINACLIP:$adapterDB:2:40:15 LEADING:$TRIM_PHRED TRAILING:$TRIM_PHRED SLIDINGWINDOW:$SLIDINGWINDOW MINLEN:$MIN_LEN $AVG_QUAL TOPHRED$PHRED_OUT 2>>$logfile";

  open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
  print LOGFILE "# ----------------------------------------------------------------------\n";
  print LOGFILE "# Trimming and cleaning reads with Trimmomatic\n";
  print LOGFILE "# ----------------------------------------------------------------------\n\n";
  print LOGFILE "Trimmomatic run with: \'$trimmomatic_command\'\n\n";
  print LOGFILE "BEGIN TRIMMOMATIC;\n";
  close LOGFILE;
  
  my $start_time = new Benchmark;
  
  # run Trimmomatic, redirecting run info (printed to STDERR) to $logfile
  system( $trimmomatic_command );
  
  my $end_time   = new Benchmark;
  my $difference = timediff( $end_time, $start_time );

  # re-open log file for printing
  open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
  print LOGFILE "\nTrimmomatic CPU: ", timestr( $difference ), "\n";
  print LOGFILE "Output files: $fwd_out, $rev_out\n";
  print LOGFILE "END TRIMMOMATIC;\n\n\n\n";
  close LOGFILE;
  
  # remove unpaired output files - no need for these
  system( "rm $ID\_junk_1.fq.gz $ID\_junk_2.fq.gz" );
  
  return( $fwd_out, $rev_out );
  
}
####################################################################################################
sub filter_univec{

  # this subroutine runs Bowtie2 to filter out vector contaminants

  my ( $ID, $logfile, $fwd_in, $rev_in, $bowtie2_db ) = @_;

  # establish base name for unmapped reads output
  my $unmapped_basename = "$ID\_noVec_%.fq.gz";
  ( my $fwd_out = $unmapped_basename ) =~ s/%/1/;
  ( my $rev_out = $unmapped_basename ) =~ s/%/2/;
  
  # bowtie univec command
  my $bowtie_univec_cmd = "bowtie2 -x $bowtie2_db -1 $fwd_in -2 $rev_in --phred33 --local --un-conc-gz $unmapped_basename >/dev/null 2>>$logfile";

  # open log file for printing
  open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
  print LOGFILE "# ----------------------------------------------------------------------\n";
  print LOGFILE "# Filtering UniVec reads with Bowtie2\n";
  print LOGFILE "# ----------------------------------------------------------------------\n\n";
  print LOGFILE "Bowtie2 run with: \'$bowtie_univec_cmd\'\n\n";
  print LOGFILE "BEGIN BOWTIE2 UNIVEC;\n";
  close LOGFILE;

  my $start_time = new Benchmark;

  system( $bowtie_univec_cmd );
    
  my $end_time   = new Benchmark;
  my $difference = timediff( $end_time, $start_time );

  # re-open log file for printing
  open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
  print LOGFILE "\nBowtie2 UniVec filter CPU: ", timestr( $difference ), "\n";
  print LOGFILE "Output files: $fwd_out, $rev_out\n";
  print LOGFILE "END BOWTIE2 UNIVEC;\n\n\n\n";
  close LOGFILE;

  return( $fwd_out, $rev_out );
  
}
####################################################################################################
sub filter_rRNA{

  # this subroutine runs Bowtie2 to filter out rRNA reads
  # note $rna_type is either "ssu_diatom" or "ssu_tol" [tree of life SSU] or "lsu"

  my ( $ID, $rna_type, $logfile, $fwd_in, $rev_in, $bowtie2_db ) = @_;
  
  # establish base name for unmapped reads output
  my $unmapped_basename = "$ID\_noVec_no_$rna_type\_%.fq.gz";
  ( my $fwd_unmapped = $unmapped_basename ) =~ s/%/1/;
  ( my $rev_unmapped = $unmapped_basename ) =~ s/%/2/;

  # filename for mapped reads, output by Bowtie2
  my $rrna_reads = "$ID\_$rna_type\_reads_%.fq.gz";
  ( my $fwd_rrna = $rrna_reads ) =~ s/%/1/;
  ( my $rev_rrna = $rrna_reads ) =~ s/%/2/;

  if( $PBS_ONLY ){
    my $fwd_unmapped = "$ID\_trimmed_filtered_1.fq.gz";
    my $rev_unmapped = "$ID\_trimmed_filtered_2.fq.gz";

    return( $fwd_rrna, $rev_rrna, $fwd_unmapped, $rev_unmapped );
  }else{


    # generate appropriate bowtie command
    my $bowtie_rrna_cmd;
    if( $rna_type eq "ssu_diatom" ){
      $bowtie_rrna_cmd = "bowtie2 -x $bowtie2_db -1 $fwd_in -2 $rev_in --phred33 --end-to-end --very-sensitive --al-conc-gz $rrna_reads >/dev/null 2>>$logfile";
    }else{
      $bowtie_rrna_cmd = "bowtie2 -x $bowtie2_db -1 $fwd_in -2 $rev_in --phred33 --local --very-sensitive --un-conc-gz $unmapped_basename --al-conc-gz $rrna_reads >/dev/null 2>>$logfile";
    }

      
    # open log file for printing
    open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
    print LOGFILE "# ----------------------------------------------------------------------\n";
    print LOGFILE "# Filtering ", uc( $rna_type ), " rRNA reads with Bowtie2\n";
    print LOGFILE "# ----------------------------------------------------------------------\n\n";
    print LOGFILE "Bowtie2 run with: \'$bowtie_rrna_cmd\'\n\n";
    print LOGFILE "BEGIN BOWTIE2 ", uc( $rna_type ), " RRNA;\n";
    close LOGFILE;
    
    my $start_time = new Benchmark;
    
    system( $bowtie_rrna_cmd );
    
    my $end_time   = new Benchmark;
    my $difference = timediff( $end_time, $start_time );

    # re-open log file for printing
    open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
    print LOGFILE "\nBowtie2 ", uc( $rna_type ), " rRNA filter CPU: ", timestr( $difference ), "\n";

    # after the final rRNA filter (LSU), change names of unmapped read files to be generic: $ID_trimmed_filtered_1.fq.gz and $ID_trimmed_filtered_2.fq.gz
    if( $rna_type eq "lsu" ){
      print LOGFILE "Make generic and more informative filenames for the unmapped reads, which are the cleaned & filtered nuclear+organelle reads:\n";
      print LOGFILE "     mv $fwd_unmapped $ID\_trimmed_filtered_1.fq.gz\n";
      print LOGFILE "     mv $rev_unmapped $ID\_trimmed_filtered_2.fq.gz\n";

      system( "mv $fwd_unmapped $ID\_trimmed_filtered_1.fq.gz" );
      system( "mv $rev_unmapped $ID\_trimmed_filtered_2.fq.gz" );

      $fwd_unmapped = "$ID\_trimmed_filtered_1.fq.gz";
      $rev_unmapped = "$ID\_trimmed_filtered_2.fq.gz";

      print LOGFILE "Output files: $fwd_rrna, $rev_rrna, $fwd_unmapped, $rev_unmapped\n";
    }elsif( $rna_type eq "ssu_tol" ){
      print LOGFILE "Output files: $fwd_rrna, $rev_rrna, $fwd_unmapped, $rev_unmapped\n";
    }else{
      print LOGFILE "Output files: $fwd_rrna, $rev_rrna\n";
    }

    print LOGFILE "END BOWTIE2 ", uc( $rna_type ), " RRNA;\n\n\n\n";
    close LOGFILE;
  }

  return( $fwd_rrna, $rev_rrna, $fwd_unmapped, $rev_unmapped );

}
####################################################################################################
sub merge_reads_with_BBMerge{

  # this subroutine runs BBMerge to merge paired-end reads

  my ( $ID, $fwd_in, $rev_in, $logfile, $type ) = @_;

  if( $PBS_ONLY ){
    return( "$ID\_$type\_merged_fwd.fq.gz", "$ID\_$type\_merged_rev.fq.gz" );
  }else{

    # BBMerge command
    my $bbmerge = "bbmerge.sh vstrict=t in1=$fwd_in in2=$rev_in out=$ID\_merged.fq outu1=$ID\_unmerged1.fq outu2=$ID\_unmerged2.fq.gz 2>>$logfile";
    
    # open log file for printing
    open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
    print LOGFILE "# ----------------------------------------------------------------------\n";
    print LOGFILE "# Merging $type reads with BBMerge\n";
    print LOGFILE "# ----------------------------------------------------------------------\n\n";
    print LOGFILE "BBMERGE run with: \'$bbmerge\'\n\n";
    print LOGFILE "BEGIN ", uc $type, " BBMERGE;\n";
    close LOGFILE;
    
    my $start_time = new Benchmark;
    
    system( $bbmerge );
    
    # concatenate and compress fwd unmerged and merged read files
    system( "cat $ID\_unmerged1.fq $ID\_merged.fq | gzip > $ID\_$type\_merged_fwd.fq.gz" );
    
    # rename rev unmerged reads
    system( "mv $ID\_unmerged2.fq.gz $ID\_$type\_merged_rev.fq.gz" );
    
    # remove intermediate files
    system( "rm $ID\_unmerged1.fq $ID\_merged.fq" );
    
    my $end_time   = new Benchmark;
    my $difference = timediff( $end_time, $start_time );
    
    # re-open log file for printing
    open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
    print LOGFILE "\nBBMerge $type CPU: ", timestr( $difference ), "\n";
    print LOGFILE "File actions:\n";
    print LOGFILE "     cat $ID\_unmerged1.fq $ID\_merged.fq | gzip > $ID\_$type\_merged_fwd.fq.gz\n";
    print LOGFILE "     mv $ID\_unmerged2.fq.gz $ID\_$type\_merged_rev.fq.gz\n";
    print LOGFILE "     rm $ID\_unmerged1.fq $ID\_merged.fq\n";
    print LOGFILE "END ", uc $type, " BBMERGE;\n\n\n\n";
    close LOGFILE;
  }

  return( "$ID\_$type\_merged_fwd.fq.gz", "$ID\_$type\_merged_rev.fq.gz" );
  
}
####################################################################################################
sub writeTrinity_Transdecoder_PBS{

  my( $logfile, $fwd_trimmed, $rev_trimmed, $fwd_nuc_org, $rev_nuc_org, $fwd_nuc_org_merged, $rev_nuc_org_merged, $fwd_ssu_diatom, $rev_ssu_diatom, $fwd_ssu_tol, $rev_ssu_tol, $fwd_lsu, $rev_lsu, $ID, $genus, $species, $cultureID, $parent, $directory_id, $trinity_plugins, $trinity_utils ) = @_;

  my $ppn; # number processors
  my $mem; # available memory
  my $trinity_pbs      = "trinity_$ID.pbs";
  my $transdecoder_pbs = "transdecoder_$ID.pbs";

  # files names for the unzipped nuclear reads
  ( my $fwd_fq = $fwd_nuc_org ) =~ s/\.gz$//;
  ( my $rev_fq = $rev_nuc_org ) =~ s/\.gz$//;
  
  # name of the Trinity assembly FASTA file when '--output' is not specified
  my $trinity_assembly_file = 'trinity_out_dir.Trinity.fasta';

  # for selecting a random queue
  my @queue_list = ( 'mem512GB64core', 'mem768GB32core' );

  # choose the queue if random
  $QUEUE eq 'random' and $QUEUE = $queue_list[rand @queue_list];

  if( $QUEUE eq 'aja' ){
    $ppn = 16;
    $mem = 32;
  }elsif( $QUEUE eq 'mem96GB12core' ){
    $ppn = 12;
    $mem = 96;
  }elsif( $QUEUE eq 'mem512GB64core' ){
    $ppn = 64;
    $mem = 512;
  }elsif( $QUEUE eq 'mem768GB32core' ){
    $ppn = 32;
    $mem = 768;
  }else{
    die "Queue \'$QUEUE\' not recognized\n";
  }

  # format the normalization parameter
  if( $NORMALIZE ){
    $NORMALIZE = '';
  }else{
    $NORMALIZE = '--no_normalize_reads';
  }

  
  ##################### generate Trinity/RSEM/Transrate/BUSCO PBS job script #####################
  open( TRINITY, '>', $trinity_pbs ) || die "Can't open $trinity_pbs: $!\n";
  
  # print PBS job parameters
  print TRINITY "#PBS -N assemble_$ID\n";
  print TRINITY '#PBS -j oe', "\n";
  print TRINITY '#PBS -m abe', "\n";
  print TRINITY '#PBS -M aja@uark.edu', "\n";
  print TRINITY '#PBS -q ' . $QUEUE . "\n";
  print TRINITY '#PBS -l nodes=1:ppn='. $ppn . "\n";
  print TRINITY '#PBS -l walltime=' . $WALLTIME . ':00:00', "\n";
  print TRINITY "#PBS -o assemble_$ID\." . '$PBS_JOBID' . "\n\n";
  print TRINITY 'cd $PBS_O_WORKDIR', "\n\n";

  # load modules
  print TRINITY "# load Trinity modules\n";
  print TRINITY 'module purge', "\n";
  print TRINITY 'module load gcc/4.8.2', "\n";
  print TRINITY 'module load gcc/4.9.1', "\n";
  print TRINITY 'module load boost/1.57.0', "\n";
  print TRINITY 'module load perl/5.10.1', "\n";
  print TRINITY 'module load python/3.3.2', "\n";
  print TRINITY 'module load java/sunjdk_1.8.0', "\n";
  print TRINITY 'module load bowtie', "\n";
  print TRINITY 'module load samtools/0.1.19', "\n";
  print TRINITY 'module load trinity/2.3.2', "\n";
  print TRINITY 'module load blast/2.2.29+', "\n";
  print TRINITY 'module load hmmer/3.1b2', "\n";
  print TRINITY 'module load augustus/3.2.2', "\n";
  print TRINITY 'module load busco/1.1b1', "\n";
  print TRINITY 'module load cdhit/4.6.5', "\n\n";


  # assemble diatom SSU rRNA reads
  print TRINITY "# assemble diatom SSU rRNA reads\n";
  if( $STRANDED ){
    print TRINITY "Trinity --seqType fq --min_kmer_cov $MIN_KMER_RNA --max_memory ", $mem, 'G ', '--CPU ', $ppn, " --left $fwd_ssu_diatom --right $rev_ssu_diatom --SS_lib_type RF --no_normalize_reads --full_cleanup\n";
  }else{
    print TRINITY "Trinity --seqType fq --min_kmer_cov $MIN_KMER_RNA --max_memory ", $mem, 'G ', '--CPU ', $ppn, " --left $fwd_ssu_diatom --right $rev_ssu_diatom --no_normalize_reads --full_cleanup\n";
  }
  print TRINITY "mv $trinity_assembly_file $ID\_ssu_diatom_rrna.fa\n\n";

  # assemble TOL SSU rRNA reads
  print TRINITY "# assemble TOL SSU rRNA reads\n";
  if( $STRANDED ){
    print TRINITY "Trinity --seqType fq --min_kmer_cov $MIN_KMER_RNA --max_memory ", $mem, 'G ', '--CPU ', $ppn, " --left $fwd_ssu_tol --right $rev_ssu_tol --SS_lib_type RF --no_normalize_reads --full_cleanup\n";
  }else{
    print TRINITY "Trinity --seqType fq --min_kmer_cov $MIN_KMER_RNA --max_memory ", $mem, 'G ', '--CPU ', $ppn, " --left $fwd_ssu_tol --right $rev_ssu_tol --no_normalize_reads --full_cleanup\n";
  }
  print TRINITY "mv $trinity_assembly_file $ID\_ssu_tol_rrna.fa\n\n";

  # assemble LSU rRNA reads
  print TRINITY "# assemble LSU rRNA reads\n";
  if( $STRANDED ){
    print TRINITY "Trinity --seqType fq --min_kmer_cov $MIN_KMER_RNA --max_memory ", $mem, 'G ', '--CPU ', $ppn, " --left $fwd_lsu --right $rev_lsu --SS_lib_type RF --no_normalize_reads --full_cleanup\n";
  }else{
    print TRINITY "Trinity --seqType fq --min_kmer_cov $MIN_KMER_NUC --max_memory ", $mem, 'G ', '--CPU ', $ppn, " --left $fwd_lsu --right $rev_lsu --no_normalize_reads --full_cleanup\n";
  }
  print TRINITY "mv $trinity_assembly_file $ID\_lsu_rrna.fa\n\n";

  # assemble nuclear reads
  print TRINITY "# assemble nuclear reads\n";
  if( $STRANDED ){
    print TRINITY "Trinity --seqType fq --min_kmer_cov $MIN_KMER_NUC --max_memory ", $mem, 'G ', '--CPU ', $ppn, " --left $fwd_nuc_org_merged --right $rev_nuc_org_merged --SS_lib_type RF --full_cleanup $NORMALIZE\n\n";
  }else{
    print TRINITY "Trinity --seqType fq --min_kmer_cov $MIN_KMER_NUC --max_memory ", $mem, 'G ', '--CPU ', $ppn, " --left $fwd_nuc_org_merged --right $rev_nuc_org_merged --full_cleanup $NORMALIZE\n\n";
  }
  print TRINITY "mv $trinity_assembly_file $ID\_nuclear.fa\n\n";
  
  # BASH loop to simplify Trinity FASTA headers
  print TRINITY "# BASH loop to simplify Trinity FASTA headers\n";
  print TRINITY "for fasta in *.fa; do
     sed 's/\|/_/' \$fasta > m\n";
  if( $genus ){
    if( $species ){
      print TRINITY "     perl -anwe \'if(/^>/){print \"\$F[0] $ID\_$genus\_$species\\n\"}else{print}\' m > \$fasta\n";
    }else{
      print TRINITY "     perl -anwe \'if(/^>/){print \"\$F[0] $ID\_$genus\\n\"}else{print}\' m > \$fasta\n";
    }
  }else{
    print TRINITY "     perl -anwe \'if(/^>/){print \"\$F[0] $ID\\n\"}else{print}\' m > \$fasta\n";
  }
  print TRINITY "done\n\n";
  
  # BLAST SSU rRNA assemblies to diatom or TOL SSU databases and record the top hits
  print TRINITY "# BLAST diatom SSU rRNA contigs to diatom SSU database\n";
  print TRINITY "blastn -query $ID\_ssu_diatom_rrna.fa -db $blast_diatom_rrna_db -max_target_seqs 1 -max_hsps 1 -num_threads ", $ppn, " -outfmt 6 > $ID\_ssu_diatom_hits.blastn\n\n";

  print TRINITY "# BLAST TOL SSU rRNA contigs to TOL SSU database\n";
  print TRINITY "blastn -query $ID\_ssu_tol_rrna.fa -db $blast_tol_rrna_db -max_target_seqs 1 -max_hsps 1 -num_threads ", $ppn, " -outfmt 6 > $ID\_ssu_tol_hits.blastn\n\n";

  # print report for nuclear assembly
  print TRINITY "# print report for nuclear+organelle assembly\n";
  print TRINITY "$trinity_utils/TrinityStats.pl $ID\_nuclear.fa > $ID\_nuclear.stats\n\n";

  # estimate LSU transcript abundance with RSEM
  print TRINITY "# estimate LSU transcript abundance with RSEM\n";
  print TRINITY "$trinity_utils/align_and_estimate_abundance.pl --thread_count ", $ppn, " --transcripts $ID\_lsu_rrna.fa --seqType fq --left $fwd_trimmed --right $rev_trimmed --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --output_dir $ID\_lsu_RSEM\n\n";

  # estimate diatom SSU transcript abundance with RSEM
  print TRINITY "# estimate diatom SSU transcript abundance with RSEM\n";
  print TRINITY "$trinity_utils/align_and_estimate_abundance.pl --thread_count ", $ppn, " --transcripts $ID\_ssu_diatom_rrna.fa --seqType fq --left $fwd_trimmed --right $rev_trimmed --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --output_dir $ID\_ssu_diatom_RSEM\n\n";
  
  # estimate TOL SSU transcript abundance with RSEM
  print TRINITY "# estimate TOL SSU transcript abundance with RSEM\n";
  print TRINITY "$trinity_utils/align_and_estimate_abundance.pl --thread_count ", $ppn, " --transcripts $ID\_ssu_tol_rrna.fa --seqType fq --left $fwd_trimmed --right $rev_trimmed --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --output_dir $ID\_ssu_tol_RSEM\n\n";

  # estimate nuclear transcript abundance with RSEM
  print TRINITY "# estimate nuclear transcript abundance with RSEM\n";
  print TRINITY "$trinity_utils/align_and_estimate_abundance.pl --thread_count ", $ppn, " --transcripts $ID\_nuclear.fa --seqType fq --left $fwd_trimmed --right $rev_trimmed --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --output_dir $ID\_nuclear_RSEM\n\n";

  # assess nuclear assembly quality with transrate
  print TRINITY "# assess assembly quality with transrate\n";
  print TRINITY "$transrate_executable --threads ", $ppn, " --assembly $ID\_nuclear.fa --left $fwd_nuc_org --right $rev_nuc_org --reference $REFERENCE\n";
  print TRINITY "mv transrate_results/assemblies.csv transrate.csv\n\n";

  # identify number of conserved orthologs with BUSCO
  print TRINITY "# identify number of conserved orthologs with BUSCO\n";
  print TRINITY "python3 busco -o $ID\_BUSCO -in $ID\_nuclear.fa -l /storage/aja/busco_dbs/eukaryota -m trans -c ", $ppn, " \n\n";

  # move RSEM output
  print TRINITY "# move RSEM output to appropriate directories\n";
  print TRINITY "mv $ID\_lsu_rrna.fa.bowtie.* $ID\_lsu_rrna.fa.RSEM.* $ID\_lsu_RSEM/\n";
  print TRINITY "mv $ID\_ssu_diatom_rrna.fa.bowtie.* $ID\_ssu_diatom_rrna.fa.RSEM.* $ID\_ssu_diatom_RSEM/\n";
  print TRINITY "mv $ID\_ssu_tol_rrna.fa.bowtie.* $ID\_ssu_tol_rrna.fa.RSEM.* $ID\_ssu_tol_RSEM/\n";
  print TRINITY "mv $ID\_nuclear.fa.bowtie.* $ID\_nuclear.fa.RSEM.* $ID\_nuclear_RSEM/\n\n";

  #remove intermediate output files
  print TRINITY "# remove intermediate output files\n";
  print TRINITY "rm -r temp m transrate_results run_$ID\_BUSCO/hmmer_output run_$ID\_BUSCO/translated_proteins *.gene_trans_map\n";

  close TRINITY;


  ##################### generate Transdecoder/CD-HIT PBS job script #####################
  open( TRANSDECODER, '>', $transdecoder_pbs ) || die "Can't open $transdecoder_pbs: $!\n";

  # print PBS job parameters
  print TRANSDECODER "#PBS -N transdecoder_$ID\n";
  print TRANSDECODER '#PBS -j oe', "\n";
  print TRANSDECODER '#PBS -m abe', "\n";
  print TRANSDECODER '#PBS -M aja@uark.edu', "\n";
  print TRANSDECODER '#PBS -q med12core', "\n";
  print TRANSDECODER '#PBS -l nodes=1:ppn=12', "\n";
  print TRANSDECODER '#PBS -l walltime=72:00:00', "\n";
  print TRANSDECODER "#PBS -o transdecoder_$ID\." . '$PBS_JOBID' . "\n\n";
  print TRANSDECODER 'cd $PBS_O_WORKDIR', "\n\n";

  # load modules
  print TRANSDECODER 'module purge', "\n";
  print TRANSDECODER 'module load perl/5.10.1', "\n";
  print TRANSDECODER 'module load hmmer/3.1b2', "\n";
  print TRANSDECODER 'module load trinity/2.0.2', "\n";
  print TRANSDECODER 'module load transdecoder/2.0.1', "\n";
  print TRANSDECODER 'module load blast/2.2.29+', "\n";
  print TRANSDECODER 'module load cdhit/4.6.5', "\n\n";
  
  # use TransDecoder to translate nuclear contigs
  print TRANSDECODER "# TransDecoder to translate nuclear contigs\n";
  print TRANSDECODER "# TransDecoder step 1: Predict all long open reading frames\n";
  print TRANSDECODER "TransDecoder.LongOrfs -t $ID\_nuclear.fa\n\n";
  
  my $long_orfs = "$ID\_nuclear.fa\.transdecoder_dir/longest_orfs.pep";

  print TRANSDECODER "# TransDecoder step 2: Homology searches to SwissProt and Pfam\n";
  print TRANSDECODER "# TransDecoder step 2.1: Parallel BLAST search long ORFs from Step 1 to SwissProt db\n";
  print TRANSDECODER "cat $long_orfs | $parallel_executable --will-cite --block 100k --progress --recstart \'>\' --pipe /share/apps/blast/ncbi-blast-2.2.29+/bin/blastp -outfmt 6 -max_target_seqs 1 -evalue 1e-3 -db $swissprot_db -num_threads 4 -query - > uniref.blastp\n\n";
  print TRANSDECODER "# TransDecoder step 2.2: HMMR search search long ORFs to Pfam database\n";
  print TRANSDECODER "hmmscan --cpu 12 --domtblout pfam.domtblout $pfam_db $long_orfs\n\n";
  
  print TRANSDECODER "# TransDecoder step 3: Use Transdecoder to translate nuclear contigs using homology search information\n";
  print TRANSDECODER "TransDecoder.Predict -t $ID\_nuclear.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits uniref.blastp\n\n";

  # run CD-HIT to remove redundant contigs
  print TRANSDECODER "# Run CD-HIT to remove redundant/overlapping contigs\n";
  print TRANSDECODER "cd-hit -i $ID\_nuclear.fa.transdecoder.pep -o $ID\_cd-hit.pep -c 0.99 -n 5 -d 1000 -M 64000 -T 0\n";
  print TRANSDECODER 'rm -r *.clstr ', "$ID\_nuclear.fa.transdecoder_dir $ID\_nuclear.fa.transdecoder.gff3 $ID\_nuclear.fa.transdecoder.bed", ' *.bowtie.* pfam.domtblout uniref.blastp', "\n\n";

  # rsync output from razor to /storage/aja/rnaseq_assemblies ($copy_output_to) and storage07:/data/data/rnaseq_assemblies ($copy_output_to2)
  # Reminder: $parent = '/scratch/aja/rnaseq';
  print TRANSDECODER "Copying output directory to \'$copy_output_to\' and \'$copy_output_to2\'\n";
  print TRANSDECODER "    \'rsync -azv $parent/$directory_id $copy_output_to/\'\n";
  print TRANSDECODER "    \'rsync -azv $parent/$directory_id $copy_output_to2/\'\n\n";
  print TRANSDECODER "rsync -azv $parent/$directory_id $copy_output_to/alverson_lab\n";
  print TRANSDECODER "rsync -azv $parent/$directory_id $copy_output_to2/\n";

  close TRANSDECODER;

  # if PBS only, skip writing out to log file
  $PBS_ONLY and return ($trinity_pbs, $transdecoder_pbs);

  # get Trinity version
  open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
  print LOGFILE "# ----------------------------------------------------------------------\n";
  print LOGFILE "# Trinity version and PBS job script\n";
  print LOGFILE "# ----------------------------------------------------------------------\n\n";
  close LOGFILE;
  system ( "Trinity --version | grep '^Trinity' >> $logfile" );
  
  # identify Trinity PBS job script
  open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
  print LOGFILE "Trinity PBS job script written to: \'$trinity_pbs\'\n";
  print LOGFILE "Transdecoder PBS job script written to: \'$transdecoder_pbs\'\n\n\n";
  close LOGFILE;

  return ($trinity_pbs, $transdecoder_pbs);
  
}
####################################################################################################
sub parseArgs{

  my $usage = "\nUsage: $0 <ID> [options]

   Required: (1) Alverson lab RNA ID (e.g. 'R12' or 'R12a'), Library ID (e.g. L155), *OR*
             (2) base name for the read files (e.g. 'index-AGTTCC' for 'index-AGTTCC_1.fq.gz' and 'index-AGTTCC_2.fq.gz')


   General options
          --trinity_queue  - job queue ('mem96GB12core', 'mem512GB64core' [default], 'mem768GB32core', 'aja', 'random')
          --db             - MySQL database (default: 'transcriptome')
          --table          - MySQL table to read (default: 'RNA_Book')
          --walltime       - walltime for Trinity run (integer [hrs], default: 48)
          --pbs_only       - generate a Trinity job script only, skipping Trimmomatic and Bowtie2 (default: false)


   Trimmomatic options
          --window    - Trimmomatic SLIDINGWINDOW parameter specified as <windowSize><requiredQuality> (default: '4:2')
          --min_len   - filter reads shorter than this length threshold (default: 30)
          --avg_qual  - filter reads with average quality below the specified level (default: 28 [select --avg_qual=none for no avg_qual filter])
          --phred_in  - input  Phred scale (default: 33)
          --phred_out - output Phred scale (default: 64)


   Trinity options
          --min_kmer_rna  - Trinity --min_kmer_cov parameter for rRNA assembly (default: 1)
          --min_kmer_nuc  - Trinity --min_kmer_cov parameter for nuclear assembly (default: 1)
          --stranded      - Trinity --SS_lib_type parameter, set to RF when true (boolean, default: true [library is stranded])
          --normalize     - Trinity --no_normalize_reads parameter; Trinity normalizes by default (boolean, default: do not normalize)

   Transrate options
          --reference  - specify whether to use Phaeodactylum (phaeo) or T. pseudonana (thaps)
                         as the reference proteome (default: phaeo)

\n";
  
  my $result = GetOptions	(
	 'trinity_queue=s' => \$QUEUE,
	 'db=s'            => \$DB,
	 'table=s'         => \$TABLE,
	 'walltime=s'      => \$WALLTIME,
	 'pbs_only'        => \$PBS_ONLY,

	 'window=s'     => \$SLIDINGWINDOW,
	 'trim_phred=s' => \$TRIM_PHRED,
	 'min_len=s'    => \$MIN_LEN,
	 'avg_qual=s'   => \$AVG_QUAL,
	 'phred_out=s'  => \$PHRED_OUT,
#	 'phred_in=s'   => \$PHRED_IN,
				 
	 'min_kmer_rna=i'  => \$MIN_KMER_RNA,
	 'min_kmer_nuc=i'  => \$MIN_KMER_NUC,
	 'stranded!'       => \$STRANDED,
	 'normalize!'      => \$NORMALIZE,
				 
	 'reference=s'     => \$REFERENCE,
				);
  
  $ARGV[0] or die $usage;

}
####################################################################################################
