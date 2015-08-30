#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Benchmark;
use DBI;
use DBD::mysql;

#----------------------------------------------------------------------
# set global variables
#----------------------------------------------------------------------
my $QUEUE     = 'mem512GB64core'; # job queue for Trinity PBS script
my $USE_SQL;                      # boolean: 1 = query the MySQL database for genus, species, etc.
my $HOST      = 'tgv';            # the MySQL server on UARK RAZOR cluster
my $DB        = 'transcriptome';  # MySQL datbase
my $TABLE     = 'RNA_Book';       # table within MySQL database
my $NORMALIZE;        # whether or not to perform kmer normalization of nuclear reads with bbnorm
my $PBS_ONLY;         # skip analyses and just generate a Trinity job script
my $WALLTIME   = 18;  # wall time for Trinity run
my $HEADCROP   = 10;  # for Trimmomatic, number of 5' bp to trim
my $TRIM_PHRED = 5;   # for Trimmomatic, bases below this Phred score will be trimmed from 5' and 3' ends of reads
my $AVG_PHRED  = 32;  # for Trimmomatic, remove reads below this average Phred score
my $MIN_LEN    = 72;  # for Trimmomatic, remove reads shorter than this
my $MIN_KMER_RNA = 1; # for Trinity, min_kmer_cov parameter (rRNA assembly)
my $MIN_KMER_ORG = 1; # for Trinity, min_kmer_cov parameter (organelle assembly)
my $MIN_KMER_NUC = 2; # for Trinity, min_kmer_cov parameter (nuclear assembly)

#----------------------------------------------------------------------
# set paths to directories, bowtie2 databases, and executables
#----------------------------------------------------------------------
my $parent              = '/scratch/aja/rnaseq';
my $trinity_executable  = '/share/apps/trinity/trinityrnaseq-2.0.2/Trinity';
my $trinity_plugins     = '/share/apps/trinity/trinityrnaseq-2.0.2/trinity-plugins';
my $trinity_utils       = '/share/apps/trinity/trinityrnaseq-2.0.2/util';
my $bowtie2_executable  = '/share/apps/bowtie2/bowtie2-2.2.3/bowtie2';
my $bowtie_univec_db    = '/scratch/aja/bowtie_univec_database/UniVec_Core';
my $bowtie_rrna_db      = '/scratch/aja/bowtie_rrna_database/diatom_rRNA';
my $bowtie_organelle_db = '/scratch/aja/bowtie_organelle_database/diatom_cp-mt';
my $hmmer3_executable   = '/share/apps/hmmer-3.0/src/hmmscan';
my $pfam_db             = '/scratch/aja/pfam/Pfam-A.hmm';
my $diatom_aa_db        = '/scratch/aja/diatom_protein_dbs/Phatr_Thaps_aa.fa';
my $copy_output_to      = 'storage07:/data/data/rnaseq';
my $trinity_PBS_script; # name of Trinity PBS script


# read command line options
parseArgs();

# get the RNA ID and declare related variables
my $rnaID  = shift @ARGV;
my $suffix = "";
my( $logfile, $directory_id, $genus, $species, $cultureID, $fwd_reads, $rev_reads );


#----------------------------------------------------------------------
# parse the RNA ID
#----------------------------------------------------------------------
# determine whether this an Alverson lab RNA ID
if( $rnaID =~ /(\d+)([a-zA-Z])?/ ){ # Alverson lab RNA ID
  my $rnaNUM = $1;
  $2 and $suffix = $2;

  # extract metadata from MySQL server
  ( $genus, $species, $cultureID, $logfile, $directory_id ) = sql_fetch( $rnaNUM, );

  # establish file names for forward and reverse raw reads
  $fwd_reads = "$rnaID$suffix\_1.fq.gz";
  $rev_reads = "$rnaID$suffix\_2.fq.gz";

}else{                              # generic/other RNA ID
  # establish names for log file and output directory
  $logfile = "clean_filter_$rnaID$suffix.log";
  $directory_id = $rnaID;

  # establish file names for forward and reverse raw reads
  $fwd_reads = "$rnaID\_1.fq.gz";
  $rev_reads = "$rnaID\_2.fq.gz";
}


#----------------------------------------------------------------------
# verify read files exist; make output directory if it doesn't already exist
#----------------------------------------------------------------------
-e $fwd_reads or die "Expecting to find $fwd_reads, but it doesn't exist\n";
-e $rev_reads or die "Expecting to find $rev_reads, but it doesn't exist\n";
-d "$parent/$directory_id" or system( "mkdir $parent/$directory_id" );


#----------------------------------------------------------------------
# PBS file only? If so, make it and exit.
#----------------------------------------------------------------------
if( $PBS_ONLY ){
  # accessing these subroutines only to get file names for rRNA, organelle, and nuclear reads
  my( $fwd_rrna, $rev_rrna ) = filter_rRNA( $rnaID, $suffix );
  my( $fwd_nuc, $rev_nuc, $fwd_organelle, $rev_organelle ) = filter_organelle( $rnaID, $suffix );
  my( $fwd_nuc_norm, $rev_nuc_norm ) = normalize_with_BBNorm( $fwd_nuc, $rev_nuc, $logfile );

  # write Trinity PBS job script
  if( $NORMALIZE ){
    $trinity_PBS_script = writeTrinityPBS( $logfile, $fwd_nuc_norm, $rev_nuc_norm, $fwd_nuc, $rev_nuc, $fwd_rrna, $rev_rrna, $fwd_organelle, $rev_organelle, $rnaID, $suffix, $parent, $directory_id, $trinity_executable, $trinity_plugins, $trinity_utils );
  }else{
    $trinity_PBS_script = writeTrinityPBS( $logfile, $fwd_nuc, $rev_nuc, $fwd_rrna, $rev_rrna, $fwd_organelle, $rev_organelle, $rnaID, $suffix, $parent, $directory_id, $trinity_executable, $trinity_plugins, $trinity_utils );
  }
  
  # write Trinity PBS job script to the output directory for this RNA ID
  system( "mv $trinity_PBS_script $parent/$directory_id/" );

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
  print LOGFILE "$genus ", lc $species, ", RNA ID: $rnaID$suffix\n\n";
  close LOGFILE;


  #----------------------------------------------------------------------
  # use Trimmomatic to remove adapters, trim reads, and remove low-quality read pairs
  #----------------------------------------------------------------------
  my( $fwd_cleaned, $rev_cleaned ) = run_Trimmomatic( $rnaID, $suffix, $logfile, $fwd_reads, $rev_reads );

  
  #----------------------------------------------------------------------
  # filter UniVec, rRNA, and organelle sequences with Bowtie2
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
  system ( "$bowtie2_executable --version >> $logfile" );
  system ( "echo >> $logfile" );


  #---------------------------------------------------------------------
  # run Bowtie2 to filter out vector contaminants
  #----------------------------------------------------------------------
  my( $fwd_noVec, $rev_noVec ) = filter_univec( $rnaID, $suffix, $logfile, $fwd_cleaned, $rev_cleaned, $bowtie2_executable, $bowtie_univec_db );


  #----------------------------------------------------------------------
  # run Bowtie2 to filter out rRNA sequences
  #----------------------------------------------------------------------
  my( $fwd_noVec_norRNA, $rev_noVec_norRNA, $fwd_rrna, $rev_rrna ) = filter_rRNA( $rnaID, $suffix, $logfile, $fwd_noVec, $rev_noVec, $bowtie2_executable, $bowtie_rrna_db );


  #----------------------------------------------------------------------
  # run Bowtie2 to filter out organelle reads
  #----------------------------------------------------------------------
  my( $fwd_cleaned_filtered, $rev_cleaned_filtered, $fwd_organelle, $rev_organelle ) = filter_organelle( $rnaID, $suffix, $logfile, $fwd_noVec_norRNA, $rev_noVec_norRNA, $bowtie2_executable, $bowtie_organelle_db );


  #----------------------------------------------------------------------
  # if specified, use BBNorm to perform kmer normalization of nuclear reads
  #----------------------------------------------------------------------
  $NORMALIZE and my( $fwd_cleaned_filtered_normalized, $rev_cleaned_filtered_normalized ) = normalize_with_BBNorm( $fwd_cleaned_filtered, $rev_cleaned_filtered, $logfile );
  

  #----------------------------------------------------------------------
  # use BBMerge to merge rRNA, organelle, and nuclear reads
  #----------------------------------------------------------------------
  # merge rRNA reads wth BBMerge
  merge_reads_with_BBMerge( $fwd_rrna, $rev_rrna, $logfile, "rRNA" );

  # merge organelle reads wth BBMerge
  merge_reads_with_BBMerge( $fwd_organelle, $rev_organelle, $logfile, "organelle" );

  # merge nuclear reads wth BBMerge
  if( $NORMALIZE ){
    merge_reads_with_BBMerge( $fwd_cleaned_filtered_normalized, $rev_cleaned_filtered_normalized, $logfile, "nuclear" );
  }else{
    merge_reads_with_BBMerge( $fwd_cleaned_filtered, $rev_cleaned_filtered, $logfile, "nuclear" );
  }

  
  #----------------------------------------------------------------------
  # generate Trinity PBS job script
  #----------------------------------------------------------------------
  if( $NORMALIZE ){
    $trinity_PBS_script = writeTrinityPBS( $logfile, $fwd_cleaned_filtered_normalized, $rev_cleaned_filtered_normalized, $fwd_cleaned_filtered, $rev_cleaned_filtered, $fwd_rrna, $rev_rrna, $fwd_organelle, $rev_organelle, $rnaID, $suffix, $parent, $directory_id, $trinity_executable, $trinity_plugins, $trinity_utils, $hmmer3_executable, $pfam_db );
  }else{
    $trinity_PBS_script = writeTrinityPBS( $logfile, $fwd_cleaned_filtered, $rev_cleaned_filtered, $fwd_rrna, $rev_rrna, $fwd_organelle, $rev_organelle, $rnaID, $suffix, $parent, $directory_id, $trinity_executable, $trinity_plugins, $trinity_utils, $hmmer3_executable, $pfam_db );
  }


  #----------------------------------------------------------------------
  # move or delete output files
  #----------------------------------------------------------------------
  open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
  print LOGFILE "\n# ----------------------------------------------------------------------\n";
  print LOGFILE "# Cleaning up and moving output files\n";
  print LOGFILE "# ----------------------------------------------------------------------\n\n";
  
  # remove Trimmomatic and intermediate Bowtie2 UniVec output files - don't need these
  print LOGFILE "Deleting intermediate files:\n";
  print LOGFILE "    $fwd_cleaned\, $rev_cleaned\, $fwd_noVec\, $rev_noVec\, $fwd_noVec_norRNA\, $rev_noVec_norRNA\n\n";
  system( "rm $fwd_cleaned $rev_cleaned $fwd_noVec $rev_noVec $fwd_noVec_norRNA $rev_noVec_norRNA" );

  # move files to output directory created for this species/RNA_ID
  print LOGFILE "Moving log file; Trinity PBS script; nuclear, organelle, and rRNA reads to:\n";
  print LOGFILE "    $parent/$directory_id/\n\n";

  if( $NORMALIZE ){
    system( "mv $logfile $trinity_PBS_script $fwd_cleaned_filtered_normalized $rev_cleaned_filtered_normalized $fwd_cleaned_filtered $rev_cleaned_filtered $fwd_rrna $rev_rrna $fwd_organelle $rev_organelle $parent/$directory_id/" );
  }else{
    system( "mv $logfile $trinity_PBS_script $fwd_cleaned_filtered $rev_cleaned_filtered $fwd_rrna $rev_rrna $fwd_organelle $rev_organelle $parent/$directory_id/" );
  }
  
  # create shell script for deleting raw reads from /scratch/aja
  my $rm_reads = "rm_$rnaID$suffix\_raw_reads.sh";
  open( RMR, ">", $rm_reads  ) || die "Can't open $rm_reads: $!\n";
  print RMR "rm $fwd_reads $rev_reads\n";
  close RMR;
  
  # create shell script to:
  # (1) remove processed read files (can recreate these), (2) rsync output to storage07, (3) delete output directory from scratch
  my $rm_outdir = "rm_$rnaID$suffix\_outdir.sh";
  open( RMD, ">", "$parent/$rm_outdir"  ) || die "Can't open $rm_outdir: $!\n";
  print RMD "rm $parent/$directory_id/*.fq.gz\n";
  print RMD "rsync -a $parent/$directory_id storage07:/data/data/rnaseq/\n";
  print RMD "rm -r $parent/$directory_id/\n";
  close RMD;
  
  # make these rm shell scripts executable
  system( "chmod u+x $rm_reads" );
  system( "chmod u+x $parent/$rm_outdir" );

}


# exit the program
exit;


############################################SUBROUTINES#############################################
sub sql_fetch{

  my $rnaNUM = shift;
  
  # for connecting to the MySQL server
  my $ds     = "DBI:mysql:$DB:$HOST";
  my $user   = 'aja';
  my $passwd = 'aja123';
  
  # connect to the MySQL server
  my $dbh = DBI->connect( $ds, $user, $passwd, {
	  AutoCommit => 0,
	  RaiseError => 1,
	  }) or die $DBI::errstr;

  my $sth;  # holds the SQL prepare statement

  # prepare and execute the SQL SELECT statement for this RNA_ID
  $sth = $dbh->prepare( "SELECT * FROM $TABLE WHERE RNA_ID = $rnaNUM" );
  $sth->execute();
  
  # fetch the SQL entry for this RNA_ID into a hash
  my $hash = $sth->fetchrow_hashref();
  
  # finish the query
  $sth->finish;
  
  # build directory name (e.g., 'gyrosigma_fasciola_ecr2aja-02_R19' )
  my $dir_id = lcfirst $hash->{Genus} . "_" . lcfirst $hash->{Species} . "_" .  $hash->{Culture_ID} . "_R" . $rnaNUM . $suffix;

  # build name for log file name
  my $logname = "clean_filter_$rnaID$suffix.log";

  # disconnect from the MySQL database
  $dbh->disconnect;

  return( $hash->{Genus}, $hash->{Species}, $hash->{Culture_ID}, $logname, $dir_id );
  
}
####################################################################################################
sub run_Trimmomatic{

  # this subroutine uses Trimmomatic to trim reads and remove adapter sequences

  my ( $rnaID, $suffix, $logfile, $fwd_in, $rev_in ) = @_;

  # establish name and location of Illumina adapter database
  my $adapterDB = "/scratch/aja/illumina_adapter_database/TruSeq_adapters.fa";
  
  # 1 - remove Illumina adapters (Trimmomatic: ILLUMINACLIP)
  # 2 - trim 5' $HEADCROP (10) bp from each read (Trimmomatic: HEADCROP)
  # 3 - trim 5' and 3' bases with Phred <= $TRIM_PHRED (Trimmomatic: LEADING, TRAILING)
  # 4 - remove reads with avg Phred < 32 (Trimmomatic: AVGQUAL)
  # 5 - remove reads <= 72 bp in length (Trimmomatic: MINLEN)

  my $fwd_out = "$rnaID$suffix\_cleaned_1.fq.gz"; # (fwd)
  my $rev_out = "$rnaID$suffix\_cleaned_2.fq.gz"; # (rev)

  my $trimmomatic_command = "java -jar /share/apps/trimmomatic/Trimmomatic-0.32/trimmomatic-0.32.jar PE -phred64 $fwd_in $rev_in $fwd_out junk_$rnaID\_1.fq.gz $rev_out junk_$rnaID\_2.fq.gz ILLUMINACLIP:$adapterDB:2:40:15 HEADCROP:$HEADCROP LEADING:$TRIM_PHRED TRAILING:$TRIM_PHRED AVGQUAL:$AVG_PHRED MINLEN:$MIN_LEN TOPHRED33 2>>$logfile";


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
  system( "rm junk_$rnaID\_1.fq.gz junk_$rnaID\_2.fq.gz" );
  
  return( $fwd_out, $rev_out );
  
}
####################################################################################################
sub filter_univec{

  # this subroutine runs Bowtie2 to filter out vector contaminants

  my ( $rnaID, $suffix, $logfile, $fwd_in, $rev_in, $bowtie2_executable, $bowtie2_db ) = @_;

  # establish base name for unmapped reads output
  my $unmapped_basename = "$rnaID$suffix\_noVec_%.fq.gz";
  ( my $fwd_out = $unmapped_basename ) =~ s/%/1/;
  ( my $rev_out = $unmapped_basename ) =~ s/%/2/;
  
  # bowtie univec command
  my $bowtie_univec_cmd = "$bowtie2_executable -x $bowtie2_db -1 $fwd_in -2 $rev_in --phred33 --local --un-conc-gz $unmapped_basename >/dev/null 2>>$logfile";

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

  my ( $rnaID, $suffix, $logfile, $fwd_in, $rev_in, $bowtie2_executable, $bowtie2_db ) = @_;

  # establish base name for unmapped reads output
  my $unmapped_basename = "$rnaID$suffix\_noVec_noRRNA_%.fq.gz";
  ( my $fwd_out = $unmapped_basename ) =~ s/%/1/;
  ( my $rev_out = $unmapped_basename ) =~ s/%/2/;
  
  # filename for mapped rRNA reads, output by Bowtie2
  my $rrna_reads = "$rnaID$suffix\_rrna_reads_%.fq.gz";
  ( my $fwd_rrna = $rrna_reads ) =~ s/%/1/;
  ( my $rev_rrna = $rrna_reads ) =~ s/%/2/;

  if( $PBS_ONLY ){
    return( $fwd_rrna, $rev_rrna );
  }else{

    # bowtie univec command
    my $bowtie_rrna_cmd = "$bowtie2_executable -x $bowtie2_db -1 $fwd_in -2 $rev_in --phred33 --local --un-conc-gz $unmapped_basename --al-conc-gz $rrna_reads >/dev/null 2>>$logfile";

    # open log file for printing
    open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
    print LOGFILE "# ----------------------------------------------------------------------\n";
    print LOGFILE "# Filtering rRNA reads with Bowtie2\n";
    print LOGFILE "# ----------------------------------------------------------------------\n\n";
    print LOGFILE "Bowtie2 run with: \'$bowtie_rrna_cmd\'\n\n";
    print LOGFILE "BEGIN BOWTIE2 RRNA;\n";
    close LOGFILE;
    
    my $start_time = new Benchmark;
    
    system( $bowtie_rrna_cmd );
    
    my $end_time   = new Benchmark;
    my $difference = timediff( $end_time, $start_time );
    
    # re-open log file for printing
    open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
    print LOGFILE "\nBowtie2 rRNA filter CPU: ", timestr( $difference ), "\n";
    print LOGFILE "Output files: $fwd_out, $rev_out, $fwd_rrna, $rev_rrna\n";
    print LOGFILE "END BOWTIE2 RRNA;\n\n\n\n";
    close LOGFILE;
  }

  return( $fwd_out, $rev_out, $fwd_rrna, $rev_rrna );

}
####################################################################################################
sub filter_organelle{

  # this subroutine runs Bowtie2 to filter out organelle reads

  my ( $rnaID, $suffix, $logfile, $fwd_in, $rev_in, $bowtie2_executable, $bowtie2_db ) = @_;

  # establish names for unmapped reads output
  my $unmapped_basename = "$rnaID$suffix\_cleaned_filtered_%.fq.gz";
  ( my $fwd_out = $unmapped_basename ) =~ s/%/1/;
  ( my $rev_out = $unmapped_basename ) =~ s/%/2/;
  
  # filename for mapped organelle reads, output by Bowtie2
  my $organelle_reads = "$rnaID$suffix\_organelle_reads_%.fq.gz";
  ( my $fwd_organelle = $organelle_reads ) =~ s/%/1/;
  ( my $rev_organelle = $organelle_reads ) =~ s/%/2/;

  if( $PBS_ONLY ){
    return( $fwd_out, $rev_out, $fwd_organelle, $rev_organelle );
  }else{

    # bowtie univec command
    my $bowtie_organelle_cmd = "$bowtie2_executable -x $bowtie2_db -1 $fwd_in -2 $rev_in --phred33 --local --un-conc-gz $unmapped_basename --al-conc-gz $organelle_reads >/dev/null 2>>$logfile";

    # open log file for printing
    open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
    print LOGFILE "# ----------------------------------------------------------------------\n";
    print LOGFILE "# Filtering organelle reads with Bowtie2\n";
    print LOGFILE "# ----------------------------------------------------------------------\n\n";
    print LOGFILE "Bowtie2 run with: \'$bowtie_organelle_cmd'\n\n";
    print LOGFILE "BEGIN BOWTIE2 ORGANELLE;\n";
    close LOGFILE;
    
    my $start_time = new Benchmark;
    
    system( $bowtie_organelle_cmd );
    
    my $end_time   = new Benchmark;
    my $difference = timediff( $end_time, $start_time );
    
    # re-open log file for printing
    open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
    print LOGFILE "\nBowtie2 organelle filter CPU: ", timestr( $difference ), "\n";
    print LOGFILE "Output files: $fwd_out, $rev_out, $fwd_organelle, $rev_organelle\n";
    print LOGFILE "END BOWTIE2 ORGANELLE;\n\n\n\n";
    close LOGFILE;
  }
  
    return( $fwd_out, $rev_out, $fwd_organelle, $rev_organelle );
  
}
####################################################################################################
sub merge_reads_with_BBMerge{

  # this subroutine runs BBMerge to merge paired-end reads

  my ( $fwd_in, $rev_in, $logfile, $type ) = @_;

  # BBMerge command
  my $bbmerge = "bbmerge.sh strict=t in1=$fwd_in in2=$rev_in out=merged.fq outu1=unmerged1.fq outu2=unmerged2.fq.gz 2>>$logfile";
  
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
  system( "cat unmerged1.fq merged.fq | gzip > $fwd_in" );

  # rename rev unmerged reads
  system( "mv unmerged2.fq.gz $rev_in" );

  # remove intermediate files
  system( "rm unmerged1.fq merged.fq" );

  my $end_time   = new Benchmark;
  my $difference = timediff( $end_time, $start_time );
  
  # re-open log file for printing
  open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
  print LOGFILE "\nBBMerge $type CPU: ", timestr( $difference ), "\n";
  print LOGFILE "File actions:\n";
  print LOGFILE "     cat unmerged1.fq merged.fq | gzip > $fwd_in\n";
  print LOGFILE "     mv unmerged2.fq.gz $rev_in\n";
  print LOGFILE "     rm unmerged1.fq merged.fq\n";
  print LOGFILE "END ", uc $type, " BBMERGE;\n\n\n\n";
  close LOGFILE;

}
####################################################################################################
sub normalize_with_BBNorm{

  # this subroutine runs BBNorm to perform kmer normalization of the reads

  my ( $fwd_in, $rev_in, $logfile ) = @_;
  my $version;
  ( my $fwd_out = $fwd_in ) =~ s/_1/_normalized_1/;
  ( my $rev_out = $rev_in ) =~ s/_2/_normalized_2/;
  
  if( $PBS_ONLY ){
    return( $fwd_out, $rev_out );
  }else{

    # get BBNorm "version"
    open( V, "bbnorm.sh | " );
    while( <V> ){ 
      chomp;
      /Last modified/ and $version .= $_;
    }
    close V;

    # BBNorm command
    # my $bbnorm = "bbnorm.sh target=50 in=$fwd_in in2=$rev_in out=$fwd_out out2=$rev_out 2>>/dev/null";
    my $bbnorm = "bbnorm.sh target=50 in=$fwd_in in2=$rev_in out=$fwd_out out2=$rev_out 2>>$logfile";
    
    # open log file for printing
    open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
    print LOGFILE "# ----------------------------------------------------------------------\n";
    print LOGFILE "# Normalizing nuclear reads with BBNorm (ver. \"$version\")\n";
    print LOGFILE "# ----------------------------------------------------------------------\n\n";
    print LOGFILE "BBNORM run with: \'$bbnorm\'\n\n";
    print LOGFILE "BEGIN BBNORM;\n";
    close LOGFILE;
    
    my $start_time = new Benchmark;
    
    system( $bbnorm );
    
    my $end_time   = new Benchmark;
    my $difference = timediff( $end_time, $start_time );
    
    # re-open log file for printing
    open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
    print LOGFILE "\nBBNorm CPU: ", timestr( $difference ), "\n";
    print LOGFILE "Output files: $fwd_out, $rev_out\n";
    print LOGFILE "END BBNORM;\n\n\n\n";
    close LOGFILE;
  }
  
  return( $fwd_out, $rev_out );

}
####################################################################################################
sub writeTrinityPBS{

  my( $logfile, $fwd_nuc_norm, $rev_nuc_norm, $fwd_nuc_all, $rev_nuc_all, $fwd_rna, $rev_rna, $fwd_org, $rev_org, $rnaID, $suffix, $parent, $directory_id, $trinity_executable, $trinity_plugins, $trinity_utils, $hmmer3_executable, $pfam_db );
  
  if( $NORMALIZE ){
    ( $logfile, $fwd_nuc_norm, $rev_nuc_norm, $fwd_nuc_all, $rev_nuc_all, $fwd_rna, $rev_rna, $fwd_org, $rev_org, $rnaID, $suffix, $parent, $directory_id, $trinity_executable, $trinity_plugins, $trinity_utils, $hmmer3_executable, $pfam_db ) = @_;
  }else{
    ( $logfile, $fwd_nuc_all, $rev_nuc_all, $fwd_rna, $rev_rna, $fwd_org, $rev_org, $rnaID, $suffix, $parent, $directory_id, $trinity_executable, $trinity_plugins, $trinity_utils, $hmmer3_executable, $pfam_db ) = @_;
  }
  
  my $ppn; # number processors
  my $mem; # available memory
  my $outfile = "trinity_assemble_$rnaID$suffix.pbs";

  # files names for the unzipped nuclear reads
  ( my $fwd_fq = $fwd_nuc_all ) =~ s/\.gz$//;
  ( my $rev_fq = $rev_nuc_all ) =~ s/\.gz$//;
  
  # name of the Trinity assembly FASTA file when '--output' is not specified
  my $trinity_assembly_file = 'trinity_out_dir.Trinity.fasta';

  # for selecting a random queue
  my @queue_list = ( 'mem512GB64core', 'mem768GB32core');

  # choose the queue if random
  $QUEUE eq 'random' and $QUEUE = $queue_list[rand @queue_list];

  if( $QUEUE eq 'aja' ){
    $ppn = 16;
    $mem = 32;
  }elsif( $QUEUE eq 'mem512GB64core' ){
    $ppn = 64;
    $mem = 512;
  }elsif( $QUEUE eq 'mem768GB32core' ){
    $ppn = 32;
    $mem = 768;
  }else{
    die "Queue \'$QUEUE\' not recognized\n";
  }

  open( PBS, '>', $outfile ) || die "Can't open $outfile: $!\n";
  
  # print PBS job parameters
  print PBS "#PBS -N assemble_$rnaID$suffix\n";
  print PBS '#PBS -j oe', "\n";
  print PBS '#PBS -m abe', "\n";
  print PBS '#PBS -M aja@uark.edu', "\n";
  print PBS '#PBS -q ' . $QUEUE . "\n";
  print PBS '#PBS -l nodes=1:ppn='. $ppn . "\n";
  print PBS '#PBS -l walltime=' . $WALLTIME . ':00:00', "\n";
  print PBS "#PBS -o trinity_assemble_$rnaID$suffix\." . '$PBS_JOBID' . "\n\n";
  print PBS 'cd $PBS_O_WORKDIR', "\n\n";

  # load Trinity modules
  print PBS "# load Trinity modules\n";
  print PBS 'module purge', "\n";
  print PBS 'module load bowtie', "\n";
  print PBS 'module load samtools/0.1.19', "\n";
  print PBS 'module load trinity/2.0.2', "\n";
  print PBS 'module load blast/2.2.29+', "\n";
  print PBS 'module load transdecoder', "\n\n";

  # assemble rRNA reads
  print PBS "# assemble rRNA reads\n";
  print PBS "$trinity_executable --seqType fq --min_kmer_cov $MIN_KMER_RNA --max_memory ", $mem-6, 'G ', '--CPU ', $ppn-2, " --left $fwd_rna --right $rev_rna --full_cleanup\n";
  print PBS "mv $trinity_assembly_file $rnaID$suffix\_rrna.fa\n\n";

  # assemble organelle reads
  print PBS "# assemble organelle reads\n";
  print PBS "$trinity_executable --seqType fq --min_kmer_cov $MIN_KMER_ORG --max_memory ", $mem-6, 'G ', '--CPU ', $ppn-2, " --left $fwd_org --right $rev_org --full_cleanup\n";
  print PBS "mv $trinity_assembly_file $rnaID$suffix\_organelle.fa\n\n";
  
  # assemble nuclear reads
  print PBS "# assemble nuclear reads\n";
  if( $NORMALIZE ){
    print PBS "$trinity_executable --seqType fq --min_kmer_cov $MIN_KMER_NUC --max_memory ", $mem-6, 'G ', '--CPU ', $ppn-2, " --left $fwd_nuc_norm --right $rev_nuc_norm --full_cleanup\n";
  }else{
    print PBS "$trinity_executable --seqType fq --min_kmer_cov $MIN_KMER_NUC --max_memory ", $mem-6, 'G ', '--CPU ', $ppn-2, " --left $fwd_nuc_all --right $rev_nuc_all --full_cleanup\n";
  }
  print PBS "mv $trinity_assembly_file $rnaID$suffix\_nuclear.fa\n\n";
  
  # print report for nuclear assembly
  print PBS "# print assembly report\n";
  print PBS "$trinity_utils/TrinityStats.pl $rnaID$suffix\_nuclear.fa > $rnaID$suffix\_nuclear.stats\n\n";

  # estimate transcript abundance with RSEM
  print PBS "# estimate nuclear transcript abundance with RSEM\n";
  print PBS "$trinity_utils/align_and_estimate_abundance.pl --thread_count ", $ppn-2, " --transcripts $rnaID$suffix\_nuclear.fa --seqType fq --left $fwd_fq --right $rev_fq --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference\n\n";
  
  # BLAST and hmmer3 searches to diatom proteins and Pfam databases, respectively
  # NEED TO INCORPORATE GNU PARALLEL OR GRIDRUNNER
  # NEED TO REMOVE BLAST AND PFAM OUTPUT
  # print PBS "blastp -query transdecoder_dir/longest_orfs.pep -db uniprot_sprot.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 > blastp.outfmt6\n";
  # print PBS "$hmmer3_executable --cpu 8 --domtblout pfam.domtblout $pfam_db transdecoder_dir/longest_orfs.pep\n";
  
  # use TransDecoder to translate nuclear reads
  #print PBS "# use TransDecoder to translate nuclear reads\n";
  #print PBS "TransDecoder.LongOrfs -t $rnaID$suffix\_nuclear.fa\n\n";
  #print PBS "TransDecoder.Predict -t target_transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6\n\n";
  
  # remove intermediate output files
  print PBS "# remove intermediate output files\n";
  print PBS "rm -r *.fq *.readcount $rnaID$suffix\_nuclear.fa.* RSEM.stat RSEM.isoforms.results.ok transdecoder.tmp.* bowtie.*\n";

  # rsync output from razor to storage07:/data/data/rnaseq
  # Reminder: $parent = '/scratch/aja/rnaseq';
  # print LOGFILE "Copying output directory to \'$copy_output_to\'q:\n";
  # print LOGFILE "    \'rsync -av $parent/$directory_id/ $copy_output_to/$directory_id\'\n\n";
  # system( "rsync -av $parent/$directory_id/ $copy_output_to/$directory_id\n" );
  
  close PBS;

  # if PBS only, skip writing out to log file
  $PBS_ONLY and return $outfile;

  # get Trinity version
  open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
  print LOGFILE "# ----------------------------------------------------------------------\n";
  print LOGFILE "# Trinity version and PBS job script\n";
  print LOGFILE "# ----------------------------------------------------------------------\n\n";
  close LOGFILE;
  system ( "$trinity_executable --version >> $logfile" );
  
  # identify Trinity PBS job script
  open( LOGFILE, ">>", $logfile  ) || die "Can't open $logfile: $!\n";
  print LOGFILE "Trinity PBS job script written to: \'$outfile\'\n\n\n";
  close LOGFILE;

  return $outfile;
  
}
####################################################################################################
sub parseArgs{

  my $usage = "\nUsage: $0 <ID> [options]

   Required: (1) Alverson lab RNA ID (e.g. '12' or '12a'), *OR*
             (2) base name for the read files (e.g. 'index-AGTTCC' for 'index-AGTTCC_1.fq.gz' and 'index-AGTTCC_2.fq.gz')


   General options
          --queue     - job queue ('mem256GB48core', 'mem512GB64core' [default], 'mem768GB32core', 'aja', 'random')
          --db        - MySQL database (default: 'transcriptome')
          --table     - MySQL table to read (default: 'RNA_Book')
          --walltime  - walltime for Trinity run (integer [hrs], default: 12)
          --normalize - perform kmer normalization of *nuclear* reads with BBNorm (default: false)
          --pbs_only  - generate a Trinity job script only, skipping Trimmomatic and Bowtie2 (default: false)


   Trimmomatic options
          --headcrop   - number of 5' bases to trim (default: 10)
          --trim_phred - bases below this minimum Phred quality will be trimmed from 5' and 3' ends (default: 5)
          --avg_phred  - filter reads if average quality is below this threshold (default: 32)
          --min_len    - filter reads shorter than this length threshold (default: 72)


   Trinity options
          --min_kmer_rna  - Trinity min_kmer_cov parameter for rRNA assembly (default: 1)
          --min_kmer_org  - Trinity min_kmer_cov parameter for organelle assembly (default: 1)
          --min_kmer_nuc  - Trinity min_kmer_cov parameter for nuclear assembly (default: 2)

\n";
  
  my $result = GetOptions	(
	 'queue=s'      => \$QUEUE,
	 'db=s'         => \$DB,
	 'table=s'      => \$TABLE,
	 'walltime=s'   => \$WALLTIME,
	 'normalize!'   => \$NORMALIZE,
	 'pbs_only'     => \$PBS_ONLY,

	 'headcrop=s'   => \$HEADCROP,
	 'trim_phred=s' => \$TRIM_PHRED,
	 'avg_phred=s'  => \$AVG_PHRED,
	 'min_len=s'    => \$MIN_LEN,

	 'min_kmer_rna=i'  => \$MIN_KMER_RNA,
	 'min_kmer_org=i'  => \$MIN_KMER_ORG,
	 'min_kmer_nuc=i'  => \$MIN_KMER_NUC,
				);
  
  $ARGV[0] or die $usage;

}
####################################################################################################
