### Perl scripts for processing RNAseq data

- `write_rnaseq_clean_filter_jobscript.pl`: creates PBS job script for running `rnaseq_clean_filter.pl`
- `rnaseq_clean_filter.pl`: pipeline used by Alverson lab for the following:
 - correcting reads with rcorrector
 - trimming and cleaning reads with Trimmomatic
 - filtering vector contaminants
 - filtering rRNA reads
 - merging overlapping reads
 - generating a PBS job script for assembling nuclear reads with Trinity  
 - generating a PBS job script for translating the assembled contigs with Transdecoder

>Note: although parts of it can be borrowed, this script is set up to run on the [AHPCC](http://hpc.uark.edu/hpc/) and interface with the MySQL databases that store our sample information ; a typical set of commands is listed in `write_rnaseq_clean_filter_jobscript.md`

