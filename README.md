### Perl scripts for processing RNAseq data

- `write_rnaseq_clean_filter_jobscript.pl`: creates PBS job script for running `rnaseq_clean_filter.pl`
- `rnaseq_clean_filter.pl`: pipeline used by Alverson lab for
 - trimming and cleaning reads
 - filtering rRNA and organelle reads
 - normalizing kmers (optional)
 - generating a PBS job script for assembling nuclear reads with Trinity  

>Note: although parts of it certainly can be borrowed, this script is highly customized to interface with the MySQL databases that store our sample information as well as run on the [AHPCC](http://hpc.uark.edu/hpc/); a typical set of commands is listed in `write_rnaseq_clean_filter_jobscript.md`

- `trinity_assembly_counts.pl`: count the number of unique genes and transcripts in Trinity assembly
