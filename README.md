### Perl and Python scripts for processing RNAseq data

`compare_md5sum.pl`: compares two files with pre- and post-download md5sum values

`write_rnaseq_clean_filter_jobscript.pl`: creates PBS job script for running `rnaseq_clean_filter.pl`

`rnaseq_clean_filter.pl`: pipeline used by Alverson lab for the following:
- correcting reads with rcorrector
- trimming and cleaning reads with Trimmomatic
- filtering vector contaminants
- filtering diatom rRNA reads
- filtering bacterial rRNA reads
- merging overlapping reads
- generating a PBS job script for assembling nuclear reads with Trinity  
- generating a PBS job script for translating the assembled contigs with Transdecoder
>This script is set up to run on the [AHPCC](http://hpc.uark.edu/hpc/) and interface with the MySQL databases that store our sample information ; a typical set of commands is listed in `write_rnaseq_clean_filter_jobscript.md`

`filter_crosstalk.py`: use CD-HIT to remove obvious index crosstalk in samples that were multiplexed and sequenced together

`parse_rnaseq_clean_filter_log.pl`: parse `rnaseq_clean_filter.pl` log file and various assembly files for read characteristics and assembly statistics

`pickH_from_rsem.pl`: pick isoform with highest expression level
