## Perl and Python scripts for processing RNAseq data

### Data QC
`compare_md5sum.pl`: compares two files with pre- and post-download md5sum values


### Transcriptome assembly
`write_rnaseq_clean_filter_jobscript.pl`: creates PBS job script for running `rnaseq_clean_filter.pl`

`rnaseq_clean_filter.pl`: pipeline used by Alverson lab for the following:
- correcting reads with rcorrector
- trimming and cleaning reads with Trimmomatic
- filtering vector contaminants
- filtering diatom rRNA reads
- filtering bacterial rRNA reads
- merging overlapping reads
- generating a PBS job script for assembling nuclear reads with Trinity  
- generating a PBS job script for translating assembled contigs with Transdecoder
- assembly QC with BUSCO
>Note: this script is set up to run on the [AHPCC](http://hpc.uark.edu/hpc/) and interface with the MySQL databases that store our sample information ; a typical set of commands is listed in `write_rnaseq_clean_filter_jobscript.md`

`parse_rnaseq_clean_filter_log.pl`: parses log files of `rnaseq_clean_filter.pl` for data and assembly statistics


### Filtering index crosstalk
`remove_croco-filtered_from_transdecoder.py`: removes CroCo-identified crosstalk (FASTA format) from Transdecoder FASTA file

`cat_croco_dubious_contam.py`: concatenate sequences flagged as contaminant and dubious into a single file

`filter_crosstalk.py`: use CD-HIT to remove obvious index crosstalk in samples that were multiplexed and sequenced on the same Illumina lane


### Orthofinder scripts
`summarize_orthogroup_membership.py`: parses `Orthogroups.GeneCount.csv` file to summarize ingroup/outgroup taxon occupancy in OrthoFinder-based orthogroups

`link_to_orthogroups_by_taxon_occupancy.py`: parse output of `summarize_orthogroup_membership.py` to create symbolic links to FASTA files of orthogroups with a minimum specified taxon occupancy

`orthogroups_to_fasta.py`: creates FASTA files from OrthoFinder's 'Orthogroups.csv' output file


### Transdecoder scripts
`select_longest_orf_from_transdecoder.py`

`select_longest_orf_from_transdecoder_longest_orfs.py`


### Miscellaneous
`clean_trinity_fasta_header.pl`: remove special character, etc. from Trinity FASTA headers
embl_to_cds.py

`pickH_from_rsem.pl`: pick isoform with highest expression level
