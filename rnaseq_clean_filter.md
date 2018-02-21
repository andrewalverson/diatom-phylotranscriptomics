#Cleaning, filtering, and assembling Illumina RNAseq data with ``rnaseq_clean_filter.pl``

##Cleaning and filtering reads

## Notes and Prerequisites
- The script `rnaseq_clean_filter.pl` uses _rcorrector_ to perform error correction, then it uses _Trimmomatic_ and _Bowtie2_ to trim reads, remove Illumina adapters, and filter out reads matching known sequencing vectors, diatom rRNA sequences, and contaminant rRNA sequences
- The script uses _BBMerge_ to merge overlapping reads prior to assembly
- File names for raw reads **must** be formatted as, for example, `sampleID_1.fq.gz` (fwd) and `sampleID_2.fq.gz` (rev)
>Note: for Alverson lab samples, the sample ID will correspond to the RNA or library ID
>Note: if an RNA is re-sequenced, the later read files are named as `R8a_1.fq.gz`, `R8b_1.fq.gz`

- Make sure that the MySQL databases on *Razor* are up to date, so the sample ID can be used to pull information from lab databases
>**INSTRUCTIONS FORTHCOMING**

- All analyses are run on the AHPCC's *Razor* cluster

####Create PBS job script for running `rnaseq_clean_filter.pl`

- Use `write_rnaseq_clean_filter_jobscript.pl` to create a PBS job script for running the cleaning and filtering job, which can take 18--20 hrs for large datasets

- Usage:

	`write_rnaseq_clean_filter_jobscript.pl <RNA ID> [options]`

- Output:

    `run_rnaseq2clean_filter_[RL]\d+.pbs`

- Double-check the PBS script to make sure `rnaseq_clean_filter.pl` has all the options you want

##Clean and filter reads
- The PBS job script just created runs `rnaseq_clean_filter.pl`, which does the following, in order:

1. Trimmomatic

    java -jar trimmomatic-0.32.jar PE -phred64 <R#_1.fq.gz> <R#_2.fq.gz> <R#_cleaned_1.fq.gz> junk1.fq.gz <R#_cleaned_2.fq.gz> junk2.fq.gz ILLUMINACLIP:TruSeq_adapters.fa:2:40:15 HEADCROP:10 LEADING:5 TRAILING:5 AVGQUAL:32 MINLEN:72 TOPHRED33 2>>clean_filter_R#.log

Trimmomatic options:

 - `PE -phred64`: inputting paired-end reads with offset-64 Phred scores
 - `ILLUMINACLIP:TruSeq_adapters.fa:2:40:15`: clip Illumina adapters, as specified in 'Truseq_adapters.fa'
 - `HEADCROP:10`: trim 10 bp from 5' end of each read
 - `LEADING:5`: trim 5' bases with quality <= 5
 - `TRAILING:5`: trim 5' bases with quality <= 5
 - `AVGQUAL:32`: remove reads with avg Phred < 32
 - `MINLEN:72`: remove reads <= 72 bp in length
 - `TOPHRED33`: output offset-33 Phred scores

Input files are raw sequencing reads, e.g.:

 - `$fwd_in`: R8_1.fq.gz
 - `$rev_in`: R8_2.fq.gz

Output files are clipped and cleaned reads:

 - **R8_cleaned_1.fq.gz**  and  **R8_cleaned_2.fq.gz**: clipped and cleaned reads

1. Remove vector contaminants
`bowtie2 -x $bowtie_univec_db -1 $fwd_in -2 $rev_in --phred33 --local --un-conc-gz $unmapped_basename >/dev/null 2>>$logfile`

1. Remove rRNA sequences
`bowtie2 -x $bowtie_rrna_db -1 $fwd_in -2 $rev_in --phred33 --local --un-conc-gz $unmapped_basename --al-conc-gz $rrna_reads >/dev/null 2>>$logfile`

1. Merge overlapping reads with BBMerge
`bbmerge.sh strict=t in1=$fwd_in in2=$rev_in out=merged.fq outu1=unmerged1.fq outu2=unmerged2.fq.gz 2>>$logfile`

1. Write Trinity PBS job script


###Assembling the cleaned and filtered reads
    Run the *.pbs Trinity job script output by ``rnaseq_clean_filter.pl``

> Written with [StackEdit](https://stackedit.io/).
