
## Scripts for analyzing diatom phylotranscriptomic dataset

### R scripts

- `check-monophyly.R`: check gene trees for monophyly of
  Thalassiosirales and raphid pennates

- `Densi.mod.R`: implementation of DensiTree-style function to draw a
  cloud of diversification histories across a set of phylogenetic trees

- `Part2_plotMultiMisse.R`: plot the DensiTree-style results

- `plot-busco-omark.R`: plot BUSCO and OMArk completeness

- `trees-ordination-and-consensus.R`: NMDS ordination of species trees

### Data QC

`compare_md5sum.pl`: compares two files with pre- and post-download
md5sum values

### Transcriptome assembly

`write_rnaseq_clean_filter_jobscript.pl`: creates PBS job script for
running `rnaseq_clean_filter.pl`

`rnaseq_clean_filter.pl`: pipeline for:

- correcting reads with rcorrector
- trimming and cleaning reads with Trimmomatic
- filtering vector contaminants
- filtering diatom rRNA reads
- filtering bacterial rRNA reads
- merging overlapping reads
- generating a PBS job script for assembling nuclear reads with
  Trinity  
- generating a PBS job script for translating assembled contigs with
  Transdecoder
- assembly QC with BUSCO and transrate

**Note**: this script is set up to run on the
[AHPCC](http://hpc.uark.edu/hpc/) and interface with lab-specific SQL
databases that store our sample information

### Filtering index crosstalk identified by [CroCo](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-018-0486-7)

`remove_croco-filtered_from_transdecoder.py`: removes CroCo-identified
crosstalk (FASTA format) from Transdecoder FASTA file

`cat_croco_dubious_contam.py`: concatenate sequences flagged as
contaminant and dubious into a single file

`filter_crosstalk.py`: use CD-HIT to remove obvious index crosstalk in
samples that were multiplexed and sequenced on the same Illumina lane

### Orthofinder scripts

`summarize_orthogroup_membership.py`: parses `Orthogroups.GeneCount.csv`
file to summarize ingroup/outgroup taxon occupancy in OrthoFinder-based
orthogroups

`link_to_orthogroups_by_taxon_occupancy.py`: parse output of

`summarize_orthogroup_membership.py` to create symbolic links to FASTA
files of orthogroups with a minimum specified taxon occupancy

`orthogroups_to_fasta.py`: creates FASTA files from OrthoFinder’s
‘Orthogroups.csv’ output file

### Miscellaneous

`alignment_stats_to_file_bins.py`: bin alignments by various criteria

`check_outgroups.py`: check whether multi-species FASTA file includes
specified outgroup taxa

`clean_trinity_fasta_header.pl`: remove special characters, etc. from
Trinity FASTA headers

`count_ortholog_types.py`

`extract_busco_summary_from_batch_outfile.py`: extract key values from
batch BUSCO output into a csv file for plotting in R

`link_to_orthogroups_by_taxon_occupancy.py`

`standardize_fasta_headers_and_filenames.pl`: standardize FASTA headers
among different output files

`standardize_mmetsp_headers_and_filenames.pl`: standardize FASTA headers
among different output files for MMETSP data

`summarize_treeshrink.py`: summarize tips removed by
[TreeShrink](https://github.com/uym2/TreeShrink)

`tip-label2taxon-label.py`: parses the taxon name from FASTA headers,
following the standardized format of this dataset

`trimal_AA_to_CDS_coords.py`: convert amino acid coordinates from
trimAl’s \#ColumnsMap output to CDS coordinates for identical trimming
of corresponding CDS alignments

### Write SLURM and PBS scripts for different analyses

- `write_IQ-TREE_AA_tree_PBS_scripts.py`

- `write_mafft_pbs.py`

- `write_razor_transrate.py`

- `write_IQ-TREE_AA_tree_slurm_scripts.py`

- `write_mafft_slurm.py`

- `write_rnaseq_clean_filter_jobscript.pl`

- `write_IQ-TREE_CDS_tree_PBS_scripts.py`

- `write_pinnacle_corset_filtering.py`

- `write_treesearch_homolog_slurm.py`

- `write_IQ-TREE_CDS_tree_slurm_scripts.py`

- `write_pinnacle_hmmer_slurm.py`

- `write_treesearch_pbs.py`

- `write_IQ-TREE_ortholog_AA.py`

- `write_pinnacle_transdecoder_slurm.py`

- `write_treesearch_slurm.py`

- `write_UPP_PBS_scripts.py`

- `write_raxml_AA_prelim_tree_PBS.py`

- `write_trimal_gt_commands.py`

- `write_UPP_bash_scripts.py`

- `write_razor_hmmer.py`

- `write_trimal_with_selectcols.py`

- `write_UPP_slurm.py`

- `write_razor_transdecoder.py`
