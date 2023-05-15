# Candida_genome_visualization
Simple genome-level visualization of CNV and SNP/LOH data from alignment files

Bash and R scripts that can be run on a cluster or locally. Intended workflow is:
candida_gc_correct.sh (optional) -> candida_ymap.sh (which calls berman_count_snps_v5.py and genome_vis.R)
genome_vis.R can easily be run interactively in RStudio.

Input files are BAM alignments. Intermediate outputs are depth and putative SNP tab-delimited files used in R to produce a linear chromosome view of copy number and LOH density. Final output from R includes image files and underlying dataframes as excel spreadsheets.
