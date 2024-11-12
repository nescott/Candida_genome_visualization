# Candida_genome_visualization

Genome-level visualization of CNV and SNP/LOH data from BAM alignment files.

To plot CNV and LOH data for individual samples across each chromosome, the intended workflow is: candida_gc_correct.sh (optional) -\> candida_ymap.sh (which calls berman_count_snps_v5.py and genome_vis.R).

A heatmap-style plot of LOH patterns can be generated for multiple samples using the intermediate files generated from plotting individual samples.