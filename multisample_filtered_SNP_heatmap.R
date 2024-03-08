## ---------------------------
## Purpose: Generate a heatmap of all SNPs for a group of isolates, across all chromosomes
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
# Input files

genotype_table <-"align_variants_SC5314_A21_ref/clustering/Calbicans_MEC_bwa_genotypes.txt"
sample_list <- "~/umn/data/metadata/Calbicans_MEC_raxml_midpoint_tips.csv"
feature_file <- "path/to/features.txt"
label_file <- "~/umn/Candida_genome_visualization/ref_genome_files/Calbicans_SC5314_A21_chr_labels.txt"

# Load packages
library(tidyverse)
library(ggplot2)
library(writexl)
library(readxl)

# Set variables
window <- 5000 # size of window used for rolling mean and snp density

save_dir <- "~/umn/images/Calbicans/"

ref <- "sc5314" # label for file name or "" to leave out

mito <- "Ca19-mtDNA"

sample_order <- read_csv(sample_list, show_col_types = FALSE)

# SNP colors, plot function uses 2-color gradient scale
snp_low <- "white"
snp_high <- "dodgerblue3"

# Color of chromosome outlines
chrom_outline_color <- "gray26"

# Line width of chromosome outlines
chrom_line_width <- 0.2

# For overwriting scaffold names in final plot
chr_ids <- scan(label_file, what = character())

# Read in genotypes (generated from filtered VCF)
genotypes <- read.table(genotype_table,
                        header = TRUE,
                        sep="\t")

gt <- as.data.frame(genotypes)
gt <- gt %>% filter(CHROM != mito)
gt <- gt %>% mutate(snp_bin=(POS %/% window) * window +1)

snp_counts <- gt %>%
  group_by(CHROM, snp_bin) %>%
  summarise_at(vars(AMS5231:MEC374), function(x) sum(x!="0/0"))

snp_counts$x_num <- seq.int(nrow(snp_counts))

chrs <- snp_counts %>%
  group_by(CHROM) %>%
  summarise(border_start=min(x_num),
            border_stop=max(x_num),
            tick=min(x_num) + (max(x_num)-min(x_num))/2)

snp_again <- snp_counts %>%
  pivot_longer(names_to = "sample", values_to = "snp_count", cols=-c(CHROM,snp_bin, x_num))

################################################################################
# Plot heatmap
p <- snp_again %>%
  mutate(clustered_samples = fct_relevel(sample, rev(sample_order$sample))) %>%
  ggplot(aes(x=x_num, y=clustered_samples, fill=snp_count))+
  #scale_fill_gradient2(low = "white", mid="darkgoldenrod2", high = "dodgerblue3",
  #                    midpoint= 80, na.value = "white",
  #                  name="Heterozygous SNPs\n per 5kb window") +
  scale_fill_gradient(na.value = "white", high = snp_high, low = snp_low,
                      name="Total SNPs\n per 5kb window") +
  geom_tile() +
  theme_minimal() +
  geom_rect(data=chrs, aes(group=CHROM, xmin=border_start, xmax=border_stop, ymin=0.5, ymax=Inf),
            fill=NA, inherit.aes=FALSE,colour = chrom_outline_color, linejoin = "round") +
  scale_x_continuous(expand = c(0,0),
                     breaks = chrs$tick,
                     labels = chr_ids)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=6.5),
        legend.title = element_text(size=10)
  )

ggsave(paste0(save_dir,Sys.Date(),"_MEC_Calbicans_all_SNP_heatmap.png"), p, device=png, dpi=300, bg="white",
       width = 10, height = 7.5, units="in")

