## ---------------------------
## Purpose: Generate a heatmap of all SNPs for a group of isolates, across all chromosomes
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
# Input files

genotype_table <-"/home/selmecki/shared/2022_Cglabrata_Erayil/align_variants_cgd_cbs138/bwa/clustering/Cglabrata_MEC_snps.table"
sample_list <- "~/umn/data/metadata/Cglabrata_MEC_raxml_unrooted_tips.csv"
feature_file <- "~/umn/Candida_genome_visualization/ref_genome_files/Cglabrata_CGD_s05m03r02_features.txt"
label_file <- "~/umn/Candida_genome_visualization/ref_genome_files/Cglabrata_CGD_s05m03r02_chr_labels.txt"

# Load packages
library(tidyverse)
library(ggplot2)
library(writexl)
library(readxl)

# Set variables
window <- 5000 # size of window used for rolling mean and snp density

save_dir <- "~/umn/images/Cglabrata/"

ref <- "CBS138" # label for file name or "" to leave out

mito <- "mito_C_glabrata_CBS138"

ref_allele <- "0"

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

features <- read_tsv(feature_file, show_col_types = FALSE)

# Read in genotypes (generated from filtered VCF)
genotypes <- read.table(genotype_table,
                        header = TRUE,
                        sep="\t")

gt <- as.data.frame(genotypes)
gt <- gt %>% filter(CHROM != mito)
gt <- gt %>% mutate(snp_bin=(POS %/% window) * window +1)

snp_counts <- gt %>%
  group_by(CHROM, snp_bin) %>%
  summarise_at(vars(colnames(gt)[3]:colnames(gt)[ncol(gt)-1]), function(x) sum(x!=ref_allele))

snp_counts$x_pos <- seq.int(nrow(snp_counts))

chrs <- snp_counts %>%
  group_by(CHROM) %>%
  summarise(border_start=min(x_pos),
            border_stop=max(x_pos),
            tick=min(x_pos) + (max(x_pos)-min(x_pos))/2)

tidy_snp_counts <- snp_counts %>%
  pivot_longer(names_to = "sample", values_to = "snp_count", cols=-c(CHROM,snp_bin, x_pos))

################################################################################
# Plot heatmap
p <- tidy_snp_counts %>%
  mutate(clustered_samples = fct_relevel(sample, rev(sample_order$sample))) %>%
  ggplot(aes(x=x_pos, y=clustered_samples, fill=snp_count))+
  #scale_fill_gradient2(low = "white", mid="darkgoldenrod2", high = "dodgerblue3",
  #                    midpoint= 80, na.value = "white",
  #                  name="Heterozygous SNPs\n per 5kb window") +
  scale_fill_gradient(na.value = "white", high = snp_high, low = snp_low,
                      name="SNPs per\n 5kb window") +
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

ggsave(paste0(save_dir,Sys.Date(),"_MEC_Cglabrata_all_SNP_heatmap.png"), p, device=png, dpi=300, bg="white",
       width = 10, height = 7.5, units="in")

